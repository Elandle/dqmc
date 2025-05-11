module checkerboard_mod
    use iso_fortran_env, only: real64
    implicit none

    integer, parameter, private :: dp = real64
    !
    ! Outline of the checkerboard method
    !
    ! Main idea:
    !
    ! To compute updates A = exp(T) * A or A = A * exp(T)
    ! using fast (non-dense) matrix-matrix multiplication approximately
    ! using the fact that T is sparse.
    !
    ! Theory:
    !
    ! To do this, T is first split:
    !
    ! T = diag(T) + C(1) + ... + C(n)
    !
    ! Then exp(T) is approximated using the Suzuki-Trotter approximation:
    !
    ! exp(T) ~ exp(diag(T)) * exp(C(1)) * ... * exp(C(n))
    !
    ! By choosing each C(i) to be strictly sparse, multiplications by
    ! exp(C(i)) can be done fast (left multiplication: row updates, right
    ! multiplication: column updates).
    !
    ! The Suzuki-Trotter approximation gets worse the higher n is.
    ! A minimal n can be found and used by edge-colouring T.
    ! 
    ! A symmetric pair of T is a pair (a, b) along with a pair of indices
    ! i, j, such that:
    !
    ! T(i, j) = a
    ! T(j, i) = b
    !       i < j (a is above the main diagonal and b is below)
    ! (a, b) =/= (0, 0) (one can still be 0)
    !
    ! If A is a matrix with all entries 0 except for:
    !
    ! A(i, j) = a
    ! A(j, i) = b 
    !
    ! then it can be shown that multiplication by exp(A) amounts to performing either
    ! two row updates or two column updates. For a symmetric pair (a, b), call the matrix
    ! A formed this way the matrix of the symmetric pair (a, b).
    !
    ! Assuming T has no diagonal (if it did, then approximate exp(T) ~ exp(diag(T)) * exp(offdiag(T)))
    ! form a graph from T in the following way.
    !
    ! Vertices: the indices 1, ..., n (T is n x n)
    ! Edges   : i and j are connected by an edge if they form the indices of a symmetric pair of T
    !
    ! Edge colour this graph. This means for each edge (i, j), give it a colour so that no two incident
    ! edges to a vertex share the same colour. The lesser the amount of colours used, the better
    ! the approximation the checkerboard method will produce.
    !
    ! After edge colouring, each edge (i, j) corresponds to a colour (in other words, label).
    ! Each edge (i, j) corresponds to a symmetric pair of T.
    ! For each colour k, take all of these symmetric pairs (that correspond to a coloured edge), and
    ! form the matrix:
    !
    ! C(k) = A(sympair(1)) + ... A(sympair(m))
    !
    ! In other words, start with a matrix C(k) of all 0's and fill in its entries according to the rules
    ! of making a symmetric pair matrix for each symmetric pair with colour k.
    !
    ! It can be shown that multiplication by exp(C(k)) is the same as multiplication by the exponential
    ! of each symmetric pair's matrix, which is just 2 row or column updates each.
    ! In other words, to multiply by exp(C(k)), iterate through multiplication by exp(A(l)), where A(l)
    ! is the matrix of the lth symmetric pair (in colour k).
    !
    ! For best results, T should be as sparse as possible, and the largest magnitude entry
    ! of T should be less than 1 (lower is better).
    !
    ! Implementation:
    !
    ! The base data structure is a symmetric pair, which holds the information for multiplying by
    ! a symmetric matrix (done by the right_symmult and left_symmult subroutines). construct_sympair
    ! sets this information up using the pair (a, b) and their indices (i, j).
    !
    ! sympair's are stored in an array in the ckbcolour datatype, which should store all 
    ! sympair's corresponding to a single colour. The subroutines left_colourmult and right_colourmult
    ! handle multiplication by all sympair's in a single colour by iterating through their multiplication.
    !
    ! The checkerboard datatype holds an array of colours (corresponding to the matrices a matrix is split into)
    ! and handles approximate multiplication (left_ckbmult and right_ckbmult) of a matrix exponential by the
    ! checkerboard method.
    !
    ! Note: this module is coded assuming the matrix T being exponentiated by has no diagonal.
    ! This is what comes up in DQMC, and if diagonals are explicitly needed, they can either be:
    !
    ! first extracted D = diag(T), T updated T = offdiag(T), checkerboard applied to the updated T,
    ! then multiplication by the diagonal matrix D implemented
    !
    ! or this extraction process can easily be added to the checkerboard method code
    !
    type :: checkerboard
        !
        ! Stores ckbcolour's for checkerboard multiplication by iterating through
        ! multiplying by different colours
        !
        integer                      :: n
        type(ckbcolour), allocatable :: colours(:)

    endtype checkerboard


    type :: ckbcolour
        !
        ! Stores sympair's of the same colour together for multliplication
        ! at the same time
        !
        integer                    :: n
        type(sympair), allocatable :: pairs(:)

    endtype ckbcolour


    type :: sympair
        !
        ! Stores the matrix-exponential information of a symmetric pair
        ! Note that i < j for a symmetric pair
        !
        integer  :: i, j
        real(dp) :: d, ij, ji

    endtype sympair


    contains


        subroutine construct_sympair(sym, i, j, a, b)
            !
            ! Constructs a sympair
            !
            !
            ! If a * b > 0:
            ! 
            !     0  a       cosh(x)    a*sinh(x)/x     d   ij
            ! exp b  0  =  b*sinh(x)/x    cosh(x)    =  ji   d
            !
            ! x = sqrt(a * b)
            !
            ! If a * b < 0:
            !
            !     0  a       cos(x)    a*sin(x)/x     d   ij
            ! exp b  0  =  b*sin(x)/x    cos(x)    =  ji   d
            !
            ! x = sqrt(-a * b)
            !
            ! If b = 0:
            !
            !     0  a     1  a     d   ij
            ! exp 0  0  =  0  1  =  ji   d
            !
            ! If a = 0:
            !
            !     0  0     1  0     d   ij
            ! exp b  0  =  b  1  =  ji   d
            !
            ! Assumption: a and b are both not 0
            !
            ! The cases a = 0 or b = 0 are assumed to be rare,
            ! so no special multiplication subroutine or sympair
            ! storage type are created for them. If this does not
            ! happen to be rare, then creating ones might be a
            ! worthy endeavor.
            !
            type(sympair), intent(out) :: sym
            integer      , intent(in)  :: i
            integer      , intent(in)  :: j
            real(dp)     , intent(in)  :: a
            real(dp)     , intent(in)  :: b

            real(dp) :: x

            sym%i = i
            sym%j = j

            if     (a * b .gt. 0.0_dp) then
                x      = sqrt(a * b)
                sym%d  = cosh(x)
                sym%ij = a * sinh(x) / x
                sym%ji = b * sinh(x) / x
            elseif (a * b .lt. 0.0_dp) then
                x      = sqrt(-a * b)
                sym%d  = cos(x)
                sym%ij = a * sin(x) / x
                sym%ji = b * sin(x) / x
            elseif (    b .eq. 0.0_dp) then
                sym%d  = 1.0_dp
                sym%ij = a
                sym%ji = 0.0_dp
            elseif (a     .eq. 0.0_dp) then
                sym%d  = 1.0_dp
                sym%ij = 0.0_dp
                sym%ji = b
            endif


        endsubroutine construct_sympair


        subroutine right_symmult(sym, A, n, work)
            !
            ! Updates:
            !
            ! A = A * exp(sym)
            !
            ! By performing the updates:
            !
            ! A(:, i) = d * A(:, i) + ji * A(:, j)
            ! A(:, j) = d * A(:, j) + ij * A(:, i)
            !
            ! exp(sym) is the matrix exponential whose information is stored in sym
            !
            type(sympair), intent(in)    :: sym
            real(dp)     , intent(inout) :: A(n, n)
            integer      , intent(in)    :: n
            real(dp)     , intent(out)   :: work(n)
            
            ! work = A(:, i)
            call dcopy(n, A(1, sym%i), 1, work, 1)
            ! A(:, i) = d * A(:, i)
            call dscal(n, sym%d, A(1, sym%i), 1)
            ! A(:, i) = A(:, i) + ji * A(:, j)
            call daxpy(n, sym%ji, A(1, sym%j), 1, A(1, sym%i), 1)
            ! A(:, j) = d * A(:, j)
            call dscal(n, sym%d, A(1, sym%j), 1)
            ! A(:, j) = A(:, j) + ij * work
            call daxpy(n, sym%ij, work, 1, A(1, sym%j), 1)

            ! No BLAS:
            ! work        = A(:, sym%i)
            ! A(:, sym%i) = sym%d * A(:, sym%i) + sym%ji * A(:, sym%j)
            ! A(:, sym%j) = sym%d * A(:, sym%j) + sym%ij * work

        endsubroutine right_symmult


        subroutine left_symmult(sym, A, n, work)
            !
            ! Updates:
            !
            ! A = symexp * A
            !
            ! By performing the updates:
            !
            ! A(i, :) = d * A(i, :) + ij * A(j, :)
            ! A(j, :) = d * A(j, :) + ji * A(i, :)
            !
            ! symexp is the matrix exponential whose information is stored in sym
            !
            type(sympair), intent(in)    :: sym
            real(dp)     , intent(inout) :: A(n, n)
            integer      , intent(in)    :: n
            real(dp)     , intent(out)   :: work(n)
            
            ! work = A(i, :)
            call dcopy(n, A(sym%i, 1), n, work, 1)
            ! A(i, :) = d * A(i, :)
            call dscal(n, sym%d, A(sym%i, 1), n)
            ! A(i, :) = A(i, :) + ij * A(j, :)
            call daxpy(n, sym%ij, A(sym%j, 1), n, A(sym%i, 1), n)
            ! A(j, :) = d * A(j, :)
            call dscal(n, sym%d, A(sym%j, 1), n)
            ! A(j, :) = A(j, :) + ji * work
            call daxpy(n, sym%ji, work, 1, A(sym%j, 1), n)

            ! No BLAS:
            ! work        =         A(sym%i, :)
            ! A(sym%i, :) = sym%d * A(sym%i, :) + sym%ij * A(sym%j, :)
             !A(sym%j, :) = sym%d * A(sym%j, :) + sym%ji * work


        endsubroutine left_symmult
        

        subroutine initialize_ckbcolour(colour, n)
            !
            ! Gets a ckbcolour ready to store n sympairs
            !
            type(ckbcolour), intent(out) :: colour
            integer        , intent(in)  :: n

            colour%n = n
            allocate(colour%pairs(n))

        endsubroutine initialize_ckbcolour


        subroutine right_colourmult(colour, A, n, work)
            !
            ! Updates:
            !
            ! A = A * exp(colour)
            !
            ! Where exp(colour) is the matrix exponential stored in colour,
            ! done by iterating through sympair exponentials
            !
            ! In other words, if colour stores the information of sympair1, ..., sympairn,
            ! then A is updated according to:
            !
            ! A = A * exp(sympair1) * ... * exp(sympairn)
            !
            type(ckbcolour), intent(in)    :: colour
            real(dp)       , intent(inout) :: A(n, n)
            integer        , intent(in)    :: n
            real(dp)       , intent(out)   :: work(n)

            integer :: i

            do i = 1, colour%n
                call right_symmult(colour%pairs(i), A, n, work)
            enddo

        endsubroutine right_colourmult


        subroutine left_colourmult(colour, A, n, work)
            !
            ! Updates:
            !
            ! A = exp(colour) * A
            !
            ! Where exp(colour) is the matrix exponential stored in colour,
            ! done by iterating through sympair exponentials
            !
            ! In other words, if colour stores the information of sympair1, ..., sympairn,
            ! then A is updated according to:
            !
            ! A = exp(sympair1) * ... * exp(sympairn) * A
            !
            type(ckbcolour), intent(in)    :: colour
            real(dp)       , intent(inout) :: A(n, n)
            integer        , intent(in)    :: n
            real(dp)       , intent(out)   :: work(n)

            integer :: i

            do i = colour%n, 1, -1
                call left_symmult(colour%pairs(i), A, n, work)
            enddo

        endsubroutine left_colourmult


        subroutine initialize_checkerboard(ckb, n)
            !
            ! Gets a checkerboard ready to store n ckbcolours
            !
            type(checkerboard), intent(out) :: ckb
            integer           , intent(in)  :: n

            ckb%n = n
            allocate(ckb%colours(n))

        endsubroutine initialize_checkerboard


        subroutine right_ckbmult(ckb, A, n, work)
            !
            ! Updates:
            !
            ! A = A * exp(ckb)
            !
            ! Where exp(ckb) is the matrix exponential stored in ckb
            !
            ! Where exp(ckb) is the matrix exponential stored in ckb,
            ! done by iterating through colour exponentials
            !
            ! In other words, if ckb stores the information of colour1, ..., colourn,
            ! then A is updated according to:
            !
            ! A = A * exp(colour1) * ... * exp(colourn)
            !
            ! In column-major languages (such as Fortran and Julia), right multiplication by the
            ! checkerboard method is faster than left multiplication.
            !
            type(checkerboard), intent(in)    :: ckb
            real(dp)          , intent(inout) :: A(n, n)
            integer           , intent(in)    :: n
            real(dp)          , intent(out)   :: work(n)

            integer :: i

            do i = 1, ckb%n
                call right_colourmult(ckb%colours(i), A, n, work)
            enddo

        endsubroutine right_ckbmult


        subroutine left_ckbmult(ckb, A, n, work)
            !
            ! Updates:
            !
            ! A = A * exp(ckb)
            !
            ! Where exp(ckb) is the matrix exponential stored in ckb
            !
            ! Where exp(ckb) is the matrix exponential stored in ckb,
            ! done by iterating through colour exponentials
            !
            ! In other words, if ckb stores the information of colour1, ..., colourn,
            ! then A is updated according to:
            !
            ! A = exp(colour1) * ... * exp(colourn) * A
            !
            ! In row-major languages (such as C and C++), left multiplication by the checkerboard
            ! method is faster than right multiplication.
            !
            type(checkerboard), intent(in)    :: ckb
            real(dp)          , intent(inout) :: A(n, n)
            integer           , intent(in)    :: n
            real(dp)          , intent(out)   :: work(n)

            integer :: i

            do i = ckb%n, 1, -1
                call left_colourmult(ckb%colours(i), A, n, work)
            enddo

        endsubroutine left_ckbmult


        subroutine read_ckb(ckb, filename, iounit, dtau)
            !
            ! Reads in an input checkerboard file and constructs its checkerboard
            !
            ! The input file must have the following form:
            !
            ! n_colors
            !
            ! nsym1
            ! i1 j1 T(i1, j1) T(j1, i1)
            ! ...
            ! insym1 jnsym1 T(insym1, jnsym1) T(jnsym1, insym1)
            !
            ! nsym2
            ! i1 j1 T(i1, j1) T(j1, i1)
            ! ...
            ! insym2 jnsym2 T(insym2, jnsym2) T(jnsym2, insym2)
            !
            ! ...
            !
            ! nsymn_colors
            ! i1 j1 T(i1, j1) T(j1, i1)
            ! ...
            ! insymn_colors jnsymn_colors T(insymn_colors, jnsymn_colors) T(jnsymn_colors, insymn_colors)
            !
            ! In other words:
            !
            ! line 1 is the number of colors T is colored by
            ! line 2 is blank
            ! line 3 is the amount of symmetric pairs in color 1
            ! for each symmmetric pair in color 1, printed on a single line is the following:
            ! its i index, its j index, the matrix element T(i, j), the matrix element T(j, i)
            ! where each entry is separated by a single space.
            ! Once all symmetric pair information is printed for color 1, skip a line
            ! Repeat for color 2
            ! ...
            ! Repeat for the last color
            ! No blank line at the end
            !
            type(checkerboard), intent(out)           :: ckb
            real(dp)          , intent(in), optional  :: dtau
            character(len=*)  , intent(in)            :: filename
            integer           , intent(out)           :: iounit

            integer  :: k, l
            integer  :: n_colors
            integer  :: n_col
            integer  :: i, j
            real(dp) :: ij, ji
            character(len=100) :: str

            real(dp) :: scale

            ! TODO:
            ! Implement input error checking
            ! Get better at strings and rewrite (very brute force and basic here)

            if (present(dtau)) then
                scale = dtau
            else
                scale = 1.0_dp
            endif

            ! Open the input file and assign it a unit  (should have no unit beforehand)
            open(file=filename, newunit=iounit)

            ! Read line 1: the number of colors
            read(iounit, "(a100)") str
            read(str, *) n_colors
            call initialize_checkerboard(ckb, n_colors)
        
            ! Iterate reading in the different colors
            do k = 1, n_colors
                read(iounit, "(a100)") str ! This should be a blank line
                ! After each blank is the number of symmetric pairs in a color
                read(iounit, "(a100)") str
                read(str, *) n_col
                call initialize_ckbcolour(ckb%colours(k), n_col)
                ! Iterate through symmetric pairs
                do l = 1, n_col
                    read(iounit, "(a100)") str
                    read(str, *) i, j, ij, ji
                    call construct_sympair(ckb%colours(k)%pairs(l), i, j, scale * ij, scale * ji)
                enddo
            enddo


        endsubroutine read_ckb


        subroutine read_ckbT(T, N, filename, iounit, dtau)
            real(dp)          , intent(out)           :: T(N, N)
            integer           , intent(in)            :: N
            real(dp)          , intent(in) , optional :: dtau
            character(len=*)  , intent(in)            :: filename
            integer           , intent(out)           :: iounit

            integer  :: k, l
            integer  :: n_colors
            integer  :: n_col
            integer  :: i, j
            real(dp) :: ij, ji
            character(len=100) :: str

            real(dp) :: scale

            if (present(dtau)) then
                scale = dtau
            else
                scale = 1.0_dp
            endif

            T = 0

            ! Open the input file and assign it a unit  (should have no unit beforehand)
            open(file=filename, newunit=iounit)

            ! Read line 1: the number of colors
            read(iounit, "(a100)") str
            read(str, *) n_colors

            ! Iterate reading in the different colors
            do k = 1, n_colors
                read(iounit, "(a100)") str ! This should be a blank line
                ! After each blank is the number of symmetric pairs in a color
                read(iounit, "(a100)") str
                read(str, *) n_col
                ! Iterate through symmetric pairs
                do l = 1, n_col
                    read(iounit, "(a100)") str
                    read(str, *) i, j, ij, ji
                    T(i, j) = scale * ij
                    T(j, i) = scale * ji
                    write(6, "(a, i4, a, i4, a, f8.4)") "T(", i, ", ", j, ") = ", scale * ij
                    write(6, "(a, i4, a, i4, a, f8.4)") "T(", j, ", ", i, ") = ", scale * ji
                enddo
            enddo




            do i = 1, N
                do j = 1, N
                    write(6, "(f17.8)", advance="no") T(i, j)
                enddo
                write(6, "(a)") ""
            enddo
        

        endsubroutine read_ckbT



endmodule checkerboard_mod