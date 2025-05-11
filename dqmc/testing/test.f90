program main
    use iso_fortran_env, only: dp => real64, terminal => output_unit
    use checkerboard_mod
    use expm_mod
    implicit none

    interface
        function dlange(norm, m, n, A, lda, work) result(res)
            import :: dp
            character(len=*), intent(in) :: norm
            integer, intent(in) :: m, n, lda
            real(dp), intent(in) :: A(lda, *)
            real(dp), intent(out) :: work(*)
            real(dp) :: res
        end function dlange
    end interface

    

    integer, parameter :: n = 8
    type(checkerboard) :: ckb
    integer            :: iounit, i, j
    real(dp)           :: d, ij, ji
    real(dp)           :: dtau
    real(dp)           :: A(n, n), expA(n, n), ckbA(n, n), R(n, n), work(n)
    real(dp)           :: matrix(n, n)
    real(dp)           :: col1m(n, n), col2m(n, n), col3m(n, n), col4m(n, n)


    dtau = 0.08_dp
    call read_ckb(ckb, "delta02.txt", iounit, dtau)
    call read_ckbT(A, n, "delta02.txt", iounit, dtau)
    ckbA = 0.0_dp
    do i = 1, n
        ckbA(i, i) = 1.0_dp
    enddo

    call right_ckbmult(ckb, ckbA, n, work)
    call expm(A, expA)
    write(terminal, "(a)") "A = "
    call print_matrix(A, terminal)
    write(terminal, "(a)") "ckbA = "
    call print_matrix(ckbA, terminal)
    write(terminal, "(a)") "expA = "
    call print_matrix(expA, terminal)
    R = expA - ckbA
    write(terminal, "(a)") "expA - ckbA = "
    call print_matrix(R, terminal)
    write(terminal, "(a)") "expA - ckbA norms:"
    write(terminal, "(a, f17.8)") "maxabs entry   = ", dmnorm(R, 'm')
    write(terminal, "(a, f17.8)") "one       norm = ", dmnorm(R, 'o')
    write(terminal, "(a, f17.8)") "infinity  norm = ", dmnorm(R, 'i')
    write(terminal, "(a, f17.8)") "frobenius norm = ", dmnorm(R, 'f')
    write(terminal, "(a, f17.8)") "two       norm = ", dmnorm(R, 't')
    write(terminal, "(a, f17.8)") "dtau   = ", dtau
    write(terminal, "(a, f17.8)") "dtau^2 = ", dtau * dtau



    write(terminal, "(a, i4)") "ckb%n            = ", ckb%n            ! 4

    write(terminal, "(a, i4)") "ckb%colours(1)%n = ", ckb%colours(1)%n ! 4
    write(terminal, "(a, i4)") "ckb%colours(2)%n = ", ckb%colours(2)%n ! 3
    write(terminal, "(a, i4)") "ckb%colours(3)%n = ", ckb%colours(3)%n ! 3
    write(terminal, "(a, i4)") "ckb%colours(4)%n = ", ckb%colours(4)%n ! 2

    write(terminal, "(a)") "ckb%colours(1)%pairs(1) attributes:"
    i  = ckb%colours(1)%pairs(1)%i
    j  = ckb%colours(1)%pairs(1)%j
    d  = ckb%colours(1)%pairs(1)%d
    ij = ckb%colours(1)%pairs(1)%ij
    ji = ckb%colours(1)%pairs(1)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(1)%pairs(2) attributes:"
    i  = ckb%colours(1)%pairs(2)%i
    j  = ckb%colours(1)%pairs(2)%j
    d  = ckb%colours(1)%pairs(2)%d
    ij = ckb%colours(1)%pairs(2)%ij
    ji = ckb%colours(1)%pairs(2)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(1)%pairs(3) attributes:"
    i  = ckb%colours(1)%pairs(3)%i
    j  = ckb%colours(1)%pairs(3)%j
    d  = ckb%colours(1)%pairs(3)%d
    ij = ckb%colours(1)%pairs(3)%ij
    ji = ckb%colours(1)%pairs(3)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    i  = ckb%colours(1)%pairs(4)%i
    j  = ckb%colours(1)%pairs(4)%j
    d  = ckb%colours(1)%pairs(4)%d
    ij = ckb%colours(1)%pairs(4)%ij
    ji = ckb%colours(1)%pairs(4)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(1) * id = "
    matrix = id(n)
    call left_colourmult(ckb%colours(1), matrix, n, work)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "nonzero entries:"
    call print_matrix_nonzero(matrix, terminal)
    write(terminal, "(a)") "id * ckb%colours(1) = "
    matrix = id(n)
    call right_colourmult(ckb%colours(1), matrix, n, work)
    call print_matrix(matrix, terminal)
    col1m = matrix
    write(terminal, "(a)") "nonzero entries:"
    call print_matrix_nonzero(matrix, terminal)

    write(terminal, "(a)") "ckb%colours(2)%pairs(1) attributes:"
    i  = ckb%colours(2)%pairs(1)%i
    j  = ckb%colours(2)%pairs(1)%j
    d  = ckb%colours(2)%pairs(1)%d
    ij = ckb%colours(2)%pairs(1)%ij
    ji = ckb%colours(2)%pairs(1)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(2)%pairs(2) attributes:"
    i  = ckb%colours(2)%pairs(2)%i
    j  = ckb%colours(2)%pairs(2)%j
    d  = ckb%colours(2)%pairs(2)%d
    ij = ckb%colours(2)%pairs(2)%ij
    ji = ckb%colours(2)%pairs(2)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(2)%pairs(3) attributes:"
    i  = ckb%colours(2)%pairs(3)%i
    j  = ckb%colours(2)%pairs(3)%j
    d  = ckb%colours(2)%pairs(3)%d
    ij = ckb%colours(2)%pairs(3)%ij
    ji = ckb%colours(2)%pairs(3)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(2) * id = "
    matrix = id(n)
    call left_colourmult(ckb%colours(2), matrix, n, work)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "nonzero entries:"
    call print_matrix_nonzero(matrix, terminal)
    write(terminal, "(a)") "id * ckb%colours(2) = "
    matrix = id(n)
    call right_colourmult(ckb%colours(2), matrix, n, work)
    call print_matrix(matrix, terminal)
    col2m = matrix
    write(terminal, "(a)") "nonzero entries:"
    call print_matrix_nonzero(matrix, terminal)

    write(terminal, "(a)") "ckb%colours(3)%pairs(1) attributes:"
    i  = ckb%colours(3)%pairs(1)%i
    j  = ckb%colours(3)%pairs(1)%j
    d  = ckb%colours(3)%pairs(1)%d
    ij = ckb%colours(3)%pairs(1)%ij
    ji = ckb%colours(3)%pairs(1)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(3)%pairs(2) attributes:"
    i  = ckb%colours(3)%pairs(2)%i
    j  = ckb%colours(3)%pairs(2)%j
    d  = ckb%colours(3)%pairs(2)%d
    ij = ckb%colours(3)%pairs(2)%ij
    ji = ckb%colours(3)%pairs(2)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(3)%pairs(3) attributes:"
    i  = ckb%colours(3)%pairs(3)%i
    j  = ckb%colours(3)%pairs(3)%j
    d  = ckb%colours(3)%pairs(3)%d
    ij = ckb%colours(3)%pairs(3)%ij
    ji = ckb%colours(3)%pairs(3)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(3) * id = "
    matrix = id(n)
    call left_colourmult(ckb%colours(3), matrix, n, work)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "nonzero entries:"
    call print_matrix_nonzero(matrix, terminal)
    write(terminal, "(a)") "id * ckb%colours(3) = "
    matrix = id(n)
    call right_colourmult(ckb%colours(3), matrix, n, work)
    call print_matrix(matrix, terminal)
    col3m = matrix
    write(terminal, "(a)") "nonzero entries:"
    call print_matrix_nonzero(matrix, terminal)

    write(terminal, "(a)") "ckb%colours(4)%pairs(1) attributes:"
    i  = ckb%colours(4)%pairs(1)%i
    j  = ckb%colours(4)%pairs(1)%j
    d  = ckb%colours(4)%pairs(1)%d
    ij = ckb%colours(4)%pairs(1)%ij
    ji = ckb%colours(4)%pairs(1)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(4)%pairs(2) attributes:"
    i  = ckb%colours(4)%pairs(2)%i
    j  = ckb%colours(4)%pairs(2)%j
    d  = ckb%colours(4)%pairs(2)%d
    ij = ckb%colours(4)%pairs(2)%ij
    ji = ckb%colours(4)%pairs(2)%ji
    write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
    write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
    write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
    write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
    write(terminal, "(a)") "ckb%colours(4) * id = "
    matrix = id(n)
    call left_colourmult(ckb%colours(4), matrix, n, work)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "nonzero entries:"
    call print_matrix_nonzero(matrix, terminal)
    write(terminal, "(a)") "id * ckb%colours(4) = "
    matrix = id(n)
    call right_colourmult(ckb%colours(4), matrix, n, work)
    call print_matrix(matrix, terminal)
    col4m = matrix
    write(terminal, "(a)") "nonzero entries:"
    call print_matrix_nonzero(matrix, terminal)
    ! write(terminal, "(a)") "inverse:"
    ! call print_matrix(inverse(A), terminal)
    
    write(terminal, "(a)") "colour1 * colour2 * colour3 * colour4 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col1m, col2m), col3m), col4m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour1 * colour2 * colour3 * colour4 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour1 * colour2 * colour4 * colour3 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col1m, col2m), col4m), col3m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour1 * colour2 * colour4 * colour3 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour1 * colour3 * colour2 * colour4 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col1m, col3m), col2m), col4m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour1 * colour3 * colour2 * colour4 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour1 * colour3 * colour4 * colour2 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col1m, col3m), col4m), col2m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour1 * colour3 * colour4 * colour2 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour1 * colour4 * colour2 * colour3 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col1m, col4m), col2m), col3m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour1 * colour4 * colour2 * colour3 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour1 * colour4 * colour3 * colour2 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col1m, col4m), col3m), col2m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour1 * colour4 * colour3 * colour2 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour2 * colour1 * colour3 * colour4 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col2m, col1m), col3m), col4m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour2 * colour1 * colour3 * colour4 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour2 * colour1 * colour4 * colour3 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col2m, col1m), col4m), col3m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour2 * colour1 * colour4 * colour3 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour2 * colour3 * colour1 * colour4 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col2m, col3m), col1m), col4m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour2 * colour3 * colour1 * colour4 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour2 * colour3 * colour4 * colour1 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col2m, col3m), col4m), col1m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour2 * colour3 * colour4 * colour1 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour2 * colour4 * colour1 * colour3 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col2m, col4m), col1m), col3m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour2 * colour4 * colour1 * colour3 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour2 * colour4 * colour3 * colour1 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col2m, col4m), col3m), col1m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour2 * colour4 * colour3 * colour1 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour3 * colour1 * colour2 * colour4 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col3m, col1m), col2m), col4m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour3 * colour1 * colour2 * colour4 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour3 * colour1 * colour4 * colour2 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col3m, col1m), col4m), col2m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour3 * colour1 * colour4 * colour2 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour3 * colour2 * colour1 * colour4 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col3m, col2m), col1m), col4m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour3 * colour2 * colour1 * colour4 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour3 * colour2 * colour4 * colour1 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col3m, col2m), col4m), col1m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour3 * colour2 * colour4 * colour1 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour3 * colour4 * colour1 * colour2 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col3m, col4m), col1m), col2m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour3 * colour4 * colour1 * colour2 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour3 * colour4 * colour2 * colour1 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col3m, col4m), col2m), col1m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour3 * colour4 * colour2 * colour1 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour4 * colour1 * colour2 * colour3 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col4m, col1m), col2m), col3m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour4 * colour1 * colour2 * colour3 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour4 * colour1 * colour3 * colour2 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col4m, col1m), col3m), col2m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour4 * colour1 * colour3 * colour2 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour4 * colour2 * colour1 * colour3 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col4m, col2m), col1m), col3m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour4 * colour2 * colour1 * colour3 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour4 * colour2 * colour3 * colour1 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col4m, col2m), col3m), col1m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour4 * colour2 * colour3 * colour1 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour4 * colour3 * colour1 * colour2 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col4m, col3m), col1m), col2m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour4 * colour3 * colour1 * colour2 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    
    write(terminal, "(a)") "colour4 * colour3 * colour2 * colour1 (dense multiplication) = "
    matrix = matmul(matmul(matmul(col4m, col3m), col2m), col1m)
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour4 * colour3 * colour2 * colour1 (dense multiplication) = "
    call print_matrix(ckbA - matrix, terminal)
    





    write(terminal, "(a)") "colour1 * colour2 - colour2 * colour1 (dense) = "
    call print_matrix(matmul(col1m, col2m) - matmul(col2m, col1m), terminal)
    write(terminal, "(a)") "colour1 * colour3 - colour3 * colour1 (dense) = "
    call print_matrix(matmul(col1m, col3m) - matmul(col3m, col1m), terminal)
    write(terminal, "(a)") "colour1 * colour4 - colour4 * colour1 (dense) = "
    call print_matrix(matmul(col1m, col4m) - matmul(col4m, col1m), terminal)
    write(terminal, "(a)") "colour2 * colour3 - colour3 * colour2 (dense) = "
    call print_matrix(matmul(col2m, col3m) - matmul(col3m, col2m), terminal)
    write(terminal, "(a)") "colour2 * colour4 - colour4 * colour2 (dense) = "
    call print_matrix(matmul(col2m, col4m) - matmul(col4m, col2m), terminal)
    write(terminal, "(a)") "colour3 * colour4 - colour4 * colour3 (dense) = "
    call print_matrix(matmul(col3m, col4m) - matmul(col4m, col3m), terminal)


    contains

        function inverse(A) result(invA)
            real(dp), intent(in)  :: A(:, :)

            real(dp) :: invA(size(A, 1), size(A, 2))
            integer  :: P(size(A, 1))
            real(dp) :: work(64 * size(A, 1))
            integer  :: info

            call dcopy(size(A, 1)*size(A, 1), A, 1, invA, 1)
            call dgetrf(size(invA, 1), size(invA, 1), invA, size(invA, 1), P, info)
            call dgetri(size(invA, 1), invA, size(invA, 1), P, work, 4*size(invA, 1), info)


        endfunction inverse

        real(dp) function dmnorm(A, norm)
            !
            ! max(abs(A(i,j))), NORM = 'M' or 'm'
            ! norm1(A),         NORM = '1', 'O' or 'o'
            ! normI(A),         NORM = 'I' or 'i'
            ! normF(A),         NORM = 'F', 'f', 'E' or 'e'
            ! norm2(A),         NORM = 'T', 't', or '2'
            ! where  norm1  denotes the  one norm of a matrix (maximum column sum),
            ! normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
            ! normF  denotes the  Frobenius norm of a matrix (square root of sum of
            ! squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
            !
            real(dp) , intent(in) :: A(:, :)
            character, intent(in) :: norm

            real(dp), allocatable :: work(:)
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            if ((norm .eq. 'I') .or. (norm .eq. 'i')) then
                allocate(work(m))
            endif

            if ((norm .eq. '2') .or. (norm .eq. 'T') .or. (norm .eq. 't')) then
                dmnorm = twonorm(A)
            else
                dmnorm = dlange(norm, m, n, A, m, work)
            endif

            if (allocated(work)) then
                deallocate(work)
            endif

        endfunction dmnorm

        real(dp) function twonorm(A) result(norm)
            !
            ! Returns the 2 norm of the matrix A
            !
            real(dp), intent(in)  :: A(:, :)
            real(dp), allocatable :: B(:, :)
            real(dp), allocatable :: S(:)
            real(dp), allocatable :: work(:)
            integer               :: m
            integer               :: n
            integer               :: lwork
            integer               :: info
            
            m = size(A, 1)
            n = size(A, 2)
            lwork = 5 * (3 * min(m, n) + max(m, n) + 5 * min(m, n))

            allocate(B(m, n))
            allocate(S(min(m, n)))
            allocate(work(lwork))

            call dlacpy('A', m, n, A, m, B, m)
            call dgesvd('N', 'N', m, n, B, m, S, work, lwork, work, lwork, work, lwork, info)

            norm = S(1)

            deallocate(B)
            deallocate(S)
            deallocate(work)

            
        endfunction twonorm




        subroutine print_matrix(A, ounit)
            real(dp), intent(in) :: A(:, :)
            integer , intent(in) :: ounit
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            do i = 1, n
                do j = 1, m
                    write(ounit, "(f17.8)", advance="no") A(i, j)
                enddo
                write(ounit, "(a)") ""
            enddo


        endsubroutine print_matrix


        function id(n)
            integer, intent(in) :: n
            
            real(dp) :: id(n, n)
            integer  :: i

            id = 0.0_dp
            do i = 1, n
                id(i, i) = 1.0_dp
            enddo

        endfunction id


        subroutine print_matrix_nonzero(A, ounit)
            real(dp), intent(in) :: A(:, :)
            integer , intent(in) :: ounit
            integer :: i, j
            do i = 1, size(A, 1)
                do j = 1, size(A, 2)
                    if (abs(A(i, j)) .ge. 10e-4) then
                        write(ounit, "(i5, a, i5, f17.8)") i, " , ", j, A(i, j)
                    endif
                enddo
            enddo
        endsubroutine print_matrix_nonzero




endprogram main