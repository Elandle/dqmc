module customlatests_mod
    use numbertypes
    use customla_mod
    use convenientla_mod

    ! Contains procedures for testing the procedures of customla_mod
    ! Since the purpose is testing, the procedures here were not designed for
    ! efficiency (eg, they allocate memory instead of asking for a workspace)

    contains


        subroutine test_twonorm()
            real(dp) :: A(2, 2)
            real(dp) :: norm
            real(dp) :: exact

            A(1, 1) = 1
            A(2, 1) = 2
            A(1, 2) = 1.5
            A(2, 2) = 2.5

            exact = 3.67171

            print *, "Testing twonorm..."
            norm = twonorm(A)
            print *, "Exact result = ", exact
            print *, "Calculated result = ", norm
            print *, "Absolute difference = ", abs(exact - norm)

        endsubroutine test_twonorm


        subroutine test_right_diagmult()
            real(dp) :: A(3, 3)
            real(dp) :: B(3, 3)
            real(dp) :: D(3)

            A(1, 1) = 5.5  ; A(1, 2) = -2.2    ; A(1, 3) = 67.4

            A(2, 1) = -1.0 ; A(2, 2) = 1000.24 ; A(2, 3) = -3.23

            A(3, 1) = 12   ; A(3, 2) = 12.01   ; A(3, 3) = 11.89

            D(1) = 41.3 ; D(2) = -0.89 ; D(3) = 6.2

            B(1, 1) = 227.15 ; B(1, 2) = 1.958 ; B(1, 3) = 417.88
            
            B(2, 1) = -41.3 ; B(2, 2) = -890.2136 ; B(2, 3) = -20.026

            B(3, 1) = 495.6 ; B(3, 2) = -10.6889 ; B(3, 3) = 73.718

            print *, "Testing right_diagmult..."
            call right_diagmult(A, D, 3)
            print *, "2 norm of difference between calculation and exact result = ", twonorm(A - B)


        endsubroutine test_right_diagmult


        subroutine test_left_diagmult()
            real(dp) :: A(3, 3)
            real(dp) :: B(3, 3)
            real(dp) :: D(3)

            A(1, 1) = 5.5  ; A(1, 2) = -2.2    ; A(1, 3) = 67.4

            A(2, 1) = -1.0 ; A(2, 2) = 1000.24 ; A(2, 3) = -3.23

            A(3, 1) = 12   ; A(3, 2) = 12.01   ; A(3, 3) = 11.89

            D(1) = 41.3 ; D(2) = -0.89 ; D(3) = 6.2

            B(1, 1) = 2.271500e+02 ; B(1, 2) = -9.086000e+01 ; B(1, 3) = 2.783620e+03
            
            B(2, 1) = 8.900000e-01 ; B(2, 2) = -8.902136e+02 ; B(2, 3) = 2.874700e+00

            B(3, 1) = 7.440000e+01 ; B(3, 2) = 7.446200e+01  ; B(3, 3) = 7.371800e+01

            print *, "Testing left_diagmult..."
            call left_diagmult(A, D, 3)
            print *, "2 norm of difference between calculation and exact result = ", twonorm(A - B)


        endsubroutine test_left_diagmult




endmodule customlatests_mod