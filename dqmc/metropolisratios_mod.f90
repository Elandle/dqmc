module metropolisratios_mod
    use numbertypes
    use simulationsetup_mod
    implicit none

    ! TODO:
    ! Create subroutines to calculate R according to the formula:
    !
    ! R = det(id + B(beta, 0))
    !   = det(id + B(L) * ... * B(1))
    !
    ! using the ASvQRD algorithm.
    !
    ! This subroutine is not meant to be used in the final code, but rather
    ! testing purposes (ie, make sure R is calculated correctly)

    contains


        subroutine greens_Rup(S, i, l)
            !
            ! Sets:
            !
            ! S%R     = 1 + (1 - G(i, i)) * delta
            ! S%delta = exp(-2 * sigma * alpha * h(i, l)) - 1
            !
            ! In this subroutine spin up (sigma=1) is hardcoded.
            ! In other words, the update:
            !
            ! S%Rup     = 1 + (1 - Gup(i, i)) * deltaup
            ! S%deltaup = exp(-2 * alpha * h(i, l)) - 1
            !
            ! is performed.
            !
            ! Note that this subroutine assumes that h(i, l) has not been flipped (since
            ! the purpose of this subroutine is to decide whether or not to flip h(i, l))
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i
            integer         , intent(in)    :: l

            S%deltaup = exp(-2 * S%alpha * S%h(i, l)) - 1
            S%Rup     = 1 + (1 - S%Gup(i, i)) * S%deltaup
            

        endsubroutine greens_Rup


        subroutine greens_Rdn(S, i, l)
            !
            ! Sets:
            !
            ! S%R     = 1 + (1 - G(i, i)) * delta
            ! S%delta = exp(-2 * sigma * alpha * h(i, l)) - 1
            !
            ! In this subroutine spin up (sigma=1) is hardcoded.
            ! In other words, the update:
            !
            ! S%Rdn     = 1 + (1 - Gdn(i, i)) * deltadn
            ! S%deltadn = exp(2 * alpha * h(i, l)) - 1
            !
            ! is performed.
            !
            ! Note that this subroutine assumes that h(i, l) has not been flipped (since
            ! the purpose of this subroutine is to decide whether or not to flip h(i, l))
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i
            integer         , intent(in)    :: l

            S%deltadn = exp(2 * S%alpha * S%h(i, l)) - 1
            S%Rdn     = 1 + (1 - S%Gdn(i, i)) * S%deltadn
            

        endsubroutine greens_Rdn




endmodule metropolisratios_mod