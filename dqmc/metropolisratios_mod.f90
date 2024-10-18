module metropolisratios_mod
    use numbertypes
    use simulationsetup_mod
    implicit none
    !
    ! Contains procedures for calculating the Metropolis ratio:
    !
    !     R = Rup * Rdn
    !
    ! that is used for determining whether or not a flip is accepted or not.
    !
    !
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
    !
    contains


        subroutine greens_R(S, i, l, sigma)
            !
            ! Sets:
            !
            ! S%R     = 1 + (1 - G(i, i)) * delta
            ! S%delta = exp(-2 * sigma * alpha * h(i, l)) - 1
            !
            ! Where S%R is either S%Rup or S%Rdn and S%delta is either
            ! S%deltaup or S%deltadn, depending on whether sigma = 1 (spin up)
            ! or sigma = -1 (spin dn).
            !
            ! Note that this subroutine assumes that h(i, l) has not been flipped (since
            ! the purpose of this subroutine is to decide whether or not to flip h(i, l))
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            if (sigma .eq. 1) then      ! up spin (sigma = 1)
                S%deltaup = exp(-2 * sigma * S%alpha * S%h(i, l)) - 1
                S%Rup     = 1 + (1 - S%Gup(i, i)) * S%deltaup
            else                        ! dn spin (sigma = -1)
                S%deltadn = exp(-2 * sigma * S%alpha * S%h(i, l)) - 1
                S%Rdn     = 1 + (1 - S%Gdn(i, i)) * S%deltadn
            endif


        endsubroutine greens_R


endmodule metropolisratios_mod