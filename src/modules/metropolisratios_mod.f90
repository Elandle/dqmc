!> \brief Contains procedures for calculating the Metropolis ratio \f$R=R_\uparrow R_\downarrow\f$
!! that is used for determining whether or not a flip in the Hubbard-Stratonovich field \f$h\f$ is accepted.
module metropolisratios_mod
    use numbertypes
    use simulationsetup_mod
    use equalgreens_mod
    implicit none


    contains


        !> \brief Computes the Metropolis ratio \f$R_\sigma\f$ for flipping the Hubbard-Stratonivich field \f$h(i, l)\f$.
        !!
        !! The Metropolis ratio for spin \f$\sigma\f$ when proposing to flip \f$h(i, l)\f$ is:
        !! \f[R_\sigma = 1 + \delta_\sigma(1 - G_\sigma(i, i))\f]
        !! where:
        !! \f[\delta_\sigma = \exp(-2\sigma\alpha h(i, l)) - 1\f]
        !! The entire Metropolis ratio (ie, the one actually used in a Monte Carlo step) is:
        !! \f[R = |R_\uparrow R_\downarrow|\f]
        !! Note that this subroutine serves to compute the Metropolis ratio used to
        !! determine whether or not \f$h(i, l)\f$ should be flipped, so the field
        !! should not be flipped already.
        !!
        !! Since they are needed in the quick update formula, \f$R_\sigma\f$ and
        !! \f$\delta_\sigma\f$ are stored in the simulation data type `S`:
        !! `S%%Rup` and `S%%deltaup` or `S%%Rdn` and `S%%deltadn` depending on if
        !! \f$\sigma = \uparrow\f$ (`sigma = 1`) or \f$\sigma = \downarrow\f$ (`sigma = -1`).
        !!
        !! \param[inout] S (`Simulation`) Simulation data type. When done, `S%%Rup` and `S%%deltaup` (`sigma = 1`) or
        !! `S%%Rdn` and `S%%deltadn` (`sigma = -1`) will contain the Metropolis ratio for spin
        !! \f$\sigma\f$ and a factor \f$\delta_\sigma\f$ needed for the quick update formula.
        !! \param[in]   i  (`integer`) Site being proposed to be flipped at.
        !! \param[in]   l  (`integer`) Imaginary time slice being proposed to be flipped at.
        !! \param[in]   sigma (`integer`) Spin \f$\sigma\f$ (`sigma = 1` for \f$\sigma = \uparrow\f$, `sigma = -1` for \f$\sigma=\downarrow\f$).
        subroutine greens_R(S, i, l, sigma)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            if (sigma .eq. 1) then      ! up spin (sigma = 1)
                S%deltaup = exp(-2.0_dp * sigma * S%alpha * S%h(i, l)) - 1.0_dp
                S%Rup     = 1.0_dp + (1.0_dp - S%Gup(i, i)) * S%deltaup
            else                        ! dn spin (sigma = -1)
                S%deltadn = exp(-2.0_dp * sigma * S%alpha * S%h(i, l)) - 1.0_dp
                S%Rdn     = 1.0_dp + (1.0_dp - S%Gdn(i, i)) * S%deltadn
            endif
        endsubroutine greens_R


endmodule metropolisratios_mod