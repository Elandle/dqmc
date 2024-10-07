module measurements_mod
    use numbertypes
    use simulationsetup_mod
    use statistics_mod
    implicit none


    contains

        subroutine measure(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            call measure_sgn(S, i)
            call measure_den(S, i)


        endsubroutine measure

    
        subroutine measure_den(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            integer :: j

            do j = 1, S%N
                S%updenbin(i) = 1.0_dp - S%Gup(j, j)
                S%dndenbin(i) = 1.0_dp - S%Gdn(j, j)
            enddo
            S%updenbin(i) = S%sgn * (S%updenbin(i) / S%N)
            S%dndenbin(i) = S%sgn * (S%dndenbin(i) / S%N)


        endsubroutine measure_den


        subroutine measure_sgn(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            S%sgnbin(i) = S%sgn

        endsubroutine measure_sgn


        subroutine avgbin(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            S%sgnbinavgs  (i) = real(sum(S%sgnbin), dp) / S%binsize
            S%updenbinavgs(i) = vector_avg(S%updenbin, S%binsize) / S%sgnbinavgs(i)
            S%dndenbinavgs(i) = vector_avg(S%dndenbin, S%binsize) / S%sgnbinavgs(i)


        endsubroutine avgbin


        subroutine dostatistics(S)
            type(Simulation), intent(inout) :: S

            call jackknife(S%sgnbinavgs, S%nbin, S%sgnavg, S%sgnerr)
            call jackknife(S%updenbinavgs, S%nbin, S%updenavg, S%updenerr)
            call jackknife(S%dndenbinavgs, S%nbin, S%dndenavg, S%dndenerr)

        endsubroutine dostatistics


endmodule measurements_mod