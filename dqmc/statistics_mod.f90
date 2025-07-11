!> \brief Contains procedures for doing statistics on data.
module statistics_mod
    use numbertypes
    implicit none

    contains

        !> \brief Returns the average of the entries of the length \f$n\f$ real vector \f$x\f$.
        !!
        !! The average value \f$\text{avg}\f$ of a vector \f$x\f$ of length \f$n\f$ is:
        !! \f[\text{avg} = \frac{1}{n}\sum_{i=1}^nx_i\f]
        !! 
        !! \param[in] x   (`real(dp), dimension(n)`) Vector to take average of entries.
        !! \param[in] n   (`integer`)                Dimension of x
        !! \result    avg (`real(dp)`)               Average value of entries of x
        real(dp) function vector_avg(x, n) result(avg)
            real(dp), intent(in) :: x(n)
            integer , intent(in) :: n

            avg = sum(x) / n
        endfunction vector_avg

        !> Performs the Jackknife method for the average and standard deviation
        !!  of the entries of the length \f$n\f$ vector \f$n\f$.
        !!
        !! The Jackknife method computes the average and standard deviation of a
        !! collection of samples by computing the sample average and sample standard
        !! deviation of sample averages and sample standard deviations computed by
        !! leaving a single entry out.
        !!
        !! To be specific, the Jackknife mean \f$\text{avg}_{\text{jack}}\f$ of a length \f$n\f$
        !! vector \f$x\f$ is computed by first computing:
        !! \f[\text{avg}_i = \frac{1}{n-1}\sum_{j=1,j\neq i}^nx_j\f]
        !! for \f$i=1,\dots,n\f$, and then computing:
        !! \f[\begin{aligned}
        !! \text{avg}_{\text{jack}} &= \frac{1}{n}\sum_{i=1}^n\text{avg}_{i} \\
        !!                          &= \frac{1}{n}\sum_{i=1}^nx_i \\
        !!                          &= \text{avg}
        !! \end{aligned}\f]
        !! It happens that the Jackknife mean equals the sample mean computed as normal.
        !!
        !! The Jackknife variance \f$\sigma_{\text{jack}}^2\f$ (square of the
        !! Jackknife standard deviation/error) is computed by:
        !! \f[\begin{aligned}
        !! \sigma_{\text{jack}}^2 &= \frac{n-1}{n}\sum_{i=1}^n(\text{avg}_i - \text{avg}_{\text{jack}})^2 \\
        !!                        &= \frac{1}{n(n-1)}\sum_{i=1}^n(x_i - \text{avg})^2
        !! \end{aligned}\f]
        !!
        !! \param[in]  x    (`real(dp), dimension(n)`) Real vector to perform the Jackknife method on the entries of.
        !! \param[in]  n    (`integer`)                Dimension of `x`.
        !! \param[out] avg  (`real(dp)`)               Variable to hold the computed Jackknife average.
        !! \param[out] err  (`real(dp)`)               Variable to hold the computed Jackknife error.
        subroutine jackknife(x, n, avg, err)
            real(dp), intent(in)  :: x(n)
            integer , intent(in)  :: n
            real(dp), intent(out) :: avg, err

            avg = vector_avg(x, n)
            err = sqrt(sum((x - avg) ** 2) / (n*(n-1)))
        endsubroutine jackknife


        complex(dp) function zvector_avg(x, n)
            !
            ! Returns the average of the entries of the length n vector x
            !
            complex(dp), intent(in) :: x(n)
            integer    , intent(in) :: n

            zvector_avg = sum(x) / n


        endfunction zvector_avg


        subroutine zjackknife(x, n, avg, err)
            !
            ! Performs the jackknife method on the entries of the length n vector x.
            !
            ! The average obtained by the jackknife method will be stored in avg,
            ! and the error in err.
            !
            complex(dp), intent(in)  :: x(n)
            integer    , intent(in)  :: n
            complex(dp), intent(out) :: avg, err

            avg = zvector_avg(x, n)
            err = sqrt(sum((x - avg) ** 2) / (n*(n-1)))

            
        endsubroutine zjackknife


endmodule statistics_mod