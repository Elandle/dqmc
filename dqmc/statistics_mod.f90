!> \brief Contains procedures for doing statistics on data.
module statistics_mod
    use numbertypes
    implicit none

    interface vector_avg
        module procedure :: svector_avg
        module procedure :: dvector_avg
        module procedure :: cvector_avg
        module procedure :: zvector_avg
    endinterface vector_avg

    interface jackknife

        module procedure :: djackknife

        module procedure :: zjackknife
    endinterface jackknife

    contains

        ! -----------------------------------------------------------------------------------------------------
        ! vector_avg ------------------------------------------------------------------------------------------

        !> \brief Returns the average of the entries of the length \f$n\f$ `real(sp)` vector \f$x\f$.
        !!
        !! The average value \f$\text{avg}\f$ of a vector \f$x\f$ of length \f$n\f$ is:
        !! \f[\text{avg} = \frac{1}{n}\sum_{i=1}^nx_i\f]
        !! 
        !! \param[in] x   (`real(sp), dimension(n)`) Vector to take average of entries.
        !! \param[in] n   (`integer`)                Dimension of x
        !! \result    avg (`real(sp)`)               Average value of entries of x
        real(dp) function svector_avg(x, n)
            real(sp), intent(in) :: x(n)
            integer , intent(in) :: n

            svector_avg = sum(x) / n
        endfunction svector_avg

        !> \brief Returns the average of the entries of the length \f$n\f$ `real(dp)` vector \f$x\f$.
        !!
        !! The average value \f$\text{avg}\f$ of a vector \f$x\f$ of length \f$n\f$ is:
        !! \f[\text{avg} = \frac{1}{n}\sum_{i=1}^nx_i\f]
        !! 
        !! \param[in] x   (`real(dp), dimension(n)`) Vector to take average of entries.
        !! \param[in] n   (`integer`)                Dimension of x
        !! \result    avg (`real(dp)`)               Average value of entries of x
        real(dp) function dvector_avg(x, n)
            real(dp), intent(in) :: x(n)
            integer , intent(in) :: n

            dvector_avg = sum(x) / n
        endfunction dvector_avg

        !> \brief Returns the average of the entries of the length \f$n\f$ `complex(sp)` vector \f$x\f$.
        !!
        !! The average value \f$\text{avg}\f$ of a vector \f$x\f$ of length \f$n\f$ is:
        !! \f[\text{avg} = \frac{1}{n}\sum_{i=1}^nx_i\f]
        !! 
        !! \param[in] x   (`complex(sp), dimension(n)`) Vector to take average of entries.
        !! \param[in] n   (`integer`)                   Dimension of x
        !! \result    avg (`complex(sp)`)               Average value of entries of x
        complex(dp) function cvector_avg(x, n)
            complex(sp), intent(in) :: x(n)
            integer    , intent(in) :: n

            cvector_avg = sum(x) / n
        endfunction cvector_avg

        !> \brief Returns the average of the entries of the length \f$n\f$ `complex(dp)` vector \f$x\f$.
        !!
        !! The average value \f$\text{avg}\f$ of a vector \f$x\f$ of length \f$n\f$ is:
        !! \f[\text{avg} = \frac{1}{n}\sum_{i=1}^nx_i\f]
        !! 
        !! \param[in] x   (`complex(dp), dimension(n)`) Vector to take average of entries.
        !! \param[in] n   (`integer`)                   Dimension of x
        !! \result    avg (`complex(dp)`)               Average value of entries of x
        complex(dp) function zvector_avg(x, n) 
            complex(dp), intent(in) :: x(n)
            integer    , intent(in) :: n

            zvector_avg = sum(x) / n
        endfunction zvector_avg

        ! -----------------------------------------------------------------------------------------------------
        ! end vector_avg --------------------------------------------------------------------------------------

        ! -----------------------------------------------------------------------------------------------------
        ! jackknife -------------------------------------------------------------------------------------------
        
        !> Performs the Jackknife method for the average and standard deviation
        !!  of the entries of the length \f$n\f$ `real(dp)` vector \f$n\f$.
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
        subroutine djackknife(x, n, avg, err)
            real(dp), intent(in)  :: x(n)
            integer , intent(in)  :: n
            real(dp), intent(out) :: avg, err

            avg = vector_avg(x, n)
            err = sqrt(sum((x - avg) ** 2) / (n*(n-1)))
        endsubroutine djackknife

        !> Performs the Jackknife method for the average and standard deviation
        !!  of the entries of the length \f$n\f$ `complex(dp)` vector \f$n\f$.
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
        !! \param[in]  x    (`complex(dp), dimension(n)`) Real vector to perform the Jackknife method on the entries of.
        !! \param[in]  n    (`integer`)                   Dimension of `x`.
        !! \param[out] avg  (`complex(dp)`)               Variable to hold the computed Jackknife average.
        !! \param[out] err  (`complex(dp)`)               Variable to hold the computed Jackknife error.
        subroutine zjackknife(x, n, avg, err)
            complex(dp), intent(in)  :: x(n)
            integer    , intent(in)  :: n
            complex(dp), intent(out) :: avg, err

            avg = zvector_avg(x, n)
            err = sqrt(sum((x - avg) ** 2) / (n*(n-1)))
        endsubroutine zjackknife

        ! -----------------------------------------------------------------------------------------------------
        ! end jackknife ---------------------------------------------------------------------------------------

endmodule statistics_mod