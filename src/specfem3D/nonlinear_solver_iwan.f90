module nonlinear_solver_iwan

    use constants

    implicit none

    public :: compute_nonlinear_stress

contains

!=====================================================================

subroutine  compute_nonlinear_stress(mul,lambdal,lambdalplus2mul, &
                      duxdxl_NL,duydyl_NL,duzdzl_NL, &
                      duxdyl_NL,duydxl_NL,duzdxl_NL,duxdzl_NL,duzdyl_NL,duydzl_NL, &
                      dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

    implicit none

    real(kind=CUSTOM_REAL), intent(in) :: lambdal,mul,lambdalplus2mul
    real(kind=CUSTOM_REAL), intent(in) :: duxdxl_NL, duydyl_NL, duzdzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_NL,duydxl_NL,duzdxl_NL,duxdzl_NL,duzdyl_NL,duydzl_NL
    double precision, intent(out) ::  dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz


    ! local parameters
    real(kind=CUSTOM_REAL) :: Gact


    ! Elasticity test
    Gact = mul
    call MAT_IWAN_elastic (Gact,lambdal,lambdalplus2mul, &
                      duxdxl_NL,duydyl_NL,duzdzl_NL, &
                      duxdyl_NL,duydxl_NL,duzdxl_NL,duxdzl_NL,duzdyl_NL,duydzl_NL, &
                      dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)




end subroutine compute_nonlinear_stress
!=====================================================================

! Incremental elastic computation of 
! stress matrix for given strain matrix

subroutine MAT_IWAN_elastic (Gact,lambdal,lambdalplus2mul, &
                      duxdxl_NL,duydyl_NL,duzdzl_NL, &
                      duxdyl_NL,duydxl_NL,duzdxl_NL,duxdzl_NL,duzdyl_NL,duydzl_NL, &
                      dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)


    implicit none


    real(kind=CUSTOM_REAL), intent(in) :: lambdal,Gact,lambdalplus2mul
    real(kind=CUSTOM_REAL), intent(in) :: duxdxl_NL, duydyl_NL, duzdzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_NL,duydxl_NL,duzdxl_NL,duxdzl_NL,duzdyl_NL,duydzl_NL
    double precision, intent(out) ::  dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz

    real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl_NL, duxdxl_plus_duzdzl_NL, duydyl_plus_duzdzl_NL, &
                              duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL



    duxdxl_plus_duydyl_NL = duxdxl_NL + duydyl_NL
    duxdxl_plus_duzdzl_NL = duxdxl_NL + duzdzl_NL
    duydyl_plus_duzdzl_NL = duydyl_NL + duzdzl_NL
    duxdyl_plus_duydxl_NL = duxdyl_NL + duydxl_NL
    duzdxl_plus_duxdzl_NL = duzdxl_NL + duxdzl_NL
    duzdyl_plus_duydzl_NL = duzdyl_NL + duydzl_NL


    dsigma_xx = lambdalplus2mul * duxdxl_NL + lambdal * duydyl_plus_duzdzl_NL
    dsigma_yy = lambdalplus2mul * duydyl_NL + lambdal * duxdxl_plus_duzdzl_NL
    dsigma_zz = lambdalplus2mul * duzdzl_NL + lambdal * duxdxl_plus_duydyl_NL
    dsigma_xy = Gact * duxdyl_plus_duydxl_NL
    dsigma_xz = Gact * duzdxl_plus_duxdzl_NL
    dsigma_yz = Gact * duzdyl_plus_duydzl_NL


end subroutine MAT_IWAN_elastic

end module nonlinear_solver_iwan

