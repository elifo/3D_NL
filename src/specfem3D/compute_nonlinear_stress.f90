subroutine  compute_nonlinear_stress (mul,lambdal,lambdalplus2mul, &
            duxdxl_NL,duydyl_NL,duzdzl_NL, &
            duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
            duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
            dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

    use constants, only: CUSTOM_REAL

    implicit none

    real(kind=CUSTOM_REAL), intent(in) :: lambdal,mul,lambdalplus2mul
    real(kind=CUSTOM_REAL), intent(in) :: duxdxl_NL, duydyl_NL, duzdzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL
    
    real(kind=CUSTOM_REAL), intent(out) ::  dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz

    ! local parameters
    real(kind=CUSTOM_REAL), dimension(6) :: dS, de, dsigm
    ! define NSPR later
    real(kind=CUSTOM_REAL), dimension(50) :: dF
    real(kind=CUSTOM_REAL), dimension(6,6) :: Ed, Esd


!     test this : working!
    dsigma_xx = lambdalplus2mul * duxdxl_NL + lambdal * duydyl_plus_duzdzl_NL
    dsigma_yy = lambdalplus2mul * duydyl_NL + lambdal * duxdxl_plus_duzdzl_NL
    dsigma_zz = lambdalplus2mul * duzdzl_NL + lambdal * duxdxl_plus_duydyl_NL
    dsigma_xy = mul * duxdyl_plus_duydxl_NL
    dsigma_xz = mul * duzdxl_plus_duxdzl_NL
    dsigma_yz = mul * duzdyl_plus_duydzl_NL

    ! if previous test works, start writing Iwan here
    ! put elastic part in a different subroutine.


    ! notes:
    ! herseyden once elastic function yaz ve test et!
    ! get Kmod
    ! get NSPR
    ! get aktif


!     ! Nonlinearity
!     de    = 0._CUSTOM_REAL
!     dsig  = 0._CUSTOM_REAL
!     dsigm = (duxdxl_NL+ duydyl_NL+ duzdzl_NL)* Kmod
!     !dsigm = (deps(1)+deps(2)+deps(6))* Kmod

!     ! First time step
!     if (aktif == -1) then

!         dsig = MAT_IWAN_elastic(mu,Gact,lambda,deps)
!         dS    = dsig
!         dS(1) = dsig(1)- dsigm
!         dS(2) = dsig(2)- dsigm
!         dS(6) = dsig(6)- dsigm   
      
!         aktif = 0
!         S     = S+ dS
!         Sa    = 0._CUSTOM_REAL

!         return
!     endif





end subroutine compute_nonlinear_stress
