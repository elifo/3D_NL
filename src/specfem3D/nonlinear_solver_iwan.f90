module nonlinear_solver_iwan

    use constants

    implicit none

    !!!

    public :: compute_nonlinear_stress

! I allocated every variable in specfem_par;
! maybe I should create a new file just for this module's parameters
!

contains
!=====================================================================

subroutine  compute_nonlinear_stress (mul,lambdal,lambdalplus2mul, &
            duxdxl_NL,duydyl_NL,duzdzl_NL, &
            duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
            duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
            dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

    implicit none


    real(kind=CUSTOM_REAL), intent(in) :: lambdal,mul,lambdalplus2mul
    real(kind=CUSTOM_REAL), intent(in) :: duxdxl_NL, duydyl_NL, duzdzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL    
    real(kind=CUSTOM_REAL), intent(out) ::  dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz


    ! local parameters
    real(kind=CUSTOM_REAL) :: Gact, Gm0, ni !, NSPR
    integer :: n_active_surface

    integer :: surface, errorflag, start, k, j, D
    integer, dimension(6) :: INDX(6)

    real(kind=CUSTOM_REAL) :: lambda_nl, lambdalplus2mu_nl, Kmod
    real(kind=CUSTOM_REAL) :: dsigm, depsm
    real(kind=CUSTOM_REAL), dimension(6) :: dS, de
    real(kind=CUSTOM_REAL), dimension(50) :: dF
    real(kind=CUSTOM_REAL), dimension(6,6) :: Ed, Esd

! Elasticity test
    Gact = mul
    call MAT_IWAN_elastic (Gact,lambdal,lambdalplus2mul, &
            duxdxl_NL,duydyl_NL,duzdzl_NL, &
            duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
            duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
            dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)


    ! notes:
    ! get ni  in
    ! get NSPR  in
    ! get n_active_surface inout
    ! Sa(:,:) inout
    ! F(:) inout
    ! R(:) in, or compute it each time we call.
    ! CNinv(:) in, or compute it each time we call.
    ! 
    ! 2d kodda 1,2,6: x,z,y idi. dikkat!!!
    ! burda gamma ile calisicam, epsilon ile degil, aman !!!



!     Gm0   = mul !! to get as input
!     ni    = 0.2987478_CUSTOM_REAL !! to change, to get as input
!     n_active_surface = -1    


! !     ! Nonlinearity
!     de    = 0._CUSTOM_REAL
!     dsig  = 0._CUSTOM_REAL
    
!     ! update moduli
!     Kmod      = 2._CUSTOM_REAL* Gm0* (1._CUSTOM_REAL+ ni) ! E (Young modulus) 
!     lambda_nl = Kmod* ni/ (1._CUSTOM_REAL+ ni)/ (1._CUSTOM_REAL- 2._CUSTOM_REAL* ni)
!     lambdalplus2mu_nl = lambda_nl+ 2._CUSTOM_REAL* Gm0
!     Kmod   = Kmod/ (3._CUSTOM_REAL* (1._CUSTOM_REAL- 2._CUSTOM_REAL* ni)) ! K (Bulk modulus)

!     ! mean stress
!     dsigm = (duxdxl_NL+ duydyl_NL+ duzdzl_NL)* Kmod

!     ! First time step
!     if (n_active_surface == -1) then

!         call MAT_IWAN_elastic (Gact,lambdal,lambdalplus2mul, &
!                 duxdxl_NL,duydyl_NL,duzdzl_NL, &
!                 duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
!                 duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
!                 dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

!         dS(1) = dsigma_xx- dsigm
!         dS(2) = dsigma_yy- dsigm
!         dS(6) = dsigma_zz- dsigm
!         dS(3) = dsigma_xy
!         dS(4) = dsigma_xz
!         dS(5) = dsigma_yz

!         n_active_surface = 0
!         S = S+ dS
!         Sa = 0._CUSTOM_REAL

!         return

!     endif

!     depsm = (duxdxl_NL+ duydyl_NL+ duzdzl_NL)/ 3._CUSTOM_REAL
!     ! Incremental deviatoric strain
!     de(1) = duxdxl_NL- depsm
!     de(2) = duydyl_NL- depsm
!     de(6) = duzdzl_NL- depsm
!     !!! check here !!!
!     de(3) = duxdyl_plus_duydxl_NL  ! deps_xy+ deps_yx
!     de(4) = duzdxl_plus_duxdzl_NL  ! deps_xz+ deps_zx
!     de(5) = duzdyl_plus_duydzl_NL  ! deps_yz+ deps_zy


!     if (n_active_surface == 0) then
!         F (1) = MAT_IWAN_surface(S, Sa(1,:))
!         dF(1) = MAT_IWAN_dsurface(S,Sa(1,:),de)
!     endif

    
!     surface = 0
!     if (n_active_surface > 0) then

!         do j = 1,n_active_surface

!             ! New centers
!             do k = 1,6
!                 Sa(j,k) = S(k)- R(j)/ sqrt(F(j))* (S(k)-Sa(j,k))
!             enddo

!             F (j) = MAT_IWAN_surface(S, Sa(j,:))
!             dF(j) = MAT_IWAN_dsurface(S, Sa(j,:), de)            

!             if ( (dF(j) >= 0._CUSTOM_REAL)  .AND. (F(j) >= R(j)** 2) ) &
!             surface = surface+ 1  

!         enddo
!     endif



!     if ( (dF(1) >= 0._CUSTOM_REAL)  .AND. (F(1) < R(1)** 2) ) then

!         call MAT_IWAN_elastic (Gact,lambdal,lambdalplus2mul, &
!                 duxdxl_NL,duydyl_NL,duzdzl_NL, &
!                 duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
!                 duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
!                 dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

!         dS(1) = dsigma_xx- dsigm
!         dS(2) = dsigma_yy- dsigm
!         dS(6) = dsigma_zz- dsigm
!         dS(3) = dsigma_xy
!         dS(4) = dsigma_xz
!         dS(5) = dsigma_yz

!         S = S+ dS
!         return
!     endif


!     ! Ed computation
!     Ed    = 0._CUSTOM_REAL
!     start = 1

!     Ed(1,1) = 0.5_CUSTOM_REAL/ Gact
!     Ed(2,2) = 0.5_CUSTOM_REAL/ Gact
!     Ed(3,3) = 0.5_CUSTOM_REAL/ Gact
!     Ed(4,4) = 0.5_CUSTOM_REAL/ Gact
!     Ed(5,5) = 0.5_CUSTOM_REAL/ Gact
!     Ed(6,6) = 0.5_CUSTOM_REAL/ Gact

!     if (surface > 0) &
!     call MAT_IWAN_Ematris(start, Nspr, Ed, CNinv, S, Sa, F, surface)  

!     do j = surface+1,Nspr-1
!         F(j)  = MAT_IWAN_surface(S, Sa(j,:))
!         dF(j) = MAT_IWAN_dsurface(S, Sa(j,:), de)    

!         if ( (dF(j) >= 0._CUSTOM_REAL)  .AND. (F(j) >= R(j)** 2) ) then
!             surface = surface+ 1  
!             start   = surface
!             call MAT_IWAN_Ematris(start, Nspr, Ed, CNinv, S, Sa, F, surface)
!         else
!             EXIT
!         endif
!     enddo 


!     ! Alternative for inversion
!     ! solve Ed·x = de (de is used as input/ouput)
!     call LUDCMP(Ed, 6, INDX, D, errorflag)
! !     if (errorflag .ne. 0) stop 'not invertible matrix' ! here error call !!!
!     call LUBKSB(Ed, 6, INDX, de) 
!     dS = de


!     S = S+ dS    
!     dsigma_xx = dS(1)+ dsigm
!     dsigma_yy = dS(2)+ dsigm
!     dsigma_zz = dS(6)+ dsigm
!     dsigma_xy = dS(3)
!     dsigma_xz = dS(4)
!     dsigma_yz = dS(5)


!     n_active_surface = max(1, surface)



end subroutine compute_nonlinear_stress
!=====================================================================

! Incremental elastic computation of 
! stress matrix for given strain matrix

subroutine MAT_IWAN_elastic (Gact,lambdal,lambdalplus2mul, &
            duxdxl_NL,duydyl_NL,duzdzl_NL, &
            duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
            duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
            dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

    
    implicit none


    real(kind=CUSTOM_REAL), intent(in) :: lambdal,Gact,lambdalplus2mul
    real(kind=CUSTOM_REAL), intent(in) :: duxdxl_NL, duydyl_NL, duzdzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL
    real(kind=CUSTOM_REAL), intent(out) ::  dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz


    dsigma_xx = lambdalplus2mul * duxdxl_NL + lambdal * duydyl_plus_duzdzl_NL
    dsigma_yy = lambdalplus2mul * duydyl_NL + lambdal * duxdxl_plus_duzdzl_NL
    dsigma_zz = lambdalplus2mul * duzdzl_NL + lambdal * duxdxl_plus_duydyl_NL
    dsigma_xy = Gact * duxdyl_plus_duydxl_NL
    dsigma_xz = Gact * duzdxl_plus_duxdzl_NL
    dsigma_yz = Gact * duzdyl_plus_duydzl_NL


end subroutine MAT_IWAN_elastic
!=====================================================================
! ! Iwan surface(s) computation 

! function MAT_IWAN_surface(S,Sa)  result(F)

!     real(kind=CUSTOM_REAL), intent(in) :: S(:)
!     real(kind=CUSTOM_REAL), intent(in) :: Sa(:)
!     real(kind=CUSTOM_REAL) :: F

!     ! burdaki 2 katsayilarini KONTROL ET !!!

!     F = 0._CUSTOM_REAL
!     F = 0.5_CUSTOM_REAL* ( (S(1)-Sa(1))** 2+ ( S(2)- Sa(2) )** 2 &
!                   + 2._CUSTOM_REAL* ( S(3)- Sa(3) )** 2  &
!                   + 2._CUSTOM_REAL* ( S(4)- Sa(4) )** 2  &
!                   + 2._CUSTOM_REAL* ( S(5)- Sa(5) )** 2  &
!                   + ( S(6)- Sa(6) )** 2 )

! end function MAT_IWAN_surface  
! !=====================================================================
! ! Iwan surface(s) movement computation 

! function MAT_IWAN_dsurface(S,Sa,de)  result(dF)

!     real(kind=CUSTOM_REAL), intent(in) :: S(:)
!     real(kind=CUSTOM_REAL), intent(in) :: Sa(:)
!     real(kind=CUSTOM_REAL), intent(in) :: de(:)
!     real(kind=CUSTOM_REAL) :: dF

!     ! burdaki 2 katsayilarini KONTROL ET !!!


!     dF = 0._CUSTOM_REAL
!     dF = 0.5_CUSTOM_REAL* ( (S(1)- Sa(1) )* de(1)+ ( S(2)- Sa(2) )* de(2)+ &
!                   1._CUSTOM_REAL* ( S(3)- Sa(3) )* de(3) + &
!                   1._CUSTOM_REAL* ( S(4)- Sa(4) )* de(4) + &
!                   1._CUSTOM_REAL* ( S(5)- Sa(5) )* de(5) &
!                   + ( S(6)- Sa(6) )* de(6) )

! end function MAT_IWAN_dsurface  
! !=======================================================================
! ! Iwan plasticity matrix

! subroutine MAT_IWAN_Ematris (start, Nspr, Ed, CNinv, S1, Sa1, F1, aktif)

!     integer, intent(IN) :: start
!     integer, intent(IN) :: Nspr
!     real(kind=CUSTOM_REAL), INTENT(INOUT)  :: Ed(:,:)
!     real(kind=CUSTOM_REAL), INTENT(IN)     :: CNinv(:)
!     real(kind=CUSTOM_REAL), INTENT(IN)     :: S1(:)
!     real(kind=CUSTOM_REAL), INTENT(IN)     :: Sa1(:,:)
!     real(kind=CUSTOM_REAL), INTENT(INOUT)  :: F1 (:)
!     integer, intent(IN) :: aktif

!     integer :: j,m,k
!     real(kind=CUSTOM_REAL) :: ss(6)

!     ss(1) = 1._CUSTOM_REAL
!     ss(2) = 1._CUSTOM_REAL
!     ss(3) = 2._CUSTOM_REAL
!     ss(4) = 2._CUSTOM_REAL
!     ss(5) = 2._CUSTOM_REAL
!     ss(6) = 1._CUSTOM_REAL

!     j = start
!     do while ( j < aktif+1 )

!         do m = 1,6

!           do k = 1,6    
!             Ed(m,k) = Ed(m,k)+ CNinv(j)* ss(k)* (S1(m)- Sa1(j,m)) &
!                       * (S1(k)- Sa1(j,k))/ (2.0* F1(j))
!           enddo

!         enddo
!         j = j+1

!     enddo

! end subroutine MAT_IWAN_Ematris




!=====================================================================

end module nonlinear_solver_iwan
