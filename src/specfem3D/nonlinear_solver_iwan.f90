module nonlinear_solver_iwan

    use constants

    implicit none

    !!!

    public :: MAT_IWAN_init, compute_nonlinear_stress

! I allocated every variable in specfem_par;
! maybe I should create a new file just for this module's parameters
!

contains
!=====================================================================

subroutine MAT_IWAN_init()

    use specfem_par

    use specfem_par_elastic, only: ispec_is_elastic, &
                    R_NL, CNinv_NL, S_NL, F_NL, n_active_surface, &
                    Sa_NL_xx, Sa_NL_yy, Sa_NL_zz, &
                    Sa_NL_xy, Sa_NL_yz, Sa_NL_xz, &
                    NSPR 

    implicit none

    integer :: i,j,k,ispec
    logical :: is_OverBurden, testing_sil
    real(kind=CUSTOM_REAL) :: Gm0, gamma_ref

    ! Force the code 
    is_OverBurden = .False.



    ! shifts shear moduli
    do ispec = 1,NSPEC_AB

        ! skips non elastic elements
        if (ispec_is_elastic(ispec) .eqv. .false.) cycle

        ! determines XXX for each GLL point
        do k = 1,NGLLZ
            do j = 1,NGLLY
                do i = 1,NGLLX        

                    R_NL(:,i,j,k,ispec) = 0e0_CUSTOM_REAL
                    CNinv_NL(:,i,j,k,ispec) = 0e0_CUSTOM_REAL
                    !S_NL(:,i,j,k,ispec) = 1e-3_CUSTOM_REAL
                    S_NL(:,i,j,k,ispec) = 0e0_CUSTOM_REAL
                    F_NL(:,i,j,k,ispec) = 0e0_CUSTOM_REAL

                    n_active_surface(i,j,k,ispec) = -1

                    Sa_NL_xx(:,i,j,k,ispec)= 0e0_CUSTOM_REAL
                    Sa_NL_yy(:,i,j,k,ispec)= 0e0_CUSTOM_REAL
                    Sa_NL_zz(:,i,j,k,ispec)= 0e0_CUSTOM_REAL
                    Sa_NL_xy(:,i,j,k,ispec)= 0e0_CUSTOM_REAL
                    Sa_NL_xz(:,i,j,k,ispec)= 0e0_CUSTOM_REAL
                    Sa_NL_yz(:,i,j,k,ispec)= 0e0_CUSTOM_REAL

                    ! construct the backbone curve
                    if (.not. is_OverBurden) then
                        gamma_ref = 0.000365e0_CUSTOM_REAL
                        Gm0 = mustore(i,j,k,ispec)


                        call MAT_IWAN_backbone_elem (gamma_ref, NSPR, Gm0, &
                                    R_NL(:,i,j,k,ispec),CNinv_NL(:,i,j,k,ispec) ) 

!                        if (ispec == 50 .and. i==1 .and. j==1 .and. k ==1) then
!                            print*, 'Found shear modulus: ', Gm0
!                            print*, R_NL(:,i,j,k,ispec), CNinv_NL(:,i,j,k,ispec)
!                        endif


! not working
                         call MAT_IWAN_backbone_elem_testing (gamma_ref, NSPR, Gm0, &
                                     R_NL(:,i,j,k,ispec),CNinv_NL(:,i,j,k,ispec) ) 


!                        if (ispec == 50 .and. i==1 .and. j==1 .and. k ==1) then
!                            print*, 'Found shear modulus: ', Gm0
!                            print*, R_NL(:,i,j,k,ispec), CNinv_NL(:,i,j,k,ispec)
!                        endif


                    endif


                enddo
            enddo
        enddo

    enddo

end subroutine MAT_IWAN_init
!=====================================================================

subroutine MAT_IWAN_backbone_elem_testing(gref,Nspr,mu,R,CNinv)

    implicit none

    integer, intent(in) :: Nspr
    real(kind=CUSTOM_REAL), intent(in) :: gref, mu
    real(kind=CUSTOM_REAL), intent(inout) :: R(Nspr), CNinv(Nspr-1)

    real(kind=CUSTOM_REAL) :: x0, xu, dx, gamma(Nspr), G(Nspr), summy
    integer :: i

    x0 = -6e0_CUSTOM_REAL
    xu = log10(0.1e0_CUSTOM_REAL)
    dx = (xu- x0)/ (Nspr- 1)

    gamma  = 0e0_CUSTOM_REAL
    G      = 0e0_CUSTOM_REAL
    R      = 0e0_CUSTOM_REAL
    CNinv  = 0e0_CUSTOM_REAL

    ! G/Gmax=1 to test elasticity
    do i = 1,Nspr
        gamma (i) = 10e0_CUSTOM_REAL** ( x0+ dx* (i-1) )
        G(i)    = 1e0_CUSTOM_REAL
        R(i)    = G(i)* mu* gamma(i)
    enddo
    
    summy = 0e0_CUSTOM_REAL
    do i = 1,Nspr-1
        CNinv(i) = ( gamma(i+1)/2e0_CUSTOM_REAL - gamma(i)/2e0_CUSTOM_REAL )/ &
                        ( R(i+1)- R(i) )- 0.5e0_CUSTOM_REAL/mu- summy
        summy = summy+ CNinv(i)
    enddo   

end subroutine MAT_IWAN_backbone_elem_testing
!=====================================================================

subroutine MAT_IWAN_backbone_elem(gref,Nspr,mu,R,CNinv)

    implicit none

    integer, intent(in) :: Nspr
    real(kind=CUSTOM_REAL), intent(in) :: gref, mu
    real(kind=CUSTOM_REAL), intent(inout) :: R(Nspr), CNinv(Nspr-1)

    real(kind=CUSTOM_REAL) :: x0, xu, dx, gamma(Nspr), G(Nspr), summy
    integer :: i

    x0 = -6e0_CUSTOM_REAL
    xu = log10(0.1e0_CUSTOM_REAL)
    dx = (xu- x0)/ (Nspr- 1)

    gamma  = 0e0_CUSTOM_REAL
    G      = 0e0_CUSTOM_REAL
    R      = 0e0_CUSTOM_REAL
    CNinv  = 0e0_CUSTOM_REAL

    do i = 1,Nspr
        gamma (i) = 10e0_CUSTOM_REAL** ( x0+ dx* (i-1) )
        G(i)    = 1e0_CUSTOM_REAL/ ( 1e0_CUSTOM_REAL+ abs( gamma(i)/ gref) )
        R(i)    = G(i)* mu* gamma(i)
    enddo
    
    summy = 0e0_CUSTOM_REAL
    do i = 1,Nspr-1
        CNinv(i) = ( gamma(i+1)/2e0_CUSTOM_REAL - gamma(i)/2e0_CUSTOM_REAL )/ &
                        ( R(i+1)- R(i) )- 0.5e0_CUSTOM_REAL/mu- summy
        summy = summy+ CNinv(i)
    enddo   

end subroutine MAT_IWAN_backbone_elem
!=====================================================================

subroutine  compute_nonlinear_stress (ii,jj,kk,ispec, Nspr, n_active_surface, &
            S, F, R, CNinv, Sa_xx, Sa_yy, Sa_zz, Sa_xy, Sa_xz, Sa_yz, &
            mul,lambdal,lambdalplus2mul, &
            duxdxl_NL,duydyl_NL,duzdzl_NL, &
            duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
            duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
             duxdyl_NL, duxdzl_NL, duydzl_NL, &
            dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

    implicit none

    integer, intent(in) :: Nspr, ispec, ii,jj,kk
    integer, intent(inout) :: n_active_surface
    real(kind=CUSTOM_REAL), intent(inout) :: S(6)
    real(kind=CUSTOM_REAL), intent(inout) :: F(Nspr)
    real(kind=CUSTOM_REAL), intent(inout) :: R(Nspr)
    real(kind=CUSTOM_REAL), intent(inout) :: CNinv(Nspr-1)
    real(kind=CUSTOM_REAL), intent(inout) :: Sa_xx(Nspr), Sa_yy(Nspr), Sa_zz(Nspr)
    real(kind=CUSTOM_REAL), intent(inout) :: Sa_xy(Nspr), Sa_xz(Nspr), Sa_yz(Nspr)

    real(kind=CUSTOM_REAL), intent(in) :: lambdal,mul,lambdalplus2mul
    real(kind=CUSTOM_REAL), intent(in) :: duxdxl_NL, duydyl_NL, duzdzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_NL, duxdzl_NL, duydzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL    
    real(kind=CUSTOM_REAL), intent(out) :: dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz


    ! local parameters
    integer :: surface, errorflag, start, k, j, D
    integer :: INDX(6)
    real(kind=CUSTOM_REAL) :: Gm0, ni    
    real(kind=CUSTOM_REAL) :: lambda_nl, lambdalplus2mu_nl, Kmod
    real(kind=CUSTOM_REAL) :: dsigm, depsm
    real(kind=CUSTOM_REAL) :: dS(6), de(6)
    real(kind=CUSTOM_REAL) :: Ed(6,6), Esd(6,6)
    real(kind=CUSTOM_REAL) :: dF(Nspr)
    real(kind=CUSTOM_REAL) :: dummy, dummy2, infinity


    

! ELASTICITY TEST I
!    Gm0   = mul 
!    call MAT_IWAN_elastic (Gm0,lambdal,lambdalplus2mul, &
!            duxdxl_NL,duydyl_NL,duzdzl_NL, &
!            duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
!            duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
!            dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)
!    return




! ELASTICITY TEST II
!    ! this part can be moved to prep before time loop !!!
!    ! update moduli
!    Gm0   = mul 
!    ni    = 0.3875e0_CUSTOM_REAL 
!    Kmod = 2e0_CUSTOM_REAL* Gm0* (1e0_CUSTOM_REAL+ ni) ! E (Young modulus) 
!    lambda_nl = Kmod* ni/ (1e0_CUSTOM_REAL+ ni)/ (1e0_CUSTOM_REAL- 2e0_CUSTOM_REAL* ni)
!    lambdalplus2mu_nl = lambda_nl+ 2e0_CUSTOM_REAL* Gm0
!    Kmod = Kmod/ (3e0_CUSTOM_REAL* (1e0_CUSTOM_REAL- 2e0_CUSTOM_REAL* ni)) ! K (Bulk modulus)

!    call MAT_IWAN_elastic (Gm0,lambda_nl,lambdalplus2mu_nl, &
!            duxdxl_NL,duydyl_NL,duzdzl_NL, &
!            duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
!            duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
!             dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)
!    return




! ELASTICITY TEST III

    Gm0   = mul 
    ni    = 0.3875e0_CUSTOM_REAL 
    Kmod = 2e0_CUSTOM_REAL* Gm0* (1e0_CUSTOM_REAL+ ni) ! E (Young modulus) 
    lambda_nl = Kmod* ni/ (1e0_CUSTOM_REAL+ ni)/ (1e0_CUSTOM_REAL- 2e0_CUSTOM_REAL* ni)
    lambdalplus2mu_nl = lambda_nl+ 2e0_CUSTOM_REAL* Gm0
    Kmod = Kmod/ (3e0_CUSTOM_REAL* (1e0_CUSTOM_REAL- 2e0_CUSTOM_REAL* ni)) ! K (Bulk modulus)


    ! mean stress
    dsigm = (duxdxl_NL+ duydyl_NL+ duzdzl_NL)* Kmod

    ! First time step
    if (n_active_surface == -1) then

        call MAT_IWAN_elastic (Gm0,lambda_nl,lambdalplus2mu_nl, &
                duxdxl_NL,duydyl_NL,duzdzl_NL, &
                duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
                duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
                dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

        S(1) = S(1)+ dsigma_xx- dsigm
        S(2) = S(2)+ dsigma_yy- dsigm
        S(6) = S(6)+ dsigma_zz- dsigm
        S(3) = S(3)+ dsigma_xy
        S(4) = S(4)+ dsigma_xz
        S(5) = S(5)+ dsigma_yz

        n_active_surface = 0

        Sa_xx = 0e0_CUSTOM_REAL
        Sa_yy = 0e0_CUSTOM_REAL
        Sa_zz = 0e0_CUSTOM_REAL
        Sa_xy = 0e0_CUSTOM_REAL
        Sa_xz = 0e0_CUSTOM_REAL
        Sa_yz = 0e0_CUSTOM_REAL

        return
    endif




    de    = 0e0_CUSTOM_REAL
    depsm = (duxdxl_NL+ duydyl_NL+ duzdzl_NL)/ 3e0_CUSTOM_REAL
    ! Incremental deviatoric strain
    de(1) = duxdxl_NL- depsm
    de(2) = duydyl_NL- depsm
    de(6) = duzdzl_NL- depsm
    de(3) = duxdyl_NL
    de(4) = duxdzl_NL
    de(5) = duydzl_NL
    !!! check here !!!
!    de(3) = duxdyl_plus_duydxl_NL  ! deps_xy+ deps_yx
!    de(4) = duzdxl_plus_duxdzl_NL  ! deps_xz+ deps_zx
!    de(5) = duzdyl_plus_duydzl_NL  ! deps_yz+ deps_zy


    if (n_active_surface == 0) then
        F(1)  = MAT_IWAN_surface( S, Sa_xx(1), Sa_yy(1), Sa_zz(1), Sa_xy(1), Sa_xz(1), Sa_yz(1) )
        dF(1) = MAT_IWAN_dsurface(S, Sa_xx(1), Sa_yy(1), Sa_zz(1), Sa_xy(1), Sa_xz(1), Sa_yz(1), de)
    endif


   
    surface = 0
    if (n_active_surface > 0) then

        do j = 1,n_active_surface



            ! New centers (traceback: here is the bug)  
            Sa_xx(j) = S(1)- R(j)/ sqrt(F(j))* (S(1)- Sa_xx(j))
            Sa_yy(j) = S(2)- R(j)/ sqrt(F(j))* (S(2)- Sa_yy(j))
            Sa_zz(j) = S(6)- R(j)/ sqrt(F(j))* (S(6)- Sa_zz(j))
            Sa_xy(j) = S(3)- R(j)/ sqrt(F(j))* (S(3)- Sa_xy(j))
            Sa_xz(j) = S(4)- R(j)/ sqrt(F(j))* (S(4)- Sa_xz(j))
            Sa_yz(j) = S(5)- R(j)/ sqrt(F(j))* (S(5)- Sa_yz(j))

            F(j)  = MAT_IWAN_surface( S, Sa_xx(j), Sa_yy(j), Sa_zz(j), Sa_xy(j), Sa_xz(j), Sa_yz(j) )
            dF(j) = MAT_IWAN_dsurface ( S, Sa_xx(j), Sa_yy(j), Sa_zz(j), Sa_xy(j), Sa_xz(j), Sa_yz(j), de )

            if ( (dF(j) >= 0e0_CUSTOM_REAL)  .AND. (F(j) >= R(j)** 2) ) &
            surface = surface+ 1  

        enddo

    endif


    if ( (dF(1) >= 0e0_CUSTOM_REAL)  .AND. (F(1) < R(1)** 2) ) then

        call MAT_IWAN_elastic (Gm0,lambda_nl,lambdalplus2mu_nl, &
                duxdxl_NL,duydyl_NL,duzdzl_NL, &
                duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
                duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
                dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

        S(1) = S(1)+ dsigma_xx- dsigm
        S(2) = S(2)+ dsigma_yy- dsigm
        S(6) = S(6)+ dsigma_zz- dsigm
        S(3) = S(3)+ dsigma_xy
        S(4) = S(4)+ dsigma_xz
        S(5) = S(5)+ dsigma_yz

        return
    endif

    ! Ed computation: the matrix that relates deviatoric stress and strain
    Ed    = 0e0_CUSTOM_REAL
    start = 1


    ! elastic part 
    dummy = 0.5_CUSTOM_REAL/ Gm0
    Ed(1,1) = dummy
    Ed(2,2) = dummy
    Ed(3,3) = dummy
    Ed(4,4) = dummy
    Ed(5,5) = dummy
    Ed(6,6) = dummy

    if (surface > 0) &
    call MAT_IWAN_Ematris(start, Nspr, Ed, CNinv, S, Sa_xx, Sa_yy, Sa_zz, Sa_xy, Sa_xz, Sa_yz, F, surface)  



    do j = surface+1,Nspr-1

        F(j)  = MAT_IWAN_surface( S, Sa_xx(j), Sa_yy(j), Sa_zz(j), Sa_xy(j), Sa_xz(j), Sa_yz(j) )
        dF(j) = MAT_IWAN_dsurface ( S, Sa_xx(j), Sa_yy(j), Sa_zz(j), Sa_xy(j), Sa_xz(j), Sa_yz(j), de )

        if ( (dF(j) >= 0e0_CUSTOM_REAL)  .AND. (F(j) >= R(j)** 2) ) then
            surface = surface+ 1  
            start   = surface
            call MAT_IWAN_Ematris(start, Nspr, Ed, CNinv, S, Sa_xx, Sa_yy, Sa_zz, Sa_xy, Sa_xz, Sa_yz, F, surface)  
        else
            EXIT
        endif
    enddo 


    !print*, 'STEP IV'
    ! Alternative for inversion
    ! solve Ed·x = de (de is used as input/ouput)
    call LUDCMP(Ed, 6, INDX, D, errorflag)
!    if (errorflag .ne. 0) stop 'not invertible matrix' ! here error call !!!
    if (errorflag .ne. 0) print*, '*** not invertible matrix' 
    call LUBKSB(Ed, 6, INDX, de) 
    


    dS = de
    S = S+ dS    

    dsigma_xx = dS(1)+ dsigm
    dsigma_yy = dS(2)+ dsigm
    dsigma_zz = dS(6)+ dsigm
    dsigma_xy = dS(3)
    dsigma_xz = dS(4)
    dsigma_yz = dS(5)


    n_active_surface = max(1, surface)


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
! Iwan surface(s) computation 

function MAT_IWAN_surface(S, Sa_xx, Sa_yy, Sa_zz, Sa_xy, Sa_xz, Sa_yz)  result(F)

    implicit none

    real(kind=CUSTOM_REAL), intent(in) :: S(6)
    real(kind=CUSTOM_REAL), intent(inout) :: Sa_xx, Sa_yy, Sa_zz
    real(kind=CUSTOM_REAL), intent(inout) :: Sa_xy, Sa_xz, Sa_yz
    real(kind=CUSTOM_REAL) :: F


    F = 0.e0_CUSTOM_REAL

    F = 0.5e0_CUSTOM_REAL* ( (S(1)-Sa_xx )** 2+ ( S(2)- Sa_yy )** 2 &
                  + 2e0_CUSTOM_REAL* ( S(3)- Sa_xy )** 2  &
                  + 2e0_CUSTOM_REAL* ( S(4)- Sa_xz )** 2  &
                  + 2e0_CUSTOM_REAL* ( S(5)- Sa_yz )** 2  &
                  + ( S(6)- Sa_zz )** 2 )

    !!!
    F = max()


end function MAT_IWAN_surface  
!=====================================================================
! Iwan surface(s) movement computation 

function MAT_IWAN_dsurface(S,Sa_xx,Sa_yy,Sa_zz,Sa_xy,Sa_xz,Sa_yz,de) result(dF)

    implicit none

    real(kind=CUSTOM_REAL), intent(in) :: S(6)
    real(kind=CUSTOM_REAL), intent(inout) :: Sa_xx, Sa_yy, Sa_zz
    real(kind=CUSTOM_REAL), intent(inout) :: Sa_xy, Sa_xz, Sa_yz
    real(kind=CUSTOM_REAL), intent(in) :: de(6)
    real(kind=CUSTOM_REAL) :: dF

    dF = 0e0_CUSTOM_REAL

    dF = 0.5e0_CUSTOM_REAL* ( ( S(1)- Sa_xx )* de(1)+ ( S(2)- Sa_yy )* de(2)+ &
                  2e0_CUSTOM_REAL* ( S(3)- Sa_xy )* de(3) + &
                  2e0_CUSTOM_REAL* ( S(4)- Sa_xz )* de(4) + &
                  2e0_CUSTOM_REAL* ( S(5)- Sa_yz )* de(5) + &
                  ( S(6)- Sa_zz )* de(6) )

end function MAT_IWAN_dsurface  
!=======================================================================
! Iwan plasticity matrix

subroutine MAT_IWAN_Ematris (start, Nspr, Ed, CNinv, S1, Sa_xx, Sa_yy, Sa_zz, Sa_xy, Sa_xz, Sa_yz, F1, aktif)

    implicit none

    integer, intent(IN) :: start
    integer, intent(IN) :: Nspr
    real(kind=CUSTOM_REAL), INTENT(INOUT)  :: Ed(6,6)
    real(kind=CUSTOM_REAL), INTENT(IN)     :: CNinv(Nspr-1)
    real(kind=CUSTOM_REAL), INTENT(IN)     :: S1(6)
    real(kind=CUSTOM_REAL), INTENT(IN)     :: Sa_xx(Nspr), Sa_yy(Nspr), Sa_zz(Nspr)
    real(kind=CUSTOM_REAL), INTENT(IN)     :: Sa_xy(Nspr), Sa_xz(Nspr), Sa_yz(Nspr)
    real(kind=CUSTOM_REAL), INTENT(INOUT)  :: F1 (Nspr)
    integer, intent(IN) :: aktif

    integer :: j,m,k
    real(kind=CUSTOM_REAL) :: ss(6), Sa(Nspr,6)

    ss(1) = 1e0_CUSTOM_REAL
    ss(2) = 1e0_CUSTOM_REAL
    ss(3) = 2e0_CUSTOM_REAL
    ss(4) = 2e0_CUSTOM_REAL
    ss(5) = 2e0_CUSTOM_REAL
    ss(6) = 1e0_CUSTOM_REAL

    Sa(:,1) = Sa_xx
    Sa(:,2) = Sa_yy
    Sa(:,3) = Sa_xy
    Sa(:,4) = Sa_xz
    Sa(:,5) = Sa_yz
    Sa(:,6) = Sa_zz


    j = start
    do while ( j < aktif+1 )

        do m = 1,6
            do k = 1,6    
                Ed(m,k) = Ed(m,k)+ CNinv(j)* ss(k)* (S1(m)- Sa(j,m)) &
                          * (S1(k)- Sa(j,k))/ (2.0* F1(j))
            enddo
        enddo
        j = j+1

    enddo

end subroutine MAT_IWAN_Ematris
!=====================================================================

subroutine LUDCMP(A,N,INDX,D,CODE)

    IMPLICIT NONE

    integer, parameter :: nmax = 100
    real(kind=CUSTOM_REAL), parameter :: tiny = 1.5e-16_CUSTOM_REAL

    real(kind=CUSTOM_REAL), intent(inout), dimension(N,N) :: A
    integer, intent(in) :: N
    integer, intent(out) :: D, CODE
    integer, intent(out), dimension(N) :: INDX
    REAL(kind=CUSTOM_REAL) :: AMAX, DUM, SUMM, VV(NMAX)
    INTEGER :: I,J,K,IMAX


    D=1; CODE=0

    DO I=1,N
        AMAX=0.d0
        DO J=1,N
            IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
        END DO ! j loop
        IF(AMAX.LT.TINY) THEN
            CODE = 1
            RETURN
        END IF
        VV(I) = 1.d0 / AMAX
    END DO ! i loop

    DO J=1,N
        DO I=1,J-1
            SUMM = A(I,J)
            DO K=1,I-1
                SUMM = SUMM - A(I,K)*A(K,J) 
            END DO ! k loop
            A(I,J) = SUMM
        END DO ! i loop
        AMAX = 0.d0
        DO I=J,N
            SUMM = A(I,J)
            DO K=1,J-1
                SUMM = SUMM - A(I,K)*A(K,J) 
            END DO ! k loop
            A(I,J) = SUMM
            DUM = VV(I)*ABS(SUMM)
            IF(DUM.GE.AMAX) THEN
                IMAX = I
                AMAX = DUM
            END IF
        END DO ! i loop  

        IF(J.NE.IMAX) THEN
            DO K=1,N
                DUM = A(IMAX,K)
                A(IMAX,K) = A(J,K)
                A(J,K) = DUM
            END DO ! k loop
            D = -D
            VV(IMAX) = VV(J)
        END IF

        INDX(J) = IMAX
        IF(ABS(A(J,J)) < TINY) A(J,J) = TINY

        IF(J.NE.N) THEN
            DUM = 1.d0 / A(J,J)
            DO I=J+1,N
                A(I,J) = A(I,J)*DUM
            END DO ! i loop
        END IF 
    END DO ! j loop

    RETURN

END subroutine LUDCMP
!=====================================================================
!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************

subroutine LUBKSB(A, N, INDX, B)

!     implicit none 

    integer, intent(in) :: N 
    real(kind=CUSTOM_REAL), intent(in), dimension(N,N) :: A
    integer, intent(in), dimension(N) :: INDX
    real(kind=CUSTOM_REAL), intent(inout), dimension(N) :: B

    REAL(kind=CUSTOM_REAL)::  SUMM
    integer :: II, LL,J,I

    II = 0

    DO I=1,N
        LL = INDX(I)
        SUMM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
            DO J=II,I-1
            SUMM = SUMM - A(I,J)*B(J)
            END DO ! j loop
        ELSE IF(SUMM.NE.0.d0) THEN
            II = I
        END IF
        B(I) = SUMM
    END DO ! i loop

    DO I=N,1,-1
        SUMM = B(I)
        IF(I < N) THEN
            DO J=I+1,N
                SUMM = SUMM - A(I,J)*B(J)
            END DO ! j loop
        END IF
        B(I) = SUMM / A(I,I)
    END DO ! i loop

    RETURN

END subroutine LUBKSB
!=====================================================================

end module nonlinear_solver_iwan
