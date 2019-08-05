module nonlinear_solver_iwan

    implicit none

    public :: MAT_IWAN_init, compute_nonlinear_stress

! I allocated every variable in specfem_par;
! maybe I should create a new file just for this module's parameters
!

contains
!=====================================================================

subroutine MAT_IWAN_init()

    use constants, only: NGLLZ,NGLLX,NGLLY,CUSTOM_REAL

    use specfem_par, only: NSPEC_AB, mustore

    use specfem_par_elastic, only: ispec_is_elastic, &
                    R_NL, CNinv_NL, S_NL, F_NL, n_active_surface, &
                    Sa_NL_xx, Sa_NL_yy, Sa_NL_zz, &
                    Sa_NL_xy, Sa_NL_yz, Sa_NL_xz, &
                    NSPR 

    implicit none

    integer :: i,j,k,ispec
    double precision :: Gm0, gamma_ref


    ! shifts shear moduli
    do ispec = 1,NSPEC_AB

        ! skips non elastic elements
        if (ispec_is_elastic(ispec) .eqv. .false.) cycle

        ! determines XXX for each GLL point
        do k = 1,NGLLZ
            do j = 1,NGLLY
                do i = 1,NGLLX        

                    R_NL(:,i,j,k,ispec) = 0d0
                    CNinv_NL(:,i,j,k,ispec) = 0d0
                    S_NL(:,i,j,k,ispec) = 1d-3
                    F_NL(:,i,j,k,ispec) = 0d0

                    n_active_surface(i,j,k,ispec) = -1

                    Sa_NL_xx(:,i,j,k,ispec)= 0d0
                    Sa_NL_yy(:,i,j,k,ispec)= 0d0
                    Sa_NL_zz(:,i,j,k,ispec)= 0d0
                    Sa_NL_xy(:,i,j,k,ispec)= 0d0
                    Sa_NL_xz(:,i,j,k,ispec)= 0d0
                    Sa_NL_yz(:,i,j,k,ispec)= 0d0

                    ! construct the backbone curve
                    gamma_ref = 0.000365d0
                    Gm0 = dble(mustore(i,j,k,ispec))

                    call MAT_IWAN_backbone_elem (gamma_ref, NSPR, Gm0, &
                                R_NL(:,i,j,k,ispec),CNinv_NL(:,i,j,k,ispec) ) 


                enddo
            enddo
        enddo

    enddo

end subroutine MAT_IWAN_init
!=====================================================================

subroutine MAT_IWAN_backbone_elem(gref,Nspr,mu,R,CNinv)

    implicit none

    integer, intent(in) :: Nspr
    double precision, intent(in) :: gref, mu
    double precision, intent(inout) :: R(Nspr), CNinv(Nspr-1)

    double precision :: x0, xu, dx, gamma(Nspr), G(Nspr), summy
    integer :: i

    x0 = -6d0
    xu = log10(0.1d0)
    dx = (xu- x0)/ (Nspr- 1)

    gamma  = 0d0
    G      = 0d0
    R      = 0d0
    CNinv  = 0d0

    do i = 1,Nspr
        gamma (i) = 10d0** ( x0+ dx* (i-1) )
        G(i)    = 1d0/ ( 1d0+ abs( gamma(i)/ gref) )
        R(i)    = G(i)* mu* gamma(i)
    enddo
    
    summy = 0d0
    do i = 1,Nspr-1
        CNinv(i) = ( gamma(i+1)/2d0 - gamma(i)/2d0 )/ &
                        ( R(i+1)- R(i) )- 0.5d0/mu- summy
        summy = summy+ CNinv(i)
    enddo   

end subroutine MAT_IWAN_backbone_elem
!=====================================================================

subroutine  compute_nonlinear_stress (Nspr, n_active_surface, &
            S, F, R, CNinv, Sa_xx, Sa_yy, Sa_zz, Sa_xy, Sa_xz, Sa_yz, &
            mul,lambdal,lambdalplus2mul, &
            duxdxl_NL,duydyl_NL,duzdzl_NL, &
            duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL, &
            duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL, &
             duxdyl_NL, duxdzl_NL, duydzl_NL, &
            dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz)

    use constants, only: CUSTOM_REAL

    implicit none

    integer, intent(in) :: Nspr
    integer, intent(inout) :: n_active_surface
    double precision, intent(inout) :: S(6)
    double precision, intent(inout) :: F(Nspr)
    double precision, intent(inout) :: R(Nspr)
    double precision, intent(inout) :: CNinv(Nspr-1)
    double precision, intent(inout) :: Sa_xx(Nspr), Sa_yy(Nspr), Sa_zz(Nspr)
    double precision, intent(inout) :: Sa_xy(Nspr), Sa_xz(Nspr), Sa_yz(Nspr)

    real(kind=CUSTOM_REAL), intent(in) :: lambdal,mul,lambdalplus2mul
    real(kind=CUSTOM_REAL), intent(in) :: duxdxl_NL, duydyl_NL, duzdzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_NL, duxdzl_NL, duydzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL    
    double precision, intent(out) :: dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz


    ! local parameters
    integer :: surface, errorflag, start, k, j, D
    integer :: INDX(6)
    double precision :: Gm0, ni    
    double precision :: lambda_nl, lambdalplus2mu_nl, Kmod
    double precision :: dsigm, depsm
    double precision:: dS(6), de(6)
    double precision:: Ed(6,6), Esd(6,6)
    double precision :: dF(Nspr)
    double precision :: dummy, dummy2, infinity


    
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

    Gm0   = dble(mul)
    ni    = 0.3875d0 
    Kmod = 2d0* Gm0* (1d0+ ni) ! E (Young modulus) 
    lambda_nl = Kmod* ni/ (1d0+ ni)/ (1d0- 2d0* ni)
    lambdalplus2mu_nl = lambda_nl+ 2d0* Gm0
    Kmod = Kmod/ (3d0* (1d0- 2d0* ni)) ! K (Bulk modulus)


    ! mean stress
    dsigm = dble(duxdxl_NL+ duydyl_NL+ duzdzl_NL)* Kmod

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

        Sa_xx = 0d0
        Sa_yy = 0d0
        Sa_zz = 0d0
        Sa_xy = 0d0
        Sa_xz = 0d0
        Sa_yz = 0d0

        return
    endif


    de    = 0d0
    depsm = dble(duxdxl_NL+ duydyl_NL+ duzdzl_NL)/ 3d0
    ! Incremental deviatoric strain
    de(1) = duxdxl_NL- depsm
    de(2) = duydyl_NL- depsm
    de(6) = duzdzl_NL- depsm
    !!! check here !!!   
    de(3) = duxdyl_NL
    de(4) = duxdzl_NL
    de(5) = duydzl_NL

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

            if ( (dF(j) >= 0d0)  .AND. (F(j) >= R(j)** 2) ) &
            surface = surface+ 1  

        enddo
    endif


    if ( (dF(1) >= 0d0)  .AND. (F(1) < R(1)** 2) ) then

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
    Ed    = 0d0
    start = 1


    ! elastic part 
    dummy = 0.5d0/ Gm0
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

        if ( (dF(j) >= 0d0)  .AND. (F(j) >= R(j)** 2) ) then
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
    
    use constants, only: CUSTOM_REAL

    implicit none

    double precision, intent(in) :: lambdal,Gact,lambdalplus2mul
    real(kind=CUSTOM_REAL), intent(in) :: duxdxl_NL, duydyl_NL, duzdzl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duydyl_plus_duzdzl_NL, duxdxl_plus_duzdzl_NL, duxdxl_plus_duydyl_NL
    real(kind=CUSTOM_REAL), intent(in) :: duxdyl_plus_duydxl_NL, duzdxl_plus_duxdzl_NL, duzdyl_plus_duydzl_NL
    double precision, intent(out) ::  dsigma_xx, dsigma_yy, dsigma_zz, dsigma_xy, dsigma_xz, dsigma_yz


    dsigma_xx = lambdalplus2mul * dble(duxdxl_NL) + lambdal * dble(duydyl_plus_duzdzl_NL)
    dsigma_yy = lambdalplus2mul * dble(duydyl_NL) + lambdal * dble(duxdxl_plus_duzdzl_NL)
    dsigma_zz = lambdalplus2mul * dble(duzdzl_NL) + lambdal * dble(duxdxl_plus_duydyl_NL)
    dsigma_xy = Gact * dble(duxdyl_plus_duydxl_NL)
    dsigma_xz = Gact * dble(duzdxl_plus_duxdzl_NL)
    dsigma_yz = Gact * dble(duzdyl_plus_duydzl_NL)

end subroutine MAT_IWAN_elastic
!=====================================================================
! Iwan surface(s) computation 

function MAT_IWAN_surface(S, Sa_xx, Sa_yy, Sa_zz, Sa_xy, Sa_xz, Sa_yz)  result(F)

    implicit none

    double precision, intent(in) :: S(6)
    double precision, intent(inout) :: Sa_xx, Sa_yy, Sa_zz
    double precision, intent(inout) :: Sa_xy, Sa_xz, Sa_yz
    double precision :: F

    F = 0d0
    F = 0.5d0* ( (S(1)-Sa_xx )** 2+ ( S(2)- Sa_yy )** 2 &
                  + 2d0* ( S(3)- Sa_xy )** 2  &
                  + 2d0* ( S(4)- Sa_xz )** 2  &
                  + 2d0* ( S(5)- Sa_yz )** 2  &
                  + ( S(6)- Sa_zz )** 2 )


end function MAT_IWAN_surface  
!=====================================================================
! Iwan surface(s) movement computation 

function MAT_IWAN_dsurface(S,Sa_xx,Sa_yy,Sa_zz,Sa_xy,Sa_xz,Sa_yz,de) result(dF)

    implicit none

    double precision, intent(in) :: S(6)
    double precision, intent(inout) :: Sa_xx, Sa_yy, Sa_zz
    double precision, intent(inout) :: Sa_xy, Sa_xz, Sa_yz
    double precision, intent(in) :: de(6)
    double precision :: dF

    dF = 0d0
    dF = 0.5d0* ( ( S(1)- Sa_xx )* de(1)+ ( S(2)- Sa_yy )* de(2)+ &
                  2d0* ( S(3)- Sa_xy )* de(3) + &
                  2d0* ( S(4)- Sa_xz )* de(4) + &
                  2d0* ( S(5)- Sa_yz )* de(5) + &
                  ( S(6)- Sa_zz )* de(6))

end function MAT_IWAN_dsurface  
!=======================================================================
! Iwan plasticity matrix

subroutine MAT_IWAN_Ematris (start, Nspr, Ed, CNinv, S1, Sa_xx, Sa_yy, Sa_zz, Sa_xy, Sa_xz, Sa_yz, F1, aktif)

    implicit none

    integer, intent(IN) :: start
    integer, intent(IN) :: Nspr
    double precision, INTENT(INOUT)  :: Ed(6,6)
    double precision, INTENT(IN)     :: CNinv(Nspr-1)
    double precision, INTENT(IN)     :: S1(6)
    double precision, INTENT(IN)     :: Sa_xx(Nspr), Sa_yy(Nspr), Sa_zz(Nspr)
    double precision, INTENT(IN)     :: Sa_xy(Nspr), Sa_xz(Nspr), Sa_yz(Nspr)
    double precision, INTENT(INOUT)  :: F1 (Nspr)
    integer, intent(IN) :: aktif

    integer :: j,m,k
    double precision :: ss(6), Sa(Nspr,6)

    ss = 1d0
    ss(3) = 2d0
    ss(4) = 2d0
    ss(5) = 2d0

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
                          * (S1(k)- Sa(j,k))/ (2d0* F1(j))
            enddo
        enddo
        j = j+1

    enddo

end subroutine MAT_IWAN_Ematris
!=====================================================================

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************

Subroutine LUDCMP(A,N,INDX,D,CODE)
    
    IMPLICIT NONE
    
    integer, parameter :: nmax = 100
    double precision, parameter :: tiny = 1.5D-16
    !real, parameter :: tiny = 1.5D-16    

    double precision, intent(inout), dimension(N,N) :: A
    integer, intent(in) :: N
    integer, intent(out) :: D, CODE
    integer, intent(out), dimension(N) :: INDX
    !f2py depend(N) A, indx

    double precision  :: AMAX, DUM, SUMM, VV(NMAX)
    INTEGER :: I,J,K,IMAX

    D=1; CODE=0

    DO I=1,N
    AMAX=0.d0
    DO J=1,N
     IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
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
     DUM = VV(I)*DABS(SUMM)
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
    IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

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

    integer, intent(in) :: N 
    double precision, intent(in), dimension(N,N) :: A
    integer, intent(in), dimension(N) :: INDX
    double precision, intent(inout), dimension(N) :: B
    !f2py depend(N) A, INDX, B

    double precision ::  SUMM
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