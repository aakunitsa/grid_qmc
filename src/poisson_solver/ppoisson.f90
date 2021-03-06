MODULE POISSON
! Note: the module should be compiled and linked togather with 
!       the standalone lebedev angular grid code && LAPACK

IMPLICIT NONE

contains

SUBROUTINE INITIALIZE_POISSON(NRADIAL, NANGULAR, ATOMIC_NUMBER) bind (C, name = 'initialize_poisson_')
! Initializes the global variables and constructs
! grids

   use iso_c_binding

   USE SINANOGLU
   USE STRUCTURE
   USE DFT

   IMPLICIT NONE 
   INTEGER, INTENT(IN) :: NRADIAL, NANGULAR, ATOMIC_NUMBER
   integer:: i, j, k

   ! The module will be able to perform atomic calculations only

   NATOM = 1
   CHARGE = 0
   ATOMX(1) = 0._8
   ATOMY(1) = 0._8
   ATOMZ(1) = 0._8
   IATOM(1) = ATOMIC_NUMBER

   ! Grid construction

   allocate(ngrid(1:natom))
   allocate(gridx(1:nradial* nangular, 1:natom))
   allocate(gridy(1:nradial * nangular, 1:natom))
   allocate(gridz(1:nradial * nangular, 1:natom))
   allocate(gridw(1:nradial * nangular, 1:natom))

   allocate(sgn(-80:80))

   do i = -80, 80
      sgn(i) = (-1.0D0)**i
   end do
   
   CALL SG_CONSTRUCT_GRID(NRADIAL, NANGULAR)

END SUBROUTINE

subroutine finalize_poisson() bind (C, name = 'finalize_poisson_')

    use iso_c_binding

    use dft
    use sinanoglu

    deallocate (ngrid, gridw, gridx, gridy, gridz, sgn)

    call sg_destroy_grid

end subroutine 

DOUBLE PRECISION FUNCTION ERI(Ci, Cj, Ck, Cl) bind (C, name = 'eri_fortran_')

   ! My very own Fortran function that calculates ERI-s
   ! Created for benchmarking purposes

   use iso_c_binding 
   USE CONSTANTS
   USE DFT
   USE STRUCTURE
   USE SINANOGLU

   implicit none

   double precision, intent(in) :: Ci(maxngrid), Cj(maxngrid), Ck(maxngrid), Cl(maxngrid)

   integer :: i, ia
   double precision :: den_ij(1:maxngrid,1), pot_ij(1:maxngrid,1), den_kl(1:maxngrid,1)

   den_ij(:,1) = Ci * Cj
   den_kl(:, 1) = Ck * Cl
   pot_ij(:, 1) = 0.0D0
   eri = 0.0D0

   call construct_potential(den_ij, pot_ij)

   DO IA=1,NATOM
    DO I=1,NGRID(IA)
     eri=eri+pot_ij(I,IA)*den_kl(I,IA)*GRIDW(I,IA)
    ENDDO
   ENDDO

END FUNCTION ERI

SUBROUTINE CONSTRUCT_potential(ORB1,ORB2) bind (C, name = 'construct_potential_')
! Construct a Coulomb potential by solving Poisson's equation

   use iso_c_binding

   USE CONSTANTS
   USE DFT
   USE STRUCTURE
   USE SINANOGLU

   IMPLICIT NONE
   DOUBLE PRECISION :: ORB1(MAXNGRID,NATOM)
   DOUBLE PRECISION :: ORB2(MAXNGRID,NATOM)
   DOUBLE PRECISION :: ORB3(MAXNGRID), ORB4(MAXNGRID)
   DOUBLE COMPLEX   :: RTEMP(0:LMAX,-LMAX:LMAX,RG)
   DOUBLE COMPLEX   :: W1(RG)
   DOUBLE PRECISION :: W2(RG,RG),W3(RG,2)
#if defined(EXPERT_LAPACK)
   DOUBLE PRECISION :: W2F(RG,RG), C(RG), R(RG), B(RG, 2)
   DOUBLE PRECISION :: RCOND, FERR(2), BERR(2)
   DOUBLE PRECISION :: WRK(4*RG)
   INTEGER          :: IWRK(RG)
   DOUBLE PRECISION, PARAMETER :: ERR_THRESH = 1E-8
   LOGICAL :: ERR_OCCURRED
   CHARACTER*1 :: EQUIL
#endif
   INTEGER :: I,J,L,M,INFO,IA,IB
   INTEGER :: INDX(RG)
!   double precision :: err, r
!   double complex :: ylm_val

   !print*, 'Maxngrid = ', MAXNGRID
   !print*, 'Lmax = ', LMAX
   !print*, NATOM

   !print*, 'Spherical harmonic table for the grid'
   !do ia = 1, natom
   !  do L = 0, LMAX
   !    do m = -L, L
   !     write(*, *)  L, m
   !     do i = 1, ngrid(ia)
   !       call RYLM(l,m,i,r,ylm_val, ia, ia)
   !       write(*, '(2f28.20)') dreal(ylm_val), dimag(ylm_val) 
   !     end do
   !    end do
   !  end do
   !end do

   ORB2=0.0D0
   DO IA=1,NATOM
    DO I=1,NGRID(IA)
     ORB3(I)=ORB1(I,IA) ! FUZZY has been removed
    ENDDO
    CALL EXPAND2YLM(ORB3,RTEMP,IA,IA)
    DO L=0,LMAX
     DO M=-L,L
      DO I=1,RG
       W1=DCMPLX(0.0D0,0.0D0)
       W1(I)=DCMPLX(1.0D0,0.0D0)
       CALL SECOND_DERIV2(W1,IA)
       !print*, 'applied second_deriv2'
       DO J=1,RG
        W2(J,I)=DREAL(W1(J))
       ENDDO
       W2(I,I)=W2(I,I)-DFLOAT(L*(L+1))/SG_GAUSS_CHEV(I,IA)**2
      ENDDO
      DO I=1,RG
       W3(I,1)=-4.0D0*PI*SG_GAUSS_CHEV(I,IA)*DREAL(RTEMP(L,M,I))
       W3(I,2)=-4.0D0*PI*SG_GAUSS_CHEV(I,IA)*DIMAG(RTEMP(L,M,I))
      ENDDO
#if defined(EXPERT_LAPACK)
      ERR_OCCURRED = .FALSE.
      B = W3 ! RHS; Solution will be saved to W3
      CALL DGESVX('E','N',RG,2,W2,RG,W2F,RG,INDX,EQUIL,R,C,B,RG,W3,RG,RCOND,FERR,BERR,WRK,IWRK,INFO)

      if ((FERR(1) .gt. ERR_THRESH) .or. (BERR(1) .gt. ERR_THRESH)) THEN
          ! Issue warning
          ERR_OCCURRED = .TRUE.
          write(*, *) '...Error threshold exceeded in poisson solver!'
          write(*, '(A, 1x, E20.10)') '...Error estimate (real part)', max(ferr(1), berr(1))
      end if

      if ((FERR(2) .gt. ERR_THRESH) .or. (BERR(2) .gt. ERR_THRESH)) THEN
          ! Issue warning
          ERR_OCCURRED = .TRUE.
          write(*, *) '...Error threshold exceeded in poisson solver!'
          write(*, '(A, 1x, E20.10)') '...Error estimate (imaginary part)', max(ferr(2), berr(2))
      end if

      if (ERR_OCCURRED) THEN
          write(*,*) 'Condition number of the linear system is ', RCOND
      end if
#else
      CALL DGESV(RG,2,W2,RG,INDX,W3,RG,INFO)
#endif
      !IF (INFO /= 0) CALL PABORT('DGESV FAILED')
      IF (INFO /= 0) stop
      DO I=1,RG
       RTEMP(L,M,I)=DCMPLX(W3(I,1),W3(I,2))/SG_GAUSS_CHEV(I,IA)
      ENDDO
     ENDDO
    ENDDO
!write(*,*) 'U(0,0,~inf)=',rtemp(0,0,1)*sg_gauss_chev(1,ia)
!write(*,*) 'U(0,0, inf)=',dsqrt(4.0d0*pi)*r1
    !print*, 'packing the density'
    DO IB=1,NATOM
     CALL PACKYLM(ORB3,RTEMP,IB,IA)
     DO I=1,NGRID(IB)
      ORB2(I,IB)=ORB2(I,IB)+ORB3(I)
     ENDDO
    ENDDO

   ENDDO
   RETURN

END SUBROUTINE


SUBROUTINE SECOND_DERIV2(R,IA)
! 7-point finite-difference second differentiation
! (1/r)(d2/dr2)r r^-1f(r) = (1/r)(d2x/dr2)(df/dx) + (1/r)(dx/dr)^2(d2f/dx2)

   USE CONSTANTS
   USE STRUCTURE
   USE SINANOGLU

   IMPLICIT NONE
   INTEGER :: IA
   DOUBLE COMPLEX   :: R(RG)
   DOUBLE PRECISION :: H
   DOUBLE COMPLEX,ALLOCATABLE :: R1D(:)
   DOUBLE COMPLEX,ALLOCATABLE :: R2D(:)
   INTEGER :: I
   DOUBLE PRECISION,PARAMETER :: BSRADI(17) = &
     (/0.472D0,0.472D0,1.370D0,0.992D0,0.803D0,0.661D0,0.614D0,0.567D0,0.472D0, &
       0.472D0,1.701D0,1.417D0,1.181D0,1.039D0,0.945D0,0.945D0,0.945D0/) ! Ne DATA WAS MISSING IN THE ORIGINAL
   DOUBLE PRECISION :: R0,X,OMEGA

   ALLOCATE(R1D(RG),R2D(RG))
   H=PI/DFLOAT(RG+1)

   ! FIRST DERIVATIVES
   R1D(1)=(-1764.0D0*R(1)+4320.0D0*R(2)-5400.0D0*R(3)+4800.0D0*R(4) &
           -2700.0D0*R(5)+864.0D0*R(6)-120.0D0*R(7))/(720.0D0*H)
   R1D(2)=(-120.0D0*R(1)-924.0D0*R(2)+1800.0D0*R(3)-1200.0D0*R(4) &
           +600.0D0*R(5)-180.0D0*R(6)+24.0D0*R(7))/(720.0D0*H)
   R1D(3)=(24.0D0*R(1)-288.0D0*R(2)-420.0D0*R(3)+960.0D0*R(4) &
           -360.0D0*R(5)+96.0D0*R(6)-12.0D0*R(7))/(720.0D0*H)
   R1D(4)=(-12.0D0*R(1)+108.0D0*R(2)-540.0D0*R(3) &
          +540.0D0*R(5)-108.0D0*R(6)+12.0D0*R(7))/(720.0D0*H)
   R1D(RG-3)=(-12.0D0*R(RG-6)+108.0D0*R(RG-5)-540.0D0*R(RG-4) &
             +540.0D0*R(RG-2)-108.0D0*R(RG-1)+12.0D0*R(RG))/(720.0D0*H)
   R1D(RG-2)=(-12.0D0*R(RG-5)+108.0D0*R(RG-4)-540.0D0*R(RG-3) &
             +540.0D0*R(RG-1)-108.0D0*R(RG))/(720.0D0*H)
   R1D(RG-1)=(+12.0D0*R(RG-5)-96.0D0*R(RG-4)+360.0D0*R(RG-3)-960.0D0*R(RG-2) &
             +420.0D0*R(RG-1)+288.0D0*R(RG))/(720.0D0*H)
   R1D(RG)=(-24.0D0*R(RG-5)+180.0D0*R(RG-4)-600.0D0*R(RG-3)+1200.0D0*R(RG-2) &
            -1800.0D0*R(RG-1)+924.0D0*R(RG))/(720.0D0*H)
   DO I=5,RG-4
    R1D(I)=(+144.0D0*R(I-4)-1536.0D0*R(I-3)+8064.0D0*R(I-2)-32256.0D0*R(I-1) &
           +32256.0D0*R(I+1)-8064.0D0*R(I+2)+1536.0D0*R(I+3)-144.0D0*R(I+4))/(40320.0D0*H)
   ENDDO

   ! SECOND DERIVATIVES
   R2D(1)=(1624.0D0*R(1)-6264.0D0*R(2)+10530.0D0*R(3)-10160.0D0*R(4) &
           +5940.0D0*R(5)-1944.0D0*R(6)+274.0D0*R(7))/(360.0D0*H**2)
   R2D(2)=(274.0D0*R(1)-294.0D0*R(2)-510.0D0*R(3)+940.0D0*R(4) &
           -570.0D0*R(5)+186.0D0*R(6)-26.0D0*R(7))/(360.0D0*H**2)
   R2D(3)=(-26.0D0*R(1)+456.0D0*R(2)-840.0D0*R(3)+400.0D0*R(4) &
           +30.0D0*R(5)-24.0D0*R(6)+4.0D0*R(7))/(360.0D0*H**2)
   R2D(4)=(4.0D0*R(1)-54.0D0*R(2)+540.0D0*R(3)-980.0D0*R(4) &
           +540.0D0*R(5)-54.0D0*R(6)+4.0D0*R(7))/(360.0D0*H**2)
   R2D(RG-3)=(4.0D0*R(RG-6)-54.0D0*R(RG-5)+540.0D0*R(RG-4)-980.0D0*R(RG-3) &
             +540.0D0*R(RG-2)-54.0D0*R(RG-1)+4.0D0*R(RG))/(360.0D0*H**2)
   R2D(RG-2)=(+4.0D0*R(RG-5)-54.0D0*R(RG-4)  +540.0D0*R(RG-3)-980.0D0*R(RG-2) &
             +540.0D0*R(RG-1)-54.0D0*R(RG))/(360.0D0*H**2)
   R2D(RG-1)=(+4.0D0*R(RG-5)-24.0D0*R(RG-4)   +30.0D0*R(RG-3)+400.0D0*R(RG-2) &
             -840.0D0*R(RG-1)+456.0D0*R(RG))/(360.0D0*H**2)
   R2D(RG)=(-26.0D0*R(RG-5)+186.0D0*R(RG-4)  -570.0D0*R(RG-3)+940.0D0*R(RG-2) &
            -510.0D0*R(RG-1)-294.0D0*R(RG))/(360.0D0*H**2)
   DO I=5,RG-4
    R2D(I)=(-36.0D0*R(I-4) +512.0D0*R(I-3)-4032.0D0*R(I-2)+32256.0D0*R(I-1)-57400.0D0*R(I) &
            +32256.0D0*R(I+1)-4032.0D0*R(I+2)+512.0D0*R(I+3)-36.0D0*R(I+4))/(20160.0D0*H**2)
   ENDDO

   R0=BSRADI(IATOM(IA))
   DO I=1,RG
     OMEGA=DFLOAT(I)*PI/DFLOAT(RG+1)
     X=DCOS(OMEGA)
     R(I)=+R1D(I)*((1.0D0-X)**3/(2.0D0*R0**2*DSIN(OMEGA)) &
          -(1.0D0-X)**4*DCOS(OMEGA)/(4.0D0*R0**2*DSIN(OMEGA)**3)) &
          +R2D(I)*((1.0D0-X)**2/(2.0D0*R0*DSIN(OMEGA)))**2
   ENDDO

   DEALLOCATE(R2D,R1D)
   RETURN
END SUBROUTINE



SUBROUTINE EXPAND2YLM(F,R,DEST,ORIG)
! Expand a given function as a linear combination of radial spherical functions

   USE CONSTANTS
   USE SINANOGLU

   IMPLICIT NONE
   INTEGER :: DEST,ORIG
!   DOUBLE PRECISION :: F(AG,RG)
   DOUBLE PRECISION :: F(AG * RG)
   DOUBLE COMPLEX   :: YLM
   DOUBLE COMPLEX   :: R(0:LMAX,-LMAX:LMAX,RG)
   DOUBLE PRECISION :: R1
   INTEGER :: I,J,K,L,M

   !F(AG,RG) = 1.0D0
   !R(LMAX,LMAX,RG) = 1.0D0
   !print*, 'Inside expand2ylm'
   !print*, RG, AG

   DO L=0,LMAX       ! LOOP OVER L
    DO M=-L,L        ! LOOP OVER M
     K=0
     DO I=1,RG       ! LOOP OVER RADIAL GRID POINTS
      R(L,M,I)=DCMPLX(0.0D0,0.0D0)
      DO J=1,AG      ! LOOP OVER ANGULAR GRID POINTS
       K=K+1
       CALL RYLM(L,M,K,R1,YLM,DEST,ORIG) ! GET R & YLM AT THE POINT
       !R(L,M,I)=R(L,M,I)+DCONJG(YLM)*F(J,I)*SG_LEBEDVW(J)
       R(L,M,I)=R(L,M,I)+DCONJG(YLM)*F(J + (I - 1) * AG )*SG_LEBEDVW(J)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   RETURN
END SUBROUTINE


SUBROUTINE PACKYLM(F,R,DEST,ORIG)
! Expand a linear combination of radial spherical functions to real space

   USE CONSTANTS
   USE SINANOGLU

   IMPLICIT NONE
   INTEGER,PARAMETER :: INTORDER=4
   INTEGER :: DEST,ORIG
   DOUBLE PRECISION :: F(AG,RG)
   DOUBLE COMPLEX   :: YLM
   DOUBLE COMPLEX   :: R(0:LMAX,-LMAX:LMAX,RG)
   DOUBLE PRECISION :: R1
   DOUBLE PRECISION :: RREAL,RIMAG,DY
   DOUBLE PRECISION :: X(RG),Y(RG)
   INTEGER :: I,J,K,L,M,N,P

   IF (DEST==ORIG) THEN
   ! SIMPLE EXPANSION
   
    K=0
    DO I=1,RG       ! LOOP OVER RADIAL GRID POINTS
     DO J=1,AG      ! LOOP OVER ANGULAR GRID POINTS
      K=K+1
      F(J,I)=0.0D0
      DO L=0,LMAX   ! LOOP OVER L
       DO M=-L,L    ! LOOP OVER M
        CALL RYLM(L,M,K,R1,YLM,DEST,ORIG) ! GET R & YLM AT THE POINT
        F(J,I)=F(J,I)+R(L,M,I)*YLM
       ENDDO
      ENDDO
     ENDDO
    ENDDO

   ELSE
   ! RATIONAL INTERPOLATION

    X=SG_GAUSS_CHEV(:,ORIG)
    K=0
    DO I=1,RG       ! LOOP OVER RADIAL GRID POINTS OF ATOM IB
     DO J=1,AG      ! LOOP OVER ANGULAR GRID POINTS OF ATOM IB
      K=K+1
      F(J,I)=0.0D0
      DO L=0,LMAX   ! LOOP OVER L CENTERED AT ATOM IA
       DO M=-L,L    ! LOOP OVER M CENTERED AT ATOM IA
        CALL RYLM(L,M,K,R1,YLM,DEST,ORIG) ! GET R & YLM AT THE POINT
        P=1
        DO N=1,RG
         IF (SG_GAUSS_CHEV(N,ORIG) > R1) P=N
        ENDDO
        P=MAX(1,MIN(P-1,RG-3))
        Y=DREAL(R(L,M,:))
        CALL RATINT(X(P),Y(P),INTORDER,R1,RREAL,DY)
        Y=DIMAG(R(L,M,:))
        CALL RATINT(X(P),Y(P),INTORDER,R1,RIMAG,DY)
        F(J,I)=F(J,I)+DCMPLX(RREAL,RIMAG)*YLM
       ENDDO
      ENDDO
     ENDDO
    ENDDO

   ENDIF

   RETURN
END SUBROUTINE

FUNCTION plgndr(l,m,x)
! Returns associated Legendre polynomial value
! (C) Copr. 1986-92 Numerical Recipes Software t.)-5i.

   IMPLICIT NONE
   INTEGER :: l,m
   DOUBLE PRECISION :: plgndr,x
   INTEGER :: i,ll
   DOUBLE PRECISION :: fact,pll,pmm,pmmp1,somx2

   !if (m < 0.or.m > l.or.dabs(x) > 1.0d0) call pabort('bad arguments in plgndr')
   pmm=1.0d0
   if(m > 0) then
     somx2=dsqrt((1.0d0-x)*(1.0d0+x))
     fact=1.0d0
     do i=1,m
       pmm=-pmm*fact*somx2
       fact=fact+2.0d0
     enddo
   endif
   if(l == m) then
     plgndr=pmm
   else
     pmmp1=x*dfloat(2*m+1)*pmm
     if(l == m+1) then
       plgndr=pmmp1
     else
       do ll=m+2,l
         pll=(x*dfloat(2*ll-1)*pmmp1-dfloat(ll+m-1)*pmm)/dfloat(ll-m)
         pmm=pmmp1
         pmmp1=pll
       enddo
       plgndr=pll
     endif
   endif
   return
END FUNCTION



SUBROUTINE RYLM(L,M,I,R,YLM,DEST,ORIG)
! Returns the value of radius and real spherical harmonics of the origin atom
! at a grid point of the destination atom

   USE CONSTANTS
   USE SINANOGLU

   IMPLICIT NONE
   DOUBLE COMPLEX :: YLM
   INTEGER :: L,M
   INTEGER :: I,J
   INTEGER :: DEST,ORIG
   DOUBLE PRECISION :: R,THETA,PHI
   DOUBLE PRECISION :: NORM

   CALL C2SPHERICAL(R,THETA,PHI,I,DEST,ORIG)
   YLM=PLGNDR(L,IABS(M),DCOS(THETA))*DCMPLX(DCOS(DFLOAT(IABS(M))*PHI),DSIN(DFLOAT(IABS(M))*PHI))
   NORM = DFLOAT(2*L+1)/4.0D0/PI
   IF (IABS(M) >= 1) THEN
    DO J=L-IABS(M)+1,L+IABS(M)
     NORM=NORM/DFLOAT(J)
    ENDDO
   ENDIF
   YLM=DSQRT(NORM)*YLM
   IF (M < 0) YLM=SGN(M)*DCONJG(YLM)
   RETURN
END SUBROUTINE


   
SUBROUTINE C2SPHERICAL(R,THETA,PHI,I,DEST,ORIG)
! Returns the spherical coordinates of the origin atom of a grid point of the destination atom

   USE CONSTANTS
   USE STRUCTURE
   USE DFT
   USE SINANOGLU
   
   IMPLICIT NONE
   INTEGER :: I
   INTEGER :: DEST,ORIG
   DOUBLE PRECISION :: X,Y,Z
   DOUBLE PRECISION :: R,THETA,PHI
   DOUBLE PRECISION :: RXY

   X=ATOMX(DEST)+GRIDX(I,DEST)-ATOMX(ORIG)
   Y=ATOMY(DEST)+GRIDY(I,DEST)-ATOMY(ORIG)
   Z=ATOMZ(DEST)+GRIDZ(I,DEST)-ATOMZ(ORIG)
   R=DSQRT(X**2+Y**2+Z**2)
   THETA=DACOS(Z/R)
   IF (X == 0.0D0) THEN
    IF (Y > 0.0D0) PHI=PI/2.0D0
    IF (Y < 0.0D0) PHI=-PI/2.0D0
   ELSE IF (Y == 0.0D0) THEN
    IF (X > 0.0D0) PHI=0.0D0
    IF (X < 0.0D0) PHI=-PI
   ELSE IF (X > 0.0D0) THEN
    PHI=DATAN(Y/X)
   ELSE
    PHI=DATAN(Y/X)+PI
   ENDIF
!   RXY = DSQRT(X**2 + Y**2) 
!   IF (RXY > TINY(RXY)) THEN
!       PHI = ACOS(X/RXY)
!   ELSE
!       PHI = 2*PI
!   ENDIF
!   IF (Y < TINY(Y)) THEN
!       PHI = 2.0D0 * PI - PHI
!   ENDIF

! DEBUG CODE ...
  IF (DABS(Z-R*DCOS(THETA))>1.0D-10) WRITE(*,*) 'Z',Z,R*DCOS(THETA)
  IF (DABS(X-R*DSIN(THETA)*DCOS(PHI))>1.0D-10) WRITE(*,*) 'X',X,R*DSIN(THETA)*DCOS(PHI), PHI/PI
  IF (DABS(Y-R*DSIN(THETA)*DSIN(PHI))>1.0D-10) WRITE(*,*) 'Y',Y,R*DSIN(THETA)*DSIN(PHI), PHI/PI
! ... END DEBUG
   RETURN
END SUBROUTINE




SUBROUTINE SG_CONSTRUCT_GRID(NRAD, NANG)
! CONSTRUCT GRID FOR NUMERICAL INTEGRATION IN SINANOGLU CALCULATIONS.
! THIS VERSION OF THE SUBROUTINE will only implement becke's grid as 
! opposed to the original Polymer code
 
   USE CONSTANTS
   USE STRUCTURE
   USE DFT
   USE SINANOGLU

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: NRAD, NANG
   DOUBLE PRECISION,PARAMETER :: BSRADI(17) = &
     (/0.472D0,0.472D0,1.370D0,0.992D0,0.803D0,0.661D0,0.614D0,0.567D0,0.472D0, &
       0.472D0,1.701D0,1.417D0,1.181D0,1.039D0,0.945D0,0.945D0,0.945D0/) ! Ne DATA WAS MISSING IN THE ORIGINAL
   INTEGER :: I,J,K,L,M,II,JJ,IA,JA
   DOUBLE PRECISION :: R1,R2,X,W
   DOUBLE PRECISION :: THETA,PHI
   REAL :: ICPUS,ICPUE

   RG = NRAD
   AG = NANG
!   IF (RG*AG == 0) CALL PABORT('SPECIFY RADIAL AND ANGULAR GRIDS')
!   IF (RG < 7) CALL PABORT('RADIAL GRID POINTS TOO FEW')
!  IF (NATOM /= 1) CALL PABORT('ATOMIC CALCULATIONS ONLY')

   ALLOCATE(SG_GAUSS_CHEV(RG,NATOM),SG_GAUSS_CHEV_W(RG,NATOM))
   ALLOCATE(SG_LEBEDVX(AG),SG_LEBEDVY(AG),SG_LEBEDVZ(AG),SG_LEBEDVW(AG))

   DO IA=1,NATOM
    R1=BSRADI(IATOM(IA))
     DO J=1,RG
      X=DCOS(DFLOAT(J)*PI/DFLOAT(RG+1))
      SG_GAUSS_CHEV(J,IA)=R1*(1.0D0+X)/(1.0D0-X)
      SG_GAUSS_CHEV_W(J,IA)=2.0D0*R1/(1.0D0-X)**2*PI/DFLOAT(RG+1)*DSIN(DFLOAT(J)*PI/DFLOAT(RG+1))*SG_GAUSS_CHEV(J,IA)**2
     ENDDO
   ENDDO

   IF (AG == 1) THEN
    LMAX=0
    SG_LEBEDVX(1)=1.0D0
    SG_LEBEDVY(1)=0.0D0
    SG_LEBEDVZ(1)=0.0D0
    SG_LEBEDVW(1)=1.0D0
   ELSE IF (AG == 6) THEN
    LMAX=1
    CALL LD0006(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 14) THEN
    LMAX=2
    CALL LD0014(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 26) THEN
    LMAX=3
    CALL LD0026(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 38) THEN
    LMAX=4
    CALL LD0038(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 50) THEN
    LMAX=5
    CALL LD0050(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 74) THEN
    LMAX=6
    CALL LD0074(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 86) THEN
    LMAX=7
    CALL LD0086(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 110) THEN
    LMAX=8
    CALL LD0110(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 146) THEN
    LMAX=9
    CALL LD0146(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 170) THEN
    LMAX=10
    CALL LD0170(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 194) THEN
    LMAX=11
    CALL LD0194(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 230) THEN
    LMAX=12
    CALL LD0230(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 266) THEN
    LMAX=13
    CALL LD0266(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 302) THEN
    LMAX=14
    CALL LD0302(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 350) THEN
    LMAX=(31-1)/2
    CALL LD0350(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 434) THEN
    LMAX=(35-1)/2
    CALL LD0434(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 590) THEN
    LMAX=(41-1)/2
    CALL LD0590(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 770) THEN
    LMAX=(47-1)/2
    CALL LD0770(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 974) THEN
    LMAX=(53-1)/2
    CALL LD0974(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 1202) THEN
    LMAX=(59-1)/2
    CALL LD1202(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 1454) THEN
    LMAX=(65-1)/2
    CALL LD1454(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 1730) THEN
    LMAX=(71-1)/2
    CALL LD1730(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 2030) THEN
    LMAX=(77-1)/2
    CALL LD2030(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 2354) THEN
    LMAX=(83-1)/2
    CALL LD2354(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 2702) THEN
    LMAX=(89-1)/2
    CALL LD2702(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 3074) THEN
    LMAX=(95-1)/2
    CALL LD3074(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 3470) THEN
    LMAX=(101-1)/2
    CALL LD3470(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 3890) THEN
    LMAX=(107-1)/2
    CALL LD3890(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 4334) THEN
    LMAX=(113-1)/2
    CALL LD4334(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 4802) THEN
    LMAX=(119-1)/2
    CALL LD4802(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 5294) THEN
    LMAX=(125-1)/2
    CALL LD5294(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE IF (AG == 5810) THEN
    LMAX=(131-1)/2
    CALL LD5810(SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW,AG)
   ELSE
!    CALL PABORT('LEBEDEV GRID NOT FOUND')
     write(*,*) 'LEBEDEV GRID NOT FOUND'
     stop
   ENDIF
   DO J=1,AG
    SG_LEBEDVW(J)=SG_LEBEDVW(J)*4.0D0*PI
   ENDDO

   ! SET GRID INFORMATION IN CONVENTIONAL DFT FASHION
   M=0
   DO IA=1,NATOM
    !WRITE(6,'(A,I0,A,A,A,I0,A,I0,A,I0,A)') 'GRID FOR ATOM ',IA,' ( ',CATOM(IATOM(IA)),' ) IS ',RG,' x ',AG,' (LMAX = ',LMAX,')'
    NGRID(IA)=RG*AG
    DO J=1,RG
     DO K=1,AG
      M=M+1
      GRIDX((J-1)*AG+K,IA)=SG_GAUSS_CHEV(J,IA)*SG_LEBEDVX(K)
      GRIDY((J-1)*AG+K,IA)=SG_GAUSS_CHEV(J,IA)*SG_LEBEDVY(K)
      GRIDZ((J-1)*AG+K,IA)=SG_GAUSS_CHEV(J,IA)*SG_LEBEDVZ(K)
      GRIDW((J-1)*AG+K,IA)=SG_GAUSS_CHEV_W(J,IA)*SG_LEBEDVW(K)
     ENDDO
    ENDDO
   ENDDO
   !WRITE(6,'(A,I0)') 'NUMBER OF GRID POINTS = ',M
   MAXNGRID = M

   RETURN
END SUBROUTINE



SUBROUTINE SG_DESTROY_GRID
! DEALLOCATE GRID RELATED ARRAYS

   USE SINANOGLU

   IMPLICIT NONE

   DEALLOCATE(SG_GAUSS_CHEV,SG_GAUSS_CHEV_W,SG_LEBEDVX,SG_LEBEDVY,SG_LEBEDVZ,SG_LEBEDVW)

   RETURN
END SUBROUTINE



DOUBLE PRECISION FUNCTION WIGNER3J(L1,L2,L3,M1,M2,M3)
! RETURNS WIGNER 3J SYMBOL NUMERICAL VALUES

   USE SINANOGLU

   IMPLICIT NONE
   INTEGER :: L1,L2,L3,M1,M2,M3
   INTEGER :: JMIN,JMAX,J
   INTEGER :: I1,I2,I3,I4,I5
   DOUBLE PRECISION :: A,B,C

   IF ((L3 > L1+L2).OR.(L3 < IABS(L1-L2)).OR.(IABS(M1) > L1).OR.(IABS(M2) > L2).OR.(IABS(M3) > L3).OR.(M1+M2+M3/=0)) THEN
!   CALL PABORT('WIGNER3J ARGUMENTS OUT OF BOUNDS')
!    CALL WARNING('WIGNER3J ARGUMENTS OUT OF BOUNDS')
    write(*,*) 'WIGNER3J ARGUMENTS OUT OF BOUNDS'
    WRITE(6,'(6I3)') L1,L2,L3,M1,M2,M3
    WIGNER3J=0.0D0
    RETURN
   ENDIF

   A=(FACTORIAL(L1+L2-L3)*FACTORIAL(L1-L2+L3)*FACTORIAL(-L1+L2+L3))/(FACTORIAL(L1+L2+L3+1))
   B=(FACTORIAL(L1+M1)*FACTORIAL(L1-M1)*FACTORIAL(L2+M2)*FACTORIAL(L2-M2)*FACTORIAL(L3+M3)*FACTORIAL(L3-M3))
   I1=L2-L3-M1
   I2=L1-L3+M2
   I3=L1-M1
   I4=L2+M2
   I5=L1+L2-L3
   JMIN=MAX(0,MAX(I1,I2))
   JMAX=MIN(I3,MIN(I4,I5))
   C=0.0D0
   DO J=JMIN,JMAX
    C=C+SGN(J)/(FACTORIAL(J)*FACTORIAL(J-I1)*FACTORIAL(J-I2)*FACTORIAL(I3-J)*FACTORIAL(I4-J)*FACTORIAL(I5-J))
   ENDDO
   WIGNER3J=(SGN(L1-L2-M3))*DSQRT(A)*DSQRT(B)*C
   RETURN
END FUNCTION

SUBROUTINE ratint(xa,ya,n,x,y,dy)
    ! Does not use any global vars
IMPLICIT NONE
INTEGER n,NMAX
DOUBLE PRECISION dy,x,y,xa(n),ya(n),TINY
PARAMETER (NMAX=10,TINY=1.d-25)
INTEGER i,m,ns
DOUBLE PRECISION dd,h,hh,t,w,c(NMAX),d(NMAX)
ns=1
hh=dabs(x-xa(1))
do i=1,n
  h=dabs(x-xa(i))
  if (h.eq.0.0d0)then
    y=ya(i)
    dy=0.0d0
    return
  else if (h.lt.hh) then
    ns=i
    hh=h
  endif
  c(i)=ya(i)
  d(i)=ya(i)+TINY
enddo
y=ya(ns)
ns=ns-1
do m=1,n-1
  do i=1,n-m
    w=c(i+1)-d(i)
    h=xa(i+m)-x
    t=(xa(i)-x)*d(i)/h
    dd=t-c(i+1)
    if(dd.eq.0.0d0) pause 'failure in ratint'
    dd=w/dd
    d(i)=c(i+1)*dd
    c(i)=t*dd
  enddo
  if (2*ns.lt.n-m)then
    dy=c(ns+1)
  else
    dy=d(ns)
    ns=ns-1
  endif
  y=y+dy
enddo
return
END

END MODULE POISSON
