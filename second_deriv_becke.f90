MODULE DERIVATIVE_AUX

    IMPLICIT NONE

INTEGER :: RG
REAL(8),PARAMETER :: PI = 3.14159265358979_8
REAL(8),PARAMETER :: BSRADI(17) = &
  (/0.472_8,0.472_8,1.370_8,0.992_8,0.803_8,0.661_8,0.614_8,0.567_8,0.472_8, &
    0.472_8,1.701_8,1.417_8,1.181_8,1.039_8,0.945_8,0.945_8,0.945_8/) ! Ne DATA WAS MISSING IN THE ORIGINAL
REAL(8),ALLOCATABLE :: SG_GAUSS_CHEV(:)      ! RADIAL GRID POINTS
REAL(8),ALLOCATABLE :: SG_GAUSS_CHEV_W(:)    ! RADIAL GRID WEIGHTS

END MODULE DERIVATIVE_AUX


MODULE DERIVATIVE

IMPLICIT NONE

CONTAINS

SUBROUTINE SECOND_DERIV(R, IA, NRAD ) bind (C, name = 'second_deriv_')
! 7-point finite-difference second differentiation
! (1/r)(d2/dr2)rf(r) = (2/r)(df/dr)+(d2f/dr2) = (2/r)(dx/dr)(df/dx) + (d2x/dr2)(df/dx) + (dx/dr)^2(d2f/dx2)
   
   use iso_c_binding
   USE DERIVATIVE_AUX

   IMPLICIT NONE
   INTEGER :: IA, NRAD
   DOUBLE PRECISION  :: R(1:NRAD)
   DOUBLE PRECISION :: H
   DOUBLE PRECISION, ALLOCATABLE :: R1D(:)
   DOUBLE PRECISION, ALLOCATABLE :: R2D(:)
   INTEGER :: I
   
   DOUBLE PRECISION :: R0,X,OMEGA

   H=PI/DFLOAT(RG+1)
   allocate(R1D(RG),R2D(RG))
   !write(*,*) 'RG inside SECOND_DERIV = ', RG 

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
   R1D(RG-2)=(12.0D0*R(RG-6)-96.0D0*R(RG-5)+360.0D0*R(RG-4)-960.0D0*R(RG-3) &
             +420.0D0*R(RG-2)+288.0D0*R(RG-1)-24.0D0*R(RG))/(720.0D0*H)
   R1D(RG-1)=(-24.0D0*R(RG-6)+180.0D0*R(RG-5)-600.0D0*R(RG-4)+1200.0D0*R(RG-3) &
             -1800.0D0*R(RG-2)+924.0D0*R(RG-1)+120.0D0*R(RG))/(720.0D0*H)
   R1D(RG)=(120.0D0*R(RG-6)-864.0D0*R(RG-5)+2700.0D0*R(RG-4)-4800.0D0*R(RG-3) &
            +5400.0D0*R(RG-2)-4320.0D0*R(RG-1)+1764.0D0*R(RG))/(720.0D0*H)
   DO I=5,RG-4
    R1D(I)=(+144.0D0*R(I-4)-1536.0D0*R(I-3)+8064.0D0*R(I-2)-32256.0D0*R(I-1) &
           +32256.0D0*R(I+1)-8064.0D0*R(I+2)+1536.0D0*R(I+3)-144.0D0*R(I+4))/(40320.0D0*H)
   ENDDO

   
   !write(*, *) 'R1D :'
   !write(*, *) 'R1D(1) = ', R1D(1)
   !write(*, *) 'R1D(2) = ', R1D(2)
   !write(*, '(10F13.6 )') (R1D(I), I = 1, RG)

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
   R2D(RG-2)=(4.0D0*R(RG-6)-24.0D0*R(RG-5)+30.0D0*R(RG-4)+400.0D0*R(RG-3) &
             -840.0D0*R(RG-2)+456.0D0*R(RG-1)-26.0D0*R(RG))/(360.0D0*H**2)
   R2D(RG-1)=(-26.0D0*R(RG-6)+186.0D0*R(RG-5)-570.0D0*R(RG-4)+940.0D0*R(RG-3) &
             -510.0D0*R(RG-2)-294.0D0*R(RG-1)+274.0D0*R(RG))/(360.0D0*H**2)
   R2D(RG)=(274.0D0*R(RG-6)-1944.0D0*R(RG-5)+5940.0D0*R(RG-4)-10160.0D0*R(RG-3) &
           +10530.0D0*R(RG-2)-6264.0D0*R(RG-1)+1624.0D0*R(RG))/(360.0D0*H**2)
   DO I=5,RG-4
    R2D(I)=(-36.0D0*R(I-4)+512.0D0*R(I-3)-4032.0D0*R(I-2)+32256.0D0*R(I-1)-57400.0D0*R(I) &
         +32256.0D0*R(I+1)-4032.0D0*R(I+2)+512.0D0*R(I+3)-36.0D0*R(I+4))/(20160.0D0*H**2)
   ENDDO

   !write(*, *) 'R2D :'
   !write(*, *) 'R2D(1) = ', R2D(1)
   !write(*, *) 'R2D(2) = ', R2D(2)
   !write(*, '(10F13.6 )') (R2D(I), I = 1, RG)

   !print*, 'Will use the grid below'
   R0=BSRADI(IA)
   DO I=1,RG
     OMEGA=DFLOAT(I)*PI/DFLOAT(RG+1)
     X=DCOS(OMEGA)
     R(I)=-2.0D0*R1D(I)/SG_GAUSS_CHEV(I)*(1.0D0-X)**2/(2.0D0*R0*DSIN(OMEGA)) &
          +R1D(I)*((1.0D0-X)**3/(2.0D0*R0**2*DSIN(OMEGA)) &
                  -(1.0D0-X)**4*DCOS(OMEGA)/(4.0D0*R0**2*DSIN(OMEGA)**3)) &
          +R2D(I)*((1.0D0-X)**2/(2.0D0*R0*DSIN(OMEGA)))**2
   ENDDO

   DEALLOCATE(R2D,R1D)


   RETURN
END SUBROUTINE SECOND_DERIV



SUBROUTINE SECOND_DERIV2(R,IA, NRAD) bind (C, name = 'second_deriv2_')
! 7-point finite-difference second differentiation
! (1/r)(d2/dr2)r r^-1f(r) = (1/r)(d2x/dr2)(df/dx) + (1/r)(dx/dr)^2(d2f/dx2)
   
   use iso_c_binding
   USE DERIVATIVE_AUX


   IMPLICIT NONE
   INTEGER :: IA, NRAD
   DOUBLE PRECISION   :: R(1:NRAD)
   DOUBLE PRECISION :: H
   DOUBLE PRECISION,ALLOCATABLE :: R1D(:)
   DOUBLE PRECISION,ALLOCATABLE :: R2D(:)
   INTEGER :: I
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

    !WRITE(*,*) R1D(1:RG)
    !WRITE(*,*) R2D(1:RG)

   R0=BSRADI(IA)
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

SUBROUTINE SECOND_DERIV2_DEBUG(R,IA, NRAD, R1D, R2D) bind (C, name = 'second_deriv2_debug_')
! 7-point finite-difference second differentiation
! (1/r)(d2/dr2)r r^-1f(r) = (1/r)(d2x/dr2)(df/dx) + (1/r)(dx/dr)^2(d2f/dx2)
   
   use iso_c_binding
   USE DERIVATIVE_AUX


   IMPLICIT NONE
   INTEGER :: IA, NRAD
   ! AUTOMATIC ARRAYS
   real(8)   :: R(1:NRAD), R1D(1:NRAD), R2D(1:NRAD)
   real(8) :: H
   INTEGER :: I
   real(8) :: R0,X,OMEGA

   H=PI/DFLOAT(RG+1)

   ! FIRST DERIVATIVES
    R1D(1)=(-1764.0_8 * R(1)+4320.0_8 *R(2)-5400.0_8*R(3)+4800.0_8*R(4) &
            -2700.0_8*R(5)+864.0_8*R(6)-120.0_8*R(7))/(720.0_8*H)
    R1D(2)=(-120.0_8*R(1)-924.0_8*R(2)+1800.0_8*R(3)-1200.0_8*R(4) &
            +600.0_8*R(5)-180.0_8*R(6)+24.0_8*R(7))/(720.0_8*H)
    R1D(3)=(24.0_8*R(1)-288.0_8*R(2)-420.0_8*R(3)+960.0_8*R(4) &
            -360.0_8*R(5)+96.0_8*R(6)-12.0_8*R(7))/(720.0_8*H)
    R1D(4)=(-12.0_8*R(1)+108.0_8*R(2)-540.0_8*R(3) &
           +540.0_8*R(5)-108.0_8*R(6)+12.0_8*R(7))/(720.0_8*H)
    R1D(RG-3)=(-12.0_8*R(RG-6)+108.0_8*R(RG-5)-540.0_8*R(RG-4) &
              +540.0_8*R(RG-2)-108.0_8*R(RG-1)+12.0_8*R(RG))/(720.0_8*H)
    R1D(RG-2)=(-12.0_8*R(RG-5)+108.0_8*R(RG-4)-540.0_8*R(RG-3) &
              +540.0_8*R(RG-1)-108.0_8*R(RG))/(720.0_8*H)
    R1D(RG-1)=(+12.0_8*R(RG-5)-96.0_8*R(RG-4)+360.0_8*R(RG-3)-960.0_8*R(RG-2) &
              +420.0_8*R(RG-1)+288.0_8*R(RG))/(720.0_8*H)
    R1D(RG)=(-24.0_8*R(RG-5)+180.0_8*R(RG-4)-600.0_8*R(RG-3)+1200.0_8*R(RG-2) &
             -1800.0_8*R(RG-1)+924.0_8*R(RG))/(720.0_8*H)
    DO I=5,RG-4
     R1D(I)=(+144.0_8*R(I-4)-1536.0_8*R(I-3)+8064.0_8*R(I-2)-32256.0_8*R(I-1) &
            +32256.0_8*R(I+1)-8064.0_8*R(I+2)+1536.0_8*R(I+3)-144.0_8*R(I+4))/(40320.0_8*H)
    ENDDO

   ! SECOND DERIVATIVES
    R2D(1)=(1624.0_8*R(1)-6264.0_8*R(2)+10530.0_8*R(3)-10160.0_8*R(4) &
            +5940.0_8*R(5)-1944.0_8*R(6)+274.0_8*R(7))/(360.0_8*H**2)
    R2D(2)=(274.0_8*R(1)-294.0_8*R(2)-510.0_8*R(3)+940.0_8*R(4) &
            -570.0_8*R(5)+186.0_8*R(6)-26.0_8*R(7))/(360.0_8*H**2)
    R2D(3)=(-26.0_8*R(1)+456.0_8*R(2)-840.0_8*R(3)+400.0_8*R(4) &
            +30.0_8*R(5)-24.0_8*R(6)+4.0_8*R(7))/(360.0_8*H**2)
    R2D(4)=(4.0_8*R(1)-54.0_8*R(2)+540.0_8*R(3)-980.0_8*R(4) &
            +540.0_8*R(5)-54.0_8*R(6)+4.0_8*R(7))/(360.0_8*H**2)
    R2D(RG-3)=(4.0_8*R(RG-6)-54.0_8*R(RG-5)+540.0_8*R(RG-4)-980.0_8*R(RG-3) &
              +540.0_8*R(RG-2)-54.0_8*R(RG-1)+4.0_8*R(RG))/(360.0_8*H**2)
    R2D(RG-2)=(+4.0_8*R(RG-5)-54.0_8*R(RG-4)  +540.0_8*R(RG-3)-980.0_8*R(RG-2) &
              +540.0_8*R(RG-1)-54.0_8*R(RG))/(360.0_8*H**2)
    R2D(RG-1)=(+4.0_8*R(RG-5)-24.0_8*R(RG-4)   +30.0_8*R(RG-3)+400.0_8*R(RG-2) &
              -840.0_8*R(RG-1)+456.0_8*R(RG))/(360.0_8*H**2)
    R2D(RG)=(-26.0_8*R(RG-5)+186.0_8*R(RG-4)  -570.0_8*R(RG-3)+940.0_8*R(RG-2) &
             -510.0_8*R(RG-1)-294.0_8*R(RG))/(360.0_8*H**2)
    DO I=5,RG-4
     R2D(I)=(-36.0_8*R(I-4) +512.0_8*R(I-3)-4032.0_8*R(I-2)+32256.0_8*R(I-1)-57400.0_8*R(I) &
             +32256.0_8*R(I+1)-4032.0_8*R(I+2)+512.0_8*R(I+3)-36.0_8*R(I+4))/(20160.0_8*H**2)
    ENDDO


   R0=BSRADI(IA)
   DO I=1,RG
     OMEGA=DFLOAT(I)*PI/DFLOAT(RG+1)
     X=DCOS(OMEGA)
     !R(I)=+R1D(I)*((1.0_8-X)**3/(2.0_8*R0**2*DSIN(OMEGA)) &
     !     -(1.0_8-X)**4*DCOS(OMEGA)/(4.0_8*R0**2*DSIN(OMEGA)**3)) &
     !     +R2D(I)*((1.0_8-X)**2/(2.0_8*R0*DSIN(OMEGA)))**2
     R(I)=+R1D(I)*( (1.0_8-X)*(1.0_8-X)*(1.0_8-X) /(2.0_8*R0*R0*DSIN(OMEGA)) &
          -(1.0_8-X)*(1.0_8-X)*(1.0_8-X)*(1.0_8-X)*DCOS(OMEGA)/(4.0_8*R0*R0*DSIN(OMEGA)*DSIN(OMEGA)*DSIN(OMEGA))) &
          +((1.0_8-X)*(1.0_8-X)*(1.0_8-X)*(1.0_8-X)/(4.0_8*R0*R0*DSIN(OMEGA)*DSIN(OMEGA))) * R2D(I)
   ENDDO

   RETURN
END SUBROUTINE

SUBROUTINE CONSTRUCT_RGRID(NRAD, IA) bind (C, name = 'construct_rgrid_')

    use iso_c_binding
    USE DERIVATIVE_AUX

    implicit none

    real(8) :: R1, X
    integer :: ia, j, nrad


    R1 = BSRADI(ia) 
    RG = NRAD


    ! Construct the radial grid

    !print*, 'Will construct the grid in second_derivative'
    ALLOCATE(SG_GAUSS_CHEV(RG),SG_GAUSS_CHEV_W(RG))

    DO J=1,RG
        X=DCOS(DFLOAT(J)*PI/DFLOAT(RG+1))
        SG_GAUSS_CHEV(J)=R1*(1.0_8+X)/(1.0_8-X)
        SG_GAUSS_CHEV_W(J)=2.0_8*R1/(1.0_8-X)**2*PI/DFLOAT(RG+1)*DSIN(DFLOAT(J)*PI/DFLOAT(RG+1))*SG_GAUSS_CHEV(J)**2
    ENDDO

    !print*, 'finished constructing the grid in second_derivative'


END SUBROUTINE


SUBROUTINE DESTROY_RGRID() bind (C, name = 'destroy_rgrid_')

    use iso_c_binding
    USE DERIVATIVE_AUX

    implicit none

    !print*, 'Will destroy becke grid below'

    DEALLOCATE(SG_GAUSS_CHEV ,SG_GAUSS_CHEV_W)


END SUBROUTINE

END MODULE DERIVATIVE


