MODULE CONSTANTS
! DEFINE PHYSICAL CONSTANTS AND ATOMIC SYMBOLS.
   IMPLICIT NONE
   DOUBLE PRECISION,PARAMETER :: PI = 3.14159265358979D0
   DOUBLE PRECISION,PARAMETER :: RTTWO = 1.41421356237309504880D0
   DOUBLE PRECISION,PARAMETER :: ANG_BOHR = 0.529177249D0
   DOUBLE PRECISION,PARAMETER :: RIEMANN_ZETA3 = 1.20205690315032D0
   DOUBLE PRECISION,PARAMETER :: BOLTZMANN = 1.0D0/3.1577464D5
   DOUBLE PRECISION,PARAMETER :: EV = 27.2113957D0
   DOUBLE PRECISION,PARAMETER :: DEBYED = -2.5418D0  ! UNITS FOR DIPOLE MOMENTS
   DOUBLE PRECISION,PARAMETER :: DEBYEQ = -1.34504D0 ! UNITS FOR QUADRUPOLE MOMENTS
   CHARACTER(2),PARAMETER :: CATOM(20) = &
      (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al',&
        'Si','P ','S ','Cl','Ar','K ','Ca'/)
   CHARACTER(5),PARAMETER :: CNUMBER(10) = &
      (/'  ONE','  TWO','THREE',' FOUR',' FIVE','  SIX','SEVEN','EIGHT',' NINE','  TEN'/)
END MODULE






MODULE STRUCTURE
! INPUT DATA FOR MOLECULAR STRUCTURES
   IMPLICIT NONE
   INTEGER :: NATOM                                     ! NUMBER OF NUCLEAI
   INTEGER :: CHARGE                                    ! CHARGE (NONZERO VALUE IS VALID ONLY FOR ATOMIC & MOLECULAR CALCULATIONS)
   INTEGER :: ICORE,IVIRTCORE,IOCC                      ! NUMBER OF CORE ORBITALS & OCCUPIED ORBITALS (INCLUDING CORE)
   INTEGER,ALLOCATABLE :: IALL(:)                       ! NUMBER OF ALL ORBITALS (NORMALLY EQUALS TO NCGS)
   INTEGER :: IALLMAX                                   ! MAXIMUM NUMBER OF ALL ORBITALS (NORMALLY EQUALS TO NCGS)
   INTEGER :: IATOM(256)                                ! ATOMIC NUMBER
   DOUBLE PRECISION :: PERIOD,HELIX                     ! TRANSLATIONAL PERIOD & HELICAL ANGLE
   DOUBLE PRECISION :: ATOMX(256),ATOMY(256),ATOMZ(256) ! CARTESIAN COORDINATES OF NUCLEAI
   DOUBLE PRECISION :: FERMI                            ! FERMI ENERGY
END MODULE





MODULE DFT
! NUMERICAL QUADRUTURE USED IN DFT CALCULATIONS
   IMPLICIT NONE
   INTEGER :: MAXNGRID,MGRID
   INTEGER,ALLOCATABLE :: NGRID(:)
   DOUBLE PRECISION,ALLOCATABLE :: GRIDX(:,:),GRIDY(:,:),GRIDZ(:,:),GRIDW(:,:)    ! CARTESIAN COORDINATES & BECKE WEIGHTS FOR GRID POINTS
!  DOUBLE PRECISION,ALLOCATABLE :: GRIDXA(:),GRIDYA(:),GRIDZA(:),GRIDWA(:)        ! CARTESIAN COORDINATES & BECKE WEIGHTS FOR GRID POINTS
   DOUBLE PRECISION,ALLOCATABLE :: FUZZY(:,:),ATOMW(:,:)                          ! A FUZZY CELL WEIGHT FUNCTION
   DOUBLE PRECISION,ALLOCATABLE :: RHO(:,:),RHO_GX(:,:),RHO_GY(:,:),RHO_GZ(:,:),GAMMA(:,:) ! ELECRTRON DENSITY, ITS GRADIENTS, & GRADIENT INVARIANTS
   DOUBLE PRECISION,ALLOCATABLE :: D1F1(:,:),D1F2(:,:),XCF(:,:)                   ! FUNCTIONAL VALUES & FUNCTIONAL FIRST DERIVATIVES
   DOUBLE PRECISION,ALLOCATABLE :: D2F1(:,:),D2F2(:,:),D2F3(:,:)                  ! FUNCTIONAL SECOND DERIVATIVES
   DOUBLE COMPLEX,ALLOCATABLE :: TRHO(:,:),TRHO_GX(:,:),TRHO_GY(:,:),TRHO_GZ(:,:) ! TRIAL ELECTRON DENSITY & ITS GRADIENTS
END MODULE





MODULE SINANOGLU
! GRID-BASED ATOMIC SINANOGLU CALCULATIONS
   IMPLICIT NONE
   INTEGER :: SALG ! 1: TRAPEZOIDAL; 2: LOGARITHM; 3: BECKE
   INTEGER :: RG, AG                                       ! NUMBER OF RADIAL / ANGULAR GRID POINTS
   INTEGER :: LMAX                                         ! PRINCIPAL ANGULAR MOMENTUM
   INTEGER :: NAO                                          ! NUMBER OF SMALL AO
   INTEGER :: NAO2                                         ! NUMBER OF LARGE AO
   DOUBLE PRECISION,ALLOCATABLE :: FACTORIAL(:)            ! FACTORIAL
   DOUBLE PRECISION,ALLOCATABLE :: SGN(:)                  ! (-1)**I
   DOUBLE PRECISION,ALLOCATABLE :: SG_GAUSS_CHEV(:,:)      ! RADIAL GRID POINTS
   DOUBLE PRECISION,ALLOCATABLE :: SG_GAUSS_CHEV_W(:,:)    ! RADIAL GRID WEIGHTS
   DOUBLE PRECISION,ALLOCATABLE :: SG_LEBEDVX(:)           ! LEBEDEV GRID POINTS
   DOUBLE PRECISION,ALLOCATABLE :: SG_LEBEDVY(:)           ! LEBEDEV GRID POINTS
   DOUBLE PRECISION,ALLOCATABLE :: SG_LEBEDVZ(:)           ! LEBEDEV GRID POINTS
   DOUBLE PRECISION,ALLOCATABLE :: SG_LEBEDVW(:)           ! LEBEDEV GRID WEIGHTS
!  DOUBLE PRECISION,ALLOCATABLE :: R12(:,:,:,:)            ! DISTANCE MATRIX
!  DOUBLE PRECISION,ALLOCATABLE :: RECIPROCALR12(:,:,:,:)  ! RECIPROCAL DISTANCE MATRIX
   DOUBLE PRECISION,ALLOCATABLE :: ORBITALS(:,:,:)         ! HF ORBITAL AMPLITUDES AT GRID POINTS 
   DOUBLE PRECISION,ALLOCATABLE :: ORBITALS_SAVE(:,:,:)    ! HF ORBITAL AMPLITUDES AT GRID POINTS
   DOUBLE PRECISION,ALLOCATABLE :: DXORBITALS(:,:,:)       ! OCCUPIED ORBITAL X-DERIVATIVES AT GRID POINTS
   DOUBLE PRECISION,ALLOCATABLE :: DYORBITALS(:,:,:)       ! OCCUPIED ORBITAL Y-DERIVATIVES AT GRID POINTS
   DOUBLE PRECISION,ALLOCATABLE :: DZORBITALS(:,:,:)       ! OCCUPIED ORBITAL Z-DERIVATIVES AT GRID POINTS
   DOUBLE COMPLEX,ALLOCATABLE   :: RORB(:,:,:,:,:)         ! RADIAL FUNCTIONS FOR ORBITALS
   DOUBLE COMPLEX,ALLOCATABLE   :: RXORB(:,:,:,:,:)        ! D/DX RADIAL FUNCTIONS FOR ORBITALS
   DOUBLE COMPLEX,ALLOCATABLE   :: RYORB(:,:,:,:,:)        ! D/DY RADIAL FUNCTIONS FOR ORBITALS
   DOUBLE COMPLEX,ALLOCATABLE   :: RZORB(:,:,:,:,:)        ! D/DZ RADIAL FUNCTIONS FOR ORBITALS
   DOUBLE COMPLEX,ALLOCATABLE   :: RN(:,:,:,:)             ! RADIAL FUNCTIONS FOR NUCLEAR ATTRACTION
   DOUBLE COMPLEX,ALLOCATABLE   :: RJ(:,:,:,:)             ! RADIAL FUNCTIONS FOR COULOMB POTENTIAL
   DOUBLE COMPLEX,ALLOCATABLE   :: RRHO(:,:,:,:)           ! RADIAL FUNCTIONS FOR ELECTRON DENSITY
   DOUBLE COMPLEX,ALLOCATABLE   :: AORB(:,:,:,:,:)         ! RADIAL FUNCTIONS FOR AO (SMALL)
   DOUBLE COMPLEX,ALLOCATABLE   :: AORB2(:,:,:,:,:)        ! RADIAL FUNCTIONS FOR AO (LARGE)
!  DOUBLE COMPLEX,ALLOCATABLE   :: AXORB(:,:,:,:,:)        ! D/DX RADIAL FUNCTIONS FOR AO
!  DOUBLE COMPLEX,ALLOCATABLE   :: AYORB(:,:,:,:,:)        ! D/DY RADIAL FUNCTIONS FOR AO
!  DOUBLE COMPLEX,ALLOCATABLE   :: AZORB(:,:,:,:,:)        ! D/DZ RADIAL FUNCTIONS FOR AO
   DOUBLE PRECISION,ALLOCATABLE :: AOBASIS(:,:,:)          ! AO (SMALL)
   DOUBLE PRECISION,ALLOCATABLE :: AOBASIS2(:,:,:)         ! AO (LARGE)
   DOUBLE PRECISION,ALLOCATABLE :: POTENTIALJ(:,:)         ! COULOMB POTENTIAL 
   DOUBLE PRECISION,ALLOCATABLE :: PREVIOUS_J(:,:)         ! COULOMB POTENTIAL 
   DOUBLE PRECISION,ALLOCATABLE :: POTENTIALN(:,:)         ! NUCLEAR POTENTIAL 
   DOUBLE PRECISION,ALLOCATABLE :: SG_EPSILON(:)           ! ORBITAL ENERGIES
END MODULE