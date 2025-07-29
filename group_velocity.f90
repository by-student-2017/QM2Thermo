!-----------------------------------------------------------------------
! Author   : H. Sato and M. Inukai
! Affiliation : [AUE / Riken]
! Contact  : [your.email@domain.edu]
! GitHub   : https://github.com/yourusername/seebeck_analysis (optional)
!-----------------------------------------------------------------------
! Program : group_velocity_fcc.f90
! Purpose : 
!   - To compute group velocity components (VX, VY, VZ) for FCC structures.
!   - To use the band structure data and tetrahedron method for DOS (Density of States) analysis.
!   - Applicable for calculating transport properties like Seebeck coefficient.
!
! Dates   : Created 2025.3.16, 5.9
!         : Revised 2025.07.15
!
! Physical constants:
!   AL = a = 5.431 [A],  Nk = 58^3
!   e   = 1.602176565 * 10^-19  [C]   (charge conversion: 1 [eV] -> [Coulomb])
!   h   = 6.62607015e-34 [J]   (1 Ha = 4.3597447222e-18 [J])
!   CEM = (e/h) * 10^10 = 0.241798935
!   CEM = 2.41798935D4 [cm/s per Hatree/Angstrom]
!
! Compilation:
!   gfortran -O2 group_velocity_fcc.f90 -o group_velocity_fcc.exe
!   ifort    -O2 group_velocity_fcc.f90 -o group_velocity_fcc.exe
!
! Required input files:
!   - wien.energy     : WIEN2k output data
!   - cf58*.dat       : differential points data
!
! Output:
!   - AKK.DATA        : Band-resolved energy, DOS, velocity^2
!-----------------------------------------------------------------------
!   src/
!   „¥„Ÿ„Ÿ modules/
!   „    „¥„Ÿ„Ÿ constants              ! Physical constants and size definitions
!   „    „¥„Ÿ„Ÿ band_data              ! AMA1 arrays and shared variable definitions (formerly COMMON)
!   „    „¥„Ÿ„Ÿ tetra_data             ! E(), FK(), CNT() etc. for DOS density
!   „    „¤„Ÿ„Ÿ index_data             ! For group velocity difference such as IP(), IM()
!   „¥„Ÿ„Ÿ subroutines/
!   „    „¥„Ÿ„Ÿ read_constants         ! 
!   „    „¥„Ÿ„Ÿ TGEN                   ! Define the tetrahedrons associated with each K point
!   „    „¥„Ÿ„Ÿ VDATA                  ! Load differential points (cf*.dat)
!   „    „¥„Ÿ„Ÿ AEIGEN                 ! Load AMA1 from wien.energy (WIEN2k)
!   „    „¥„Ÿ„Ÿ DENS                   ! PDNS applied to all tetras
!   „    „¥„Ÿ„Ÿ PDNS                   ! Calculate density of states for a single tetra
!   „    „¥„Ÿ„Ÿ SORT                   ! 4-point sort subroutine (S, F, G, O, P)
!   „    „¤„Ÿ„Ÿ GV                     ! Calculate group velocities VX, VY, and VZ using the central difference method  
!   „¤„Ÿ„Ÿ main                        ! Main Program
!-----------------------------------------------------------------------



!=============================================================================
! Module    : constants
! Purpose   : Define physical constants and common array size parameters
!             used throughout the band analysis and transport routines.
!
! Description :
!   - Contains mathematical and physical constants (pi, k_B, eV conversion)
!   - Grid-related parameters such as NPX, NPY, NPZ, NPT
!   - Safe upper bounds for array allocation
!=============================================================================
MODULE constants
  IMPLICIT NONE
  
  !-------------------------
  ! Physical constants
  !-------------------------
  REAL(KIND=8), PARAMETER :: PI   = 3.141592653589793D0   ! Pi
  REAL(KIND=8), PARAMETER :: EV   = 13.605693D0           ! Hartree to eV conversion
  REAL(KIND=8), PARAMETER :: KB   = 8.61733D-5            ! Boltzmann constant [eV/K]
  REAL(KIND=8), PARAMETER :: CEM  = 24179.8935D0          ! Coefficient e/h * 10^10
  REAL(KIND=8), PARAMETER :: B2A  = 0.52918               ! convert Bohr to Angstrom unit
  
  !-------------------------
  ! Lattice and grid size
  !-------------------------
  REAL(KIND=8) :: AL, BL, CL
  !REAL(KIND=8), PARAMETER :: AL   = 5.431                 ! Lattice constant [Angstrom] (example: Si)
  
  !-------------------------
  ! Shift reference in Rydberg (e.g., Fermi level)
  !-------------------------
  REAL(KIND=8) :: EF
  !REAL(KIND=8), PARAMETER :: EF   = 0.382737              ! Fermi energy [Ry]
  
  INTEGER :: NPX, NPY, NPZ
  INTEGER :: NPT
  !INTEGER, PARAMETER :: NPX   = 58                        ! k-mesh size in X (default: 58)
  !INTEGER, PARAMETER :: NPY   = 58                        ! same for Y
  !INTEGER, PARAMETER :: NPZ   = 58                        ! same for Z
  !INTEGER, PARAMETER :: NPT   = (NPX+1)*(NPY+1)*(NPZ+1)   ! Total number of k-points: Maximum number of rows in cf58A3.dat (4th column) : IEX
  
  !-------------------------
  ! Small numerical threshold
  !-------------------------
  REAL(KIND=8), PARAMETER :: EPS = 1.0D-7                 ! Small epsilon for comparisons
  
END MODULE constants



!=======================================================================
! Module    : band_data
! Purpose   : Provide global access to band energy eigenvalue arrays
!             used in AEIGEN, GV, and density calculation routines.
!
! Contents  :
!   - AMA1(IMAX, MML) : Raw band eigenvalue data from wien.energy
!   - MMIN            : Number of bands actually used after JCUT truncation
!   - IMAX            : Total number of k-points
!
! Usage     :
!   - Accessed via USE band_data in routines such as AEIGEN, GV, VDATA
!   - Allows centralized control of band structure input and truncation
!=======================================================================
MODULE band_data
  IMPLICIT NONE
  
  !-------------------------
  ! Array dimensions
  !-------------------------
  INTEGER :: IMAX                              ! Maximum number of k-points
  !INTEGER, PARAMETER :: IMAX = 4735            ! Maximum number of k-points
  INTEGER, PARAMETER :: MML  = 50              ! Maximum number of bands
  
  !-------------------------
  ! Band energy data
  !-------------------------
  REAL(KIND=8), ALLOCATABLE :: AMA1(:,:)
  !REAL(KIND=8), DIMENSION(IMAX, MML) :: AMA1   ! Band eigenvalues [eV]
  INTEGER :: MMIN                              ! Number of bands used (JMAX - JCUT)
  
END MODULE band_data



!===============================================================================
! Module    : tetra_data
! Purpose   : Provide shared physical arrays and variables used in
!             tetrahedron-based density of states and carrier calculations.
!
! Contents:
!   - Energy-dependent arrays:
!       E     : interpolated energy values per k-point
!       EM    : energy grid for DOS integration
!       FK,FKB,FKC,FKD : group velocity-related values
!       CNT,SNT : density/occupation state counters
!       FNT,FNTB,FNTC,FNTD : transport-weighted DOS integrals
!
!   - Scalar parameters:
!       DE    : energy grid spacing
!       V     : tetrahedron volume weight
!       BOA,COA : optional physical constants
!       NE    : number of energy grid points
!
!   - Shared sizes (used across loops):
!       NPT   : total number of k-points
!       NEMAX : number of EM points
!       NL,NT : auxiliary counters
!
! Usage:
!   - Access via USE tetra_data in PDNS, DENS, BS, etc.
!===============================================================================
MODULE tetra_data
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  IMPLICIT NONE
  
  !-------------------------
  ! Energy and velocity arrays
  !-------------------------
  REAL(KIND=8), ALLOCATABLE :: E(:), FK(:), FKB(:), FKC(:), FKD(:), EM(:)
  !REAL(KIND=8), DIMENSION(NPT) :: E    ! Energy per k-point (interpolated energy values per k-point)
  !REAL(KIND=8), DIMENSION(NPT) :: FK   ! velocity^2 or transport function (group velocity-related values)
  !REAL(KIND=8), DIMENSION(NPT) :: FKB  ! auxiliary velocity channel B (group velocity-related values)
  !REAL(KIND=8), DIMENSION(NPT) :: FKC  ! auxiliary velocity channel C (group velocity-related values)
  !REAL(KIND=8), DIMENSION(NPT) :: FKD  ! auxiliary velocity channel D (group velocity-related values)
  !REAL(KIND=8), DIMENSION(NPT) :: EM   ! energy mesh grid [eV] (energy grid for DOS integration)
  
  !-------------------------
  ! Density counters and integrals
  !-------------------------
  REAL(KIND=8), ALLOCATABLE :: CNT(:), SNT(:), FNT(:), FNTB(:), FNTC(:), FNTD(:)
  !REAL(KIND=8), DIMENSION(NPT) :: CNT  ! total DOS (density/occupation state counters)
  !REAL(KIND=8), DIMENSION(NPT) :: SNT  ! accumulated state number (density/occupation state counters)
  !REAL(KIND=8), DIMENSION(NPT) :: FNT  ! weighted DOS with FK  (transport-weighted DOS integrals)
  !REAL(KIND=8), DIMENSION(NPT) :: FNTB ! weighted DOS with FKB (transport-weighted DOS integrals)
  !REAL(KIND=8), DIMENSION(NPT) :: FNTC ! weighted DOS with FKC (transport-weighted DOS integrals)
  !REAL(KIND=8), DIMENSION(NPT) :: FNTD ! weighted DOS with FKD (transport-weighted DOS integrals)
  
  !-------------------------
  ! Scalar parameters
  !-------------------------
  REAL(KIND=8) :: DE   ! Energy grid spacing
  REAL(KIND=8) :: V    ! Tetrahedron volume weight
  REAL(KIND=8) :: BOA  ! Optional physical constants (Optional coefficient A)
  REAL(KIND=8) :: COA  ! Optional physical constants (Optional coefficient B)
  
  !-------------------------
  ! Grid counters
  !-------------------------
  INTEGER, PARAMETER :: NEMAX = 60001 ! Number of energy grid points
  INTEGER :: NL        ! Auxiliary index (e.g. band loop)
  INTEGER :: NT        ! Auxiliary index (e.g. temperature points)
  
END MODULE tetra_data



!===============================================================================
! Module    : index_data
! Purpose   : Provide k-point neighbor indices and mapping arrays used in
!             tetrahedron-based transport and group velocity (GV) calculations.
!
! Contents:
!   - IC(0:NPX+1,0:NPY+1,0:NPZ+1) : 3D k-mesh -> linear index mapping table
!   - N1-N7(NPT)            : Tetrahedron neighbor vertex indices for each k-point
!
!   - IP(3,15,NPT)          : Positive-side neighbor index for 15-point stencil (CK1-CK15) (x,y,z)
!   - IM(3,15,NPT)          : Negative-side neighbor index for 15-point stencil (CK1-CK15) (x,y,z)
!
!   - IEX(NPT)              : Remapping index (used to reorder AMA1 values for symmetry)
!
! Usage:
!   - Used in:
!       > TGEN              : Generates IC, N1-N7 for tetrahedral mesh
!       > DENS, PDNS        : Uses N1-N7 for DOS integration
!       > GV                : Uses IP/IM for group velocity evaluation
!       > VDATA             : Reads IP/IM and IEX from external stencil files
!
! Notes:
!   - All arrays indexed by NPT, the total number of k-points ((NPX+1)*(NPY+1)*(NPZ+1))
!   - Designed for modular access and avoids COMMON blocks
!===============================================================================
MODULE index_data
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  IMPLICIT NONE
  
  !-------------------------
  ! Tetrahedron mesh mapping
  !-------------------------
  INTEGER, ALLOCATABLE :: IC(:,:,:)
  !INTEGER, DIMENSION(0:NPX+1, 0:NPY+1, 0:NPZ+1) :: IC
  INTEGER, ALLOCATABLE :: N1(:), N2(:), N3(:), N4(:), N5(:), N6(:), N7(:)
  !INTEGER, DIMENSION(NPT) :: N1, N2, N3, N4, N5, N6, N7
  
  !-------------------------
  ! Neighbor index arrays
  !-------------------------
  INTEGER, ALLOCATABLE :: IP(:,:,:), IM(:,:,:)
  !INTEGER, DIMENSION(3,15,NPT) :: IP    ! Forward (positive) k-point stencil indices
  !INTEGER, DIMENSION(3,15,NPT) :: IM    ! Backward (negative) k-point stencil indices
  
  !-------------------------
  ! Symmetry mapping array
  !-------------------------
  INTEGER, ALLOCATABLE :: IEX(:)
  !INTEGER, DIMENSION(NPT)     :: IEX    ! Index mapping for AMA1 reordering
  
END MODULE index_data



SUBROUTINE read_constants
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  USE band_data  ! Band energies: AMA1(IMAX,MML), total bands MMIN
  IMPLICIT NONE
  
  CHARACTER(LEN=100) :: line
  REAL(KIND=8) :: LA, LB, LC, alpha, beta, gamma
  INTEGER :: iostat

  OPEN(UNIT=90, FILE='wien.dos1', STATUS='OLD', IOSTAT=iostat)
  READ(90, '(A)')
  READ(90, '(5X, f9.5)') EF
  WRITE(6,'(A, f9.5)') "Fermi energy [Ry]: ", EF
  CLOSE(90)

  OPEN(UNIT=91, FILE='wien.struct', STATUS='OLD', IOSTAT=iostat)
  READ(91, '(A)')
  READ(91, '(A)')
  READ(91, '(A)')
  READ(91, *) LA, LB, LC, alpha, beta, gamma
  AL = LA * B2A   ! B2A = 0.52918
  BL = LB * B2A   ! B2A = 0.52918
  CL = LC * B2A   ! B2A = 0.52918
  WRITE(6,'(A, f9.5)') "Lattice constant, a [A]: ", AL
  WRITE(6,'(A, f9.5)') "Lattice constant, b [A]: ", BL
  WRITE(6,'(A, f9.5)') "Lattice constant, c [A]: ", CL
  CLOSE(91)

  OPEN(UNIT=92, FILE='wien.kgen', STATUS='OLD', IOSTAT=iostat)
  READ(92, '(1X, I9)') IMAX
  WRITE(6,'(A,I9)') "Maximum number of k-points:", IMAX
  CLOSE(92)

  OPEN(UNIT=93, FILE='wien.klist', STATUS='OLD', IOSTAT=iostat)
  READ(93, '(85X, 3(1X,I2))') NPX, NPY, NPZ
  WRITE(6,'(A,I2)') "k-mesh size X, NPX:", NPX
  WRITE(6,'(A,I2)') "k-mesh size Y, NPY:", NPY
  WRITE(6,'(A,I2)') "k-mesh size Z, NPZ:", NPZ
  NPT = (NPX+1)*(NPY+1)*(NPZ+1)
  WRITE(6,'(A,I9)') "Total number of k-points:", NPT
  CLOSE(93)

END SUBROUTINE



!===============================================================================
! Program   : MAIN
! Purpose   : Calculates energy-dependent transport properties from band structure data.
!             Includes computation of group velocities and density of states via tetrahedral integration.
!
! Output    : Results are stored in AKK.DATA for post-processing or analysis.
!===============================================================================
PROGRAM MAIN
  !=======================
  ! Module imports
  !=======================
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  USE band_data  ! Band energies: AMA1(IMAX,MML), total bands MMIN
  USE tetra_data ! Variables for DOS integration: LMO, E, FK, FKB, FKC, FKD, EM, CNT, SNT, FNT, FNTB, FNTC, FNTD, DE, V, BOA, COA, NEMAX, NL, NT
  USE index_data ! k-point mapping and stencil indices: IC, N1-N7, IP, IM, IEX
  
  IMPLICIT NONE
  
  !=======================
  ! Variable declarations
  !=======================
  
  INTEGER :: LAT
  ! Lattice type (e.g., 3 = BCC). Currently unused, placeholder for future symmetry control.
  
  ! Loop indices
  INTEGER :: I, J, IKP, LB, LE, I2
  ! I, J       : Loop indices to traverse k-points or bands.
  ! IKP        : Index for k-points used in band structure calculations.
  ! LB         : Band number under processing for energy or velocity calculations.
  ! LE         : Energy mesh index used in density of states (DOS) calculations.
  ! I2         : Re-mapped index for symmetry-adjusted eigenvalues from `IEX`.
  
  ! Energy and velocity matrices
  REAL(KIND=8), ALLOCATABLE :: AMA(:,:), AKK(:,:), AKKB(:,:), AKKC(:,:), AKKD(:,:)
  !REAL(KIND=8), DIMENSION(NPT, MML) :: AMA, AKK, AKKB, AKKC, AKKD
  ! AMA       : Reordered band energies (using IEX symmetry map)
  ! AKK       : VX^2 group velocity squared (x-direction)
  ! AKKB      : Redundant VX^2 (used for B-channel calculations)
  ! AKKC      : VZ^2 group velocity squared (z-direction)
  ! AKKD      : VX raw velocity (used in directional transport)
  
  ! Scalars for energy integration and DOS
  REAL(KIND=8) :: ED, EU, TN, EE, CA, CB, CCC, CFC, CD, CC, C2, C2B, C2C, C2D
  ! ED        : Lower energy bound for calculations.
  ! EU        : Upper energy bound for calculations.
  ! TN        : Accumulated total number of states from density calculations.
  ! EE        : Energy value at the current point in the energy mesh.
  ! CA        : Weighted density of states (DOS) for auxiliary channel A.
  ! CB        : Weighted DOS for auxiliary channel B.
  ! CCC       : Weighted DOS for auxiliary channel C.
  ! CFC       : Weighted DOS for auxiliary channel D.
  ! CD        : Weighted DOS for auxiliary channel D (overlapping role with CFC).
  ! CC        : Density of states (DOS) value for the current energy level.
  ! C2        : Total DOS weighted by velocity squared (for transport properties).
  ! C2B       : Total DOS weighted by channel B corrections.
  ! C2C       : Total DOS weighted by channel C corrections.
  ! C2D       : Total DOS weighted by channel D corrections.
  
  ! Velocity components
  REAL(KIND=8) :: VX, VY, VZ
  ! VX, VY, VZ: Group velocity components [arbitrary units]
  
  !=======================
  ! Get data from WIEN2k output files and set arrays
  !=======================
  CALL read_constants
  ALLOCATE(AMA(NPT, MML), AKK(NPT, MML), AKKB(NPT, MML), AKKC(NPT, MML), AKKD(NPT, MML))
  ALLOCATE(AMA1(IMAX, MML))
  ALLOCATE(E(NPT), FK(NPT), FKB(NPT), FKC(NPT), FKD(NPT), EM(NPT))
  ALLOCATE(CNT(NPT), SNT(NPT), FNT(NPT), FNTB(NPT), FNTC(NPT), FNTD(NPT))
  ALLOCATE(IC(0:NPX+1, 0:NPY+1, 0:NPZ+1))
  ALLOCATE(N1(NPT), N2(NPT), N3(NPT), N4(NPT), N5(NPT), N6(NPT), N7(NPT))
  ALLOCATE(IP(3, 15, NPT), IM(3, 15, NPT))
  ALLOCATE(IEX(NPT))
  
  !=======================
  ! Initialization
  !=======================
  ED = -1.0                      ! Lower energy bound (in eV)
  EU =  1.4                      ! Upper energy bound (in eV)
  DE = (EU-ED)/EV/(NEMAX-1)      ! Energy mesh spacing in Ry
  
  !=======================
  ! Step 1: Generate tetrahedral mesh and stencil indices
  !=======================
  CALL TGEN(LAT, NPX, NPY, NPZ)  ! Builds IC and N1-N7 arrays
  CALL VDATA                     ! Reads stencil index files for IP, IM, IEX
  CALL AEIGEN                    ! Reads wien.energy and fills AMA1 with eigenvalues
  
  !=======================
  ! Step 2: Apply symmetry remapping to AMA1 -> AMA
  !=======================
  DO I = 1, NPT
    I2 = IEX(I)
    DO J = 1, MMIN
      AMA(I,J) = AMA1(I2,J)
    END DO
  END DO
  
  !=======================
  ! Step 3: Calculate group velocities for each k-point and band
  !         -> store squared velocities into AKK arrays
  ! Note: dkx = 2pi/a/N, h/2pi
  !=======================
  DO I = 1, NPT
    DO J = 1, MMIN
      CALL GV(I, J, VX, VY, VZ)  ! Calculate group velocity vector (VX, VY, VZ) for band J at k-point I
       AKK(I,J) = VX**2          ! x-component squared
      AKKB(I,J) = VX**2          ! redundant B-channel
      AKKC(I,J) = VZ**2          ! z-component squared
      AKKD(I,J) = VX             ! raw velocity
    END DO
  END DO
  
  !=======================
  ! Step 4: Loop over all bands to compute DOS and transport contributions
  !         Each band -> fills energy arrays -> calls DENS (integration over tetrahedra)
  !=======================
  DO LB = 1, MMIN
    write(6,*) LB                ! print band index
    
    ! Fill tetra_data arrays with energy and velocity info for current band
    DO IKP = 1, NPT
      E(IKP)   =  AMA(IKP, LB) / EV
      FK(IKP)  =  AKK(IKP, LB)
      FKB(IKP) = AKKB(IKP, LB)
      FKC(IKP) = AKKC(IKP, LB)
      FKD(IKP) = AKKD(IKP, LB)
    END DO
    
    ! Build energy grid EM()
    DO LE = 1, NEMAX
      EM(LE) = ED/EV + (LE-1)*DE
    END DO
    
    ! Integrate DOS and transport quantities for this band
    CALL DENS
  END DO
  
  !=======================
  ! Step 5: Write results to output file AKK.DATA
  !         Each energy level -> writes total DOS and weighted contributions
  !=======================
  OPEN(UNIT=10, FILE='AKK.DATA', STATUS='UNKNOWN')
  
  TN = 0.0D0
  
  DO LE = 1, NEMAX
    EE  = EM(LE)*EV     ! Convert back to eV
    CC  = CNT(LE)       ! Total DOS
    C2  = FNT(LE)       ! Velocity^2-weighted DOS
    C2B = FNTB(LE)      ! B-channel
    C2C = FNTC(LE)      ! C-channel
    C2D = FNTD(LE)      ! D-channel
    
    
    TN  = TN + DE*CC    ! Total number of states accumulator
    
    ! Avoid division by zero: if DOS vanishes at this energy
    IF(CC <= 0.0) THEN
      CA  = 0.0
      CCC = 0.0
      CFC = 0.0
    ELSE
      CA  = C2  / CC    ! mean velocity^2 at energy
      CB  = C2B / CC
      CCC = C2C / CC
      CD  = C2D / CC
    END IF
    
    ! Write to file: Energy, velocity^2, velocity^2 * DOS, total DOS, cumulative DOS
    WRITE(10,'(10E15.8)') EE, CA, C2/EV, CNT(LE)/EV, SNT(LE)
  END DO
  CLOSE(10)
  
  DEALLOCATE(AMA, AKK, AKKB, AKKC, AKKD)
  DEALLOCATE(AMA1)
  DEALLOCATE(E, FK, FKB, FKC, FKD, EM)
  DEALLOCATE(CNT, SNT, FNT, FNTB, FNTC, FNTD)
  DEALLOCATE(IC)
  DEALLOCATE(N1, N2, N3, N4, N5, N6, N7)
  DEALLOCATE(IP,IM)
  DEALLOCATE(IEX)
  
  STOP
END PROGRAM MAIN



!-----------------------------------------------------------------------
! Subroutine : TGEN
! Purpose    : Define tetrahedra associated with each k-point on uniform grid
!
! Inputs :
!   - LAT    : Lattice type index (e.g., 3 = BCC)
!   - NPXD,NPYD,NPZD : Grid dimensions along x,y,z
!
! Outputs (via COMMON or MODULE) :
!   - IC(X,Y,Z) -> Linear index mapping
!   - N1 - N7(NP) -> Tetrahedron neighbor indices
!-----------------------------------------------------------------------
SUBROUTINE TGEN(LAT, NPXD, NPYD, NPZD)
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  USE tetra_data ! Variables for DOS integration: LMO, E, FK, FKB, FKC, FKD, EM, CNT, SNT, FNT, FNTB, FNTC, FNTD, DE, V, BOA, COA, NEMAX, NL, NT
  USE index_data ! k-point mapping and stencil indices: IC, N1-N7, IP, IM, IEX
  IMPLICIT NONE
  
  ! Input parameters
  INTEGER, INTENT(IN) :: LAT               ! Lattice type (e.g., 3 = BCC), Unused in the current code, available for expansion.
  INTEGER, INTENT(IN) :: NPXD, NPYD, NPZD  ! Mesh divisions along x, y, z
  
  ! Loop counters and local variables
  INTEGER :: I, J, K                       ! Loop indices (spatial coordinates)
  INTEGER :: NP                            ! Linear index counter for k-points
  INTEGER :: LP                            ! Linear loop index for tetrahedron setup
  REAL(KIND=8) :: DKX, DKY, DKZ            ! Mesh spacing (delta k)
  
  !-------------------------------
  ! Initialize IC: map (i,j,k) -> linear index NP
  ! Set all IC entries to zero before assigning
  !-------------------------------
  DO K = 0, NPZ
    DO J = 0, NPY
      DO I = 0, NPX
        IC(I, J, K) = 0  ! I=X, J=Y, K=Z
      END DO
    END DO
  END DO
  
  ! FCC IRREDUCIBLE ZONE DEFINED BY
  ! 0 <= KZ, KY, or, KX =< 2PI/A
  ! KX + KY + KZ <= 3PI/A
  
  ! Set up parameters
  NP = 0
  
  !-------------------------------
  ! Compute spacing in reciprocal space
  ! and tetrahedron volume V
  !-------------------------------
  DKX = 2.D0 / NPXD
  DKY = 2.D0 / NPYD
  DKZ = 2.D0 / NPZD
  
  ! Effective volume assigned to each tetrahedron (scaled with symmetry factors)
  V = (DKX*DKY*DKZ) / 3.0 / 8.0 / 2.0
  
  !-------------------------------
  ! Assign linear indices to 3D grid points via IC(I,J,K)
  !-------------------------------
  ! Define tetrahedra and assign indices
  DO K = 0, NPZD
    DO J = 0, NPYD
      DO I = 0, NPXD
        NP = NP + 1
        IC(I, J, K) = NP  ! I=X, J=Y, K=Z
      END DO
    END DO
  END DO
  
  WRITE(6, *) "Total k-points: ", NP
  
  !-------------------------------
  ! Assign neighboring tetrahedron vertices for each interior k-point
  !-------------------------------
  NP = 0
  LP = 0
  DO K = 0, NPZD
    DO J = 0, NPYD
      DO I = 0, NPXD
        LP = LP + 1
        NP = IC(I, J, K)  ! I=X, J=Y, K=Z: Linear index of current point
        
        ! Boundary points (cannot form full tetrahedra)
        IF (I == NPXD .OR. J == NPYD .OR. K == NPZD) THEN
          N1(NP) = 0
          N2(NP) = 0
          N3(NP) = 0
          N4(NP) = 0
          N5(NP) = 0
          N6(NP) = 0
          N7(NP) = 0
        ELSE
          ! Interior points -> assign 7 adjacent vertex indices
          ! GENERAL POINT (I=X, J=Y, K=Z)
          ! Calculate the adjacent index around the lattice point.
          N1(NP) = IC(I+1, J,   K  )   ! +x
          N2(NP) = IC(I+1, J+1, K  )   ! +x,+y
          N3(NP) = IC(I+1, J+1, K+1)   ! +x,+y,+z
          N4(NP) = IC(I+1, J,   K+1)   ! +x,+z
          N5(NP) = IC(I  , J,   K+1)   ! +z
          N6(NP) = IC(I  , J+1, K+1)   ! +y,+z
          N7(NP) = IC(I  , J+1, K  )   ! +y
        END IF
      END DO
    END DO
  END DO
  
  RETURN
END SUBROUTINE TGEN



!-----------------------------------------------------------------------
! Subroutine : PDNS
! Purpose    : Evaluate DOS and transport quantities for one tetrahedron
!              given 4 vertex indices and their associated weights.
!
! Inputs  :
!   - I1, I2, I3, I4 : Indices of the four vertices of the tetrahedron
!
! Outputs (indirect) :
!   - CNT, SNT       : Accumulated DOS and state count
!   - FNT, FNTB, FNTC, FNTD : Transport-weighted integrals (velocity channels)
!
! Dependencies:
!   - Uses global arrays E(), FK(), FKB(), FKC(), FKD(), EM()
!   - Writes into CNT(), SNT(), FNT(), etc. at grid points LE
!-----------------------------------------------------------------------
SUBROUTINE PDNS(I1, I2, I3, I4)
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  USE tetra_data ! Variables for DOS integration: LMO, E, FK, FKB, FKC, FKD, EM, CNT, SNT, FNT, FNTB, FNTC, FNTD, DE, V, BOA, COA, NEMAX, NL, NT
  IMPLICIT NONE
  
  INTEGER :: I1, I2, I3, I4              ! indices of tetrahedron vertices
  INTEGER :: LE, NE1, NE2                ! energy mesh index variables
  REAL(KIND=8) :: S1, EB, EE             ! energies and current energy value
  REAL(KIND=8) :: DELTA, AM, D1, D2, D4  ! interpolation parameters
  REAL(KIND=8) :: CN, SN                 ! DOS contribution and state count
  REAL(KIND=8) :: FN, FNB, FNC, FND      ! transport-weighted DOS values
  REAL(KIND=8) :: FRAC, V1, DEE          ! intermediate calculation values
  REAL(KIND=8) :: CN1, CN2               ! partial DOS terms (used in triangle decomposition)
  REAL(KIND=8) :: DF, DF1, DF2           ! derivative-related weights
  REAL(KIND=8) :: SSS                    ! shared interpolation scaling
  
  REAL(KIND=8), DIMENSION(4) :: S, F, FB, FC, FD  ! energies and weighted velocities at vertices
  
  !-------------------------------------------------------------
  ! Load energies and transport values for each tetrahedron vertex
  !-------------------------------------------------------------
  S(1)=E(I1); F(1)=FK(I1); FB(1)=FKB(I1); FC(1)=FKC(I1); FD(1)=FKD(I1)
  S(2)=E(I2); F(2)=FK(I2); FB(2)=FKB(I2); FC(2)=FKC(I2); FD(2)=FKD(I2)
  S(3)=E(I3); F(3)=FK(I3); FB(3)=FKB(I3); FC(3)=FKC(I3); FD(3)=FKD(I3)
  S(4)=E(I4); F(4)=FK(I4); FB(4)=FKB(I4); FC(4)=FKC(I4); FD(4)=FKD(I4)
  
  CALL SORT(S, F, FB, FC, FD)  ! sort all arrays by ascending energy
  
  !-------------------------------------------------------------
  ! Determine energy mesh range for which this tetrahedron contributes
  !-------------------------------------------------------------
  S1 = S(1)
  IF(S1 >= EM(NEMAX)) RETURN   ! no contribution if below energy mesh
  
  EB = EM(1)                   ! base mesh energy
  NE1 = INT((S1-EB)/DE)        ! starting index for energy integration
  IF(S1 < EB) NE1 = -1         ! 
  NE1 = NE1 + 2                ! Correct Index
  NE2 = NEMAX                  ! maximum mesh index
  
  !-------------------------------------------------------------
  ! Loop through relevant energy mesh points and accumulate contribution
  !------------- 
  ! Note:
  ! In the tetrahedron interpolation method, a simple linear approximation may not be able to capture the sudden changes and 
  ! curvature of the energy-dependent physical quantities F (e.g., velocity and conductivity). Therefore, 
  ! by combining linear and nonlinear approximations, it is possible to more faithfully reproduce the changes within the energy range.
  ! - Linear component (e.g. F(1)*CN)
  ! - Nonlinear component (correction term using the interpolation gradients DF1 and DF2)
  !-------------------------------------------------------------
  DO 25 LE = NE1, NE2
    EE = EM(LE)
    
    IF(EE <= S(1)) GO TO 25  ! skip energy below tetrahedron range
    IF(EE >= S(4)) GO TO 20  ! skip energy above range -> only SN=V later
    
    !------------- CASE 1: S(1) < EE <= S(2) ----------------
    IF(EE <= S(2)) GO TO 21
    
    !------------- CASE 3: S(3) <= EE < S(4) ----------------
    IF(EE >= S(3)) GO TO 22
    
  !------------- CASE 2: S(2) < EE < S(3) ----------------
  ! Interpolate in central range of tetrahedron
  !------------- ----------
  DELTA = S(4)+S(3)-S(2)-S(1)         ! Energy range spanned by tetrahedron vertices
  AM = (S(4)*S(3)-S(2)*S(1)) / DELTA  ! Weighted average energy in the tetrahedron
  D1 = S(1) - AM                      ! Deviation of lowest vertex from average energy
  D2 = S(2) - AM                      ! Deviation of second-lowest vertex from average energy
  V1 = V / DELTA                      ! Normalized volume per energy unit
  DEE = EE - AM                       ! Energy difference from average to evaluation point
  FRAC = DEE * DEE / D2 / D1          ! Curvature term for parabolic weight
  CN = 3.D0 * V1 * (1.D0-FRAC)        ! Density of states contribution
  SN = 3.D0 * V1 * DEE-V1 * (D1+D2+DEE*FRAC)  ! Transport quantity contribution
  !------------- ----------
  
  !------------- ----------
  ! Handle near-degenerate tetrahedron geometry
  !------------- ----------
  IF (ABS(S(2)-S(1)) <= 0.0) THEN
    SSS = 3.D0 * V * (EE-S(1))/(S(3)-S(1))/(S(4)-S(1))  ! Interpolation of density of states  (SSS)
    CN1 = SSS * (S(3)-EE)/(S(3)-S(1))                   ! Contribution of small triangle part (CN1)
    CN2 = SSS * (S(4)-EE)/(S(4)-S(1))                   ! Contribution of the Great Triangle  (CN2)
    
    !------------- ----------
    ! weighted interpolation for FN and channels
    ! ----------
    ! Note (DF1): This is a first-order correction term (DF1) used to find the interpolated value near EE using
    !             the change (gradient) of the F value within the tetrahedron.
    ! EE is the energy point to be evaluated (the energy to be calculated)
    ! S(n) is the energy of the four vertices of the tetrahedron
    ! F(n) is the physical quantity at that vertex (e.g. velocity or conductivity).
    ! ----------
    ! Note (DF2): DF2 is adjusted by adjusting the weighting of the coefficients of DF1 in a different way,
    !             so that the contribution of F(4) is larger. This is a method to improve the accuracy of
    !             the interpolation when the position of EE is close to S(4).
    ! ----------
    ! Note (FN): The overall contribution of F is calculated by interpolating the distribution of 
    !            F in the energy direction within the tetrahedron. The resulting FN is a weighted average of
    !            the physical quantities corresponding to EE.
    ! F(1)*CN: Direct contribution of the base point (vertex 1) (linear)
    ! DF1*CN1 + DF2*CN2: Contribution of the interpolation gradient term (nonlinear)
    ! /3.D0: Averaging across tetras (statistical balance)
    !------------- ----------
! case F:
    DF1 =  (F(2)-F(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
      & +2*(F(3)-F(1))*(EE-S(1))/(S(3)-S(1)) &
      & +  (F(4)-F(1))*(EE-S(1))/(S(4)-S(1))
    DF2 =  (F(2)-F(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
      & +  (F(3)-F(1))*(EE-S(1))/(S(3)-S(1)) &
      & +2*(F(4)-F(1))*(EE-S(1))/(S(4)-S(1))
    FN  =   F(1)*CN + (DF1*CN1+DF2*CN2)/3.D0
! case FB:
    DF1 =  (FB(2)-FB(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
      & +2*(FB(3)-FB(1))*(EE-S(1))/(S(3)-S(1)) &
      & +  (FB(4)-FB(1))*(EE-S(1))/(S(4)-S(1))
    DF2 =  (FB(2)-FB(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
      & +  (FB(3)-FB(1))*(EE-S(1))/(S(3)-S(1)) &
      & +2*(FB(4)-FB(1))*(EE-S(1))/(S(4)-S(1))
    FNB =   FB(1)*CN + (DF1*CN1+DF2*CN2)/3.D0
! case FC:
    DF1 =  (FC(2)-FC(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
      & +2*(FC(3)-FC(1))*(EE-S(1))/(S(3)-S(1)) &
      & +  (FC(4)-FC(1))*(EE-S(1))/(S(4)-S(1))
    DF2 =  (FC(2)-FC(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
      & +  (FC(3)-FC(1))*(EE-S(1))/(S(3)-S(1)) &
      & +2*(FC(4)-FC(1))*(EE-S(1))/(S(4)-S(1))
    FNC =   FC(1)*CN + (DF1*CN1+DF2*CN2)/3.D0
! case FD:
    DF1 =  (FD(2)-FD(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
      & +2*(FD(3)-FD(1))*(EE-S(1))/(S(3)-S(1)) &
      & +  (FD(4)-FD(1))*(EE-S(1))/(S(4)-S(1))
    DF2 =  (FD(2)-FD(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
      & +  (FD(3)-FD(1))*(EE-S(1))/(S(3)-S(1)) &
      & +2*(FD(4)-FD(1))*(EE-S(1))/(S(4)-S(1))
    FND =   FD(1)*CN+ (DF1*CN1+DF2*CN2)/3.D0
    ! -----
    GOTO 23
  END IF
  
  !------------- ----------
  ! When the energy difference between the vertices of a tetrahedron is very small
  ! Since normal interpolation can become unstable in this case, special logic is used to stabilize the process.
  ! ----------
  ! Note (DF): Approximate the energy-direction interpolation difference (gradient) of the physical quantity F (e.g. velocity, conductivity, etc.)
  ! F(2)-F(1) and F(3)-F(1) are the changes from vertex 1 to the adjacent vertex.
  ! (F(4)-F(1))*(EE-S(1))/(S(4)-S(1)) is the interpolation term (linear interpolation) for the target energy EE.
  ! ----------
  ! Note (FN): Calculate the physical quantity FN as an interpolated value within the tetrahedron at the energy EE
  ! F(1)*CN: Contribution of vertex 1 (constant component)
  ! DF*CN/3: Gradient-based interpolation component (first order approximation)
  !------------- ----------
  IF(ABS(S(3)-S(2)) <= 0.0) THEN
    ! simplified linear formula for degenerate configuration
! case F:
    DF = F(2)-F(1)+F(3)-F(1)+(F(4)-F(1))*(EE-S(1))/(S(4)-S(1))
    FN =  F(1)*CN + DF*CN/3.D0
! case FB:
    DF = FB(2)-FB(1)+FB(3)-FB(1)+(FB(4)-FB(1))*(EE-S(1))/(S(4)-S(1))
    FNB= FB(1)*CN + DF*CN/3.D0
! case FD:
    DF = FD(2)-FD(1)+FD(3)-FD(1)+(FD(4)-FD(1))*(EE-S(1))/(S(4)-S(1))
    FND= FD(1)*CN + DF*CN/3.D0
! case FC:
    DF = FC(2)-FC(1)+FC(3)-FC(1)+(FC(4)-FC(1))*(EE-S(1))/(S(4)-S(1))
    FNC= FC(1)*CN + DF*CN/3.D0
    ! -----
    GOTO 23
  ENDIF
  
  !------------- ----------
  ! general interpolation with CN2 and CN1 decomposition
  ! CN1=small tri, CN2=big tri
  !------------- ----------
  CN2=3.D0*V *(EE-S(1))**2/(S(2)-S(1))/(S(3)-S(1))/(S(4)-S(1))
  CN1=CN2-CN
  
  !------------- ----------
  ! Note (DF2): Using vertex 1 as the base point, we calculate the rate of change (gradient) of the physical quantity F at the other vertices (2 to 4).
  ! Each term is the slope normalized by the energy difference between the vertices.
  ! We sum them together and multiply them by the distance from the evaluation energy E to vertex 1 (EE - S(1)) to
  !   get the contribution from linear interpolation (a straight-line approximation).
  !-------------
  ! Note (DF1): A more complex interpolation structure allows for asymmetric variations over the energy range S(1) to S(4).
  ! The weighted part of F(2)-F(1) is designed so that the influence of EE becomes stronger when it approaches vertex 2 and weakens when it moves away.
  ! In the second half, the contributions of F(3) and F(4) are complemented toward EE by adjusting the slope from interval S(2).
  !-------------
  ! Note (FN): This is a formula for highly accurate interpolation and weighted averaging of physical quantity F at energy E within a tetrahedron.
  ! F(1)*CN: Density contribution (basis component) obtained directly from the physical quantities of vertex 1.
  ! (DF2*CN2 - DF1*CN1)/3: This is the addition/subtraction term due to the interpolation gradient, which is used to improve the approximation accuracy.
  ! CN, CN1, and CN2 indicate the volume coefficients and contribution rates in the respective energy regions, which are averaged over the whole.
  !------------- ----------
! case F:
  DF2 = ( (F(2)-F(1))/(S(2)-S(1)) +(F(3)-F(1))/(S(3)-S(1)) &
      & + (F(4)-F(1))/(S(4)-S(1)) )*(EE-S(1))
  DF1 = (F(2)-F(1))*( (EE-S(1))/(S(2)-S(1))+(S(3)-EE)/(S(3)-S(2)) &
      & + (S(4)-EE)/(S(4)-S(2)) ) &
      & + ( (F(3)-F(1))/(S(3)-S(2))+(F(4)-F(1))/(S(4)-S(2)) )*(EE-S(2))
  FN  =  F(1)*CN + (DF2*CN2-DF1*CN1)/3.D0
! case FB:
  DF2 = ( (FB(2)-FB(1))/(S(2)-S(1)) +(FB(3)-FB(1))/(S(3)-S(1)) &
      & + (FB(4)-FB(1))/(S(4)-S(1)) )*(EE-S(1))
  DF1 = (FB(2)-FB(1))*( (EE-S(1))/(S(2)-S(1))+(S(3)-EE)/(S(3)-S(2)) &
      & + (S(4)-EE)/(S(4)-S(2)) ) &
      & + ( (FB(3)-FB(1))/(S(3)-S(2))+(FB(4)-FB(1))/(S(4)-S(2)) )*(EE-S(2))
  FNB = FB(1)*CN + (DF2*CN2-DF1*CN1)/3.D0
! case FD:
  DF2 = ( (FD(2)-FD(1))/(S(2)-S(1)) +(FD(3)-FD(1))/(S(3)-S(1)) &
      & + (FD(4)-FD(1))/(S(4)-S(1)) )*(EE-S(1))
  DF1 = (FD(2)-FD(1))*( (EE-S(1))/(S(2)-S(1))+(S(3)-EE)/(S(3)-S(2)) &
      & + (S(4)-EE)/(S(4)-S(2)) ) &
      & + ( (FD(3)-FD(1))/(S(3)-S(2))+(FD(4)-FD(1))/(S(4)-S(2)) )*(EE-S(2))
  FND = FD(1)*CN + (DF2*CN2-DF1*CN1)/3.D0
! case FC:
  DF2 = ( (FC(2)-FC(1))/(S(2)-S(1)) +(FC(3)-FC(1))/(S(3)-S(1)) &
      & + (FC(4)-FC(1))/(S(4)-S(1)) )*(EE-S(1))
  DF1 = (FC(2)-FC(1))*( (EE-S(1))/(S(2)-S(1))+(S(3)-EE)/(S(3)-S(2)) &
      & + (S(4)-EE)/(S(4)-S(2)) ) &
      & + ( (FC(3)-FC(1))/(S(3)-S(2))+(FC(4)-FC(1))/(S(4)-S(2)) )*(EE-S(2))
  FNC = FC(1)*CN + (DF2*CN2-DF1*CN1)/3.D0
  ! -----
  GOTO 23

  !------------- CASE 3: EE > S(3) ----------------
  ! When the energy difference between the vertices of a tetrahedron is very small
  ! Since normal interpolation can become unstable in this case, special logic is used to stabilize the process.
  !------------- 
  22 IF(ABS(S(4)-S(3)) <= 0.0) THEN
    CN=0.0; SN=0.0
    FN=0.0; FNB=0.0; FNC=0.0; FND=0.0
    GOTO 23
  ENDIF
  
  !------------- ----------
  ! standard integration formulas for upper segment
  ! Accurately interpolate the energy distribution within a tetrahedron using the Energy Interpolation Criterion (AM) 
  ! to efficiently calculate the density of states (CN) and number of states (SN)
  !------------- ----------
  DELTA = S(4)+S(3)-S(2)-S(1)         ! Energy span across the tetrahedron vertices
  AM = (S(4)*S(3)-S(2)*S(1)) / DELTA  ! Weighted mean energy of the tetrahedron
  D4 = S(4) - AM                      ! Deviation of top vertex from average energy
  V1 = V / DELTA                      ! Normalized tetrahedron volume by energy span
  DEE = EE - S(4)                     ! Energy offset from top vertex to current energy
  FRAC = DEE * DEE / D4 / (S(4)-S(3)) ! Quadratic curvature weight from vertex gradient
  CN = 3.D0 * V1 * FRAC               ! DOS contribution from this tetrahedron
  SN = V1 * (DELTA+DEE*FRAC)          ! Transport contribution (e.g. conductivity term)
  !------------- ----------
  
  !------------- ----------
  ! Note (DF): Interpolation amount with gradient correction of physical quantity F based on vertex 1
  ! F(n)-F(1): Change in physical quantity from vertex 1 to other vertices
  ! S(4)-EE and EE-S(n): Relative distance between the energy position EE and each vertex
  ! (S(4)-S(n)): Energy difference between vertices (interpolation range)
  !-------------
  ! Note (FN): By summing, the value of the physical quantity E at the evaluation energy
  !            F is smoothly reproduced based on the tetrahedron interpolation method.
  ! F(1)*CN: Direct contribution from the physical quantity of vertex 1 (linear component)
  ! DF*CN/3: Interpolation component with gradient correction (nonlinear approximation)
  !------------- ----------
! case F:
  DF =  (F(2)-F(1))*(S(4)-EE)/(S(4)-S(2))  &
     & +(F(3)-F(1))*(S(4)-EE)/(S(4)-S(3))  &
     & +(F(4)-F(1))*((EE-S(1))/(S(4)-S(1)) &
     & +(EE-S(2))/(S(4)-S(2))+(EE-S(3))/(S(4)-S(3)))
  FN = F(1)*CN + DF*CN/3.D0
! case FB:
  DF =  (FB(2)-FB(1))*(S(4)-EE)/(S(4)-S(2))  &
     & +(FB(3)-FB(1))*(S(4)-EE)/(S(4)-S(3))  &
     & +(FB(4)-FB(1))*((EE-S(1))/(S(4)-S(1)) &
     & +(EE-S(2))/(S(4)-S(2))+(EE-S(3))/(S(4)-S(3)))
  FNB = FB(1)*CN + DF*CN/3.D0
! case FD:
  DF =  (FD(2)-FD(1))*(S(4)-EE)/(S(4)-S(2))  &
     & +(FD(3)-FD(1))*(S(4)-EE)/(S(4)-S(3))  &
     & +(FD(4)-FD(1))*((EE-S(1))/(S(4)-S(1)) &
     & +(EE-S(2))/(S(4)-S(2))+(EE-S(3))/(S(4)-S(3)))
  FND = FD(1)*CN + DF*CN/3.D0
! case FC:
  DF =  (FC(2)-FC(1))*(S(4)-EE)/(S(4)-S(2))  &
     & +(FC(3)-FC(1))*(S(4)-EE)/(S(4)-S(3))  &
     & +(FC(4)-FC(1))*((EE-S(1))/(S(4)-S(1)) &
     & +(EE-S(2))/(S(4)-S(2))+(EE-S(3))/(S(4)-S(3)))
  FNC=FC(1)*CN +DF*CN/3.D0
  GO TO 23

  !------------- CASE 1: EE < S(2) ----------------
  ! When the energy difference between the vertices of a tetrahedron is very small
  ! Since normal interpolation can become unstable in this case, special logic is used to stabilize the process.
  !------------- 
  21 IF(ABS(S(2)-S(1)) <= 0.0) THEN
    CN=0.0; SN=0.0
    FN=0.0; FNB=0.0; FNC=0.0; FND=0.0
    GOTO 23
  ENDIF
  
  !------------- ----------
  ! compute interpolation for low-energy region
  ! This procedure is used to calculate the contributions near the minimum energy of the tetrahedrons, 
  ! and efficiently interpolates between specific energy points to allow stable calculation of DOS and SN values.
  !------------- ----------
  DELTA = S(4)+S(3)-S(2)-S(1)          ! Total energy spread across tetrahedron vertices
  AM = (S(4)*S(3)-S(2)*S(1)) / DELTA   ! Weighted average energy used for integration center
  D1 = S(1) - AM                       ! Deviation of lowest vertex from the average energy
  V1 = V / DELTA                       ! Volume normalization factor per energy span
  DEE = EE - S(1)                      ! Energy offset from lowest vertex to evaluation energy
  FRAC = DEE * DEE / D1 / (S(1)-S(2))  ! Quadratic weighting factor for integration
  CN = 3.D0 * V1 * FRAC                ! Density of states contribution from the tetrahedron
  SN = V1 * FRAC * DEE                 ! Transport contribution scaled by energy offset
  !------------- ----------
  
  !------------- ----------
  ! Note (DF): The rate of change (gradient) of the physical quantity F at each vertex is calculated based on vertex 1.
  ! (F(n)-F(1))/(S(n)-S(1)) is the slope from vertex 1 to the nth vertex (the rate of change of the physical quantity with respect to energy).
  ! By adding them together, we obtain the average change trend around EE.
  ! Finally, by multiplying by DEE = EE - S(1), we obtain the interpolation amount (first order approximation) to the evaluation energy position.
  !-------------
  ! Note (FN): The overall contribution value FN of the physical quantities in the evaluated energy EE is obtained.
  ! F(1)*CN: Directly reflects the contribution of the physical quantity at vertex 1 (stationary component).
  ! DF*CN/3: Add an interpolated adjustment value based on the slope (DF) (change component).
  !------------- ----------
! case F:
  DF = ( (F(2)-F(1))/(S(2)-S(1))+ (F(3)-F(1))/(S(3)-S(1)) &
     &  +(F(4)-F(1))/(S(4)-S(1)) )*DEE
  FN = F(1)*CN + DF*CN/3.D0
! case FB:
  DF = ( (FB(2)-FB(1))/(S(2)-S(1))+ (FB(3)-FB(1))/(S(3)-S(1)) &
     &  +(FB(4)-FB(1))/(S(4)-S(1)) )*DEE
  FNB = FB(1)*CN + DF*CN/3.D0
! case FD:
  DF = ( (FD(2)-FD(1))/(S(2)-S(1))+ (FD(3)-FD(1))/(S(3)-S(1)) &
     &  +(FD(4)-FD(1))/(S(4)-S(1)) )*DEE
  FND = FD(1)*CN + DF*CN/3.D0
! case FC:
  DF = ( (FC(2)-FC(1))/(S(2)-S(1))+ (FC(3)-FC(1))/(S(3)-S(1)) &
     &  +(FC(4)-FC(1))/(S(4)-S(1)) )*DEE
  FNC = FC(1)*CN + DF*CN/3.D0
  GO TO 23
  
   !------------- CASE 4: EE > S(4) ----------------
20 SN = V  ! assign full volume to state count
   CN = 0.D0
   FN = 0.D0; FNB = 0.D0; FNC = 0.D0; FND = 0.D0
    
   !------------- Final accumulation ----------------
23 CONTINUE
   SNT(LE) = SNT(LE)  + SN   ! Sum of Number of States
   CNT(LE) = CNT(LE)  + CN   ! Density of States, DOS
   FNT(LE) = FNT(LE)  + FN   ! Corrected DOS
  FNTB(LE) = FNTB(LE) + FNB  ! Corrected DOS for Channel B
  FNTC(LE) = FNTC(LE) + FNC  ! Corrected DOS for Channel C
  FNTD(LE) = FNTD(LE) + FND  ! Corrected DOS for Channel D
  
25 CONTINUE
  
  RETURN
END SUBROUTINE PDNS



!================================================================================
! Subroutine : DENS
! Purpose    : Compute density of states (DOS) and transport-related quantities
!              by integrating contributions from all tetrahedra over the 3D k-grid.
!
! Description:
!   - Iterates over all interior k-points defined on a uniform mesh.
!   - For each cube (8 surrounding points), constructs six tetrahedra using standard split.
!   - Calls PDNS for each tetrahedron to evaluate its energy-resolved contribution.
!
! Inputs :
!   - IC(I,J,K) : Maps 3D k-grid coordinates to linear index (NPT-based)
!   - N1-N7(NPT): Neighbor point indices forming tetrahedron vertices
!
! Dependencies:
!   - constants     : Provides NPX, NPT, etc.
!   - index_data    : Provides IC(), N1-N7()
!   - PDNS          : Called to evaluate DOS contributions per tetrahedron
!
! Notes:
!   - Assumes k-mesh cube size of (NPX+1)3 points
!   - Boundary points (I=NPX etc.) are skipped to avoid out-of-range access
!================================================================================
SUBROUTINE DENS()
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  USE index_data ! k-point mapping and stencil indices: IC, N1-N7, IP, IM, IEX
  IMPLICIT NONE
  
  ! Loop indices for spatial traversal
  INTEGER :: I, J, K
  
  ! Linear index variable
  INTEGER :: LP
  
  ! Vertex indices for tetrahedra
  INTEGER :: I0, I1, I2, I3, I4, I5, I6, I7
  
  !---------------------------------------------------------------
  ! Iterate through all mesh points in 3D k-space
  !---------------------------------------------------------------
  DO K = 0, NPX
    DO J = 0, NPX
      DO I = 0, NPX
        
        ! Convert grid (I,J,K) to linear index LP
        LP = IC(I, J, K) ! I=X, J=Y, K=Z
        I0 = LP          ! Central point of current cube
        I1 = N1(LP)
        I2 = N2(LP)
        I3 = N3(LP)
        I4 = N4(LP)
        I5 = N5(LP)
        I6 = N6(LP)
        I7 = N7(LP)
        ! I1-I7 define the surrounding cube corners used to construct tetrahedra
        
        ! Skip edge/boundary points where neighbor indices may go out-of-range
        IF(I == NPX .OR. J == NPX .OR. K == NPX) CYCLE
        
        !-----------------------------------------------------------
        ! Evaluate DOS and transport integrals via 6 tetrahedra split
        ! Each PDNS call processes one tetrahedron defined by 4 vertices
        !-----------------------------------------------------------
        CALL PDNS(I0, I1, I2, I3)  ! Lower left tetrahedron
        CALL PDNS(I0, I1, I3, I4)  ! Lower front
        CALL PDNS(I0, I3, I4, I5)  ! Right front
        CALL PDNS(I0, I3, I5, I6)  ! Upper front
        CALL PDNS(I0, I3, I6, I7)  ! Upper back
        CALL PDNS(I0, I2, I3, I7)  ! Upper left
        
        ! Each tetrahedron contributes to CNT(), SNT(), FNT(), etc.
        ! accumulators defined in tetra_data, binned by energy mesh EM()
      END DO
    END DO
  END DO
  
  RETURN
END SUBROUTINE DENS



!-----------------------------------------------------------------------
! Subroutine : SORT
! Purpose    : Sort the array S(1:4) in increasing order, and reorder
!              associated physical arrays F, G, O, P accordingly.
!
! Inputs/Outputs:
!   - S(4) : Energy values to sort (ascending)
!   - F(4), G(4), O(4), P(4) : Associated transport or weighting arrays
!
! SORT ELEMENTS OF S IN INCREASING ORDER
!-----------------------------------------------------------------------
SUBROUTINE SORT(S, F, G, O, P)
  IMPLICIT NONE
  
  ! Input/output arrays:
  ! S : Energy values at 4 tetrahedron vertices
  ! F, G, O, P : Associated quantities (e.g. velocities, transport weights)
  REAL(KIND=8), DIMENSION(4) :: S, F, G, O, P
  
  ! Loop variables
  INTEGER :: J, I, NP
  
  ! Temporary variables used during element swapping
  REAL(KIND=8) :: ST, FT, GT, OT, PT
  
  !-----------------------------------------------
  ! Bubble sort implementation for 4 elements.
  ! Although bubble sort is not efficient for large datasets,
  ! it's sufficient here due to fixed small size (n=4).
  ! Time complexity: O(n^2), but negligible for n=4.
  !-----------------------------------------------
  DO J = 1, 3
    NP = 4 - J  ! On each pass, the largest element bubbles to the end,
                ! so we reduce comparison range by J
    DO I=1, NP
      ! Compare adjacent energy values S(I) and S(I+1)
      IF(S(I+1) >= S(I)) CYCLE ! If already ordered, skip to next pair
      
      ST = S(I)
      FT = F(I)
      GT = G(I)
      OT = O(I)
      PT = P(I)
      
      S(I) = S(I+1)
      F(I) = F(I+1)
      G(I) = G(I+1)
      O(I) = O(I+1)
      P(I) = P(I+1)
      
      S(I+1) = ST
      F(I+1) = FT
      G(I+1) = GT
      O(I+1) = OT
      P(I+1) = PT
      
    END DO
  END DO
  
  RETURN
END SUBROUTINE SORT



!-----------------------------------------------------------------------
! Subroutine : AEIGEN
! Purpose    : Reads band eigenvalues from 'wien.energy' and stores
!              them in AMA1(IMAX, MML) after unit conversion and energy alignment.
!
! Description:
!   - Applies a reference shift EF to all eigenvalues (aligning bands)
!   - Converts values from Rydberg to eV using EV
!   - Tracks the maximum number of bands (JMAX)
!   - Final count of usable bands stored in MMIN = JMAX - JCUT
!
! Assumptions:
!   - Each block corresponds to a k-point with JC bands listed
!   - Only bands above cutoff JCUT are retained
!-----------------------------------------------------------------------
SUBROUTINE AEIGEN()
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  USE band_data  ! Band energies: AMA1(IMAX,MML), total bands MMIN
  IMPLICIT NONE
  
  ! Band cutoff index: bands below JCUT will be ignored
  INTEGER, PARAMETER :: JCUT = 0
  
  ! Local variables
  REAL(KIND=8) :: TI     ! Unused: could be used for symmetry checks
  REAL(KIND=8) :: EIGA   ! Raw energy value read from file [Ry]
  REAL(KIND=8) :: EIGAA  ! Converted energy [eV]
  INTEGER :: IA, IB, JC  ! IA: k-point index, IB: unused, JC: bands per k-point
  INTEGER :: JJ          ! Original band index (from file)
  INTEGER :: JA          ! Shifted band index (JJ - JCUT)
  INTEGER :: JMAX        ! Minimum number of bands found across all k-points
  INTEGER :: J           ! Loop index over bands
  
  ! Initialize
  TI = 2.0 * PI
  JMAX = 9000            ! Set upper bound for bands (overwritten below)
  
  OPEN(UNIT=4, FILE='wien.energy', STATUS='OLD')
  READ(4, '(1/)')
  
20 CONTINUE
  ! Read header line for each k-point: IA, IB = indices; JC = band count
  READ(4,'(62X, I6, I6, I6)') IA, IB, JC
  
  ! Update minimum band count seen so far
  IF(JC < JMAX) JMAX = JC
  
  ! Read band energies for this k-point
  DO J = 1, JC
    ! JJ: band number, EIGA: energy [Ry]: (old version: "READ(4,'(I12, E28.15)') JJ, EIGA")
    READ(4,*) JJ, EIGA
    
    ! Ignore bands below cutoff
    IF(JJ > JCUT) THEN
      EIGAA = (EIGA - EF) * EV  ! Shift and convert to eV
      JA = JJ - JCUT            ! Re-index to store into AMA1
      
      AMA1(IA, JA) = EIGAA      ! Save band energy
      
      ! Optional debug output: print first k-point data
      IF(IA == 1) WRITE(6, '(2I4, 2F15.8, I4)') IA, JA, EIGAA, EIGA, J
    END IF
  END DO
  
  ! If last k-point reached, exit loop
  IF(IA == IMAX) GOTO 30
  
  ! Otherwise continue reading next k-point
  GOTO 20
  
  ! Finalize
  CLOSE(4)
  
30 CONTINUE
  ! Set number of usable bands
  MMIN = JMAX - JCUT
  
  RETURN
END SUBROUTINE AEIGEN



!-----------------------------------------------------------------------
! Subroutine : VDATA
! Purpose    : Load precomputed finite-difference stencil indices used
!              for group velocity (GV) calculation via high-order differencing.
!
! Description:
!   - Reads 9 external data files: cf{A,B,C}{1,2,3}.dat
!   - Each file contains 6 pairs of forward/backward k-point indices
!     for a given direction (x, y, z).
!   - Final file in each group contains last 3 pairs + symmetry mapping index.
!
! Result:
!   - Populates:
!       IP(3,15,NPT) : forward stencil indices for [x,y,z]
!       IM(3,15,NPT) : backward stencil indices for [x,y,z]
!       IEX(NPT)     : k-point remapping index
!
! Usage:
!   - Must be called before GV to initialize stencil layout
!-----------------------------------------------------------------------
SUBROUTINE VDATA()
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  USE index_data ! k-point mapping and stencil indices: IC, N1-N7, IP, IM, IEX
  IMPLICIT NONE
  
  INTEGER :: L   ! Loop over all k-points
  
  !----------------------------
  ! Open stencil files for each direction
  ! A -> x, B -> y, C -> z
  ! Each direction has 3 files (15-point stencil split into blocks)
  ! ----------
  ! IP(dir, k): The positive neighbor index in the dir direction from the kth point.
  ! IM(dir, k): The negative neighbor index in the dir direction from the kth point.
  ! Corresponding to the grid point number (integer), it describes the relationship between each point and the surrounding points
  ! A dictionary for extracting nearby points that are the subject of difference operations when calculating group velocities and energy gradients.
  !----------------------------
  OPEN(UNIT=01, FILE='cfA1.dat', STATUS='OLD')  ! x-direction: points  1-6 : Note: exist -1 (Line 54, Column 11)
  OPEN(UNIT=02, FILE='cfA2.dat', STATUS='OLD')  ! x-direction: points  7-12
  OPEN(UNIT=03, FILE='cfA3.dat', STATUS='OLD')  ! x-direction: points 13-15 + IEX
  
  OPEN(UNIT=11, FILE='cfB1.dat', STATUS='OLD')  ! y-direction: points  1-6 : Note: exist 0 - -16
  OPEN(UNIT=12, FILE='cfB2.dat', STATUS='OLD')  ! y-direction: points  7-12
  OPEN(UNIT=13, FILE='cfB3.dat', STATUS='OLD')  ! y-direction: points 13-15
  
  OPEN(UNIT=21, FILE='cfC1.dat', STATUS='OLD')  ! z-direction: points  1-6 : Note: exist 0 - -56
  OPEN(UNIT=22, FILE='cfC2.dat', STATUS='OLD')  ! z-direction: points  7-12
  OPEN(UNIT=23, FILE='cfC3.dat', STATUS='OLD')  ! z-direction: points 13-15
  
  !----------------------------
  ! Loop over all k-points and read their stencil indices
  ! For each k-point:
  !   - Read 15 pairs of forward/backward indices in x, y, z
  !   - IEX: remapping index (from symmetry)
  ! ----------
  ! IP(direction, stencil_index, k_point): Forward  (positive) k-point index used in central difference method
  ! IM(direction, stencil_index, k_point): Backward (negative) k-point index used in central difference method
  !  - direction: 1 = x-axis, 2 = y-axis, 3 = z-axis
  !  - stencil_index: Index of the stencil point (1 to 15) used for high-order finite difference
  !  - k_point: Index of the current k-point (1 to NPT)
  !----------------------------
  DO L = 1, NPT
    ! x-direction stencil: IP(1,1-15,L), IM(1,1-15,L)
    read(01,103) IP(1, 1,L),IM(1, 1,L), IP(1, 2,L),IM(1, 2,L), IP(1, 3,L),IM(1, 3,L), &
               & IP(1, 4,L),IM(1, 4,L), IP(1, 5,L),IM(1, 5,L), IP(1, 6,L),IM(1, 6,L)
    read(02,103) IP(1, 7,L),IM(1, 7,L), IP(1, 8,L),IM(1, 8,L), IP(1, 9,L),IM(1, 9,L), &
               & IP(1,10,L),IM(1,10,L), IP(1,11,L),IM(1,11,L), IP(1,12,L),IM(1,12,L)
    read(03,104) IP(1,13,L),IM(1,13,L), IP(1,14,L),IM(1,14,L), IP(1,15,L),IM(1,15,L), &
               & IEX(L)  ! symmetry-based index remapping for AMA1
    
    ! y-direction stencil: IP(2,1-15,L), IM(2,1-15,L)
    read(11,103) IP(2, 1,L),IM(2, 1,L), IP(2, 2,L),IM(2, 2,L), IP(2, 3,L),IM(2, 3,L), &
               & IP(2, 4,L),IM(2, 4,L), IP(2, 5,L),IM(2, 5,L), IP(2, 6,L),IM(2, 6,L)
    read(12,103) IP(2, 7,L),IM(2, 7,L), IP(2, 8,L),IM(2, 8,L), IP(2, 9,L),IM(2, 9,L), &
               & IP(2,10,L),IM(2,10,L), IP(2,11,L),IM(2,11,L), IP(2,12,L),IM(2,12,L)
    read(13,105) IP(2,13,L),IM(2,13,L), IP(2,14,L),IM(2,14,L), IP(2,15,L),IM(2,15,L)
    
    ! z-direction stencil: IP(3,1-15,L), IM(3,1-15,L)
    read(21,103) IP(3, 1,L),IM(3, 1,L), IP(3, 2,L),IM(3, 2,L), IP(3, 3,L),IM(3, 3,L), &
               & IP(3, 4,L),IM(3, 4,L), IP(3, 5,L),IM(3, 5,L), IP(3, 6,L),IM(3, 6,L)
    read(22,103) IP(3, 7,L),IM(3, 7,L), IP(3, 8,L),IM(3, 8,L), IP(3, 9,L),IM(3, 9,L), &
               & IP(3,10,L),IM(3,10,L), IP(3,11,L),IM(3,11,L), IP(3,12,L),IM(3,12,L)
    read(23,105) IP(3,13,L),IM(3,13,L), IP(3,14,L),IM(3,14,L), IP(3,15,L),IM(3,15,L)
  END DO
103 FORMAT(12I6)  ! Each line contains 12 integer values, read as 6 index pairs
104 FORMAT( 7I6)  ! Each line contains  7 integer values, read as 3 index pairs + IEX
105 FORMAT( 6I6)  ! Each line contains  6 integer values, read as 3 index pairs
  
  !----------------------------
  ! Close all files after reading
  !----------------------------
  CLOSE(UNIT=01); CLOSE(UNIT=02); CLOSE(UNIT=03)
  CLOSE(UNIT=11); CLOSE(UNIT=12); CLOSE(UNIT=13)
  CLOSE(UNIT=21); CLOSE(UNIT=22); CLOSE(UNIT=23)
  
  RETURN
END SUBROUTINE VDATA



!================================================================================
! Subroutine : GV
! Purpose    : Compute group velocity components (VX, VY, VZ) for a given k-point I
!              and band index J using a high-order finite difference method.
!
! Method     :
!   - Applies 15-point central differences along each spatial direction (x, y, z)
!   - Evaluates dE/dk using band eigenvalues from AMA1
!   - Scales results by lattice constant and unit conversion factor CEM (e/h * 10^10)
!     group velocity: [Angstrom/s] -> [m/s]
!
! Inputs :
!   - I : k-point index
!   - J : band index
!
! Outputs :
!   - VX, VY, VZ : Group velocity components
!
! Dependencies :
!   - Uses AMA1 from band_data (eigenvalues)
!   - Uses IP and IM from index_data (stencil indices)
!================================================================================
SUBROUTINE GV(I, J, VX, VY, VZ)
  USE constants  ! Physical and numerical constants: PI, EV, KB, CEM, AL, BL, CL, EF, NPX, NPY, NPZ, NPT, MML, IMAX, NEMAX, EPS
  USE band_data  ! Band energies: AMA1(IMAX,MML), total bands MMIN
  USE index_data ! k-point mapping and stencil indices: IC, N1-N7, IP, IM, IEX
  IMPLICIT NONE
  
  ! Inputs
  INTEGER :: I, J  ! I = k-point index, J = band index
  
  ! Coefficients for 15-point central difference
  REAL(KIND=8) :: CK1, CK2, CK3, CK4, CK5, CK6, CK7, CK8, CK9
  REAL(KIND=8) :: CK10, CK11, CK12, CK13, CK14, CK15
  
  ! Intermediate derivatives along x (AAA), y (AAB), z (AAC)
  REAL(KIND=8) :: AAA, AAB, AAC
  
  ! Output velocities
  REAL(KIND=8) :: VX, VY, VZ
  
  !--------------------------------------------
  ! Define finite difference weights for 15-point stencil
  ! These coefficients improve accuracy of numerical differentiation
  !--------------------------------------------
  CK1  =  1.875000      /  2.0
  CK2  = -1.544118      /  4.0
  CK3  =  1.115196      /  6.0
  CK4  = -0.7043343     /  8.0
  CK5  =  0.3873839     / 10.0
  CK6  = -0.1844685     / 12.0
  CK7  =  7.5464390E-02 / 14.0
  CK8  = -2.6248485E-02 / 16.0
  CK9  =  7.6558078E-03 / 18.0
  CK10 = -1.8373939E-03 / 20.0
  CK11 =  3.5334498E-04 / 22.0
  CK12 = -5.2347405E-05 / 24.0
  CK13 =  5.6086506E-06 / 26.0
  CK14 = -3.8680349E-07 / 28.0
  CK15 =  1.2893449E-08 / 30.0
  
  !--------------------------------------------
  ! Compute derivative dE/dkx using symmetric 15-point difference
  ! AAA stores gradient along x
  !--------------------------------------------
  AAA = (CK1*(AMA1(IP(1, 1,I),J)-AMA1(IM(1, 1,I),J))+ &
      &  CK2*(AMA1(IP(1, 2,I),J)-AMA1(IM(1, 2,I),J))+ &
      &  CK3*(AMA1(IP(1, 3,I),J)-AMA1(IM(1, 3,I),J))+ &
      &  CK4*(AMA1(IP(1, 4,I),J)-AMA1(IM(1, 4,I),J))+ &
      &  CK5*(AMA1(IP(1, 5,I),J)-AMA1(IM(1, 5,I),J))+ &
      &  CK6*(AMA1(IP(1, 6,I),J)-AMA1(IM(1, 6,I),J))+ &
      &  CK7*(AMA1(IP(1, 7,I),J)-AMA1(IM(1, 7,I),J))+ &
      &  CK8*(AMA1(IP(1, 8,I),J)-AMA1(IM(1, 8,I),J))+ &
      &  CK9*(AMA1(IP(1, 9,I),J)-AMA1(IM(1, 9,I),J))+ &
      & CK10*(AMA1(IP(1,10,I),J)-AMA1(IM(1,10,I),J))+ &
      & CK11*(AMA1(IP(1,11,I),J)-AMA1(IM(1,11,I),J))+ &
      & CK12*(AMA1(IP(1,12,I),J)-AMA1(IM(1,12,I),J))+ &
      & CK13*(AMA1(IP(1,13,I),J)-AMA1(IM(1,13,I),J))+ &
      & CK14*(AMA1(IP(1,14,I),J)-AMA1(IM(1,14,I),J))+ &
      & CK15*(AMA1(IP(1,15,I),J)-AMA1(IM(1,15,I),J))  &
      & )*(AL*1.0*NPX)*CEM
  
  !--------------------------------------------
  ! Compute derivative dE/dky -> AAB (y-direction)
  !--------------------------------------------
  AAB = (CK1*(AMA1(IP(2, 1,I),J)-AMA1(IM(2, 1,I),J))+ &
      &  CK2*(AMA1(IP(2, 2,I),J)-AMA1(IM(2, 2,I),J))+ &
      &  CK3*(AMA1(IP(2, 3,I),J)-AMA1(IM(2, 3,I),J))+ &
      &  CK4*(AMA1(IP(2, 4,I),J)-AMA1(IM(2, 4,I),J))+ &
      &  CK5*(AMA1(IP(2, 5,I),J)-AMA1(IM(2, 5,I),J))+ &
      &  CK6*(AMA1(IP(2, 6,I),J)-AMA1(IM(2, 6,I),J))+ &
      &  CK7*(AMA1(IP(2, 7,I),J)-AMA1(IM(2, 7,I),J))+ &
      &  CK8*(AMA1(IP(2, 8,I),J)-AMA1(IM(2, 8,I),J))+ &
      &  CK9*(AMA1(IP(2, 9,I),J)-AMA1(IM(2, 9,I),J))+ &
      & CK10*(AMA1(IP(2,10,I),J)-AMA1(IM(2,10,I),J))+ &
      & CK11*(AMA1(IP(2,11,I),J)-AMA1(IM(2,11,I),J))+ &
      & CK12*(AMA1(IP(2,12,I),J)-AMA1(IM(2,12,I),J))+ &
      & CK13*(AMA1(IP(2,13,I),J)-AMA1(IM(2,13,I),J))+ &
      & CK14*(AMA1(IP(2,14,I),J)-AMA1(IM(2,14,I),J))+ &
      & CK15*(AMA1(IP(2,15,I),J)-AMA1(IM(2,15,I),J))  &
      & )*(BL*1.0*NPY)*CEM
  
  !--------------------------------------------
  ! Compute derivative dE/dkz -> AAC (z-direction)
  !--------------------------------------------
  AAC = (CK1*(AMA1(IP(3, 1,I),J)-AMA1(IM(3, 1,I),J))+ &
      &  CK2*(AMA1(IP(3, 2,I),J)-AMA1(IM(3, 2,I),J))+ &
      &  CK3*(AMA1(IP(3, 3,I),J)-AMA1(IM(3, 3,I),J))+ &
      &  CK4*(AMA1(IP(3, 4,I),J)-AMA1(IM(3, 4,I),J))+ &
      &  CK5*(AMA1(IP(3, 5,I),J)-AMA1(IM(3, 5,I),J))+ &
      &  CK6*(AMA1(IP(3, 6,I),J)-AMA1(IM(3, 6,I),J))+ &
      &  CK7*(AMA1(IP(3, 7,I),J)-AMA1(IM(3, 7,I),J))+ &
      &  CK8*(AMA1(IP(3, 8,I),J)-AMA1(IM(3, 8,I),J))+ &
      &  CK9*(AMA1(IP(3, 9,I),J)-AMA1(IM(3, 9,I),J))+ &
      & CK10*(AMA1(IP(3,10,I),J)-AMA1(IM(3,10,I),J))+ &
      & CK11*(AMA1(IP(3,11,I),J)-AMA1(IM(3,11,I),J))+ &
      & CK12*(AMA1(IP(3,12,I),J)-AMA1(IM(3,12,I),J))+ &
      & CK13*(AMA1(IP(3,13,I),J)-AMA1(IM(3,13,I),J))+ &
      & CK14*(AMA1(IP(3,14,I),J)-AMA1(IM(3,14,I),J))+ &
      & CK15*(AMA1(IP(3,15,I),J)-AMA1(IM(3,15,I),J))  &
      & )*(CL*1.0*NPZ)*CEM
  
  !--------------------------------------------
  ! Combine directional derivatives to get velocity vector
  ! VX, VY, VZ: group velocity components at point I, band J
  ! Averaging helps reduce noise and improve symmetry
  !--------------------------------------------
  VX = (AAB + AAC) / 2.0
  VY = (AAC + AAA) / 2.0
  VZ = (AAA + AAB) / 2.0
  
  RETURN
END SUBROUTINE GV