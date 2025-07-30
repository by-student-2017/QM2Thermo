!-----------------------------------------------------------------------
! Author   : H. Sato and M. Inukai
! Affiliation : [AUE / Riken]
! Contact  : [uichiro25@gmail.com] (U. Mizutani: https://doi.org/10.1007/s11669-024-01103-0)
! GitHub   : https://github.com/by-student-2017/LBT-TETRA (optional)
!-----------------------------------------------------------------------
! Program: GenerateStencil
! Purpose: Constructs finite-difference stencils for k-space derivatives
!          using symmetry operations from WIEN2k output files.
!
! Overview:
!   This program reads symmetry operations and k-point data from WIEN2k
!   output files, applies symmetry transformations to identify equivalent
!   k-points, and builds central finite-difference stencils for use in
!   band structure interpolation or derivative calculations.
!
! Components:
!
! 1. Module: constants
!    - Defines global parameters and allocatable arrays:
!        * Grid dimensions: NPX, NPY, NPZ
!        * Total k-points: NPT
!        * Symmetry count: Nsym
!        * SymOp type: rotation matrix + translation vector
!        * Arrays: IXX, IYY, IZZ (k-point coordinates), WTT (weights),
!                  IC (grid-to-kpoint mapping), LEX (stencil lookup)
!    - Provides Wrap function for periodic boundary conditions
!
! 2. Subroutines:
!    - read_constants:
!        Reads grid size and symmetry count from wien.kgen, wien.klist, wien.struct
!
!    - ReadSymmetry:
!        Loads symmetry operations from wien.struct
!        Each operation includes a 3 * 3 integer rotation matrix and a translation vector
!
!    - ApplySymmetry:
!        Applies a symmetry operation to a k-point (I1, J1, K1)
!        NOTE: In reciprocal space, translation vectors are ignored
!              because k-point equivalence is determined by rotation only.
!
!    - ReadKList:
!        Reads k-point coordinates and weights from wien.klist
!
!    - BuildIC:
!        For each grid point (LX, LY, LZ), finds the nearest symmetry-equivalent
!        k-point from the input list using all symmetry operations.
!        Constructs:
!          * IC(LX,LY,LZ): maps grid point to k-point index
!          * LEX(NP): lookup list of k-point IDs for stencil generation
!
!    - WriteStencils:
!        For each grid point, extracts +/-15 neighbors in x, y, z directions
!        using periodic wrapping, and writes stencil indices to output files:
!          * cfA1.dat - cfA3.dat : x-direction stencils
!          * cfB1.dat - cfB3.dat : y-direction stencils
!          * cfC1.dat - cfC3.dat : z-direction stencils
!
! 3. Main Program: GenerateStencil
!    - Orchestrates the full workflow:
!        a) Load constants and allocate arrays
!        b) Read symmetry operations
!        c) Read k-point list
!        d) Build symmetry-based interpolation map
!        e) Generate and write stencil indices
!
! Input Files:
!   - wien.struct : crystal structure and symmetry operations
!   - wien.klist  : list of k-points and weights
!   - wien.kgen   : number of k-points (IMAX)
!
! Output Files:
!   - cfA1.dat - cfA3.dat : x-direction stencil indices
!   - cfB1.dat - cfB3.dat : y-direction stencil indices
!   - cfC1.dat - cfC3.dat : z-direction stencil indices
!
! Notes:
!   - Translation vectors in symmetry operations are ignored in ApplySymmetry
!     because they do not affect k-point equivalence in reciprocal space.
!   - This design ensures consistency with legacy Fortran77 implementations
!     and avoids rounding errors in stencil matching.
!-----------------------------------------------------------------------
!generate_stencil.f90            <- Main Fortran source file
!+----- constants module         <- Defines grid parameters, symmetry type, and utility functions
!|----+----- SymOp type: rotation matrix + translation vector
!|    +----- Wrap function for periodic index handling
!| 
!+----- read_constants           <- Reads NPX, NPY, NPZ, IMAX, Nsym from WIEN2k files
!+----- ReadSymmetry             <- Loads symmetry operations from wien.struct
!+----- ApplySymmetry            <- Applies rotation (ignoring translation) to k-point coordinates
!+----- ReadKList                <- Loads k-point coordinates and weights from wien.klist
!+----- BuildIC                  <- Maps grid points to symmetry-equivalent k-points
!+----- WriteStencils            <- Generates and writes +/-15-point stencil indices
!| 
!+----- wien.struct              <- Input: crystal structure and symmetry operations
!+----- wien.klist               <- Input: list of k-points and weights
!+----- wien.kgen                <- Input: number of k-points (IMAX)
!| 
!+----- cfA1.dat - cfA3.dat  <- Output: x-direction stencil indices
!+----- cfB1.dat - cfB3.dat  <- Output: y-direction stencil indices
!+----- cfC1.dat - cfC3.dat  <- Output: z-direction stencil indices
!-----------------------------------------------------------------------



!---------------------------------------------------------------------
! Module: constants
! Purpose:
!   Defines global parameters, data structures, and utility functions
!   used throughout stencil generation and symmetry operations.
!
! Contents:
!   - Grid parameters: NPX, NPY, NPZ, NPT
!   - Symmetry count: Nsym
!   - SymOp type: 3 * 3 rotation matrix + translation vector
!   - Arrays:
!       * IXX, IYY, IZZ : k-point coordinates
!       * WTT           : k-point weights
!       * IC            : grid-to-kpoint mapping
!       * LEX           : stencil lookup table
!   - Function:
!       * Wrap : applies periodic boundary conditions to grid indices
!---------------------------------------------------------------------
MODULE constants
  IMPLICIT NONE
  
  INTEGER :: IMAX
  
  INTEGER :: NPX, NPY, NPZ
  INTEGER :: NPT
  !INTEGER, PARAMETER :: NPX   = 58                        ! k-mesh size in X (default: 58)
  !INTEGER, PARAMETER :: NPY   = 58                        ! same for Y
  !INTEGER, PARAMETER :: NPZ   = 58                        ! same for Z
  !INTEGER, PARAMETER :: NPT   = (NPX+1)*(NPY+1)*(NPZ+1)   ! Total number of k-points: Maximum number of rows in cf58A3.dat (4th column) : IEX
  
  INTEGER :: Nsym
  INTEGER :: num_sym_ops
  CHARACTER(LEN=15) :: sym_group_name
  
  INTEGER :: ALX, ALY, ALZ
  INTEGER :: BLX, BLY, BLZ
  INTEGER :: CLX, CLY, CLZ
  
  !-----------------------------------------------------------
  ! User-defined type representing a symmetry operation
  ! - R: 3*3 integer rotation matrix
  ! - T: translation vector (usually fractional)
  !-----------------------------------------------------------
  TYPE :: SymOp
    INTEGER :: R(3,3)
    REAL(8) :: T(3)
  END TYPE
  
  TYPE(SymOp), ALLOCATABLE :: sym(:)
  
  INTEGER, ALLOCATABLE :: IXX(:), IYY(:), IZZ(:)
  REAL(8), ALLOCATABLE :: WTT(:)
  
  INTEGER, ALLOCATABLE :: IC(:,:,:)
  INTEGER, ALLOCATABLE :: LEX(:)
  
CONTAINS
  
  !---------------------------------------------------------------------
  ! Function: Wrap
  ! Purpose: Applies periodic boundary conditions to an integer index
  !---------------------------------------------------------------------
  INTEGER FUNCTION Wrap(idx, KMAX)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx, KMAX
    
    Wrap = idx
    IF (Wrap > KMAX) Wrap = Wrap - KMAX
    IF (Wrap < 0)    Wrap = Wrap + KMAX
  END FUNCTION Wrap
  
END MODULE constants



!---------------------------------------------------------------------
! Subroutine: read_constants
! Purpose:
!   Reads grid dimensions and symmetry count from WIEN2k input files:
!     - wien.kgen   -> IMAX (number of k-points)
!     - wien.klist  -> NPX, NPY, NPZ (grid size)
!     - wien.struct -> Nsym (number of symmetry operations)
!---------------------------------------------------------------------
SUBROUTINE read_constants()
  USE constants   ! IMAX, NPX, NPY, NPZ, NPT, Nsym, sym, sym_ops, ALX-CLZ, IXX, IYY, IZZ, WTT, IC, LEX, Wrap function
  IMPLICIT NONE
  
  CHARACTER(LEN=100) :: line
  INTEGER :: iostat
  INTEGER :: check_Nsym

  OPEN(UNIT=92, FILE='wien.kgen', STATUS='OLD', IOSTAT=iostat)
  READ(92, '(1X, I9)') IMAX
  WRITE(6,'(A,I9)') "   Maximum number of k-points:", IMAX
  CLOSE(92)

  OPEN(UNIT=93, FILE='wien.klist', STATUS='OLD', IOSTAT=iostat)
  READ(93, '(85X, 3(1X,I2))') NPX, NPY, NPZ
  WRITE(6,'(A,I2)') "   k-mesh size X, NPX: ", NPX
  WRITE(6,'(A,I2)') "   k-mesh size Y, NPY: ", NPY
  WRITE(6,'(A,I2)') "   k-mesh size Z, NPZ: ", NPZ
  NPT = (NPX+1)*(NPY+1)*(NPZ+1)
  WRITE(6,'(A,I9)') "   Total number of k-points:", NPT
  CLOSE(93)
  
  ! Search for "NUMBER OF SYMMETRY OPERATIONS" line
  OPEN(UNIT=94, FILE='wien.struct', STATUS='OLD', IOSTAT=iostat)
  READ(94,*)
  READ(94,'(30X, I3, 1X, A)') num_sym_ops, sym_group_name
  WRITE(*,*) "Space group: ", num_sym_ops, " , Symbol: ", sym_group_name
  DO
    READ(94,'(A)',IOSTAT=iostat) line
    IF (iostat /= 0) EXIT
    IF (INDEX(line, 'NUMBER OF SYMMETRY OPERATIONS') > 0) THEN
      READ(line(1:5),*) Nsym
      EXIT
    END IF
  END DO
  WRITE(6,'(A,I3)') "   Number of symmetries:", Nsym
  CLOSE(94)
  
  IF (num_sym_ops == 1 .or. num_sym_ops == 2 .or. &
      num_sym_ops == 47 .or. num_sym_ops == 123 .or. &
      num_sym_ops == 221) THEN
    ! Space group: P (1), P-1 (2), Pmmm (47), P4/mmm (123), SC (221)
    ALX =  1; ALY =  0; ALZ =  0
    BLX =  0; BLY =  1; BLZ =  0
    CLX =  0; CLY =  0; CLZ =  1
  ELSE IF (num_sym_ops == 227 .or. num_sym_ops == 139) THEN
    ! Space group: FCC (high symmetry)
    ALX = -1; ALY =  1; ALZ =  1
    BLX =  1; BLY = -1; BLZ =  1
    CLX =  1; CLY =  1; CLZ = -1
  ELSE IF (num_sym_ops == 229 .or. num_sym_ops == 225) THEN
    ! Space group: BCC (high symmetry)
    ALX =  1; ALY =  1; ALZ =  0
    BLX =  0; BLY =  1; BLZ =  1
    CLX =  1; CLY =  0; CLZ =  1
  ELSE
    WRITE(*,*) "Error: Unsupported space group. Supported groups: 1 (P), 2 (P-1), 47 (Pmmm), &"
    WRITE(*,*) "       123 (P4/mmm), 221 (SC), 227 (FCC), 139 (FCC), 229 (BCC), and 225 (BCC)."
    STOP
  END IF
  WRITE(*,*) ALX, ALY, ALZ
  WRITE(*,*) BLX, BLY, BLZ
  WRITE(*,*) CLX, CLY, CLZ

  RETURN
END SUBROUTINE



!---------------------------------------------------------------------
! Subroutine: ReadSymmetry
! Purpose:
!   Reads symmetry operations from wien.struct.
!   Each operation includes:
!     - 3 * 3 integer rotation matrix
!     - translation vector (usually fractional)
!
! Note:
!   Translation vectors are stored but ignored in ApplySymmetry
!   because they do not affect k-point equivalence in reciprocal space.
!---------------------------------------------------------------------
SUBROUTINE ReadSymmetry(filename)
  USE constants   ! IMAX, NPX, NPY, NPZ, NPT, Nsym, sym, sym_ops, ALX-CLZ, IXX, IYY, IZZ, WTT, IC, LEX, Wrap function
  IMPLICIT NONE
  
  CHARACTER(LEN=*), INTENT(IN) :: filename          ! Path to wien.struct
  INTEGER :: i, j, op_id
  INTEGER :: Rmat(3,3)
  REAL(8) :: Tvec(3)
  INTEGER :: check_Nsym
  
  CHARACTER(LEN=200) :: line
  INTEGER :: ios
  
  OPEN(UNIT=99, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ios)
  IF (ios /= 0) THEN
    PRINT *, 'Error opening file:', filename
    STOP
  END IF
  
  ! Skip until symmetry section (typically after "NUMBER OF SYMMETRY OPERATIONS")
  DO
    READ(99,'(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (INDEX(line, 'NUMBER OF SYMMETRY OPERATIONS') > 0) EXIT
  END DO
  
  ! Read symmetry operations
  DO op_id = 1, Nsym
    DO i = 1,3
      READ(99,'(3I2,1X,F10.8)') rmat(i,1), rmat(i,2), rmat(i,3), tvec(i)
      !WRITE(*,'(3(1X,I2),2X, F10.8)') rmat(i,1), rmat(i,2), rmat(i,3), tvec(i)
    END DO
    READ(99,*) check_Nsym
    sym(op_id)%R = rmat
    sym(op_id)%T = tvec
  END DO
  WRITE(6,'(A,I3)') "   Check number of symmetries:", check_Nsym
  CLOSE(99)
  
END SUBROUTINE ReadSymmetry



!---------------------------------------------------------------------
! Subroutine: ApplySymmetry
! Purpose:
!   Applies a symmetry transformation to a k-point (I1, J1, K1)
!   using the rotation matrix from sym(III).
!
! Output:
!   Transformed coordinates (I4, J4, K4), rounded to nearest integer.
!
! Note:
!   Translation vectors are intentionally ignored because in reciprocal
!   space (k-space), k-point equivalence is determined solely by rotation.
!   This ensures consistency with Fortran77 implementations.
!---------------------------------------------------------------------
SUBROUTINE ApplySymmetry(III, I1, J1, K1, I4, J4, K4)
  USE constants   ! IMAX, NPX, NPY, NPZ, NPT, Nsym, sym, sym_ops, ALX-CLZ, IXX, IYY, IZZ, WTT, IC, LEX, Wrap function
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: III, I1, J1, K1
  INTEGER, INTENT(OUT):: I4, J4, K4

  INTEGER :: i
  REAL(8) :: v_in(3), v_out(3)

  ! Convert input to real-valued vector
  v_in = (/ REAL(I1), REAL(J1), REAL(K1) /)

  ! Apply rotation and translation
  DO i = 1, 3
    !v_out(i) = DOT_PRODUCT(sym(III)%R(i,:), v_in) + sym(III)%T(i)
    v_out(i) = DOT_PRODUCT(sym(III)%R(i,:), v_in)
  END DO

  ! Round transformed coordinates back to integer lattice
  I4 = INT(v_out(1))
  J4 = INT(v_out(2))
  K4 = INT(v_out(3))
END SUBROUTINE ApplySymmetry



!---------------------------------------------------------------------
! Subroutine: ReadKList
! Purpose:
!   Reads k-point coordinates and weights from wien.klist.
!   Populates arrays: IXX, IYY, IZZ, WTT
!---------------------------------------------------------------------
SUBROUTINE ReadKList(kfile)
  USE constants   ! IMAX, NPX, NPY, NPZ, NPT, Nsym, sym, sym_ops, ALX-CLZ, IXX, IYY, IZZ, WTT, IC, LEX, Wrap function
  IMPLICIT NONE
  
  CHARACTER(LEN=*), INTENT(IN) :: kfile
  
  CHARACTER(3) :: ENDC
  CHARACTER(10) :: num
  INTEGER :: i, iw, KB, ios
  
  OPEN(1, FILE=kfile, STATUS='OLD', IOSTAT=ios)
  IF (ios /= 0) STOP 'Error opening wien.klist'
  
  DO KB = 1, IMAX
    READ(1,'(A10,4I10,F5.2)') num, IXX(KB), IYY(KB), IZZ(KB), iw, WTT(KB)
    !WRITE(*,'(A10,4I10,F5.2)') num, IXX(KB), IYY(KB), IZZ(KB), iw, WTT(KB)
  END DO
  CLOSE(1)
END SUBROUTINE ReadKList



!---------------------------------------------------------------------
! Subroutine: BuildIC
! Purpose:
!   For each grid point (LX, LY, LZ), finds the nearest symmetry-equivalent
!   k-point from the input list using all symmetry operations.
!
! Constructs:
!   - IC(LX,LY,LZ): maps grid point to k-point index
!   - LEX(NP): lookup list of k-point IDs for stencil generation
!
! Note:
!   Uses a 5 * 5 * 5 search window around transformed coordinates to
!   account for symmetry-induced shifts and periodicity.
!
! Limitations and Extension Considerations:
!   - This implementation assumes that symmetry-transformed coordinates
!     remain integer-valued, which is valid only for specific space groups
!     (e.g., SC, FCC, BCC) with purely rotational operations.
!
!   - Space groups involving fractional translation vectors (e.g., 1/2, 1/4)
!     or non-orthogonal axes may produce non-integer coordinates after
!     symmetry application, making integer-based matching unreliable.
!
!   - To extend this method to general space groups:
!       * Consider rounding transformed coordinates (e.g., NINT) with care.
!       * Introduce tolerance-based matching (e.g., within +/- 0.5) to absorb
!         floating-point errors.
!       * Redesign IC() to support real-valued indexing or use hash-based
!         lookup structures.
!       * Extract and apply full symmetry operations (rotation + shift)
!         from wien.struct, and validate equivalence in real or reciprocal space.
!
!   - Without these extensions, the current method is limited to high-symmetry
!     lattices where transformed points fall back onto the integer grid.
!
! Note:
!     221, Pm-3m,  SC , no shift
!     225, Fm-3m,  FCC, High symmetry and can be handled by coordinate transformation
!     229, Im-3m,  BCC, Consists of rotating coordinates only
!     139, I4/mmm, BCC, no shift
!     123, P4/mmm, SC , Stable at integer coordinates
!      47, Pmmm,   SC , no shift
!       2, P-1,    Triclinic, no shift
!       1, P,      ALL
! Note:

!      SC(  1): LX = L2X, LY = L2Y, LZ = L2Z
!     FCC(227): LX = -L2X + L2Y + L2Z, LY =  L2X - L2Y + L2Z, LZ =  L2X + L2Y - L2Z
!     BCC(229): LX = L2X + L2Y, LY = L2Y + L2Z, LZ = L2Z + L2X
!---------------------------------------------------------------------
SUBROUTINE BuildIC()
  USE constants   ! IMAX, NPX, NPY, NPZ, NPT, Nsym, sym, sym_ops, ALX-CLZ, IXX, IYY, IZZ, WTT, IC, LEX, Wrap function
  IMPLICIT NONE
  
  INTEGER :: L2X, L2Y, L2Z, LX, LY, LZ, I3, J3, K3
  INTEGER :: I1, J1, K1, I4, J4, K4
  INTEGER :: IG1, IG2, IG3, III, L1, NP
  
  NP = 0
  
  DO L2Z = 0, NPZ
    WRITE(*,*) "STEP:", L2Z, "/", NPZ
    DO L2Y = 0, NPY
      DO L2X = 0, NPX
        NP = NP + 1
        
        !---------------------------------------------------------------------
        ! Convert grid coordinates (L2X, L2Y, L2Z) to Fortran77-style lattice coordinates (LX, LY, LZ)
        !
        ! Purpose:
        !   This transformation redefines the grid point coordinates to align with
        !   a symmetry-adapted lattice basis used in the original Fortran77 implementation.
        !
        !   By applying this coordinate change, symmetry operations (especially rotations)
        !   can be performed without requiring explicit translation vectors (shifts),
        !   and the transformed coordinates remain integer-valued.
        !
        !   This approach avoids floating-point rounding errors and simplifies
        !   the matching of symmetry-equivalent k-points.
        !
        ! Transformation (e.g., FCC):
        !   - LX = -L2X + L2Y + L2Z
        !   - LY =  L2X - L2Y + L2Z
        !   - LZ =  L2X + L2Y - L2Z
        !
        !   These expressions correspond to a change of basis that reflects
        !   the symmetry of face-centered cubic (FCC) or similar lattices.
        !   They ensure that symmetry-transformed points fall back onto
        !   the integer grid, enabling exact matching.
        !---------------------------------------------------------------------
        ! Convert grid coordinates to Fortran77-style lattice coordinates
        !LX = -L2X + L2Y + L2Z
        !LY =  L2X - L2Y + L2Z
        !LZ =  L2X + L2Y - L2Z
        ! Calculate with a matrix corresponding to the space group
        LX = ALX*L2X + ALY*L2Y + ALZ*L2Z
        LY = BLX*L2X + BLY*L2Y + BLZ*L2Z
        LZ = CLX*L2X + CLY*L2Y + CLZ*L2Z
        
        !---------------------------------------------------------------
        ! Loop over all k-points defined in the input list (1 to IMAX)
        !
        ! For each index L1:
        !   - Retrieve the corresponding k-point coordinates (I1, J1, K1)
        !     from the arrays IXX, IYY, and IZZ.
        !   - These coordinates represent a point in reciprocal space
        !     (typically on a regular grid), and will be evaluated
        !     under symmetry operations in later steps.
        !
        ! Purpose:
        !   - To access each original k-point from the list for symmetry
        !     transformation and stencil mapping.
        !   - These values are passed to ApplySymmetry() in the next block.
        !
        ! Note:
        !   - This loop forms the outer sweep that attempts to find
        !     a symmetry-transformed match to the target grid coordinates
        !     defined earlier by (LX, LY, LZ).
        !---------------------------------------------------------------
        DO L1 = 1, IMAX
          I1 = IXX(L1)
          J1 = IYY(L1)
          K1 = IZZ(L1)
          
          DO III = 1, Nsym
            !--------------------------------------------------------------------------
            ! Applies the III-th symmetry operation from the array 'sym'
            ! to the lattice point (I1, J1, K1), and returns the transformed
            ! coordinates (I4, J4, K4) as integers.
            ! This transformation uses the rotation matrix and translation vector
            ! defined in sym(III), consistent with the symmetry rules from wien.struct.
            !--------------------------------------------------------------------------
            CALL ApplySymmetry(III, I1, J1, K1, I4, J4, K4)
            
            !---------------------------------------------------------------------------
            ! Attempt to match the symmetry-transformed point (I4, J4, K4) to the current
            ! grid coordinate (LX, LY, LZ) using a localized 3D search window.
            !
            ! The nested DO loops iterate over IG1, IG2, IG3 in [-2, 2], generating a set
            ! of candidate grid points (I3, J3, K3) around the transformed point.
            !
            ! These candidate points are computed by applying weighted combinations of
            ! IG1-IG3, scaled by NPX, NPY, NPZ respectively. This effectively spans
            ! symmetry-related lattice translations without explicitly using shift vectors.
            !
            ! Conceptually:
            !   - The transformation (I1 -> I4) applies a symmetry rotation.
            !   - The search (I4 -> I3) explores nearby grid points that may match the
            !     current grid location (LX, LY, LZ) under symmetry.
            !
            ! If a candidate point (I3, J3, K3) matches the current grid location,
            ! it means that the grid point (LX, LY, LZ) is symmetry-equivalent to
            ! the original k-point (I1, J1, K1), and the mapping is recorded:
            !   - IC(LX, LY, LZ) <- L1 : index of matching k-point
            !   - LEX(NP)        <- L1 : stencil reference for this grid point
            !
            ! Upon successful match, EXIT terminates the innermost loop early.
            ! This avoids redundant checks once a valid symmetry match is found.
            !---------------------------------------------------------------------------
            DO IG1 = -2, 2
              DO IG2 = -2, 2
                DO IG3 = -2, 2
                  !I3 = I4 - (-IG1 + IG2 + IG3) * NPX
                  !J3 = J4 - ( IG1 - IG2 + IG3) * NPY
                  !K3 = K4 - ( IG1 + IG2 - IG3) * NPZ
                  ! Calculate with a matrix corresponding to the space group
                  I3 = I4 - (ALX*IG1 + ALY*IG2 + ALZ*IG3) * NPX
                  J3 = J4 - (BLX*IG1 + BLY*IG2 + BLZ*IG3) * NPY
                  K3 = K4 - (CLX*IG1 + CLY*IG2 + CLZ*IG3) * NPZ
                  
                  IF (LX == I3 .AND. LY == J3 .AND. LZ == K3) THEN
                    ! Optional debug: WRITE(6,'(12I6)') NP, L2X, L2Y, L2Z, LX, LY, LZ, I1, J1, K1, L1, III
                    IC(L2X,L2Y,L2Z) = L1     ! Link this grid point to k-point L1
                    LEX(NP) = L1             ! Record stencil reference for NP-th grid node
                    GOTO 200                 ! Since the target was found, proceed to the next LX, LY, LZ loop.
                  END IF
                END DO
              END DO
            END DO
            !
          END DO
        END DO
        !
  200 CONTINUE
      END DO
    END DO
  END DO
END SUBROUTINE BuildIC



!---------------------------------------------------------------------
! Subroutine: WriteStencils
! Purpose:
!   For each grid point, extracts +/-15 neighbors in x, y, z directions
!   using periodic wrapping, and writes stencil indices to output files.
!
! Output Files:
!   - cfA1.dat - cfA3.dat : x-direction stencils
!   - cfB1.dat - cfB3.dat : y-direction stencils
!   - cfC1.dat - cfC3.dat : z-direction stencils
!
! Application:
!   Used for central finite-difference approximations of band derivatives.
!---------------------------------------------------------------------
SUBROUTINE WriteStencils()
  USE constants   ! IMAX, NPX, NPY, NPZ, NPT, Nsym, sym, sym_ops, ALX-CLZ, IXX, IYY, IZZ, WTT, IC, LEX, Wrap function
  IMPLICIT NONE
  
  INTEGER :: L3X, L3Y, L3Z, I, J, K, NP, I2, ID
  INTEGER :: LA(30), LB(30), LC(30)
  
  ! Open stencil output files for each direction
  OPEN(01,FILE='cfA1.dat'); OPEN(02,FILE='cfA2.dat'); OPEN(03,FILE='cfA3.dat')
  OPEN(11,FILE='cfB1.dat'); OPEN(12,FILE='cfB2.dat'); OPEN(13,FILE='cfB3.dat')
  OPEN(21,FILE='cfC1.dat'); OPEN(22,FILE='cfC2.dat'); OPEN(23,FILE='cfC3.dat')
  
  NP = 0
  
  DO L3Z = 0, NPZ
    DO L3Y = 0, NPY
      DO L3X = 0, NPX
        I = L3X
        J = L3Y
        K = L3Z
        NP = NP + 1
        
        I2 = LEX(NP)
        
        !---------------------------------------------------------------------------
        ! For each grid point (I, J, K) in reciprocal space:
        !   - Collect 30 stencil indices in each direction (x, y, z)
        !   - The stencil uses 15 points forward and 15 points backward
        !     relative to the current location, forming a centered window
        !
        ! This layout is ideal for constructing central finite difference schemes
        ! where symmetric neighbors are required on both sides of the evaluation point.
        !
        ! Specifically:
        !   - LA(:) stores neighbors along the x-direction
        !   - LB(:) stores neighbors along the y-direction
        !   - LC(:) stores neighbors along the z-direction
        !
        ! Each point is located by offsetting the central coordinate
        !   - +/-ID steps along one direction (ID = 1 to 15)
        !   - Wrap() ensures proper periodic boundary conditions
        !
        ! For example:
        !   LA(ID*2-1) <- forward  point at (L3X + ID, J, K)
        !   LA(ID*2)   <- backward point at (L3X - ID, J, K)
        !
        ! These indices are used later to evaluate dE/dkx, dE/dky, dE/dkz
        ! using central difference formulas such as:
        !   (E(k + dk) - E(k - dk)) / (2dk)
        !---------------------------------------------------------------------------
        ! Collect +/-15 indices in each direction using periodic wrapping
        DO ID = 1, 15
          LA(ID*2-1) = IC(Wrap(L3X + ID, NPX), J, K)   ! IP series
          LA(ID*2)   = IC(Wrap(L3X - ID, NPX), J, K)   ! IM series
          LB(ID*2-1) = IC(I, Wrap(L3Y + ID, NPY), K)   ! IP series
          LB(ID*2)   = IC(I, Wrap(L3Y - ID, NPY), K)   ! IM series
          LC(ID*2-1) = IC(I, J, Wrap(L3Z + ID, NPZ))   ! IP series
          LC(ID*2)   = IC(I, J, Wrap(L3Z - ID, NPZ))   ! IM series
        END DO
        
        ! Write stencils to respective files (each contains 12/12/6+1 entries)
        WRITE(01,'(12I6)') LA( 1:12)
        WRITE(02,'(12I6)') LA(13:24)
        WRITE(03,'( 7I6)') LA(25:30), I2
        
        WRITE(11,'(12I6)') LB( 1:12)
        WRITE(12,'(12I6)') LB(13:24)
        WRITE(13,'( 6I6)') LB(25:30)
        
        WRITE(21,'(12I6)') LC( 1:12)
        WRITE(22,'(12I6)') LC(13:24)
        WRITE(23,'( 6I6)') LC(25:30)
      END DO
    END DO
  END DO
  
  ! Close all output units
  CLOSE(01); CLOSE(02); CLOSE(03)
  CLOSE(11); CLOSE(12); CLOSE(13)
  CLOSE(21); CLOSE(22); CLOSE(23)
END SUBROUTINE WriteStencils



!---------------------------------------------------------------------
! Program: GenerateStencil
! Purpose:
!   Orchestrates stencil construction using symmetry operations
!   and k-point data from WIEN2k output files.
!
! Workflow:
!   1. Read grid and symmetry parameters
!   2. Load symmetry operations
!   3. Load k-point list
!   4. Build symmetry-based interpolation map
!   5. Generate and write stencil indices
!---------------------------------------------------------------------
PROGRAM GenerateStencil
  USE constants   ! IMAX, NPX, NPY, NPZ, NPT, Nsym, sym, sym_ops, ALX-CLZ, IXX, IYY, IZZ, WTT, IC, LEX, Wrap function
  IMPLICIT NONE
  
  !---------------------------------------------
  ! Step 1: Get data from WIEN2k output files and set arrays
  !---------------------------------------------
  WRITE(*,*) "Step 1: Get data from WIEN2k output files and set arrays"
  CALL read_constants
  ALLOCATE(sym(Nsym))  ! Allocate array according to actual number
  ALLOCATE(IXX(IMAX), IYY(IMAX), IZZ(IMAX), WTT(IMAX))
  ALLOCATE(IC(0:NPX,0:NPY,0:NPZ))
  ALLOCATE(LEX(NPT))
  
  !---------------------------------------------
  ! Step 2: Load symmetry operations from struct file
  !---------------------------------------------
  WRITE(*,*) "Step 2: Load symmetry operations from struct file"
  CALL ReadSymmetry('wien.struct')
  
  !---------------------------------------------
  ! Step 3: Load k-point list and weights
  !---------------------------------------------
  WRITE(*,*) "Step 3: Load k-point list and weights"
  CALL ReadKList('wien.klist')
  
  !---------------------------------------------
  ! Step 4: Build interpolation map using symmetry
  !---------------------------------------------
  WRITE(*,*) "Step 4: Build interpolation map using symmetry"
  CALL BuildIC()
  
  !---------------------------------------------
  ! Step 5: Write stencil indices to output files
  !---------------------------------------------
  WRITE(*,*) "Step 5: Write stencil indices to output files"
  CALL WriteStencils()
  
  WRITE(*,*) "Stencil generation successfully completed."
  
  DEALLOCATE(sym)
  DEALLOCATE(IXX, IYY, IZZ, WTT)
  DEALLOCATE(IC)
  DEALLOCATE(LEX)
  
END PROGRAM GenerateStencil
