!-----------------------------------------------------------------------
! Author   : H. Sato and M. Inukai
! Affiliation : [AUE / Riken]
! Contact  : [your.email@domain.edu]
! GitHub   : https://github.com/yourusername/seebeck_analysis (optional)
!-----------------------------------------------------------------------
! Program      : chemical_potential.f90
! Purpose      : Calculate temperature-dependent chemical potential mu(T)
!                by solving the carrier conservation condition using
!                Fermi-Dirac distribution and numerical DOS integration.
!
! Description  : For each temperature, applies the bisection method
!                to find chemical potential mu such that total carrier
!                density matches the target value derived from DOS(E).
!
! Method       :
!   - Input: DOS(E) from 'AKK.DATA'
!   - Solve: N(T) = integral DOS(E) * f(E,mu,T) dE = constant
!   - Algorithm: Bisection root-finding for mu at each T
!
! Physical constants:
!   TEM  = Thermal energy at 300 [K] = k_B * 300 = 0.0258199727 [eV]
!   TEM1 = Temperature [K] (swept from high to low)
!
! Key dates    : Initial version - 2022.09.16, 2022.09.27
!                Final update    - 2024.12.31
!                Comment revision - 2025.07.13
!
! Compilation:
!   gfortran -O2 chemical_potential.f90 -o chemical_potential.exe
!   ifort    -O2 chemical_potential.f90 -o chemical_potential.exe
!
! Required input files:
!   - AKK.DATA     : Energy grid, density of states, velocity^2 etc.
!                    Energy [eV], DOS [states/eV], velocity [m/s]
!
! Output:
!   - apot.data    : Table of mu(T): temperature vs. chemical potential [eV]
!
! Extension-ready:
!   - Can be linked with transport models (Seebeck, tau(E))
!   - DOS(E) can be replaced with ab initio data
!-----------------------------------------------------------------------
! Possible extensions:
!   - Integrate with transport model codes to compute Seebeck coefficient S(T)
!     using dmu/dT from the computed mu(T) values.
!
!   - Use mu(T) in Boltzmann transport calculations or tau(E,T) modeling to evaluate
!     carrier mobility and conductivity.
!
!   - Replace AKK.DATA with ab initio generated data (e.g. from DFT or Wannier-based methods)
!     for material-specific analysis.
!
!   - Expand to handle spin or band-index resolved chemical potentials (e.g. mu_up(T), mu_down(T))
!     for magnetic or multiband systems.
!
!   - Visualize mu(T) with an external plotting tool (Python, gnuplot) to identify trends
!     and anomalies across temperatures.
!
!   - Perform automatic fitting or spline interpolation to enable derivative analysis
!     such as dmu/dT for thermoelectric response functions.
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Module    : carrier_data
! Purpose   : Define and share energy-resolved physical data arrays
!             across subroutines used for chemical potential evaluation.
!
! Contents  :
!   - EE(N)      : Energy grid [eV]
!   - DOS(N)     : Density of states (DOS) [states/eV]
!   - V2(N)      : Velocity squared [arbitrary units]
!   - VD(N)      : Scattering-related quantity (used in transport)
!   - EN(N)      : Additional energy-resolved parameter
!   - Auxiliary arrays (V2A, VDA, DOSA) for alternate models
!   - TEM        : Thermal energy (k_B T) at current temperature [eV]
!   - DEF        : Energy shift offset (band alignment etc.)
!
! Usage    :
!   - Used by BS and DAT subroutines
!   - Shared by chemical_potential main program
!
! Input provenance:
!   Arrays EE(:), DOS(:), V2(:), VD(:), EN(:) are derived from the 4th-column eigenvalues in wien(Si58).energy
!   AKK.DATA reflects post-processing from raw band structure across 4735 k-points and multiple bands
!
! Extension-ready:
!   - Suitable for linking with transport calculations (Seebeck, tau(E))
!   - Can incorporate band-structure derived quantities or interpolated grids
!-----------------------------------------------------------------------
MODULE carrier_data
  IMPLICIT NONE
  
  ! Parameters
  !-----------------------------------------------------------------------
  ! Read energy-resolved data from AKK.DATA
  ! Number of energy points (IMAX) is set to 60000 in band generation routine
  ! This corresponds to NE=60001 used in the original FORTRAN code for DOS output
  ! Energy range: ED = -1.0 [eV] to EU = +1.0 [eV], converted from atomic units
  ! Energy step: DE = (EU - ED) / (NE - 1) / EV
  !-----------------------------------------------------------------------
  INTEGER, PARAMETER :: N     = 261000     ! Maximum array size (safe upper bound)
  INTEGER, PARAMETER :: IMAX  = 60000      ! Actual number of data points to load from AKK.DATA

  ! Arrays for energy grid and transport properties
  REAL(8) :: EE(N)      ! Energy grid (eV)
  REAL(8) :: V2(N)      ! Velocity squared
  REAL(8) :: VD(N)      ! Scattering rate component
  REAL(8) :: DOS(N)     ! Density of states
  REAL(8) :: EN(N)      ! Energy-dependent transport parameter

  ! Auxiliary arrays for alternative transport models
  REAL(8) :: V2A(N)     ! Temporary velocity^2 array (alternative model / workspace)
  REAL(8) :: VDA(N)     ! Temporary scattering array (for alternate analysis)
  REAL(8) :: DOSA(N)    ! Temporary DOS array (backup or modified DOS)


  REAL(8), PARAMETER :: KB = 8.61733D-5   ! Boltzmann constant [eV/K]
  REAL(8), PARAMETER :: EV = 13.605698D0  ! Conversion constant: Hartree -> eV
  REAL(8), PARAMETER :: PI = 3.141592654D0     ! Mathematical constant Pi
  
  REAL(8) :: TEM        ! Thermal energy [eV], scaled from TEM1
  REAL(8) :: DEF        ! Energy offset (shift applied to EE)
  
  REAL(8) :: VEC        ! The valence electron concentration (VEC)
END MODULE carrier_data



!-----------------------------------------------------------------------
! Program      : chemical_potential
! Purpose      : Find mu(T) via carrier density matching from DOS
!
! Description  :
!   This program computes the chemical potential mu(T)
!   as a function of temperature by solving the carrier
!   conservation condition numerically.
!
!   The total carrier density is calculated via integration
!   of the density of states (DOS(E)) weighted by the
!   Fermi-Dirac distribution at each temperature.
!
!   A bisection algorithm is applied to find the chemical
!   potential mu that matches the target carrier density.
!
!   Input : AKK.DATA -> Contains energy grid, DOS(E), and related quantities
!   Output: apot.data -> Table of mu(T) [eV] versus temperature T [K]
!
!   This code is modularized with carrier_data module and is
!   extendable toward transport coefficient evaluations.
!-----------------------------------------------------------------------
PROGRAM chemical_potential
  USE carrier_data  ! Get physical data arrays and constants from a module
  IMPLICIT NONE
  
  ! ----------------------
  ! Fundamental constants
  ! ----------------------
  INTEGER, PARAMETER :: NT = 25          ! Number of temperature values

  ! ----------------------
  ! Variable declarations
  ! ----------------------
  REAL(8) :: TT(NT)     ! Array of temperature values [K]
  REAL(8) :: TEM1       ! Temperature at current loop [K]
  REAL(8) :: DEFF       ! Shift in boundary due to offset (placeholder)
  
  REAL(8) :: X1, X2     ! Bisection bracket values for chemical potential
  REAL(8) :: X10, X20   ! Initial guess boundaries
  REAL(8) :: CP         ! Chemical potential (mu) at each temperature
  REAL(8) :: F1, F2     ! Function values at boundaries
  REAL(8) :: FC         ! Current function value in bisection
  REAL(8) :: XC         ! Current chemical potential candidate
  REAL(8) :: YY         ! Product test for bisection convergence
  INTEGER :: I          ! Iteration counter (used inside bisection loop)
  INTEGER :: L          ! Temperature loop index

  ! ----------------------
  ! Check input file existence
  ! ----------------------
  LOGICAL :: EX  ! For file existence check using INQUIRE

  ! ----------------------
  ! Temperature points from high to low
  ! (Temperature scan from 1000 K to 5 K)
  ! ----------------------
  DATA TT /1000D0,950D0,900D0,850D0,800D0,750D0,700D0,650D0,600D0,550D0, &
           500D0,450D0,400D0,350D0,300D0,250D0,200D0,150D0,100D0,50D0, &
           40D0,30D0,20D0,10D0,5D0/

  ! ===========================================================
  ! Step 1: Read input DOS data from AKK.DATA
  !         -> EE(:), DOS(:), V2(:), etc. are populated
  ! ===========================================================
  CALL DAT

  ! ===========================================================
  ! Step 2: Open output file for mu(T) results
  ! ===========================================================
  INQUIRE(FILE='AKK.DATA', EXIST=EX)
  IF (.NOT. EX) THEN
    WRITE(*,*) 'ERROR: AKK.DATA not found.'
    STOP
  END IF
  OPEN(UNIT=4, FILE='apot.data')   ! Format: T [K], mu(T) [eV]

  ! ===========================================================
  ! Step 3: Main loop over temperature values
  !         -> For each T, solve for mu via bisection
  ! ===========================================================
  DO L = 1, NT
     DEFF = DEF * 1.0
     !IF (L >= 11) THEN
     !   DEFF = 0.0D0
     !END IF

     TEM1 = TT(L)      ! Set current temperature [K]
     TEM  = KB * TEM1  ! Thermal energy [eV], scaled from TEM1: (Set as a global variable)
     CP   = 0.0D0

     ! Set initial bracket for chemical potential [eV]
     X10 = -8.00D0
     X20 =  8.00D0 + DEFF
     X1  = X10
     X2  = X20

     WRITE(*,'(A,F8.1)') "temp(K) = ", TEM1

     ! -------------------------------------------------------
     ! Step 3.1: Evaluate function at bracket edges
     !          -> Function computes carrier density deviation
     ! -------------------------------------------------------
     CALL BS(1, F1, X1)
     CALL BS(2, F2, X2)
     WRITE(*,'(A,F9.5,2X,A,F9.5)') "F1 = ", F1, "F2 = ", F2

     ! -------------------------------------------------------
     ! Step 3.2: Begin bisection iterations
     ! -------------------------------------------------------
     I = 0
     DO
        XC = 0.5D0 * (X1 + X2)  ! Midpoint for new test chemical potential
        CALL BS(I, FC, XC)      ! Evaluate error function at midpoint
        !WRITE(*,*) 'TEMP:', TEM, 'F1:', F1, 'F2:', F2, 'FC:', FC, 'XC:', XC, 'X1:', X1, 'X2:', X2
        I = I + 1

        ! Check for convergence: small enough error AND sufficient iterations
        IF (I >= 100 .AND. ABS(FC) < 1.0D-15) THEN
        !IF ((I >= 1000 .AND. ABS(FC) < 1.0D-2) .OR. (I >= 10 .AND. ABS(FC) < 1.0D-15)) THEN
           WRITE(*,'(A)') "finale"
           WRITE(*,'(F10.1,1X,F14.8)') TEM1, XC  ! Write result to console
           WRITE(4,'(F10.1,1X,F14.8)') TEM1, XC  ! Write result to output file
           EXIT
        END IF

        ! Bisection rule: if sign changes between F1 and FC, move upper bound
        YY = F1 * FC
        IF (YY < 0.0D0) THEN
           WRITE(*,'(I4,1X,E16.8,1X,F14.8)') I, FC, XC
           X2 = XC
           CYCLE
        END IF

        ! Otherwise, if FC and F2 change sign, move lower bound
        YY = FC * F2
        IF (YY <= 0.0D0) THEN
           WRITE(*,'(I4,1X,E16.8,1X,F14.8)') I, FC, XC
           X1 = XC
           CYCLE
        END IF

        ! If no sign change, convergence may be ambiguous but we accept
        !WRITE(*,'(A,I3,2X,E16.8)') "Converged at I=", I, YY
        !EXIT
     END DO
  END DO
  WRITE(*,*) "mu(T) results written to: apot.data"
  WRITE(*,*) "Output saved in apot.data"

  ! ===========================================================
  ! Step 4: Close output file
  ! ===========================================================
  WRITE(4, *) "! Energy shift offset (EF = 0 -> EF = DEF) [eV]:", DEF
  WRITE(4, *) "! The valence electron concentration (VEC):", VEC
  CLOSE(4)
END PROGRAM chemical_potential



!-----------------------------------------------------------------------
! Subroutine : BS
! Purpose    : Evaluate carrier occupation difference based on
!              Fermi-Dirac distribution and DOS(E), for given mu.
!
! Inputs  :
!   L   : Identifier or label (optional; currently unused for branching)
!   CP  : Candidate chemical potential [eV]
!
! Outputs :
!   T   : Total carrier density above valence band edge (E > 0)
!         minus valence baseline (E =< 0)
!
! Method :
!   - Integrate DOS(E) * f_FD(E, mu, T) over E
!   - Subtract baseline carrier population below E = 0
!-----------------------------------------------------------------------
SUBROUTINE BS(L, T, CP)
  USE carrier_data      ! Get physical data arrays and constants from a module
  IMPLICIT NONE
  INTEGER :: I, L       ! I = energy loop index, L = label (not used here)
  REAL(8) :: E, DE      ! Energy and energy interval
  REAL(8) :: FD1        ! Fermi-Dirac distribution value
  REAL(8) :: F0, FF     ! Temporary variables for DOS contributions
  REAL(8) :: T, T0      ! Integrated carrier densities
  REAL(8) :: CP         ! Candidate chemical potential [eV]
  
  ! Trapezoidal rule or Simpson's Rule
  REAL(8) :: FD1_prev   ! Fermi-Dirac distribution value at the previous step
                        ! Used in Trapezoidal Rule to compute the area based on two points.
  
  ! Variables specifically for Simpson's Rule
  REAL(8) :: E0          ! Energy value at the previous point (left side of interval)
  REAL(8) :: E1          ! Energy value at the current point (middle of interval)
  REAL(8) :: E2          ! Energy value at the next point (right side of interval)
  REAL(8) :: F1          ! Fermi-Dirac weighted density of states at the middle point (E1)
  REAL(8) :: F2          ! Fermi-Dirac weighted density of states at the next point (E2)
  REAL(8) :: FD0         ! Fermi-Dirac distribution value at E0 (previous point)
  REAL(8) :: FD2         ! Fermi-Dirac distribution value at E2 (next point)
  REAL(8) :: FD1_next    ! Precomputed Fermi-Dirac value for the next step (used to optimize Simpson's Rule)
  INTEGER :: DIMAX       ! Adjustment variable for the maximum index (IMAX)
                         ! Ensures the loop processes only even-numbered intervals,
                         ! which is necessary for Simpson's Rule.
                         ! DIMAX = 1 if IMAX is even, or 2 if IMAX is odd.

  ! Initialization
  T  = 0.0D0
  T0 = 0.0D0
  FD1_prev = 0.0D0
  FD1_next = 0.0D0

  DO I = 1, IMAX
     E  = EE(I)
     DE = EE(I+1) - EE(I)

     FD1 = 1.0D0 / (1.0D0 + DEXP((E - CP) / TEM))   ! Fermi-Dirac probability

     !----------------------------------------
     ! Total weighted contribution at energy E (Riemann Sum)
     !FF = DOS(I) * FD1
     !T  = T + DE * FF
     !! Valence band contribution (E =< 0)
     !IF (E <= 0.0D0) THEN
     !   F0 = DOS(I)
     !   T0 = T0 + DE * F0
     !END IF
     !----------------------------------------
     
     !----------------------------------------
     ! Trapezoidal rule
     !IF (I > 1) THEN
     !  FD1_prev = 1.0D0 / (1.0D0 + DEXP((EE(I-1) - CP) / TEM))  ! Precomputation
     !  T = T + 0.5D0 * (DOS(I-1) * FD1_prev + DOS(I) * FD1) * DE
     !END IF
     !! Valence band contribution (E =< 0)
     !IF (E <= 0.0D0) THEN
     !   T0 = T0 + 0.5D0 * (DOS(I-1) + DOS(I)) * DE
     !END IF
     !----------------------------------------
     
     !----------------------------------------
     ! Simpson's Rule
     IF (I > 2 .AND. I < IMAX - 1 .AND. MOD(I, 2) == 0) THEN
       FD1_prev = 1.0D0 / (1.0D0 + DEXP((EE(I-1) - CP) / TEM))  ! Precomputation
       FD1_next = 1.0D0 / (1.0D0 + DEXP((EE(I+1) - CP) / TEM))  ! Precomputation
       DE = (EE(I+1) - EE(I-1)) / 2.0D0   ! Calculating section width
       T  = T  + DE / 3.0D0 * (DOS(I-1) * FD1_prev + 4.0D0 * DOS(I) * FD1 + DOS(I+1) * FD1_next)
       ! Valence band contribution (E =< 0)
       IF (E <= 0.0D0) THEN
         T0 = T0 + DE / 3.0D0 * (DOS(I-1) + 4.0D0 * DOS(I) + DOS(I+1))
       END IF
     END IF
     !----------------------------------------
  END DO

  ! Net carrier density relative to valence base
  T = T - T0
END SUBROUTINE BS



!-----------------------------------------------------------------------
! Subroutine : DAT
! Purpose    : Read energy-resolved band data from input file 'AKK.DATA'
!              including DOS(E), velocity squared V2(E), scattering VD(E)
!
! Notes      :
!   - Applies energy offset DEF (if nonzero) to EE(:)
!   - Populates shared arrays: EE(:), DOS(:), V2(:), VD(:), EN(:)
!   - Can be extended for alternate data handling (DOSA, etc.)
!   - V2A(:), VDA(:), DOSA(:) are defined in module for future use.
!   - Current routine does not populate or use them.
!
! Source note:
!   AKK.DATA is constructed from wien(Si58).energy, containing band eigenvalues.
!   Typically, the 4th column (or designated band index) across 4735 k-points is extracted.
!   These eigenvalues are post-processed into AMA1(:,:), then transformed into AMA(:,:)
!   Energy values E = AMA(I,J)/EV are converted into EE(:) with offset DEF as needed.
!   This serves as the input base for DOS(E), velocity squared, scattering, etc.
!-----------------------------------------------------------------------
SUBROUTINE DAT
  USE carrier_data        ! Get physical data arrays and constants from a module
  IMPLICIT NONE

  INTEGER :: LE           ! Line counter for data input
  REAL(8) :: AEE, DOS1    ! Raw energy & DOS read from file

  DEF = 0.0D0             ! Energy shift (can be used for alignment)
  OPEN(UNIT=90, FILE='parameter.txt', STATUS='OLD', ACTION='READ')
  READ(90, *)
  READ(90, '(25X, F12.8)') DEF
  WRITE(*,*) "DEF(Energy shift offset):", DEF
  CLOSE(90)

  !-----------------------------------------------------------------------
  ! Step : Load energy-resolved data from input file 'AKK.DATA'
  ! File format : Each line contains five floating-point values:
  !               E_raw, V2(E), VD(E), DOS(E), EN(E)
  !               - E_raw : Raw energy value [eV]
  !               - V2(E) : Velocity squared (used for transport)
  !               - VD(E) : Scattering or damping term
  !               - DOS(E): Density of states
  !               - EN(E) : Additional energy-dependent property
  !
  ! Purpose :
  !   - Read and store energy grid and associated quantities into shared arrays
  !   - Apply energy offset DEF (if necessary) to align energy axis
  !   - This data is used by BS() to compute chemical potential via carrier density
  !
  ! Implementation Notes :
  !   - STATUS='OLD' ensures file must exist; else error
  !   - ACTION='READ' specifies file will not be overwritten
  !   - IMAX = number of lines/data points to read (energy mesh size)
  !-----------------------------------------------------------------------
  OPEN(UNIT=10, FILE='AKK.DATA', STATUS='OLD', ACTION='READ')

  DO LE = 1, IMAX
     READ(10,'(5E15.8)') AEE, V2(LE), VD(LE), DOS1, EN(LE)
     EE(LE)   = AEE - DEF     ! Shift energy axis if needed
     DOS(LE)  = DOS1          ! Assign DOS
     ! Note: V2A, VDA, DOSA are unused in current routine
     IF(LE > 2 .AND. EE(LE-1) <= 0.0 .AND. EE(LE) >= 0.0) THEN
       VEC = EN(LE-1)+(EN(LE) - EN(LE-1))/(EE(LE) - EE(LE-1))
       WRITE(*,*) "The valence electron concentration (VEC):", VEC, &
                & "at DEF =", DEF, "[eV] (for EF = 0.0 -> EF = DEF) and 0 [K]"
     END IF
  END DO

  CLOSE(10)
END SUBROUTINE DAT
