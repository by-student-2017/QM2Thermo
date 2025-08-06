!-----------------------------------------------------------------------
! Author   : H. Sato and M. Inukai
! Affiliation : [AUE / ----]
! Contact  : [uichiro25@gmail.com] (U. Mizutani: https://doi.org/10.1007/s11669-024-01103-0)
! GitHub   : https://github.com/by-student-2017/LBT-TETRA (optional)
!-----------------------------------------------------------------------
! Program : seebeck_analysis.f90
! Purpose : Calculate Seebeck coefficient using band structure data
!           and multiple scattering models (phonon, impurity, DOS).
!
! Dates   : Created 2019.11.12, Revised 2020.05.16, 2024.12.31
!         : Revised 2025.07.12
!
! Physical constants:
!   TEM = k * 300 [K] = 0.000861734 * 300 = 0.0258199727 [eV] (thermal energy)
!   e   = 1.602176565 * 10^-19  [C]   (charge conversion: 1 [eV] -> [Coulomb])
!
! Compilation:
!   gfortran -O2 seebeck_analysis.f90 -o seebeck_analysis.exe
!   ifort    -O2 seebeck_analysis.f90 -o seebeck_analysis.exe
!
! Required input files:
!   - AKK.DATA        : Band-resolved energy, velocity^2, velocity^2*DOS, DOS, Cumulative (DOS)
!     (Energy [eV], group velocity^2 [(m/s)^2], group velocity^2 * Total DOS [(m/s)^2 * (states/eV)], Total DOS [states/eV], Cumulative DOS [states/eV])
!   - apot.data       : Temperature-dependent chemical potentials
!   - phononDOS.dat   : (Optional) Phonon DOS for tau(T) model
!   - lambda          : (Optional) lambda (QE) for tau(T) model
!   - a2F.dos         : (Optional) a2F.dos (QE) for tau(T) model
!
! Output:
!   - Seebeck_analysis.dat : Temperature-dependent Seebeck results
!-----------------------------------------------------------------------
! Model switches (in MODULE seebeck_data):
!   - PHONON_SCATTERING     : tau is proportional to 1 / (|E - mu| * T) <- Always enabled
!   - IMPURITY_SCATTERING   : tau is proportional to |E - mu|^n         <- User-defined exponent `n_impurity` (n = 1: Mild energy-dependence, 3.0: high energy, 4.0: more high energy)
!   - DOS_DEPENDENCE        : tau is proportional to 1 / [DOS(E) * T]   <- Enabled when DOS(E) is available
!   - PHONON_DOS_SCATTERING : tau(T) is proportional to integral D_ph(omega)/omega * [1+n_B] d(omega) <- Enabled when phononDOS.dat is present (need number of atoms)
!   ------------------------ in MODULE seebeck_data
!   Switch variables (editable in PROGRAM):
!   LOGICAL :: use_phonon     = .FALSE.            ! .FALSE. or .TRUE.
!   LOGICAL :: use_impurity   = .FALSE.            ! .FALSE. or .TRUE.
!   LOGICAL :: use_dos        = .FALSE.            ! .FALSE. or .TRUE.
!   LOGICAL :: use_a2Fdos     = .FALSE.            ! .FALSE. or .TRUE.
!   LOGICAL :: use_phononDOS  = .FALSE.            ! .FALSE. or .TRUE.
!   LOGICAL :: use_selection_filter = .FALSE.      ! Whether to use a selection rule filter
!   ------------------------ in MODULE seebeck_data
!   - tau0 = 1.0D-14                               ! Base relaxation time [s], arbitrary scaling (H. Sato original tau0 = 1.0): tau0 = 1/(1/tau_ph-ph(Alamode) + 1/tau_ph-el(EPW))
!   - REAL(KIND=8) :: n_impurity = 3.0D0           ! Impurity scattering exponent (1: Mild energy-dependence, 3.0: high energy, 4.0: more high energy)
!   - INTEGER :: N_atom = 1                        ! The number of atoms for normalization of phononDOS
!   - REAL(KIND=8), DIMENSION(NPH) :: g_eff_omega  ! Eliashberg function component: g_eff(omega_j)^2 for mode j: Electron-phonon coupling strength: Eliashberg function (default: usually 1.0)
!   - REAL(KIND=8), DIMENSION(NPH) :: w_mode       ! Selection rule weight for phonon mode j (e.g., symmetry, polarization-based filtering): Mode transition, selection-based weighting coefficient (0.0 to 1.0) (default: normally 1.0)
!
! Switch control:
!   - You may edit `n_impurity` or comment out individual tau models in `get_tau()`
!-----------------------------------------------------------------------
! Function and Subroutine Hierarchy for Seebeck_analysis.f90
!
! Main Source File:
!   Seebeck_analysis.f90
!   +----- Step 0: Load parameters
!   |       +----- Reads parameters from 'parameter.txt', and shows parameters.
!   |       +----- Reads lattice parameter from 'wien.struct', and shows lattice parameters.
!   |       +----- calculate volume and show its value.
!   |    
!   +----- Step 1: Load Electron-Phonon Coupling Data
!   |       +----- Reads phonon coupling constants from 'lambda'.
!   |       +----- Initializes arrays for broadening, lambdaArray, dosEf, and omega_ln.
!   |       +----- Uses `CALL ReadLambdaData()` for data processing.
!   |    
!   +----- Step 2: Load Chemical Potentials
!   |       +----- Reads temperature-dependent chemical potentials (AMU) from 'apot.data'.
!   |       +----- Associates each temperature (TEM) with its corresponding chemical potential.
!   |    
!   +----- Step 3: Load Band Velocity Data
!   |       +----- Reads energy mesh, density of states (DOS), squared group velocity,
!   |       |       and band energy levels from 'AKK.DATA'.
!   |       +----- Applies energy shift if necessary.
!   |       +----- Computes valence electron concentration (VEC) for given energy levels.
!   |    
!   +----- Step 4: Read Phonon DOS or a2F.dos Data
!   |       +----- Calls `read_phonon_dos()` if phononDOS is enabled.
!   |       +----- Calls `read_a2F_dos()` if a2F.dos is enabled.
!   |    
!   +----- Step 5: Initialize Output File
!   |       +----- Prepares 'Seebeck_analysis.dat' file for results storage.
!   |       +----- Includes headers with energy offsets, VEC, and Seebeck coefficient columns.
!   |    
!   +----- Step 6: Compute Seebeck Coefficient
!   |       +----- Loops through temperatures to calculate Seebeck coefficient.
!   |       +----- Computes Fermi-Dirac distribution and its derivative at each energy level.
!   |       +----- Multiplies by group velocity squared and energy-dependent scattering rate.
!   |       +----- Accumulates numerator and denominator using integration techniques:
!   |       |       +----- Options for Riemann sum, trapezoidal rule, or Simpson's rule.
!   |       +----- Outputs averaged energy deviation <E - mu>, Seebeck coefficient, and Nc.
!   |    
!   +----- Finalize:
!       +----- Closes output files.
!       +----- Deallocates arrays for phonon and scattering data if used.
!       +----- Ensures proper cleanup for each flag (use_phonon, use_phononDOS, use_a2Fdos).
!
! Input Files:
!   +----- lambda          : Electron-phonon coupling constants and DOS data.
!   +----- AKK.DATA        : Energy mesh, squared velocity, squared velocity * DOS, DOS, Cumulative DOS.
!   +----- apot.data       : Temperature (TEM) and chemical potential (AMU) values.
!   +----- phononDOS.dat   : Phonon DOS data (optional).
!   +----- a2F.dos         : Eliashberg function data (optional).
!
! Output Files:
!   +----- Seebeck_analysis.dat : Results file containing d Seebeck coefficients.
!-----------------------------------------------------------------------
! Function and Subroutine Relationships for get_tau
!
! Main Function:
!   get_tau(E1, T, I)
!   +----- Purpose:
!   |       Compute the relaxation time tau(E, T) based on multiple scattering 
!   |       mechanisms combined using Matthiessen's rule.
!   |    
!   +----- Scattering Mechanisms:
!   |       1. Phonon scattering (tau_ph)      : Dependent on |E - mu|, T, and Lambda.
!   |       2. Impurity scattering (tau_imp)   : User-defined |E - mu|^n exponent model.
!   |       3. DOS-based scattering (tau_dos)  : Dependent on DOS(E) * T.
!   |       4. Phonon DOS-based scattering (tau_phdos): Advanced scattering using phonon data.
!   |    
!   +----- Subroutines and Functions Called:
!   |       +----- GetLambda(E1, T, I)
!   |       |       +----- Purpose: Interpolates the electron-phonon coupling constant Lambda.
!   |       |       +----- Inputs:
!   |       |       |       - E1: Energy offset from chemical potential.
!   |       |       |       - T: Temperature.
!   |       |       |       - I: Energy index.
!   |       |       +----- Output:
!   |       |       |       - Lambda: Electron-phonon coupling constant.
!   |       |       +----- Dependencies:
!   |       |           - InterpolateLambda(E1, omega_ln, lambdaArray): Performs interpolation.
!   |       |           - ReadLambdaData(): Loads Lambda and phonon coupling data from 'lambda'.
!   |       |    
!   |       +----- tau_a2Fdos_ET(E1, T)
!   |       |       +----- Purpose: Computes relaxation time from Eliashberg function (a2F.dos data).
!   |       |       +----- Inputs:
!   |       |       |       - E1: Energy offset from chemical potential.
!   |       |       |       - T: Temperature.
!   |       |       +----- Output:
!   |       |           - Relaxation time based on phonon density and Bose-Einstein distribution.
!   |       |           - Dependencies: read_a2F_dos(), kernel_selection(), integration over a2F_total.
!   |       |    
!   |       +----- tau_phdos_ET(E1, T)
!   |       |       +----- Purpose: Computes relaxation time based on phonon DOS and n_B (Bose-Einstein).
!   |       |       +----- Inputs:
!   |       |       |       - E1: Energy offset from chemical potential.
!   |       |       |       - T: Temperature.
!   |       |       +----- Output:
!   |       |           - Relaxation time derived from phonon density of states (phononDOS.dat).
!   |       |           - Dependencies: read_phonon_dos(), kernel_selection(), phonon integrals, (GetLambda()).
!   |       |    
!   |       +----- read_a2F_dos()
!   |       |       +----- Purpose: Reads phonon density of states from a2F.dos.
!   |       |       +----- Output: Frequency-dependent a2F data.
!   |       |    
!   |       +----- read_phonon_dos()
!   |       |       +----- Purpose: Reads phonon density of states from phononDOS.dat.
!   |       |       +----- Output: Frequency-dependent phonon DOS data.
!   |       |    
!   |       +----- kernel_selection(E, omega)
!   |       |       +----- Purpose: Determines scattering kernel weights based on energy and phonon frequencies.
!   |       |       +----- Inputs:
!   |       |       |       - E: Energy.
!   |       |       |       - omega: Phonon frequency.
!   |       |       +----- Output: Logical weight applied to scattering integrals.
!   |       |    
!   |       +----- DOS(I)
!   |           +----- Purpose: Provides the density of states at a given energy index.
!   |           +----- Dependencies: read_constants() (for grid setup) and external DOS data.
!   |    
!   +----- Integration Techniques Used:
!   |       - Matthiessen's rule combines scattering rates (1/tau) for active mechanisms.
!   |       - Fallback to tau0 if no active scattering mechanisms are enabled.
!   |    
! Inputs to get_tau:
!   - E1: Energy offset from chemical potential (E - mu) [eV].
!   - T: Temperature [K].
!   - I: Index on energy mesh corresponding to E = EE(I).
!
! Outputs from get_tau:
!   - tau(E, T): Relaxation time [s], used in conductivity and Seebeck calculations.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! MODULE : kernel_model
! Purpose:
!   Provides kernel functions for modeling energy-dependent phonon-electron
!   scattering in Allen-type relaxation time calculations.
!
! Contents:
!   - FUNCTION: kernel_selection(E, omega)
!     Applies a basic selection rule:
!       - Scattering is active only when phonon energy omega < electron energy offset |E|
!
! Usage:
!   - Used in tau_phdos_ET(E,T) for filtering mode contributions via omega
!   - Extendable for symmetry-based selection, transition-specific kernels,
!     or mode-resolved filtering.
!
! Extension-ready Features:
!   - Implement alternative kernel forms (step, Lorentzian, delta)
!   - Enable kernel loading from external config files or runtime parameters
!-----------------------------------------------------------------------
MODULE kernel_model
IMPLICIT NONE
CONTAINS

  REAL(KIND=8) FUNCTION kernel_selection(E, omega)
    REAL(KIND=8), INTENT(IN) :: E, omega
    ! Simple selection rule: scattering only if phonon energy < |E|
    IF (omega < ABS(E)) THEN
       kernel_selection = 1.0D0
    ELSE
       kernel_selection = 0.0D0
    END IF
  END FUNCTION kernel_selection

END MODULE kernel_model



!-----------------------------------------------------------------------
! MODULE: seebeck_data
! Purpose:
!   Centralize definitions, data structures, and functions for
!   transport property analysis (Seebeck coefficient, relaxation time)
!
! Contents:
!   - Physical constants and scaling parameters
!   - Energy-dependent quantities and DOS arrays
!   - Temperature schedule and chemical potential
!   - Phonon DOS input and normalization routines
!   - tau(E,T) models combining phonon, impurity, and DOS effects
!-----------------------------------------------------------------------
MODULE seebeck_data
  IMPLICIT NONE

  ! === Constants for physical parameters and numerical settings ===
  !-----------------------------------------------
  INTEGER, PARAMETER :: MM = 60000                ! Number of energy mesh points (energy resolution)
  REAL(KIND=8), PARAMETER :: EV  = 13.605698D0    ! Conversion constant: Hartree -> eV
  REAL(KIND=8), PARAMETER :: PI  = 3.141592654D0  ! Pi
  REAL(KIND=8), PARAMETER :: CO  = 1.0D6          ! Conversion factor to muV/K for Seebeck output
  REAL(KIND=8), PARAMETER :: B2A  = 0.52918       ! convert Bohr to Angstrom unit
  REAL(KIND=8), PARAMETER :: ech = 1.602176565D-19 ! [C]  (Note: 1 [C/m^3] = 1 [Coulomb/m^3] = 1 [J])
  REAL(KIND=8), PARAMETER :: ems = 9.109383713D-31 ! [kg]
  REAL(KIND=8), PARAMETER :: hbar = 6.582119D-16  ! Planck constant [eV s]
  REAL(KIND=8), PARAMETER :: kb = 8.617333D-5     ! Boltzmann constant [eV/K]
  !-----------------------------------------------
  REAL(KIND=8) :: DEF = 0.0D0                     ! Energy offset (if needed)
  !-----------------------------------------------
  
  ! Model switching (Turn on only what you need) (.TRUE. or .FALSE.)
  !-----------------------------------------------
  LOGICAL :: use_phonon     = .FALSE.             ! .FALSE. or .TRUE.
  LOGICAL :: use_impurity   = .FALSE.             ! .FALSE. or .TRUE.
  LOGICAL :: use_dos        = .FALSE.             ! .FALSE. or .TRUE.
  LOGICAL :: use_a2Fdos     = .FALSE.             ! .FALSE. or .TRUE.
  LOGICAL :: use_phononDOS  = .FALSE.             ! .FALSE. or .TRUE.
  LOGICAL :: use_selection_filter = .FALSE.       ! see "SUBROUTINE generate_w_mode(w_mode, WPH)"
  REAL(KIND=8) :: tau0 = 1.0D-14                  ! Base relaxation time [s], arbitrary scaling, (Common values: 1.0D-14). Its absolute value cancels out for Seebeck, but remains important for conductivity
  REAL(KIND=8) :: n_impurity = 3.0D0              ! Impurity scattering exponent (1: Mild Eenergy-dependence, 3.0: high energy, 4.0: more high energy)
  INTEGER :: N_atom = 1                           ! The number of atoms for normalization of phononDOS
  REAL(KIND=8) :: L_bound                         ! Boundary scattering. (works: phononDOS = T and N_atom)
  REAL(KIND=8) :: sound_velocity                  ! sound velocity (this code uses an approximation based on the maximum frequency of phononDOS normalized to 3N) [Angstrom]
  REAL(KIND=8) :: wmax                            ! Maximum frequency in phonon DOS normalized to 3N (= 3 * Natom). [eV]
  REAL(KIND=8) :: volume                          ! Volume [A^3]
  REAL(KIND=8) :: b_para                          ! Umklapp (Klemens-Callaway type model) scattering (Ref. 1 - 2) (works: phononDOS = T)
  REAL(KIND=8) :: C_phel                          ! Phonon-Electron scattering coefficient (works: a2Fdos = F): (Ref. 2.31e10) [s^-1*eV^-2]
  REAL(KIND=8) :: B_pdef                          ! Point defect scattering coefficient: (Ref. 5.33e15) [s^-1*eV^-4]
  REAL(KIND=8) :: Bulk_modulus                    ! Bulk_modulus [GPa] (1 [eV/A^3] = 160.2 [GPa]), B = K
  REAL(KIND=8) :: dB_per_dP                       ! dB/dP from EOS
  REAL(KIND=8) :: dB_per_dV                       ! dB/dV
  REAL(KIND=8) :: V0                              ! V0 from EOS
  REAL(KIND=8) :: dG_per_dV                       ! dG/dV
  REAL(KIND=8) :: Shear_modulus                   ! Shear_modulus [GPa] = 3*Bulk_modulus*(1-2*Poisson_ratio) / (2*(1+Poisson_ratio)), G
  REAL(KIND=8) :: Young_modulus                   ! Young_modulus [GPa] = 9*B*G/(3*B+G), E
  REAL(KIND=8) :: Poisson_ratio                   ! Poisson_ratio = (3*B-2*G)/(2*(3*B+G))
  REAL(KIND=8) :: Poisson_ratio_read_value        ! Poisson_ratio_read_value <= 0.0 -> .TRUE. (use vs)
  REAL(KIND=8) :: density                         ! read [g/cm^3] unit -> density * 1000 [kg/m^3]
  REAL(KIND=8) :: vl                              ! sound velocity (longitudinal wave) [m/s] = ((Bulk_modulus*1.0D9+(4/3)*Shear Modulus*1.0D9)/(density*1000.0))**(0.5)   ! 1 [GPa] = 1.0e9 [kg/(m*s^2)
  REAL(KIND=8) :: vt                              ! sound velocity (transverse wave) [m/s] = (Bulk_modulus*1.0D9/(density*1000.0))**(0.5)                                 ! 1 [GPa] = 1.0e9 [kg/(m*s^2)
  REAL(KIND=8) :: vs                              ! Shear wave velocity
  REAL(KIND=8) :: va                              ! The average sound wave velocity (the Harmonic Mean of velocities in this code.)
  REAL(KIND=8) :: Mavg                            ! The mean atomic mass, density * 1.0D6 [g/m] * (volume/N_atom) [A^3/N] * 1.0D30/6.022D23
  REAL(KIND=8) :: MPF_phonon                      ! mean free path of phonon, l = (vl + 2*vt)/3 * tau0_phonon [s]
  REAL(KIND=8) :: tau0_phonon                     ! For phonons, it is about 10-100 times stronger than for electrons. (If it is 0.0, tau_ph * 100)
  REAL(KIND=8) :: Gruneisen_parameter             ! It corresponds to the ratio of the second-order elastic constant C to the third-order elastic constant D.
  REAL(KIND=8) :: Gruneisen_parameter_L           ! 
  REAL(KIND=8) :: Gruneisen_parameter_S           ! 
  REAL(KIND=8) :: Apara                           ! parameter A of Slack model.
  REAL(KIND=8) :: kappa_phonon_min                ! Cahill molde: (1/2)*(PI/6)**(1/3)*kb*(volume/N_atom)**(2/3)*(vl + 2*vt)
  REAL(KIND=8) :: kappa_phonon_min_Clarke         !
  REAL(KIND=8) :: kappa_phonon_min_Cahill         !
  REAL(KIND=8) :: kappa_phonon_min_Slack          !
  REAL(KIND=8) :: kappa_phonon_min_Slack_xK       !
  REAL(KIND=8) :: km                              ! Cezairliyan
  REAL(KIND=8) :: kappa_phonon                    ! kappa_phonon = (1.0D0/3.0D0) * Cv_DOS * ((vl + 2.0D0*vt)/3.0D0)**2.0D0 * tau0_phonon
  REAL(KIND=8) :: CN                              ! Coordination number
  LOGICAL :: use_Apara_gamma                      !
  !-----------------------------------------------
  REAL(KIND=8) :: Nd                              ! Doping concentration in cm^-3: n-type 1.0e14 - 1.0e18 [cm^-3]
  REAL(KIND=8) :: Nd_hole                         ! hole
  REAL(KIND=8) :: Nd_electron                     ! electron
  REAL(KIND=8) :: Nc                              ! Carrier concentration at TEM [cm^-3]
  !-----------------------------------------------
  REAL(KIND=8) :: VEC                             ! The valence electron concentration (VEC)
  REAL(KIND=8) :: VEC0                            ! The valence electron concentration (VEC) from parameter.txt
  !-----------------------------------------------
  
  !Note: The Gruneisen constant is a quantity that represents the deviation (anharmonicity) of the lattice vibration of a solid from harmonic vibration. 
  !  It corresponds to the ratio of the second-order elastic constant C to the third-order elastic constant D. 
  !  The thermal expansion coefficient alpha = Gruneisen parameter * constant-volume specific heat Cv.
  
  ! Electron-phonon coupling constant, lambda
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: broadening, lambdaArray, dosEf, omega_ln
  
  ! phononDOS
  !-----------------------------------------------
  INTEGER :: NPH                                  ! Number of phonon grid points
  !INTEGER, PARAMETER :: nNPH = 10000             ! (old vesion) Number of phonon grid points
  !-----------------------------------------------
  INTEGER :: nfiles, ndata, nmodes                ! a2F.dos
  REAL(KIND=8), ALLOCATABLE :: freq(:,:)          ! a2F.dos
  REAL(KIND=8), ALLOCATABLE :: a2F_total(:,:)     ! a2F.dos
  REAL(KIND=8), ALLOCATABLE :: a2F_mode(:,:,:)    ! a2F.dos
  !-----------------------------------------------
  !REAL(KIND=8), SAVE :: WPH(nNPH), DOSPH(nNPH)   ! (old vesion) Phonon frequency [eV], phonon DOS [arb.unit](recommend [states/eV/unitcell] or [states/eV/Angstrom^3]
  REAL(KIND=8), ALLOCATABLE :: WPH(:), DOSPH(:)   ! Phonon frequency [eV], phonon DOS
  !-----------------------------------------------
  !REAL(KIND=8), DIMENSION(nNPH) :: g_eff_omega   ! (old vesion) Mode-dependent electron-phonon coupling strength [(eV)^0.5]
  REAL(KIND=8), ALLOCATABLE :: g_eff_omega(:)     ! Mode-dependent electron-phonon coupling strength [(eV)^0.5]
  !-----------------------------------------------
  !REAL(KIND=8), DIMENSION(nNPH) :: w_mode        ! (old vesion) Mode selection filter (0.0 - 1.0): symmetry, polarization, activation criteria
  REAL(KIND=8), ALLOCATABLE :: w_mode(:)          ! Mode selection filter (0.0 - 1.0): symmetry, polarization, activation criteria
  !-----------------------------------------------
  
  ! Debye calculation using phonon dos.
  real(8) :: Theta_D, Cv_DOS, Cv_Debye
  real(8) :: Theta_D_va                           ! The Debye temperature is evaluated from the sound velocities.
  real(8) :: Theta_D_Cezairliyan_equ              ! use Cezairliyan equ.
  
  ! === Temperature list and chemical potentials ===
  REAL(KIND=8), SAVE :: TT(25), AMU(25)
  ! TT   : Temperatures [K] at which transport properties are computed
  ! AMU  : Chemical potential [eV] at corresponding temperature
  
  ! === Arrays for energy-dependent properties ===
  REAL(KIND=8), SAVE :: EE(MM), DOS(MM), EN(MM)                ! EE: shifted energy mesh, DOS: density of states, EN: raw energy values
  REAL(KIND=8), SAVE :: GV2(MM)                                ! squared group velocity
  REAL(KIND=8), SAVE :: GV2D(MM), GV2D2(MM), EW(MM), DOSW(MM)  ! GV2D: squared group velocity * DOS, other arrays reserved
  REAL(KIND=8), SAVE :: DOSA(MM), GV2DA(MM), FDOS(MM)          ! Reserved arrays for filtered or weighted outputs
  
  ! === Temperature values to iterate over ===
  DATA TT /1000D0,950D0,900D0,850D0,800D0,750D0,700D0,650D0,600D0,550D0, &
           500D0,450D0,400D0,350D0,300D0,250D0,200D0,150D0,100D0,50D0,   &
           40D0,30D0,20D0,10D0,5D0/

CONTAINS
  
  !-----------------------------------------------------------------------
  ! Function : tau_a2Fdos_ET
  !
  ! Allen-type electron-phonon scattering model using Eliashberg function.
  ! Assumes alpha^2*F(omega) = proportional to a2F_total(i,j), contributing to energy-dependent scattering rate.
  !
  ! Purpose  : Compute the energy- and temperature-dependent relaxation time tau(E,T).
  !
  ! Model    :
  !   tau(E,T)^-1 is proportional to integral a2F_total(i,j) * [1 + n_B(omega,T)] / omega d(omega)
  !   where:
  !     a2F_total(i,j) = Eliashberg function (electron-phonon interaction strength).
  !     n_B(omega,T)   = Bose-Einstein distribution = 1 / (exp(omega/kT) - 1).
  !     FD1, FD2       = Fermi-Dirac distributions (optional; for advanced scenarios).
  !     kernel_val     = Energy-matching kernel applied to scattering events.
  !
  ! Inputs:
  !   E    : Energy level for which relaxation time is computed [eV].
  !   T    : Temperature [K].
  ! Note: E <- (get_tau(): E1 = Energy offset from chemical potential, E - mu [eV])
  !
  ! Returns:
  !   tau_a2Fdos_ET : Estimated tau(E, T) in seconds [s] using scaled tau0.
  !
  ! Notes:
  ! - Requires freq(:,:), a2F_total(:,:), nfiles, and ndata arrays to be initialized.
  ! - If a2F.dos data is missing, result defaults to zero.
  ! - Used as one of the scattering contributions in get_tau().
  !-----------------------------------------------------------------------
  REAL(KIND=8) FUNCTION tau_a2Fdos_ET(E, T)
    USE kernel_model
    
    ! Input variables
    REAL(KIND=8), INTENT(IN) :: E, T ! Energy and temperature as function inputs
    
    ! Local variables
    INTEGER :: i, j                  ! Loop indices
    REAL(KIND=8) :: omega            ! Frequency [Rydberg]
    REAL(KIND=8) :: d_omega          ! delta Frequency
    REAL(KIND=8) :: n_B              ! Bose-Einstein distribution
    REAL(KIND=8) :: FD1, FD2         ! Fermi-Dirac probabilities (optional)
    REAL(KIND=8) :: CP               ! Chemical potential based on temperature (from TT and AMU arrays)
    REAL(KIND=8) :: integrand        ! Individual contribution to the relaxation time integral
    REAL(KIND=8) :: kernel_val       ! Kernel function value for scattering events
    REAL(KIND=8) :: sum              ! Accumulator for integral results
    
    !-----------------------------------------------------------------------
    ! STEP 1: Identify material-specific constant (CP) from given temperature
    !         TT(:) and AMU(:) must be initialized externally with material data.
    !DO i = 1, 25
    !  IF (TT(i) == T) THEN
    !    CP = AMU(i)  ! Assign CP based on temperature match
    !    EXIT         ! Exit the loop once matched
    !  END IF
    !END DO
    ! Note: E <- (get_tau(): E1 = Energy offset from chemical potential, E - mu [eV])
    
    !------------------------------
    ! Note: data for a2F.dos files
    !WRITE(*,'(3X,1000E16.6)') freq(i,j), a2F_total(i,j), (a2F_mode(i,k,j), k = 1, nmodes)
    
    !-----------------------------------------------------------------------
    ! STEP 2: Initialize sum to zero for integral computation
    sum = 0.0D0
    
    !-----------------------------------------------------------------------
    ! STEP 3: Loop through a2F.dos data to compute contributions to tau^-1
    DO i = 1, nfiles
      DO j = 1, ndata             ! Iterate through frequency values from a2F.dos
        omega = freq(i,j)         ! Frequency for current data point
        IF (omega > 1.0D-4) THEN  ! Ignore very small frequencies for numerical stability
          
          ! Compute Bose-Einstein distribution for the phonon mode
          n_B = 1.0D0 / (DEXP(omega / (kb * T)) - 1.0D0)     ! Bose-Einstein distribution
          
          ! Compute Fermi-Dirac probabilities for advanced calculations
          FD1 = 1.0D0 / (1.0D0 + DEXP((E - CP - omega) / (kb * T)))  ! Fermi-Dirac probability
          FD2 = 1.0D0 / (1.0D0 + DEXP((E - CP + omega) / (kb * T)))  ! Fermi-Dirac probability
          
          ! Optional: Toggle kernel interaction
          kernel_val = kernel_selection(E, omega)            ! Kernel for energy and frequency
          
          !-----------------------------------------------------------------------
          ! Compute the integrand for the relaxation time
          ! Uncomment appropriate line depending on whether FD terms are included:
          !integrand = a2F_total(i,j)*EV * (1.0D0 + n_B) * kernel_val / omega
          !-----------------------------------------------------------------------
          integrand = 2.0D0 * PI * a2F_total(i,j)*EV * (1.0D0 + 2.0D0 * n_B + FD1 + FD2) * kernel_val / omega
          !-----------------------------------------------------------------------
          ! Note: lambda = 2 * integral a2F(omega)/omega d(omega), omega [Ry], a2F(omega) [Ry]
          !-----------------------------------------------------------------------
          
          ! Add the integrand's contribution to the summation
          IF (j > 1) THEN
            d_omega = MAX(freq(i,j) - freq(i,j-1), 1.0D-12)  ! Avoid zero or negative width
          ELSE
            d_omega = freq(i,1)
          END IF
          sum = sum + integrand * d_omega
        END IF
      END DO
    END DO
    
    ! Relaxation time calculation
    ! Avoid division by zero
    tau_a2Fdos_ET = tau0 / (sum + 1.0D-12)
    
  END FUNCTION tau_a2Fdos_ET
  
  
  !-----------------------------------------------------------------------
  ! Function : tau_phdos_ET
  !
  ! Allen-type phonon scattering model using Eliashberg function approximation.
  ! Assumes alpha^2*F(omega) = near g_eff^2 * D_ph(omega), contributing to energy-dependent scattering rate.
  !                          = near lambda(omega) * Total-D_ph(omega)
  !
  ! Purpose  : Compute the energy- and temperature-dependent relaxation time tau(E,T)
  !
  !          : based on phonon density of states (DOS_phonon)
  !
  ! Model    :
  !   tau(E,T)^-1 is proportional to integral [D_ph(omega) / omega] * [1 + n_B(omega,T)] d(omega)
  !   where:
  !     D_ph(omega)  = phonon DOS at frequency omega. DOSPH(j) is nearly equal to D_ph(omega) at omega = WPH(j)
  !     n_B(omega,T) = Bose-Einstein distribution = 1 / (exp(omega/kT) - 1)
  !     FD1, FD2       = Fermi-Dirac distributions (optional; for advanced scenarios).
  !     kernel_val     = Energy-matching kernel applied to scattering events.
  !
  ! Inputs:
  !   - E    : Energy level for which relaxation time is computed [eV]
  !   - T    : Temperature [K]
  ! Note: E <- (get_tau(): E1 = Energy offset from chemical potential, E - mu [eV])
  !
  ! Returns:
  !   tau_phdos_ET : Estimated tau(E, T) in seconds [s] using scaled tau0
  !
  ! Notes:
  ! - Requires the following arrays to be properly initialized:
  !   - WPH(:) : Phonon frequencies [eV]
  !   - DOSPH(:) : Phonon density of states (DOS) [normalized to E]
  ! - If `phononDOS.dat` is missing or contains zeros, tau_phdos_ET will default to zero.
  ! - Used as one of the scattering contributions in get_tau()
  !
  ! Extension-ready for Allen model:
  ! - Eliashberg function approximation via g_eff_omega(j)^2 * D_ph(omega)
  ! - Mode filtering and kernel K(E, omega) for scattering selectivity
  ! - External configuration for material-specific parameters
  ! - g_eff_omega(j)^2 * w_mode(j) * D_ph(omega_j) represents alpha^2*F(omega) with selection filter
  ! - kernel_selection(E, omega_j) applies energy-matching condition for scattering events
  !-----------------------------------------------------------------------
  REAL(KIND=8) FUNCTION tau_phdos_ET(E, T)
    USE kernel_model
    
    ! Input variables
    REAL(KIND=8), INTENT(IN) :: E, T  ! Energy and temperature as function inputs
    
    ! Local variables
    INTEGER :: i, j                   ! Loop indices
    REAL(KIND=8) :: omega             ! Phonon frequency [Rydberg] -> [eV]
    REAL(KIND=8) :: d_omega           ! delta Frequency
    REAL(KIND=8) :: n_B               ! Bose-Einstein distribution value
    REAL(KIND=8) :: FD1, FD2          ! Fermi-Dirac probabilities for absorption/emission
    REAL(KIND=8) :: CP                ! Chemical potential or material constant (read from TT and AMU)
    REAL(KIND=8) :: integrand         ! Individual contribution to the relaxation time integral
    REAL(KIND=8) :: kernel_val        ! Kernel function value for energy/frequency filtering
    REAL(KIND=8) :: sum               ! Accumulator for the integral results
    
    ! Constant for k-point selection (currently fixed to 1)
    INTEGER :: I1 = 1                 ! Index in energy mesh corresponding to E = EE(I)
    
    !-----------------------------------------------------------------------
    ! STEP 1: Identify material-specific constant (CP) from given temperature
    !         TT(:) and AMU(:) must be initialized externally with material data.
    !DO i = 1, 25
    !  IF (TT(i) == T) THEN
    !    CP = AMU(i)  ! Assign CP based on temperature match
    !    EXIT         ! Exit the loop once matched
    !  END IF
    !END DO
    ! Note: E <- (get_tau(): E1 = Energy offset from chemical potential, E - mu [eV])
    
    !-----------------------------------------------------------------------
    ! STEP 2: Initialize sum to zero for integral computation
    sum = 0.0D0
    
    !-----------------------------------------------------------------------
    ! STEP 3: Loop through all phonon modes to compute contributions to tau^-1
    DO j = 1, NPH
       omega = WPH(j)            ! Phonon frequency from array WPH [eV]
       IF (omega > 1.0D-4) THEN  ! Ignore very small frequencies for numerical stability
          
          ! Compute Bose-Einstein distribution for the phonon mode
          n_B = 1.0D0 / (DEXP(omega / (kb*T)) - 1.0D0)  ! Bose-Einstein distribution
          
          ! Compute Fermi-Dirac probabilities for absorption/emission processes
          FD1 = 1.0D0 / (1.0D0 + DEXP((E - CP - omega) / (kb * T)))
          FD2 = 1.0D0 / (1.0D0 + DEXP((E - CP + omega) / (kb * T)))
          
          ! Optionally, you can weight terms with omega < E as an approximation of K(E, omega), K(E, omega) = kernel_selection(E, omega)
          ! Note: g_eff_omega = lambda / omega
          !g_eff_omega(j) = GetLambda(omega, T, I1) / omega  ! Lambda based on energy, temperature, and mode index
          
          ! Optionally, compute kernel function based on energy and frequency
          kernel_val = kernel_selection(E, omega)           ! Kernel for energy and frequency
          
          !-----------------------------------------------------------------------
          !integrand = w_mode(j) * DOSPH(j) * (1.0D0 + n_B) * kernel_val
          integrand = w_mode(j) * DOSPH(j) * (1.0D0 + 2.0D0 * n_B + FD1 + FD2) * kernel_val
          !-----------------------------------------------------------------------
          ! Allen-type contribution: Compute the integrand for the relaxation rate
          ! g_eff_omega(j)^2 * w_mode(j) * DOSPH(j) * omega is an approximation for alpha^2*F(omega) [eV]
          ! Uncomment the appropriate line depending on whether FD terms are included:
          !integrand = g_eff_omega(j)**2.0D0 * w_mode(j) * DOSPH(j) * (1.0D0 + n_B) * kernel_val
          !-----------------------------------------------------------------------
          !integrand = g_eff_omega(j)**2.0D0 * w_mode(j) * DOSPH(j) * (1.0D0 + 2.0D0 * n_B + FD1 + FD2) * kernel_val
          !-----------------------------------------------------------------------
          
          ! Add the integrand's contribution to the summation
          IF (j > 1) THEN
            d_omega = MAX(WPH(j) - WPH(j-1), 1.0D-12)  ! Avoid zero or negative width
          ELSE
            d_omega = WPH(1)
          END IF
          sum = sum + integrand * d_omega
       END IF
    END DO
    
    tau_phdos_ET = tau0 / (sum + 1.0D-12)  ! Avoid division by zero
  END FUNCTION tau_phdos_ET
  
  
  !-----------------------------------------------------------------------
  ! Function : kernel_selection
  !
  ! Purpose  : This function determines whether a given frequency (omega) 
  !            contributes to the energy (E) during scattering processes.
  !
  ! Model    :
  !   The function returns:
  !     - 1.0D0 : If omega is less than the absolute value of E (|E|).
  !     - 0.0D0 : Otherwise.
  !
  ! Inputs:
  !   - E     : Energy value [Rydberg].
  !   - omega : Frequency value [Rydberg].
  !
  ! Returns:
  !   kernel_selection : Logical weight (1 or 0) applied to scattering integrals.
  !
  ! Notes:
  ! - This function is used in the calculation of relaxation times and 
  !   determines energy-matching conditions for scattering events.
  ! - The logic can be extended or modified for more complex models involving 
  !   energy and frequency matching.
  !-----------------------------------------------------------------------
  REAL(KIND=8) FUNCTION kernel_selection(E, omega)
    REAL(KIND=8), INTENT(IN) :: E, omega
    
    ! Check if omega contributes to the scattering process
    IF (omega < ABS(E)) THEN
       kernel_selection = 1.0D0  ! Match condition satisfied
    ELSE
       kernel_selection = 0.0D0  ! Match condition not satisfied
    END IF
    
  END FUNCTION kernel_selection
  
  
  !-----------------------------------------------------------------------
  ! Function : GetLambda
  !
  ! Purpose  : Compute the electron-phonon coupling constant (lambda) 
  !            for a given energy (E1), temperature (T), and mode index (I).
  !
  ! Model    :
  !   This function interpolates the coupling constant using:
  !     - E1: Energy input [eV].
  !     - omega_ln: Logarithmic average phonon frequency [eV].
  !     - lambdaArray: Precomputed array of coupling constants.
  !
  ! Notes:
  ! - The interpolation is based on energy (E1) and omega_ln, ensuring 
  !   smooth transitions between data points.
  ! - Requires `omega_ln` and `lambdaArray` to be properly initialized 
  !   before calling this function.
  ! - The function can be extended for material-specific adjustments 
  !   or higher-dimensional interpolation methods.
  !
  ! Inputs:
  !   E1 : Energy level for which lambda is computed [eV].
  !   T  : Temperature [K].
  !   I  : Mode index or k-point (integer index).
  !
  ! Returns:
  !   lambda : Interpolated electron-phonon coupling constant.
  !
  ! Usage:
  ! - Designed to support electron-phonon scattering calculations
  !   within the Eliashberg or Allen-type models.
  !-----------------------------------------------------------------------
  FUNCTION GetLambda(E1, T, I) RESULT(lambda)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: E1  ! Energy input
    REAL(KIND=8), INTENT(IN) :: T   ! Temperature
    INTEGER, INTENT(IN) :: I        ! 
    REAL(KIND=8) :: lambda          ! Electron-phonon coupling constant
    
    ! Interpolate the lambda value based on input energy and other factors
    lambda = InterpolateLambda(E1, omega_ln, lambdaArray)  ! Interpolation uses omega_ln ([eV] unit)
    
  END FUNCTION GetLambda
  
  
  !-----------------------------------------------------------------------
  ! Subroutine : ReadLambdaData
  !
  ! Purpose    : This subroutine reads and processes the electron-phonon 
  !              coupling data (lambda) from a file named "lambda".
  !
  ! Model      :
  !   - Reads the "lambda" file, which contains:
  !       1. Broadening values.
  !       2. Electron-phonon coupling constants (lambdaArray).
  !       3. DOS at Fermi level (dosEf).
  !       4. Logarithmic phonon frequency (omega_ln).
  !   - Converts omega_ln from [K] to [eV].
  !   - Allocates memory dynamically for the data arrays.
  !   - Outputs the read data to the console for verification.
  !
  ! Inputs:
  !   - "lambda" file : The input file must follow a specific format:
  !       * Contains a header row (3 lines).
  !       * Columns: Broadening, lambdaArray, dosEf, omega_ln.
  !
  ! Outputs:
  !   - broadening(:) : Array of broadening values.
  !   - lambdaArray(:): Array of electron-phonon coupling constants.
  !   - dosEf(:)      : DOS at the Fermi level.
  !   - omega_ln(:)   : Logarithmic phonon frequency in [eV].
  !   - nfiles        : Number of data entries processed.
  !
  ! Notes:
  ! - If the "lambda" file is missing or improperly formatted, the subroutine will fail.
  ! - Ensure consistent units for input and output (e.g., omega_ln converted to eV).
  !-----------------------------------------------------------------------
  SUBROUTINE ReadLambdaData()
    IMPLICIT NONE
    
    ! Declare variables
    INTEGER :: nLines, i             ! Number of lines and loop index
    CHARACTER(LEN=200) :: line       ! Line buffer for reading file
    
    !-----------------------------------------------------------------------
    ! Step 1: Count lines (excluding header)
    nLines = 0                       ! Initialize line counter
    OPEN(UNIT=10, FILE='lambda', STATUS='OLD', ACTION='READ')  ! Open the file
    DO
      READ(10, '(A)', END=100) line  ! Read lines until EOF
      nLines = nLines + 1
    END DO
100 CONTINUE
    REWIND(10)                       ! Reset file pointer to the beginning
    nLines = nLines - 3              ! Adjust for header rows (3 lines)
    
    !-----------------------------------------------------------------------
    ! Step 2: Allocate memory for data arrays
    ALLOCATE(broadening(nLines), lambdaArray(nLines), dosEf(nLines), omega_ln(nLines))
    
    ! Skip header rows
    DO i = 1, 3
      READ(10, *)                    ! Skip the header rows
    END DO
    
    !-----------------------------------------------------------------------
    ! Step 3: Read and process data lines
    DO i = 1, nLines
      READ(10, '(12X, F8.4, 12X, F8.4, 8X, F8.4, 13X, F12.4)', END=200) broadening(i), &
        & lambdaArray(i), dosEf(i), omega_ln(i)
      omega_ln(i) = omega_ln(i) * kb     ! Convert omega_ln from [K] to [eV]
      
      ! Output read data for verification
      WRITE(*,'(3(F8.4, 2X), F12.4, 2X, F8.4)') broadening(i), lambdaArray(i), dosEf(i), omega_ln(i)/kb, omega_ln(i)
    END DO
200 CONTINUE
    CLOSE(10)                        ! Close the file
    
    !-----------------------------------------------------------------------
    ! Step 4: Store the number of processed entries
    nfiles = SIZE(omega_ln)          ! Total number of entries
    WRITE(*,*) "Number of a2Fdos files: ", nfiles
    
  END SUBROUTINE ReadLambdaData
  
  
  !-----------------------------------------------------------------------
  ! Function : InterpolateLambda
  !
  ! Purpose  : Perform linear interpolation to compute the electron-phonon 
  !            coupling constant (lambda) for a given energy (E1).
  !
  ! Model    :
  !   - Interpolates the value of lambda between two points in the input
  !     arrays: omega_ln (logarithmic phonon frequencies) and lambdaArray
  !     (precomputed coupling constants).
  !   - If the input energy (E1) lies outside the range of omega_ln, the 
  !     function extrapolates using the nearest boundary values.
  !
  ! Inputs:
  !   E1         : Energy input for which lambda is calculated [eV].
  !   omega_ln   : Array of logarithmic phonon frequencies [eV].
  !   lambdaArray: Array of electron-phonon coupling constants corresponding
  !                to omega_ln values.
  !
  ! Returns:
  !   interpolatedLambda : Interpolated or extrapolated coupling constant.
  !
  ! Notes:
  ! - Requires omega_ln and lambdaArray to be preinitialized and sorted in 
  !   ascending order for accurate interpolation.
  ! - If E1 lies within omega_ln(i) and omega_ln(i+1), linear interpolation 
  !   is used to calculate lambda.
  ! - If E1 is outside the range of omega_ln, extrapolation is performed 
  !   using the nearest values.
  !-----------------------------------------------------------------------
  REAL(KIND=8) FUNCTION InterpolateLambda(E1, omega_ln, lambdaArray)
    IMPLICIT NONE
    
    ! Input variables
    REAL(KIND=8), INTENT(IN) :: E1      ! Energy input
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: omega_ln, lambdaArray
    
    ! Local variables
    REAL(KIND=8) :: interpolatedLambda  ! Resulting lambda value
    INTEGER :: i, nPoints               ! Loop index and size of data arrays
    
    ! Determine the number of points in the input arrays
    nPoints = SIZE(omega_ln)  ! Should equal the number of files/entries
    
    !-----------------------------------------------------------------------
    ! Perform linear interpolation
    DO i = 1, nPoints - 1
      IF (E1 >= omega_ln(i) .AND. E1 <= omega_ln(i+1)) THEN
        ! Linear interpolation formula
        interpolatedLambda = lambdaArray(i) + (lambdaArray(i+1) - lambdaArray(i)) * &
                           & (E1 - omega_ln(i)) / (omega_ln(i+1) - omega_ln(i))
        RETURN
      END IF
    END DO
    
    !-----------------------------------------------------------------------
    ! Handle extrapolation for values outside the range
    IF (E1 < omega_ln(1)) THEN
      interpolatedLambda = lambdaArray(1)        ! Use first value for low extrapolation
    ELSE IF (E1 > omega_ln(nPoints)) THEN
      interpolatedLambda = lambdaArray(nPoints)  ! Use last value for high extrapolation
    END IF
    
  END FUNCTION InterpolateLambda
  
  
  !-----------------------------------------------------------------------
  ! Function : get_tau
  ! Purpose  : Compute energy- and temperature-dependent relaxation time tau(E,T)
  !          : based on Matthiessen’s rule combining multiple scattering mechanisms:
  !            - Phonon scattering      :  tau is proportional to  1 / (|E - mu| * T)
  !            - Impurity scattering    :  tau is proportional to  |E - mu|^n   (user-defined exponent)
  !            - Electronic DOS-based   :  tau is proportional to  1 / [DOS(E) * T]
  !            - Phonon DOS-based       :  tau is proportional to  1 / integral D_ph(omega) * [1 + n(omega,T)] / omega d(omega)
  !
  ! Inputs:
  !   E1  : Energy offset from chemical potential, E - mu [eV]
  !   T   : Temperature [K]
  !   I   : Index in energy mesh corresponding to E = EE(I)
  !
  ! Returns:
  !   get_tau : Relaxation time tau(E,T) [s]
  !
  ! Notes:
  ! - All active scattering channels are combined via 1/tau_total = Sum( 1/tau_i )
  ! - Each mechanism contributes independently and is numerically stabilized
  ! - tau(E,T) affects conductivity and Seebeck integrals
  !-----------------------------------------------------------------------
  REAL(KIND=8) FUNCTION get_tau(E1, T, I)
    REAL(KIND=8), INTENT(IN) :: E1 ! Energy offset from chemical potential: E - mu [eV]
    REAL(KIND=8), INTENT(IN) :: T  ! Temperature [K]
    INTEGER, INTENT(IN) :: I       ! Energy mesh index corresponding to E1 = EE(I) - mu
    !------------------------------
    REAL(KIND=8) :: tau_ph         ! Relaxation time due to phonon scattering   : Assumes more scattering at higher temperature and away from chemical potential
    REAL(KIND=8) :: tau_imp        ! Relaxation time due to impurity scattering : Typically used in heavily doped semiconductors
    REAL(KIND=8) :: tau_dos        ! Relaxation time based on density of states : Accounts for scattering probability increasing with available final states
    REAL(KIND=8) :: tau_phdos      ! Relaxation time based on phonon density of states
    REAL(KIND=8) :: dosE           ! D(E) at E = EE(I) : Used for DOS-dependent scattering calculations
    REAL(KIND=8) :: inv_sum        ! Accumulator for inverse tau according to Matthiessen's rule
    REAL(KIND=8) :: lambda         ! Electron-phonon coupling constant (Lambda)
    !------------------------------
    REAL(KIND=8) :: tau_Umk        ! Umklapp (Klemens-Callaway type model) scattering. (works: phononDOS = T)
    REAL(KIND=8) :: tau_bound      ! Boundary scattering. (works: phononDOS = T and N_atom) (L_bound corresponds to the film thickness and grain size.)
    !------------------------------
    REAL(KIND=8) :: tau_phel       ! Phonon-Electron scattering coefficient. (works: a2Fdos = F)
    REAL(KIND=8) :: tau_pdef       ! Point defect scattering coefficient.
    !------------------------------
    
    ! === Electron-Phonon Coupling Constant (Lambda) ===
    ! - Lambda represents the strength of the electron-phonon interaction.
    ! - Quantum ESPRESSO (QE) can compute Lambda using Density Functional Perturbation Theory (DFPT).
    !   * Outputs from QE include Lambda, Eliashberg spectral function alpha^2*F(omega), and related phonon data.
    ! - ABINIT can also calculate Lambda using DFPT. 
    !   * It provides information on electron-phonon matrix elements, phonon modes, and Lambda for various q-points.
    
    ! --- Phonon scattering: tau is proportional to 1 / (|E - mu| * T)
    IF (use_phonon) THEN
       IF (ABS(E1) > 1.0D-3 .AND. T > 1.0D-3) THEN
          lambda = GetLambda(E1, T, I)  ! Interpolation function for Lambda data (ReadLambdaData available)
          IF (lambda < 0.0D0 .OR. lambda > 10.0D0) THEN
            WRITE(*, *) "Warning: Lambda out of expected range !"
            lambda = 0.5D0  ! Fallback to default value
          END IF
          tau_ph = tau0 / (ABS(E1) * T * lambda) ! Apply Lambda to the relaxation time
       ELSE
          tau_ph = 0.0D0
       END IF
    ELSE
       tau_ph = 0.0D0
    END IF
    
    ! --- Impurity scattering (disabled here): tau is proportional to energy exponent
    IF (use_impurity) THEN
       IF (ABS(E1) > 1.0D-6) THEN
          tau_imp = tau0 * ABS(E1)**n_impurity
       ELSE
          tau_imp = 0.0D0
       END IF
    ELSE
       tau_imp = 0.0D0
    END IF
    
    ! --- DOS-based scattering: tau is proportional to 1 / (DOS * T)
    IF (use_dos) THEN
       dosE = DOS(I)
       IF (dosE > 1.0D-12 .AND. T > 1.0D-3) THEN
          tau_dos = tau0 / (dosE * T)
       ELSE
          tau_dos = 0.0D0
       END IF
    ELSE
       tau_dos = 0.0D0
    END IF
    
    ! Phonon DOS-based relaxation time (global temperature-based)
    IF (use_a2Fdos) THEN
       tau_phdos = tau_a2Fdos_ET(E1, T)
    ELSE IF (use_phononDOS) THEN
       tau_phdos = tau_phdos_ET(E1, T)
    ELSE
       tau_phdos = 0.0D0
    END IF
    
    ! Umklapp (Klemens-Callaway type model) scattering (T >> R.T.)
    IF ((b_para > 0.0D0) .and. use_phononDOS) THEN
       ! Debay temperature is approximate (hbar*wmax/kb) for nomarilzed phonon dos with 3N.
       !tau_Umk = tau0 / (E1**2.0D0 * T * exp(-(hbar*wmax/kb)/(b_para*T)) + 1.0D-12)
       tau_Umk = tau0 / (E1**2.0D0 * T * exp(-Theta_D/(b_para*T)) + 1.0D-12)
    ELSE
       tau_Umk = 0.0D0
    END IF
    
    ! Boundary scattering
    IF ((L_bound > 0.0D0) .and. use_phononDOS) THEN
       !sound_velocity = wmax / (6 * PI**2 * N_atom / volume)**(1.0D0/3.0D0)
       sound_velocity = Theta_D / (6 * PI**2 * N_atom / volume)**(1.0D0/3.0D0)
       tau_bound = L_bound / sound_velocity
    ELSE
       tau_bound = 0.0D0
    END IF
    
    ! Phonon-Electron scattering (Less accurate than use_a2Fdos)
    IF ((C_phel > 0.0D0) .and. (use_a2Fdos .eqv. .FALSE.)) THEN
       tau_phel = tau0 / (C_phel * E1**2.0D0 + 1.0D-12)
    ELSE
       tau_phel = 0.0D0
    END IF
    
    ! Point defect (Mass-difference) scattering (T <= R.T.)
    IF (B_pdef > 0.0D0) THEN
       tau_pdef = tau0 / (B_pdef * E1**4.0D0 + 1.0D-12)
    ELSE
       tau_pdef = 0.0D0
    END IF
    
    ! --- Matthiessen rule: combine non-zero components
    inv_sum = 0.0D0
    IF (tau_ph     > 0.0D0) inv_sum = inv_sum + 1.0D0 / tau_ph
    IF (tau_imp    > 0.0D0) inv_sum = inv_sum + 1.0D0 / tau_imp
    IF (tau_dos    > 0.0D0) inv_sum = inv_sum + 1.0D0 / tau_dos
    IF (tau_phdos  > 0.0D0) inv_sum = inv_sum + 1.0D0 / tau_phdos
    IF (tau_Umk    > 0.0D0) inv_sum = inv_sum + 1.0D0 / tau_Umk
    IF (tau_bound  > 0.0D0) inv_sum = inv_sum + 1.0D0 / tau_bound
    IF (tau_phel   > 0.0D0) inv_sum = inv_sum + 1.0D0 / tau_phel
    IF (tau_pdef   > 0.0D0) inv_sum = inv_sum + 1.0D0 / tau_pdef
    
    IF (inv_sum > 0.0D0) THEN
       get_tau = 1.0D0 / inv_sum
    ELSE
       get_tau = tau0
    END IF
  END FUNCTION get_tau
  
  
  !---------------------------------------------------------------
  ! Subroutine: fermi_derivative
  ! Purpose   : Compute the derivative of the Fermi-Dirac distribution:
  !           :   df/dE = exp(E/kT) / [ (1 + exp(E/kT))^2 * kT ]
  !           : Also computes (E - mu) * df/dE for transport integrals.
  !
  ! Inputs    :
  !   E1    : Energy offset from chemical potential (E - mu) [eV]
  !   TEM2  : Thermal energy kT [eV]
  !
  ! Outputs   :
  !   FDD   : Derivative of Fermi-Dirac distribution df/dE
  !   FDE   : Weighted derivative: (E - mu) * df/dE
  !
  ! Notes:
  ! - At low temperature or extreme |E|, output may be zeroed for stability.
  ! - Used to calculate Seebeck integrals and energy-dependent response.
  !---------------------------------------------------------------
  SUBROUTINE fermi_derivative(E1, TEM2, FDD, FDE)
    REAL(KIND=8), INTENT(IN)  :: E1, TEM2
    REAL(KIND=8), INTENT(OUT) :: FDD, FDE
    REAL(KIND=8) :: expo

    ! Ensure temperature is not extremely close to zero.
    ! Prevents division by zero when calculating derivatives.
    IF (ABS(TEM2) > 1.0D-12) THEN
       ! Compute exponential term: e^( (E - mu) / kT )
       ! This is part of the df/dE expression, but may overflow at large values.
       expo = DEXP(E1 / TEM2)
       ! Only proceed if exponentiation is numerically stable.
       ! This avoids floating-point overflow errors.
       IF (expo < 1.0D100) THEN
          ! Full expression of the derivative of the Fermi-Dirac distribution:
          !   df/dE = e^( (E - mu)/kT ) / [ (1 + e^( (E - mu)/kT ))^2 * kT ]
          ! This sharply peaks at E = mu and captures how energy levels contribute to transport.
          FDD = expo / (1.0D0 + expo)**2 / TEM2
          ! Multiply df/dE by (E - mu) for Seebeck numerator.
          ! The integrand of the Seebeck coefficient involves this product.
          FDE = FDD * E1
       ELSE
          ! If exponential term is unstable (too large), skip this energy point.
          ! Zeroing out avoids corrupting the integral with numerical artifacts.
          FDD = 0.0D0
          FDE = 0.0D0
       END IF
    ELSE
       ! If temperature is too low, treat derivative as negligible.
       ! Transport is suppressed in this limit, so contributions are effectively zero.
       FDD = 0.0D0
       FDE = 0.0D0
    END IF
  END SUBROUTINE fermi_derivative
  
  
  !-----------------------------------------------------------------------
  ! Subroutine : read_a2F_dos
  !
  ! Purpose    : This subroutine reads and processes data from "a2F.dos" files,
  !              which contain phonon density of states (DOS) and Eliashberg 
  !              function data. It allocates arrays dynamically and processes 
  !              multiple input files to prepare for scattering calculations.
  !
  ! Model      :
  !   - Reads "a2F.dosX" files (where X is the index of the file).
  !   - Each file contains:
  !       1. Frequencies (freq)
  !       2. Total Eliashberg function (a2F_total)
  !       3. Mode-specific contributions (a2F_mode)
  !   - Converts phonon frequencies from [Rydberg] to [eV].
  !   - Determines the number of modes (nmodes = 3N), where N is the number
  !     of atoms in the system.
  !
  ! Inputs:
  !   - "a2F.dosX" files : Input files containing phonon DOS and Eliashberg 
  !                         function data for multiple configurations.
  !   - Pre-initialized omega_ln and nfiles from ReadLambdaData().
  !
  ! Outputs:
  !   - freq(nfiles, nLines)      : Array of phonon frequencies [eV].
  !   - a2F_total(nfiles, nLines) : Total Eliashberg function data.
  !   - a2F_mode(nfiles, nmodes, nLines): Mode-specific Eliashberg contributions.
  !   - ndata                     : Number of data points per file.
  !   - nmodes                    : Total phonon modes (3*N_atoms).
  !   - N_atom                    : Number of atoms in the system.
  !
  ! Notes:
  ! - The input files must follow the format used in Quantum ESPRESSO 
  !   ("PHonon/PH/matdyn.f90"). Header rows are skipped during reading.
  ! - The subroutine dynamically allocates arrays based on the number of
  !   lines and columns in the input data.
  ! - Ensure input files are named sequentially as "a2F.dos1", "a2F.dos2", etc.
  !   and match the expected format.
  !-----------------------------------------------------------------------
  SUBROUTINE read_a2F_dos()
    IMPLICIT NONE
    
    ! Local variables
    INTEGER :: iunit, ios            ! File unit and status code
    INTEGER :: nLines, numColumns    ! Line count and column count
    INTEGER :: i, j, k               ! Loop indices
    CHARACTER(LEN=100) :: filename   ! Input file name
    CHARACTER(LEN=16003) :: line     ! Buffer for reading lines
    
    !-----------------------------------------------------------------------
    ! Step 1: Determine number of input files and lines of data
    nfiles = SIZE(omega_ln)          ! Number of files from ReadLambdaData()
    
    ! Open the first file to count lines and determine data structure
    OPEN(UNIT=80, FILE='a2F.dos1', STATUS="OLD", IOSTAT=ios)
    nLines = 0
    DO
      READ(80, '(A)', IOSTAT=ios) line   ! Read until EOF
      IF (ios /= 0) EXIT
      nLines = nLines + 1
    END DO
    CLOSE(80)
    nLines = nLines - 6                  ! Exclude header rows
    WRITE(*,*) "Information for a2F.dos data"
    WRITE(*,*) "Number of data: ", nLines
    ndata = nLines                       ! Save the data point count
    
    !-----------------------------------------------------------------------
    ! Step 2: Determine number of columns and modes
    OPEN(UNIT=80, FILE='a2F.dos1', STATUS="OLD", IOSTAT=ios)
    DO j = 1, 5
      READ(80,*)                         ! Skip header rows
    END DO
    READ(80, '(A)', IOSTAT=ios) line
    DO i = 1, 16003
      ! write (ifn, '(3X,1000E16.6)') E, dos_tot, dos_a2F(1:nmodes) in PHonon/PH/matdyn.f90 (QE)
      IF (line(4+16*(i-1):4+16*i) == "") EXIT
      !WRITE(*,*) line(4+16*(i-1):4+16*i)
    END DO
    CLOSE(80)
    numColumns = i - 1                   ! Columns minus extra entries
    nmodes = numColumns - 2              ! nmodes = columns - 2 (frequency and total)
    N_atom = nmodes/3                    ! Number of atoms
    WRITE(*,*) "Number of columns (freq, a2F_total, a2F of 3N data): ", numColumns
    WRITE(*,*) "Number of nmodes = 3N = 3*N_atom: ", nmodes
    WRITE(*,*) "Number of atoms = nat (QE): ", N_atom
    WRITE(*,*) "------------------------------------------------------"
    
    !-----------------------------------------------------------------------
    ! Step 3: Allocate memory for data arrays
    ALLOCATE(freq(nfiles,nLines), a2F_total(nfiles,nLines), a2F_mode(nfiles,nmodes,nLines))
    ! Note: write (ifn, '(3X,1000E16.6)') E, dos_tot, dos_a2F(1:nmodes) in PHonon/PH/matdyn.f90 (QE)
    
    !-----------------------------------------------------------------------
    ! Step 4: Read data from all files and process
    DO i = 1, nfiles
      iunit = 80
      WRITE(filename, "('a2F.dos', I0)") i
      !WRITE(*,*) "read file: ", filename
      OPEN(UNIT=iunit, FILE=filename, STATUS="OLD", IOSTAT=ios)
      DO j = 1, 5
        READ(iunit,*)                   ! Skip header rows
      END DO
      DO j = 1, ndata
        READ(iunit,*) freq(i,j), a2F_total(i,j), (a2F_mode(i,k,j), k = 1, nmodes)
        !WRITE(*,'(3X,1000E16.6)') freq(i,j), a2F_total(i,j), (a2F_mode(i,k,j), k = 1, nmodes)
        freq(i,j) = freq(i,j) * EV      ! Convert frequencies from [Ry] to [eV]
      END DO
      CLOSE(iunit)
      ! Optionally, print the loaded data for verification
      !WRITE(*,*) "Data loaded successfully !: ", filename
    END DO
    
  END SUBROUTINE read_a2F_dos
  
  
  !---------------------------------------------------------------
  ! Subroutine: read_phonon_dos
  ! Purpose   : Read phonon density of states from 'phononDOS.dat'
  !           : Normalize the phonon DOS to satisfy the sum rule:
  !           :   integral D_ph(omega) d(omega) = 3 * N_atom
  ! Inputs    : None
  ! Outputs   : Arrays WPH(NPH), DOSPH(NPH) are populated
  !           : DOSPH normalized so that total phonon modes = 3N
  !---------------------------------------------------------------
  SUBROUTINE read_phonon_dos()
    IMPLICIT NONE
    INTEGER :: iunit, i, k, ios, nLines
    REAL(KIND=8) :: sum, omega_avg, delta, norm_factor
    CHARACTER(LEN=100) :: filename
    CHARACTER(LEN=200) :: line

    filename = "phononDOS.dat"
    iunit = 11
    
    ! --- Count lines in the file ---
    nLines = 0
    OPEN(UNIT=iunit, FILE=filename, STATUS="OLD", IOSTAT=ios)

    ! Check if file exists
    IF (ios /= 0) THEN
      WRITE(*,*) "Warning: File '", filename, "' not found or cannot be opened."
      RETURN
    END IF

    DO
      READ(iunit, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      nLines = nLines + 1
    END DO
    CLOSE(iunit)
    nLines = nLines - 2

    ! --- Allocate arrays based on file line count ---
    NPH = nLines
    ALLOCATE(WPH(NPH), DOSPH(NPH), g_eff_omega(NPH), w_mode(NPH))

    wmax = 0.0
    ! --- Read data into allocated arrays ---
    OPEN(UNIT=iunit, FILE=filename, STATUS="OLD", IOSTAT=ios)
    READ(iunit,*)
    READ(iunit,*)
    DO i = 1, NPH
      READ(iunit,*, IOSTAT=ios) WPH(i), DOSPH(i)  ! Frequency = omega [eV], phonon DOS [arb.units]
      !WRITE(*,*) WPH(i), DOSPH(i)
      IF (WPH(i) >= wmax) THEN
        wmax = WPH(i)
      END IF
      wmax = WPH(i)
    END DO
    CLOSE(iunit)

    ! Optionally, print out the loaded data
    WRITE(*,*) "Phonon DOS data loaded successfully!"
    WRITE(*,*) "Number of points: ", NPH

    ! --- normalization ---
    sum = 0.0D0
    !--------------------------------------------------
    ! Trapezoidal rule
    !DO i = 1, NPH - 1
    !   delta     = WPH(i+1) - WPH(i)
    !   omega_avg = 0.5D0 * (WPH(i+1) + WPH(i))
    !   IF (delta > 0.0D0) sum = sum + omega_avg * delta
    !END DO
    !--------------------------------------------------
    ! Simpson's Rule
    DO i = 2, NPH - 2, 2  ! Loop through pairs, ensuring the range works for odd NPH
      delta = (WPH(i+1) - WPH(i-1)) / 2.0D0  ! Calculate width (average of front and back)
      IF (delta > 0.0D0) THEN                ! Execute only if delta is positive
        sum = sum + delta / 3.0D0 * (DOSPH(i-1) + 4.0D0 * DOSPH(i) + DOSPH(i+1))
      END IF
    END DO
    !-----------------------
    ! Handle the last point separately if NPH is odd
    IF (MOD(NPH, 2) /= 0) THEN
      delta = (WPH(NPH) - WPH(NPH-1)) / 2.0D0
      IF (delta > 0.0D0) THEN
        sum = sum + delta / 3.0D0 * (DOSPH(NPH-1) + DOSPH(NPH))
      END IF
    END IF
    !--------------------------------------------------

    norm_factor = 3.0D0 * N_atom / (sum + 1.0D-12)   ! 3N = 3 * N_atom

    DO i = 1, NPH
       DOSPH(i) = DOSPH(i) * norm_factor
    END DO

    ! Initialize mode-dependent electron-phonon coupling strength and selection filter
    ! - g_eff_omega(i): default 1.0 for Eliashberg function approximation
    ! - w_mode(i)     : default 1.0 for uniform mode activation (used in kernel filtering)
    DO i = 1, NPH
       g_eff_omega(i) = 1.0D0  ! For Transition to Eliashberg functions
       w_mode(i) = 1.0D0       ! For Selection rule filtering (symmetry and modal analysis)
    END DO

    ! Optional selection rule filtering
    IF (use_selection_filter) THEN
       CALL generate_w_mode()  ! use w_mode, WPH
    END IF

    CLOSE(iunit)
  END SUBROUTINE read_phonon_dos
  
  
  !-----------------------------------------------------------------------
  ! Subroutine : generate_w_mode
  !
  ! Purpose    : Generate a weighting factor (w_mode) for phonon modes based 
  !              on their frequencies (WPH). This helps categorize modes 
  !              as acoustic, intermediate, or optical for scattering or 
  !              dynamical calculations.
  !
  ! Model      :
  !   - Phonon frequencies (omega) are divided into categories based on 
  !     thresholds to determine their contribution weight:
  !       1. Acoustic modes: Strong contribution (w_mode = 1.0D0).
  !       2. Intermediate modes: Moderate contribution (w_mode = 0.5D0).
  !       3. Optical modes: Suppressed contribution (w_mode = 0.1D0).
  !
  ! Filtering Logic:
  !   - Additional refined filtering logic can be applied for more 
  !     precise control of weights based on frequency ranges.
  !
  ! Inputs:
  !   WPH  : Array of phonon frequencies [Rydberg].
  !
  ! Outputs:
  !   w_mode : Array of weights corresponding to phonon frequencies.
  !
  ! Notes:
  ! - The subroutine is designed to handle varying size arrays dynamically.
  ! - Frequency thresholds (e.g., 0.01D0, 0.03D0) can be adjusted depending 
  !   on the system or application.
  ! - Refined logic is included as optional code for more detailed filtering.
  !-----------------------------------------------------------------------
  SUBROUTINE generate_w_mode() ! use w_mode, WPH
    IMPLICIT NONE
    
    ! Local variables
    INTEGER :: j               ! Loop index
    REAL(KIND=8) :: omega      ! Frequency for current mode
    !REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: WPH
    !REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: w_mode

    !-----------------------------------------------------------------------
    ! Loop through all phonon frequencies and assign weights
    DO j = 1, SIZE(WPH)
       omega = WPH(j)

       ! --- Basic filtering logic for phonon frequencies ---
       IF (omega < 0.01D0*EV) THEN
          w_mode(j) = 1.0D0          ! Strong weight for acoustic modes
       ELSE IF (omega < 0.03D0*EV) THEN
          w_mode(j) = 0.5D0          ! Moderate weight for intermediate modes
       ELSE
          w_mode(j) = 0.1D0          ! Suppressed weight for optical modes
       END IF
       
       ! --- Refined filtering logic for precise control ---
       ! Uncomment the block below for more granular filtering:
       !IF (omega < 0.005D0*EV) THEN
       !  w_mode(j) = 1.0D0             ! Strong weight (low-frequency acoustic modes)
       !ELSE IF (omega < 0.01D0*EV) THEN
       !  w_mode(j) = 0.8D0             ! Slightly reduced weight
       !ELSE IF (omega < 0.02D0*EV) THEN
       !  w_mode(j) = 0.5D0             ! Moderate weight for mid-frequency modes
       !ELSE IF (omega < 0.04D0*EV) THEN
       !  w_mode(j) = 0.3D0             ! Suppressed weight for higher frequencies
       !ELSE
       !  w_mode(j) = 0.1D0             ! Minimal contribution for optical modes
       !END IF
    END DO
  END SUBROUTINE generate_w_mode

  !---------------------------------------------------------------
  ! Subroutine: compute_Cv_DOS
  ! Purpose: Compute heat capacity from phonon DOS at temperature T
  !---------------------------------------------------------------
  function compute_Cv_DOS(T) result(Cv)
    implicit none
    real(8), intent(in) :: T
    real(8) :: Cv, x, integrand, d_omega
    integer :: i

    Cv = 0.0d0
    d_omega = WPH(2) - WPH(1)  ! [eV]

    do i = 1, NPH
       ! x = hbar * WPH(i) / (kb * T)
       x = WPH(i) / (kb * T)
       ! Immediately after reading with read_phonon_dos(), the sum of DOSPH(i) is 3N = 3*N_atom [dimensionless].
       integrand = x**2 * exp(x) / ((exp(x) - 1.0d0)**2  + 1.0e-12) * DOSPH(i)
       IF (i > 1) THEN
         d_omega = MAX(WPH(i) - WPH(i-1), 1.0D-12)
       ELSE
         d_omega = WPH(1)
       END IF
       Cv = Cv + integrand * d_omega
    end do

    ! 1 [eV/K] * 1.60218e-19 [C] * 6.02214e23 [1/mol] = 9.6485e4 [J/(mol K)]
    Cv = Cv * kb * 9.6485e4
  end function compute_Cv_DOS

  !---------------------------------------------------------------
  ! Function: compute_Cv_Debye
  ! Purpose: Compute Debye heat capacity at temperature T and Debye temp Theta_D
  !---------------------------------------------------------------
  function compute_Cv_Debye(T, Theta_D) result(Cv)
    implicit none
    real(8), intent(in) :: T, Theta_D
    real(8) :: Cv, x, integrand, dx
    integer :: i, n_int

    n_int = 2500
    dx = Theta_D / T / n_int
    Cv = 0.0d0

    do i = 1, n_int
       x = i * dx
       integrand = x**4 * exp(x) / ((exp(x) - 1.0d0)**2 + 1.0e-12)
       Cv = Cv + integrand * dx
    end do
    
    ! 1 [eV/K] * 1.60218e-19 [C] * 6.02214e23 [1/mol] = 9.6485e4 [J/(mol K)]
    Cv = 9.0D0 * kb * (T / Theta_D)**3.0D0 * Cv * 9.6485e4
  end function compute_Cv_Debye

  !---------------------------------------------------------------
  ! Subroutine: find_matching_Theta_D
  ! Purpose: Find Debye temperature that matches Cv_DOS using bisection
  !---------------------------------------------------------------
  subroutine find_matching_Theta_D(T, Theta_D_match, Cv_DOS_out, Cv_Debye_out)
    implicit none
    real(8), intent(in) :: T
    real(8), intent(out) :: Theta_D_match, Cv_DOS_out, Cv_Debye_out
    real(8) :: Theta_D_low, Theta_D_high, Theta_D_mid, Cv_Debye_mid, Cv_DOS
    real(8), parameter :: tol = 1.0d-3
    integer :: i, n_divisions
    real(8) :: step_size, closest_diff

    Theta_D_low = wmax/kb*0.25D0
    Theta_D_high = wmax/kb*1.0D0
    Cv_DOS = compute_Cv_DOS(T)
  
    n_divisions = 2500
    step_size = (Theta_D_high - Theta_D_low) / real(n_divisions)
    closest_diff = HUGE(1.0D0) ! Initialize with a large number

    do i = 0, n_divisions
      Theta_D_mid = Theta_D_low + i * step_size
      Cv_Debye_mid = compute_Cv_Debye(T, Theta_D_mid)
      !WRITE(*,*) Theta_D_mid, Cv_Debye_mid, Cv_DOS
      if (abs(Cv_Debye_mid - Cv_DOS) < closest_diff) then
         closest_diff = abs(Cv_Debye_mid - Cv_DOS)
         Theta_D_match = Theta_D_mid
      end if
    end do

    Cv_DOS_out = Cv_DOS
    Cv_Debye_out = compute_Cv_Debye(T, Theta_D_match)
  end subroutine find_matching_Theta_D

END MODULE seebeck_data



!---------------------------------------------------------------
! PROGRAM: seebeck_analysis
! Purpose:
!   Analyze Seebeck coefficient as a function of temperature
!   using energy-resolved transport data and relaxation time models.
! Inputs:
!   - AKK.DATA       : Band-resolved energy, velocity^2, velocity^2*DOS, DOS, Cumulative (DOS)
!   - apot.data      : Chemical potentials (mu) vs. temperature
!   - phononDOS.dat  : Optional phonon DOS for tau(T) model
! Outputs:
!   - Seebeck_analysis.dat : Columns of T, mu, <E - mu>, dmu/dT * 10^6
!---------------------------------------------------------------
PROGRAM seebeck_analysis
  USE seebeck_data
  IMPLICIT NONE

  INTEGER :: I    ! Energy mesh index (used to scan EE(I) or DOS(I))
  INTEGER :: IT   ! Temperature mesh index (for looping T-dependent quantities)
  INTEGER :: LE   ! Loop index for energy-related integration or averaging (e.g., LE = 1 to MM)
  INTEGER :: L    ! Generic loop index (can be reused locally in nested structures)
  
  REAL(KIND=8) :: EE1, GV21, GV2D1, DOS1  ! Temporary variables for reading input files
  REAL(KIND=8) :: TEM, TEM2, CP, T, T1  ! TEM: temperature [K], TEM2: kT [eV], CP: chemical potential
  REAL(KIND=8) :: T_hole, T_electron    ! Conductivity (hole, electron)
  REAL(KIND=8) :: E, E1, DE             ! E: energy value, E1: E - mu, DE: differential dE
  REAL(KIND=8) :: FD, FX, FDD, FDE      ! FD: Fermi function, FDD: derivative, FDE: E-weighted
  REAL(KIND=8) :: FFOS, FDEE            ! Numerators/denominators for Seebeck integrals
  
  ! Trapezoidal rule or Simpson's Rule
  REAL(KIND=8) :: FFOS_PREV, FFOS_CURR, FFOS_NEXT             ! Conductivity contributions
  REAL(KIND=8) :: FDEE_PREV, FDEE_CURR, FDEE_NEXT
  REAL(8) :: FDD_PREV, FDD_NEXT
  REAL(8) :: FDE_PREV, FDE_NEXT
  REAL(8) :: m_eff_PREV
  REAL(8) :: m_eff_CURR
  REAL(8) :: m_eff_NEXT
  REAL(8) :: MFP_PREV
  REAL(8) :: MFP_CURR
  REAL(8) :: MFP_NEXT
  REAL(8) :: FD_PREV
  REAL(8) :: FD_CURR
  REAL(8) :: FD_NEXT
  REAL(8) :: DOS_PREV
  REAL(8) :: DOS_CURR
  REAL(8) :: DOS_NEXT
  REAL(8) :: kappa_electron_term1_PREV
  REAL(8) :: kappa_electron_term1_CURR
  REAL(8) :: kappa_electron_term1_NEXT
  REAL(8) :: kappa_electron_term2_PREV
  REAL(8) :: kappa_electron_term2_CURR
  REAL(8) :: kappa_electron_term2_NEXT
  
  ! --- Variable Initialization (Simpson's Rule) ---
  REAL(8) :: E0, E2                           ! Energy points
  REAL(8) :: FDE0, FDE1, FDE2                 ! Energy-dependent derivatives
  REAL(8) :: FDD0, FDD1, FDD2                 ! Derivatives of Fermi-Dirac distribution
  INTEGER :: LOOP_ODD
  
  ! Carrier concentration
  REAL(8) :: m_eff = 0.26D0          ! Effective mass (relative to m_e)
  REAL(8) :: alpha = 1.0D-20         ! Scattering coefficient (material dependent)
  REAL(8) :: beta = 0.05D0           ! Temperature spread coefficient (adjusted based on experimental values)
  REAL(8) :: T0 = 300.0D0            ! Standard temperature (e.g. room temperature 300 K)
  REAL(8) :: gamma = 0.01D0          ! Temperature dependence coefficient
  REAL(8) :: C = 1.0D0               ! Density of states constant
  REAL(8) :: E_C = 1.12D0            ! Band edge energy [eV]
  REAL(8) :: P = 10.0D0              ! Light intensity
  REAL(8) :: E_ph = 1.5D0            ! Photon energy [eV]
  REAL(8) :: N_sat = 1.0D19          ! Saturation doping concentration [cm^-3]
  REAL(8) :: Delta_E_gap = 2.5D-4    ! Bandgap narrowing coefficient [eV/K]
  REAL(8) :: pinningShift = 0.05D0   ! Fermi level pinning shift [eV]
  
  ! Lattice parameter
  REAL(8) :: LA, LB, LC, alpha_rad, beta_rad, gamma_rad
  
  REAL(8) :: Ln = 2.44E-8            ! Lorenz number [W*Ohm/K^2]
  REAL(8) :: MFP                     ! meam free path [A]
  REAL(8) :: MFP_hole                ! meam free path [A]
  REAL(8) :: MFP_electron            ! meam free path [A]
  REAL(8) :: m_eff_hole              ! m_eff = 2*E/vg^2 <- E = (1/2)*m*v^2
  REAL(8) :: m_eff_electron          ! m_eff = 2*E/vg^2 <- E = (1/2)*m*v^2
  REAL(8) :: DOS_hole
  REAL(8) :: DOS_electron
  
  REAL(8) :: temperature             ! T [K]
  REAL(8) :: chemical_potential      ! mu [eV]
  REAL(8) :: mean_energy             ! <E-mu> [eV]
  REAL(8) :: seebeck_coefficent      ! S [muV/K]
  REAL(8) :: conductivity            ! s [S/m]
  REAL(8) :: conductivity_hole       ! s [S/m]
  REAL(8) :: conductivity_electron   ! s [S/m]
  REAL(8) :: resistance              ! R [Ohm m]
  REAL(8) :: power_factor            ! PF [W/m/K^2]
  REAL(8) :: kappa_electron          ! ke [W/m/K]
  REAL(8) :: kappa_electron_term1    ! ke [W/m/K]
  REAL(8) :: kappa_electron_term2    ! ke [W/m/K]
  REAL(8) :: mean_free_path          ! MFP [A]
  REAL(8) :: mean_free_path_hole     ! MFP of hole [A]
  REAL(8) :: mean_free_path_electron ! MFP of electron [A]
  REAL(8) :: carrier_concentration   ! Nc [cm^-3]
  REAL(8) :: carrier_concentration_hole     ! Nc [cm^-3]
  REAL(8) :: carrier_concentration_electron ! Nc [cm^-3]
  REAL(8) :: hall_coefficient        ! RH [m^3/C]
  REAL(8) :: mobility_hole           ! M [cm^2/V/s]
  REAL(8) :: mobility_electron       ! M [cm^2/V/s]
  REAL(8) :: effective_mass_hole     ! Meff [kg/kg]
  REAL(8) :: effective_mass_electron ! Meff [kg/kg]
  REAL(8) :: specific_heat           ! Cv [J/(mol K)] (Specific heat at constant volume)
  REAL(8) :: ZT                      ! ZT
  REAL(8) :: A_T                     ! The numerator A(T) in Fig. 12
  REAL(8) :: B_T                     ! The denominator B(T)
  
  CHARACTER(LEN=318), PARAMETER :: hdr = "#  T [K]   mu [eV]    <E-mu> [eV] &
    & S [muV/K]    s_all [S/m]  s_hole [S/m] s_elec [S/m] R [Ohm m]    PF [W/m/K^2] ke [W/m/K]  &
    & MFP_hole [A] MFP_elec [A] Nc [cm^-3]   Nc_h [cm^-3] Nc_e [cm^-3] RH [m^3/C]  &
    & Mh[cm^2/V/s] Me[cm^2/V/s] Meffh[kg/kg] Meffe[kg/kg] Cv[J/(molK)]&
    & kp [W/m/K]   ZT           A(T)         B(T)"
  
  ! ------------------------------------------------------------------
  ! Step 0: Load "DEF(Energy shift offset)" data from 'parameter.txt'
  ! ------------------------------------------------------------------
  OPEN(UNIT=90, FILE='parameter.txt', STATUS='OLD', ACTION='READ')
  READ(90, *)
  READ(90, '(25X, F12.8)') DEF                    ! Read DEF(Energy shift offset) [eV]
  READ(90, '(25X, E12.6)') tau0                   ! Read Base relaxation time [s]
  READ(90, '(25X, L10)')   use_phonon             ! Read Phonon status
  READ(90, '(25X, L10)')   use_impurity           ! Read Impurity status
  READ(90, '(25X, F6.3)')  n_impurity             ! Read Number of impurities
  READ(90, '(25X, L10)')   use_dos                ! Read DOS status
  READ(90, '(25X, L10)')   use_a2Fdos             ! Read a2F.dos files and lambda file (automatically set phononDOS = F)
  READ(90, '(25X, L10)')   use_phononDOS          ! Read PhononDOS status
  READ(90, '(25X, I6)')    N_atom                 ! Read Number of atoms
  READ(90, '(25X, L10)')   use_selection_filter   ! Read Filter status
  READ(90, '(25X, E12.6)') L_bound                ! Boundary scattering. (works: phononDOS = T and N_atom)
  READ(90, '(25X, E12.6)') b_para                 ! Umklapp (Klemens-Callaway type model) scattering. (Ref. 1 - 2) (works: phononDOS = T)
  READ(90, '(25X, E12.6)') C_phel                 ! Phonon-Electron scattering coefficient. (works: a2Fdos = F): (Ref. 2.31e10)
  READ(90, '(25X, E12.6)') B_pdef                 ! Point defect scattering coefficient. (Ref. 5.33e15)
  READ(90, *)
  READ(90, '(25X, E12.6)') Bulk_modulus           ! Bulk_modulus [GPa] (1 [eV/A^3] = 160.2 [GPa]), B = K
  READ(90, '(25X, E12.6)') dB_per_dP              ! dB/dP from EOS
  READ(90, '(25X, E12.6)') V0                     ! Volume from EOS
  READ(90, '(25X, E12.6)') Shear_modulus          ! Shear_modulus [GPa] = 3*Bulk_modulus*(1-2*Poisson_ratio) / (2*(1+Poisson_ratio)), G
  READ(90, '(25X, E12.6)') dG_per_dV              ! dG/dV from EOS
  READ(90, '(25X, E12.6)') Young_modulus          ! Young_modulus [GPa] = 9*B*G/(3*B+G), E
  READ(90, '(25X, E12.6)') Poisson_ratio          ! Poisson_ratio = (3*B-2*G)/(2*(3*B+G))
  READ(90, '(25X, E12.6)') density                ! read [g/cm^3] unit -> density * 1000 [kg/m^3]
  READ(90, '(25X, E12.6)') tau0_phonon            ! For phonons, it is about 10-100 times stronger than for electrons. (If it is 0.0, tau_ph * 100)
  READ(90, '(25X, E12.6)') Gruneisen_parameter    ! If it is 0.0, calculate from Bulk modulus and Poisson's ratio (or Shear modulus).
  READ(90, '(25X, E12.6)') Apara                  ! A parameter of Slack model
  READ(90, '(25X, L10)')   use_Apara_gamma        ! T: use Ref:[17] or F: original
  READ(90, *)
  READ(90, '(25X, E12.6)') Nd                     ! Read Doping concentration [cm^-3]
  READ(90, '(25X, E12.6)') m_eff                  ! Effective mass (relative to m_e)
  READ(90, '(25X, E12.6)') alpha                  ! Scattering coefficient (material dependent)
  READ(90, '(25X, E12.6)') beta                   ! Temperature spread coefficient (adjusted based on experimental values)
  READ(90, '(25X, E12.6)') T0                     ! Standard temperature (e.g. room temperature 300 [K])
  READ(90, '(25X, E12.6)') gamma                  ! Temperature dependence coefficient
  READ(90, '(25X, E12.6)') C                      ! Density of states constant
  READ(90, '(25X, E12.6)') E_C                    ! Band edge energy [eV]
  READ(90, '(25X, E12.6)') P                      ! Light intensity
  READ(90, '(25X, E12.6)') E_ph                   ! Photon energy [eV]
  READ(90, '(25X, E12.6)') N_sat                  ! Saturation doping concentration [cm^-3]
  READ(90, '(25X, E12.6)') Delta_E_gap            ! Bandgap narrowing coefficient [eV/K]
  READ(90, '(25X, E12.6)') pinningShift           ! Fermi level pinning shift [eV]
  CLOSE(90)
  
  WRITE(*,*)
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "-----   List of loaded data from parameter.txt   -----"
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "DEF(Energy shift)   [eV]:", DEF
  WRITE(*,*) "Base relaxation time [s]:", tau0
  WRITE(*,*) "tau, phonon             :", use_phonon
  WRITE(*,*) "tau, impurity           :", use_impurity
  WRITE(*,*) "n_impurity              :", n_impurity
  WRITE(*,*) "tau, dos                :", use_dos
  WRITE(*,*) "tau, a2Fdos             :", use_a2Fdos
  WRITE(*,*) "tau, phononDOS          :", use_phononDOS
  IF (use_a2Fdos .and. use_phononDOS) THEN
    WRITE(*,*) "a2F.dos takes precedence over phononDOS.dat in calculating relaxation times."
  END IF
  WRITE(*,*) "N_atom                  :", N_atom
  WRITE(*,*) "tau, filter             :", use_selection_filter
  WRITE(*,*) "L_bound                 :", L_bound
  WRITE(*,*) "b_para                  :", b_para
  WRITE(*,*) "C_phel                  :", C_phel
  WRITE(*,*) "B_pdef                  :", B_pdef
  WRITE(*,*) "----- Thermal conductivity from equation of state: optional -----"
  WRITE(*,*) "Bulk modulus, B    [GPa]:", Bulk_modulus
  WRITE(*,*) "B' (=dB/dP)             :", dB_per_dP
  WRITE(*,*) "V0                 [A^3]:", V0
  WRITE(*,*) "Shear modulus, G   [GPa]:", Shear_modulus
  WRITE(*,*) "G' (=dG/dV)             :", dG_per_dV
  WRITE(*,*) "Young's modulus, E [GPa]:", Young_modulus
  WRITE(*,*) "Poisson's ratio         :", Poisson_ratio
  Poisson_ratio_read_value = Poisson_ratio
  IF (Poisson_ratio <= 0.0) THEN
    Poisson_ratio = 0.30D0
    WRITE(*,*) "(Automatically setting) Poisson_ratio:", Poisson_ratio
  END IF
  WRITE(*,*) "density         [g/cm^3]:", density
  WRITE(*,*) "relaxation time (ph) [s]:", tau0_phonon
  IF (tau0_phonon < 0.0) THEN
    tau0_phonon = tau0 * 100.0D0
    WRITE(*,*) "(Automatically setting) Base relaxation time (phonon) [s]:"
    WRITE(*,*) "  Base relaxation time (electron) [s] * 100:", tau0_phonon
  END IF
  WRITE(*,*) "Gruneisen parameter     :", Gruneisen_parameter
  WRITE(*,*) "Slack model parameter A :", Apara
  WRITE(*,*) "Apara_flag (Ref. [15])  :", use_Apara_gamma
  WRITE(*,*) "----- Carrier concentration (Nc) calclation: optional -----"
  WRITE(*,*) "Nd (Doping conc.)[cm^-3]:", Nd
  WRITE(*,*) "Electron Effective Mass :", m_eff
  WRITE(*,*) "alpha (Scattering coef.):", alpha
  WRITE(*,*) "beta (Temp. spread coef):", beta
  WRITE(*,*) "T0 (Standard temp.)  [K]:", T0
  WRITE(*,*) "gamma (Temp. depe. coef):", gamma
  WRITE(*,*) "C (DOS constant)        :", C
  WRITE(*,*) "E_C (Band edge ener)[eV]:", E_C
  WRITE(*,*) "P (Light intensity)     :", P
  WRITE(*,*) "E_ph (Photon energy)[eV]:", E_ph
  WRITE(*,*) "N_sat (Satu.)    [cm^-3]:", N_sat
  WRITE(*,*) "Delta E gap       [eV/K]:", Delta_E_gap
  WRITE(*,*) "Pinning Shift       [eV]:", pinningShift
  WRITE(*,*) "------------------------------------------------------"

  WRITE(*,*)
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "----- List of loaded data and caluclated values  -----"
  WRITE(*,*) "------------------------------------------------------"
  
  ! ------------------------------------------------------------------
  OPEN(UNIT=91, FILE='wien.struct', STATUS='OLD', ACTION='READ')
  READ(91, '(A)')
  READ(91, '(A)')
  READ(91, '(A)')
  READ(91, *) LA, LB, LC, alpha_rad, beta_rad, gamma_rad
  volume = (LA*B2A) * (LB*B2A) * (LC*B2A) * SQRT(1.0 &
       & - COS(alpha_rad*(PI/180.0))**2 - COS(beta_rad*(PI/180.0))**2 - COS(gamma_rad*(PI/180.0))**2 &
       & + 2.0 * COS(alpha_rad*(PI/180.0)) * COS(beta_rad*(PI/180.0)) * COS(gamma_rad*(PI/180.0)))
  WRITE(*,*) "Volume             [A^3]: ", volume
  !WRITE(6,'(A, f12.5)') " Volume [A^3]: ", volume
  WRITE(*,*) "------------------------------------------------------"
  CLOSE(91)
  
  ! ------------------------------------------------------------------
  ! Step 1: Load "Electron-phonon coupling" data from 'lambda'
  ! ------------------------------------------------------------------
  !IF (use_phonon .or. use_phononDOS .or. use_a2Fdos) THEN
  IF (use_phonon .or. use_a2Fdos) THEN
    WRITE(*,*) "Electron-phonon coupling constant, lambda"
    WRITE(*,*) "Broadening  lambda  dos(Ef)  omega_ln[K]  omega_ln[eV]"
    CALL ReadLambdaData()
    WRITE(*,*) "------------------------------------------------------"
  END IF
  
  ! ------------------------------------------------------------------
  ! Step 2: Load chemical potentials from 'apot.data'
  ! Each line contains temperature (TEM) and chemical potential (AMU)
  ! ------------------------------------------------------------------
  OPEN(UNIT=4, FILE='apot.data')
  DO L = 1, 25
     READ(4,'(F10.1,2E16.7)') TT(L), AMU(L)
  END DO
  CLOSE(4)
  
  ! ------------------------------------------------------------------
  ! Step 3: Load band velocity data from 'AKK.DATA'
  ! This file contains energy mesh (EE1), squared velocity, squared velocity*DOS, DOS, and Cumulative DOS
  ! ------------------------------------------------------------------
  OPEN(UNIT=10, FILE='AKK.DATA')
  DO LE = 1, MM
     READ(10,'(4E15.8,E15.8)') EE1, GV21, GV2D1, DOS1, EN(LE)
     EE(LE)   = EE1 - DEF    ! Apply energy shift if needed
     GV2(LE)  = GV21         ! Load squared velocity
     GV2D(LE) = GV2D1        ! Load squared velocity * DOS
     DOS(LE)  = DOS1         ! Load density of states
     IF(LE > 2 .AND. EE(LE-1) <= 0.0 .AND. EE(LE) > 0.0) THEN
       VEC =  EN(LE-1)+(EN(LE) - EN(LE-1))/(EE(LE) - EE(LE-1))
       WRITE(*,*) "The valence electron concentration (VEC):", VEC
       WRITE(*,*) "at DEF =", DEF, "[eV] (for EF = 0.0 -> EF = DEF) and 0 [K]"
       WRITE(*,*) "-------------------------------"
     END IF
     IF(LE > 2 .AND. EE(LE-1) + DEF <= 0.0 .AND. EE(LE) + DEF > 0.0) THEN
       VEC0 =  EN(LE-1)+(EN(LE) - EN(LE-1))/(EE(LE) - EE(LE-1))
       WRITE(*,*) "The valence electron concentration (VEC0):", VEC0
       WRITE(*,*) "at DEF =", 0.0, "[eV] and 0 [K]"
       WRITE(*,*) "------------------------------------------------------"
     END IF
  END DO
  IF (Nd == 0.0) THEN
    Nd = (VEC0 - VEC) / volume * 1.0e24
    WRITE(*,*) "Since Nd = 0.0, estimate Nd. The sign of Nd is positive = hole, negative = electron."
  END IF
  WRITE(*,*) "Carrier concentration, Nc [cm^-3] (at 0 [K]):", Nd ! Nd = Nc for non-condition
  WRITE(*,*) "Hall coefficient, RH [m^3/C]", (1.0D0 / (ech * Nd * 1.0e6))
  WRITE(*,*) "------------------------------------------------------"
  CLOSE(10)
  
  ! ------------------------------------------------------------------
  ! Calculates estimated values of various physical properties
  ! ------------------------------------------------------------------
  IF (Bulk_modulus > 0.0 .and. density > 0.0) THEN
    IF (Shear_modulus <= 0.0) THEN
      Shear_modulus = 3.0D0*Bulk_modulus*(1.0D0 - 2.0*Poisson_ratio) / (2.0D0*(1.0D0 + Poisson_ratio))
    END IF
    ! 1 [GPa] = 1.0e9 [kg/(m*s^2)], 1 [g/cm^3] = 1000 [kg/m^3], [GPa]*[kg/m^3] = 1.0e9 [(m/s)^2]
    vl = ((Bulk_modulus*1.0D9 + (4.0D0/3.0D0)*Shear_modulus*1.0D9) / (density*1.0D3))**(1.0D0/2.0D0)
    vt = (Shear_modulus*1.0D9 / (density*1.0D3))**(1.0D0/2.0D0)
    WRITE(*,*) "Shear modulus, G [GPa]", Shear_Modulus
    WRITE(*,*) "sound velocity (longitudinal wave) estimated from bulk and shear moduli, vl [m/s]:", vl
    WRITE(*,*) "sound velocity (transverse   wave) estimated from shear modulus        , vt [m/s]:", vt
    IF (Poisson_ratio_read_value <= 0.0) THEN
      WRITE(*,*) "Use vs instead of vt. (i.e., vt = vs)"
      vs = 0.87*vl
      vt = vs
      WRITE(*,*) "Shear wave velocity estimated from 0.87*vl. (i.e., vt = vs = 0.87*vl)  , vs [m/s]:", vs
    END IF
    WRITE(*,*) "-------------------------------"
    
    MPF_phonon = (vl + 2.0D0*vt)/3.0D0 * tau0_phonon * 1.0D10
    WRITE(*,*) "mean free path of phonon [A]:", MPF_phonon
    WRITE(*,*) "-------------------------------"
    
    WRITE(*,*)
    WRITE(*,*) "-------------------------------"
    WRITE(*,*) "Empirical estimation of thermal conductivity (Ref. https://github.com/houzf/empirical_thermal_conductivity)"
    WRITE(*,*) "-------------------------------"
    WRITE(*,*) "Clarke mode: Its main feature is that it estimates thermal conductivity by taking into account"
    WRITE(*,*) " the elastic modulus and average atomic volume. This model is sensitive to"
    WRITE(*,*) " the density and elastic properties of the material, and is particularly suitable for"
    WRITE(*,*) " highly crystalline materials."
    WRITE(*,*)
    !
    IF (Young_modulus <= 0.0) THEN
      IF (Poisson_ratio_read_value <= 0.0) THEN
        Young_modulus = 9.0D0 * Bulk_modulus * Shear_modulus / (3.0*Bulk_modulus + Shear_modulus)
      ELSE
        Young_modulus = 9.0D0 * Bulk_modulus * (1.0D0 - Poisson_ratio) / &
          & ( (1.0D0 + Poisson_ratio) * (1.0D0 - 2*0D0*Poisson_ratio) )
      END IF
      WRITE(*,*) "(Automatically setting) Young_modulus [GPa]:", Young_modulus
    END IF
    !
    ! 1 [GPa] = 1.0e9 [kg/(m*s^2)], 1 [g/cm^3] = 1000 [kg/m^3], [kg/(m*s^2)]/[kg/m^3] = [(m/s)^2]
    ! [eV/K] * [1/m^2] * [m/s] = [(eV/s)/(m*K)] = 1.602e-19 [W/(m*K)] = ech [W/(m*K)]
    ! kb [eV/K], ech [C], kb * ech [J/K] = [W*s/K]
    kappa_phonon_min = 0.87 * (kb * ech) * ((volume*1.0D-30)/N_atom)**(-2.0D0/3.0D0) * &
      & (Young_modulus*1.0D9/(density*1.0D3))**(1.0D0/2.0D0)
    WRITE(*,*) "Clarke model: kappa_phonon_min [W/m/K]:", kappa_phonon_min
    kappa_phonon_min_Clarke = kappa_phonon_min
    
    WRITE(*,*)
    va = (Young_modulus*1.0D9/(density*1.0D3))**(1.0D0/2.0D0)
    WRITE(*,*) "average sound velocity, va [m/s]:", va  ! The Harmonic Mean of velocities.
    !
    Theta_D_va = (2.0D0*PI*hbar/kb) * ( (3.0D0*N_atom) / (4.0D0*PI*(volume*1.0D-30)) )**(1.0D0/3.0D0) * va
    WRITE(*, *) "Debye temperature, Theta_D_va [K]:", Theta_D_va
    
    !Theta_D_Cezairliyan_equ = Theta_D_va
    
    !WRITE(*,*)
    !WRITE(*,*) "Ref.: Thermal Conductivity of the Elements: https://srd.nist.gov/jpcrdreprint/1.3253100.pdf"
    !WRITE(*,*) "Cezairliyan: k/km= [(1/3)*(T/Tm)^2 + 2/(3*(T/Tm))]^-1"
    !WRITE(*,*) "Let Tm = Debye temperature and the thermal conductivity at that time be km."
    !km = kappa_phonon_min
    !WRITE(*,*) "Clarke model, km:" , km, " at ", Theta_D_Cezairliyan_equ, " [K]"
    !TEM = 300.0
    !WRITE(*,*) "Clarke model: kappa_phonon_min [W/m/K]:",&
    !  & km * ( (1.0D0/3.0D0)*(TEM/Theta_D_Cezairliyan_equ)**2.0D0 + 2.0D0/(3.0D0*(TEM/Theta_D_Cezairliyan_equ)) )**(-1.0D0) ,&
    !  & " at ", TEM, " [K]"
    
    WRITE(*,*)
    WRITE(*,*) "---------- ----------"
    
    WRITE(*,*) "Cahill-Pohl model: The model emphasizes the speed of sound waves and calculates how fast longitudinal and"
    WRITE(*,*) " transverse waves travel. This model is flexible enough to be applied to non-crystalline materials and"
    WRITE(*,*) " nanoscale structures."
    WRITE(*,*)
    !
    ! kb [eV/K], ech [C], kb * ech [J/K] = [W*s/K]
    kappa_phonon_min = (1.0D0/2.0D0) * (PI/6.0D0)**(1.0D0/3.0D0) * (kb * ech) * &
      & (N_atom/(volume*1.0D-30))**(2.0D0/3.0D0) * (vl + 2.0D0*vt)
    WRITE(*,*) "Cahill-Pohl model: kappa_phonon_min [W/m/K]:", kappa_phonon_min
    kappa_phonon_min_Cahill = kappa_phonon_min
    
    WRITE(*,*)
    va = ( (1.0D0/3.0D0) * (1.0D0/vl**3.0D0 + 2.0D0/vt**3.0D0) )**(-1.0D0/3.0D0)
    WRITE(*,*) "average sound velocity, va [m/s]:", va  ! The Harmonic Mean of velocities.
    !
    Theta_D_va = (2.0D0*PI*hbar/kb) * ( (3.0D0*N_atom) / (4.0D0*PI*(volume*1.0D-30)) )**(1.0D0/3.0D0) * va
    WRITE(*, *) "Debye temperature, Theta_D_va [K]:", Theta_D_va
    
    !Theta_D_Cezairliyan_equ = Theta_D_va
    
    !WRITE(*,*)
    !WRITE(*,*) "Ref.: Thermal Conductivity of the Elements: https://srd.nist.gov/jpcrdreprint/1.3253100.pdf"
    !WRITE(*,*) "Cezairliyan: k/km= [(1/3)*(T/Tm)^2 + 2/(3*(T/Tm))]^-1"
    !WRITE(*,*) "Let Tm = Debye temperature and the thermal conductivity at that time be km."
    !km = kappa_phonon_min
    !WRITE(*,*) "Cahill-Pohl model, km:" , km, " at ", Theta_D_Cezairliyan_equ, " [K]"
    !TEM = 300.0
    !WRITE(*,*) "Cahill-Pohl model: kappa_phonon_min [W/m/K]:",&
    !  & km * ( (1.0D0/3.0D0)*(TEM/Theta_D_Cezairliyan_equ)**2.0D0 + 2.0D0/(3.0D0*(TEM/Theta_D_Cezairliyan_equ)) )**(-1.0D0) ,&
    !  & " at ", TEM, " [K]"
    
    WRITE(*,*)
    WRITE(*,*) "---------- ----------"
    
    WRITE(*,*) "Slack model: Its distinctive feature is that it performs detailed analysis using Debye temperature and"
    WRITE(*,*) " Gruneisen parameters. This model takes into account temperature changes and the effects of"
    WRITE(*,*) " interatomic bonds, allowing for highly accurate predictions."
    WRITE(*,*)
    !
    va = ( (1.0D0/3.0D0) * (1.0D0/vl**3.0D0 + 2.0D0/vt**3.0D0) )**(-1.0D0/3.0D0)
    WRITE(*,*) "average sound velocity, va [m/s]:", va  ! The Harmonic Mean of velocities.
    !
    Theta_D_va = (2.0D0*PI*hbar/kb) * ( (3.0D0*N_atom) / (4.0D0*PI*(volume*1.0D-30)) )**(1.0D0/3.0D0) * va
    WRITE(*, *) "Debye temperature, Theta_D_va [K]:", Theta_D_va
    !
    Theta_D_Cezairliyan_equ = Theta_D_va
    !
    IF (V0 <= 0.0) THEN
      V0 = volume
      WRITE(*,*) "V0:", V0
    END IF
    !
    IF (Gruneisen_parameter == 0.0) THEN
      IF (dB_per_dP > 0.0) THEN
        dB_per_dV = -1.0D0 * dB_per_dP * Bulk_modulus / V0
        IF (dG_per_dV <= 0.0) THEN
          dG_per_dV = 3.0D0*(1.0D0 - 2.0D0 * Poisson_ratio) / (2.0D0 * (1.0D0 + Poisson_ratio)) * dB_per_dV
          WRITE(*,*) "dG/dV:", dG_per_dV
        END IF
        Gruneisen_parameter_L = -(1.0D0/2.0D0) * V0 / (Bulk_modulus + (4.0D0/3.00)*Shear_modulus) * &
          & dB_per_dV + (4.0D0/3.0D0) * dG_per_dV - (1.0D0/6.0D0)
        Gruneisen_parameter_S = -(1.0D0/2.0D0) * V0/Shear_modulus * dG_per_dV - (1.0D0/6.0D0)
        Gruneisen_parameter = ( (Gruneisen_parameter_L**2.0D0 + 2.0D0*Gruneisen_parameter_S**2.0D0)/3.0D0 )**(1.0D0/2.0D0)
      ELSE
        Gruneisen_parameter = (9.0D0 - 12.0D0*(vt/vl)**2.0D0) / (2.0D0 + 4.0D0*(vt/vl)**2.0D0)
      END IF
      WRITE(*,*) "(Automatically setting) Gruneisen_parameter:", Gruneisen_parameter
    END IF
    !
    IF (Apara == 0.0) THEN
      IF (use_Apara_gamma) THEN
        ! Apara = 2.43D-8/(1.0-0.514/Gruneisen_parameter + 0.228/(Gruneisen_parameter**2.0D0))
        Apara = 1.0D0 / (1.0D0 + 1.0D0/Gruneisen_parameter + 8.3D5/Gruneisen_parameter**2.4D0)
      ELSE
        Apara = 3.1D-6 * ( (volume/N_atom) / (3.615**3.0D0/4.0D0) )**(1.0D0/3.0D0) * &
          & ( (LA + LB + LC) / volume**(1.0D0/3.0D0) )**(-1.0D0/2.0D0)
      END IF
      WRITE(*,*) "(Automatically setting) A estimated from Gruneisen_parameter:", Apara
    END IF
    !
    TEM = 300.0
    ! 1 [g/cm^3] = 1000 [kg/m^3], 1 [amu] = 1.660538921e-27 [kg]
    ! 1 [amu] = 1 [g/mol/NA] = 1.0/(6.02214085774e+23) [g] = 1.6605390402231174e-24 [g] = 1.660538921e-27 [kg]
    !Mavg = density * 1.0D6 / 1.660538921D-24 * volume * 1.0D-30 / N_atom  ! [amu] (atomic mass unit)
    Mavg = (density * volume / N_atom) / 1.660538921  ! [amu] (atomic mass unit)
    ! delta = (volume / N_atom)^(1/3) [A]
    kappa_phonon_min = Apara * Mavg * Theta_D_va**3.0D0 * (volume/N_atom)**(1.0D0/3.0D0) / &
      & (Gruneisen_parameter**2.0D0 * (N_atom)**(2.0D0/3.0D0) * TEM)
    WRITE(*,*) "Slack model: kappa_phonon_min [W/m/K]:", kappa_phonon_min, " at ", TEM, " [K]"
    kappa_phonon_min_Slack = kappa_phonon_min
    kappa_phonon_min_Slack_xK = Apara * Mavg * Theta_D_va**3.0D0 * (volume/N_atom)**(1.0D0/3.0D0) / &
      & (Gruneisen_parameter**2.0D0 * (N_atom)**(2.0D0/3.0D0))
    
    WRITE(*,*)
    WRITE(*,*) "Ref.: Thermal Conductivity of the Elements: https://srd.nist.gov/jpcrdreprint/1.3253100.pdf"
    WRITE(*,*) "Cezairliyan: k/km= [(1/3)*(T/Tm)^2 + 2/(3*(T/Tm))]^-1"
    WRITE(*,*) "Let Tm = Debye temperature and the thermal conductivity at that time be km."
    km = kappa_phonon_min_Slack_xK / Theta_D_Cezairliyan_equ
    WRITE(*,*) "Slack model, km:" , km, " at ", Theta_D_Cezairliyan_equ, " [K]"
    WRITE(*,*) "Slack model: kappa_phonon_min [W/m/K]:",&
      & km * ( (1.0D0/3.0D0)*(TEM/Theta_D_Cezairliyan_equ)**2.0D0 + 2.0D0/(3.0D0*(TEM/Theta_D_Cezairliyan_equ)) )**(-1.0D0) ,&
      & " at ", TEM, " [K]"
    
    WRITE(*,*)
    WRITE(*,*) "-------------------------------"
    WRITE(*, *) "Dulong-Petit approximatio (T >= Debye Temperature): Cv_DP = 3*R [J/(mol K)]:", 3*8.314
    WRITE(*, *) "Gas constant, R = 8.314 [J/(m K)]"
    WRITE(*,*) "------------------------------------------------------"
  END IF
  
  ! ------------------------------------------------------------------
  ! Step 4: Read phonon DOS (phononDOS or a2F.dos)
  ! ------------------------------------------------------------------
  IF (use_phononDOS) THEN
    CALL read_phonon_dos()
    WRITE(*,*) "------------------------------------------------------"
    TEM = 300.0
    I = 1
    DO WHILE( ABS(TEM - Theta_D) > 1.0D-2 .and. I <= 100)
      call find_matching_Theta_D(TEM, Theta_D, Cv_DOS, Cv_Debye)  ! Subroutines defined at the end of MODULE seebeck_data
      TEM = TEM * 0.7 + Theta_D * 0.3
      !WRITE(*,*) I, TEM, Theta_D
      I = I +1
    END DO
    IF (Theta_D > wmax/kb) THEN
      TEM = 300.0
      call find_matching_Theta_D(TEM, Theta_D, Cv_DOS, Cv_Debye)  ! Subroutines defined at the end of MODULE seebeck_data
    END IF
    WRITE(*, *) "Temperature, T [K]", TEM
    WRITE(*, *) "Matched Debye Temperature (Theta_D) [K]:", Theta_D
    WRITE(*, *) "Frequency (Omega) Maximum (wmax)    [K]:", wmax/kb
    WRITE(*, *) "Cv_DOS(T) [J/(mol K)]:", Cv_DOS
    WRITE(*, *) "Cv_Debye(T, Theta_D) [J/(mol K)]:", Cv_Debye
    WRITE(*, *) "Dulong-Petit approximation is assumed at Temperatre >= Debye Temperature", Theta_D, " or ",  wmax/kb
    !
    WRITE(*,*)
    WRITE(*,*) "---------- ----------"
    IF (Bulk_modulus <= 0.0 .or. density <= 0.0) THEN
      va = Theta_D / ( (2.0D0*PI*hbar/kb) * ( (3.0D0*N_atom) / (4.0D0*PI*(volume*1.0D-30)) )**(1.0D0/3.0D0) )
      WRITE(*,*) "average sound velocity, va [m/s]:", va, " from Debye temperature [K]", Theta_D
      
      ! [J/(mol*K)] * [mol/m^3] * ([m/s])^2 * [s] = [(J*s)/(m*K)] = [W/(m*K)]
      kappa_phonon = (1.0D0/3.0D0) * Cv_DOS * ((N_atom/volume)*(1.0e30/6.022e23)) * va**2.0D0 * tau0_phonon
      WRITE(*, *) "Calculation results of phonon thermal conductivity using Cv_DOS calculated from data in phononDOS.dat and"
      WRITE(*, *) " va calculated from Cv_DOS."
      WRITE(*, *) "kappa_phonon(Cv_DOS) [W/m/K]: ", kappa_phonon
      
      WRITE(*,*)
      WRITE(*,*) "---------- ----------"
      WRITE(*,*) "Ref.: Thermal Conductivity of the Elements: https://srd.nist.gov/jpcrdreprint/1.3253100.pdf"
      WRITE(*,*) "Cezairliyan: k/km= [(1/3)*(T/Tm)^2 + 2/(3*(T/Tm))]^-1"
      WRITE(*,*) "Let Tm = Debye temperature and the thermal conductivity at that time be km."
      km = kappa_phonon
      WRITE(*,*) "km:" , km, " at ", Theta_D, " [K]"
      WRITE(*,*) "kappa_phonon_min [W/m/K]:",&
        & km * ( (1.0D0/3.0D0)*(TEM/Theta_D)**2.0D0 + 2.0D0/(3.0D0*(TEM/Theta_D)) )**(-1.0D0) ,&
        & " at ", TEM, " [K]"
    ELSE
      ! [J/(mol*K)] * [mol/m^3] * ([m/s])^2 * [s] = [(J*s)/(m*K)] = [W/(m*K)]
      kappa_phonon = (1.0D0/3.0D0) * Cv_DOS * ((N_atom/volume)*(1.0e30/6.022e23)) * ((vl + 2.0D0*vt)/3.0D0)**2.0D0 * tau0_phonon
      WRITE(*, *) "Calculation results of phonon thermal conductivity using Cv_DOS calculated from data in phononDOS.dat and"
      WRITE(*, *) " vl and vt calculated from Bulk modulus and Poisson's_ratio, etc."
      WRITE(*, *) "kappa_phonon(Cv_DOS, vl, vt) [W/m/K]: ", kappa_phonon
      
      WRITE(*,*)
      WRITE(*,*) "---------- ----------"
      WRITE(*,*) "Ref.: Thermal Conductivity of the Elements: https://srd.nist.gov/jpcrdreprint/1.3253100.pdf"
      WRITE(*,*) "Cezairliyan: k/km= [(1/3)*(T/Tm)^2 + 2/(3*(T/Tm))]^-1"
      WRITE(*,*) "Let Tm = Debye temperature and the thermal conductivity at that time be km."
      km = kappa_phonon
      WRITE(*,*) "km:" , km, " at ", Theta_D, " [K]"
      WRITE(*,*) "kappa_phonon_min [W/m/K]:",&
        & km * ( (1.0D0/3.0D0)*(TEM/Theta_D)**2.0D0 + 2.0D0/(3.0D0*(TEM/Theta_D)) )**(-1.0D0) ,&
        & " at ", TEM, " [K]"
    END IF
  END IF
  WRITE(*,*) "------------------------------------------------------"
  !
  IF (use_a2Fdos) THEN
    CALL read_a2F_dos()
  END IF
  
  ! ------------------------------------------------------------------
  ! Step 5: Initialize output file and print header
  ! ------------------------------------------------------------------
  OPEN(UNIT=20, FILE='Seebeck_analysis.dat')
  WRITE(20, *) "! Energy shift offset (EF = 0 -> EF = DEF) [eV]:", DEF
  WRITE(20, *) "! The valence electron concentration (VEC)  (at 0 [K]):", VEC
  WRITE(20, *) "! The valence electron concentration (VEC0) (at 0 [K]):", VEC0
  WRITE(20, *) "! Carrier concentration, Nc [cm^-3] (at 0 [K]):", Nd  ! Nd = Nc for non-condition
  WRITE(20, *) "! tau mode: tau0", tau0, "[s], phonon(lambda)", use_phonon
  WRITE(20, *) "! tau mode: impurity", use_impurity, ", n_impurity", n_impurity
  WRITE(20, *) "! tau mode: dos(electron)", use_dos, ", a2Fdos", use_a2Fdos
  WRITE(20, *) "! tau mode: phnonDOS", use_phononDOS, ", N_atom", N_atom
  WRITE(20, *) "! tau mode: filter", use_selection_filter
  IF (use_phononDOS) WRITE(20, *) "! Matched Debye Temperature (Theta_D) [K]:", Theta_D
  IF (use_phononDOS .eqv. .FALSE.) THEN
    WRITE(20, *) "! Specific heat at constant volume is not calculated. Cv = 0.0000E+00 [J/(mol K)]"
    WRITE(20, *) "! Since phononDOS.dat is not used, the Slack model + Cezairliyan is used for the phonon thermal conductivity."
    ! -----
    WRITE(*, *) "Specific heat at constant volume is not calculated. Cv = 0.0000E+00 [J/(mol K)]"
    WRITE(*, *) "Since phononDOS.dat is not used, the Slack model + Cezairliyan is used for the phonon thermal conductivity."
  END IF
  WRITE(20,'(A)') hdr
  !
  OPEN(UNIT=21, FILE='ABGV2D.dat')
  WRITE(21,*) "! In this code, the relaxation time [s] is multiplied."
  WRITE(21,*) "! The dimensions shown below are those when the relaxation time is ignored with tau0 = 1.0 in the paper."
  WRITE(21,*) "! In reality, [s] is multiplied."
  WRITE(21,*) "! A(E,T): [(m/s)^2 * states/eV/unitcell] (In reality, [m^2/s * states/eV/unitcell])"
  WRITE(21,*) "! B(E,T): [(m/s)^2 * states/eV^2/unitcell] (In reality, [m^2/s * states/eV^2/unitcell])"
  WRITE(21,*) "! GV(E)^2*DOS(E): [(m/s)^2 * states/eV^2/unitcell] (In reality, [m^2/s * states/eV^2/unitcell]"
  WRITE(21,*) "# T [K]    E [eV]       E-mu [eV]    A(E,T)       B(E,T)       GV^2*DOS"

  ! ------------------------------------------------------------------
  ! Step 6: Loop over temperatures to compute Seebeck coefficient
  ! "I" is the energy index. If E1 is available, it is not required for 
  ! the get_tau() function, but we include it for future expansion.
  ! ------------------------------------------------------------------
  WRITE(*,*)
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "-----        Calculation and Results             -----"
  WRITE(*,*) "------------------------------------------------------"
  WRITE(6,'(A)')  hdr
  DO IT = 1, 25
     TEM  = TT(IT)                               ! Temperature [K]
     TEM2 = kb * TEM                             ! kT in eV (scaled)
     CP   = AMU(IT)                              ! Chemical potential mu at T

     T  = 0.0D0                                  ! Denominator accumulator
     T1 = 0.0D0                                  ! Numerator accumulator
     
     T_hole     = 0.0
     T_electron = 0.0
     
     MFP = 0.0D0                                 ! meam free path [A]
     MFP_hole = 0.0D0                            ! meam free path [A]
     MFP_electron = 0.0D0                        ! meam free path [A]
     
     m_eff_hole = 0.0D0
     m_eff_electron = 0.0D0
     
     kappa_electron = 0.0D0
     kappa_electron_term1 = 0.0D0
     kappa_electron_term2 = 0.0D0
     
     DOS_hole = 0.0D0
     DOS_electron = 0.0D0
     
     Nd = 0.0
     Nd_hole  = 0.0
     Nd_electron = 0.0

     ! Initialize previous values outside the loop
     FFOS_PREV = 0.0D0
     FFOS_CURR = 0.0D0
     FFOS_NEXT = 0.0D0
     
     FDEE_PREV = 0.0D0
     FDEE_CURR = 0.0D0
     FDEE_NEXT = 0.0D0
     
     FDD_PREV  = 0.0D0
     FDE_PREV  = 0.0D0
     
     FDD_NEXT  = 0.0D0
     FDE_NEXT  = 0.0D0
     
     MFP_PREV  = 0.0D0
     MFP_CURR  = 0.0D0
     MFP_NEXT  = 0.0D0

     m_eff_PREV = 0.0D0
     m_eff_CURR = 0.0D0
     m_eff_NEXT = 0.0D0
     
     kappa_electron_term1_PREV = 0.0D0
     kappa_electron_term1_CURR = 0.0D0
     kappa_electron_term1_NEXT = 0.0D0
     
     kappa_electron_term2_PREV = 0.0D0
     kappa_electron_term2_CURR = 0.0D0
     kappa_electron_term2_NEXT = 0.0D0
     
     FD_PREV = 0.0D0
     FD_CURR = 0.0D0
     FD_NEXT = 0.0D0
     
     DOS_PREV = 0.0D0
     DOS_CURR = 0.0D0
     DOS_NEXT = 0.0D0
     
     WRITE(21,*) 
     WRITE(21,*) "# Chemical potential [eV]:", CP

     IF (MOD(((MM - 2) - 1), 2) /= 0) LOOP_ODD = (MM - 2) - 1
     DO I = 2, MM - 2
        E  = EE(I)
        E1 = E - CP                              ! Shifted energy w.r.t. mu
        DE = EE(I+1) - EE(I)                     ! dE for integration

        ! --- Compute Fermi-Dirac distribution and its derivative ---
        FD = 1.0D0 / (1.0D0 + DEXP(E1 / TEM2))   ! f(E)
        FX = 1.0D0 - FD
        IF (E1 >= 0.0) FX = FD

        IF ((VEC0 - VEC) >= 0.0 .and. E1 <= 0.0) THEN
          Nd_hole = Nd_hole + DOS(I) * (1.0D0 - FD) * DE
        ELSE IF ((VEC0 - VEC) <= 0.0 .and. E1 >= 0.0) THEN
           Nd_electron =  Nd_electron - DOS(I) * FD * DE
        END IF

        CALL fermi_derivative(E1, TEM2, FDD, FDE)

        ! --- Multiply by velocity squared (from band data) ---
        ! simga(E) = D(E) * v(E)^2 * ta(E): if tau(E) = 1 -> simga(E) = D(E) * v(E)^2
        ! FDD = -df/dE, GV2D(E) = v(E)^2 * DOS(E) (Group velocity squared: [(m/s)^2]) (v(E) [m/s])
        ! Modified conductivity integrand using energy-dependent tau(E)
        !
        ! H. Sato et al., J. Phase Equilib. Diffus. 45, 397-415 (2024).: https://doi.org/10.1007/s11669-024-01086-y
        ! Note: The appendices are free to read, and the basics are the same, so it's a good idea to use them as reference.
        ! Theory: Seebeck coefficient, S(T) = -1/T * A(T)/(|e|*B(T)) [Equ. (2)]
        ! (Theory) A [(m/2)^2*states/unitcell], B [(m/2)^2*states/eV/unitcell], E [eV], and T [K] case: S(T) = -1/(|e|*T) * A(T)/B(T) [V/K]
        ! (On this code) A and B [(m/2)^2*states/eV/unitcell], E [eV], and T [K] case: S(T) = -1/T * A(T)/B(T)
        ! (Theory) <|v|^2 * DOS> [Equ. (2)] = (On this code) GV2D [(m/2)^2*states/eV/unitcell]
        ! In this code: S(T) = -1/TEM * T1/T
        FFOS = FDD * GV2D(I) * get_tau(E1, TEM, I)   ! Denominator: sigma(E) * df/dE with phonon and impurity scattering: B(E,T)
        FDEE = FDE * GV2D(I) * get_tau(E1, TEM, I)   ! Numerator:  (E - mu) * simga(E) * df/dE: A(E,T)

        ! Output (ABGV2D.dat): T [K], E [eV], E1 [eV], A(E,T), B(E,T), GV2D including get_tau(E1, TEM, I)
        IF (MOD(I,50) == 0) THEN
          WRITE(21,'(F8.1,1X,5(1X,E12.4))') TEM, E, E1, FDEE, FFOS, (GV2D(I) * get_tau(E1, TEM, I))
        END IF

        !---------------------------------------------
        ! Riemann Sum
        !T  = T  + DE * FFOS  ! Accumulate denominator, B(E1,TEM)
        !T1 = T1 + DE * FDEE  ! Accumulate numerator, A(E1,TEM)
        !MFP = MFP + DE * GV2(I)**(0.5) * DOS(I) * FX * get_tau(E1, TEM, I) * FX  ! [m/s]*[s] = [m]
        !kappa_electron_term1 = kappa_electron_term1 + DE * FDEE * E1
        !kappa_electron_term2 = kappa_electron_term2 + DE * FDEE**2.0D0
        !
        !IF (E1 <= 0.0) THEN
        !  T_hole   = T_hole   + DE * FFOS  ! Accumulate denominator
        !  MFP_hole = MFP_hole + DE * GV2(I)**(0.5) * DOS(I) * get_tau(E1, TEM, I)
        !  DOS_hole = DOS_hole + DOS(I) * FX
        !  IF (GV2(I) > 1.0E-12) THEN
        !    m_eff_hole = m_eff_hole + ABS(DE * 2.0D0 * E1 / GV2(I)) * DOS(I) * FX
        !  END IF
        !ELSE
        !  T_electron   = T_electron   + DE * FFOS  ! Accumulate denominator
        !  MFP_electron = MFP_electron + DE * GV2(I)**(0.5) * DOS(I) * FX * get_tau(E1, TEM, I)
        !  DOS_electron = DOS_electron + DOS(I) * FX
        !  IF (GV2(I) > 1.0E-12) THEN
        !    m_eff_electron = m_eff_electron + ABS(DE * 2.0D0 * E1 / GV2(I)) * DOS(I) * FX
        !  END IF
        !END IF
        !---------------------------------------------
        
        !---------------------------------------------
        ! Trapezoidal rule
        !DE = EE(I) - EE(I-1)  ! dE for integration
        !
        !MFP_CURR = GV2(I)**(0.5) * get_tau(E1, TEM, I) * DOS(I) * FX
        !DOS_CURR = DOS(I) * FX
        !
        !kappa_electron_term1_CURR = FDEE * E1
        !kappa_electron_term2_CURR = FDEE**2.0D0
        !
        !IF (GV2(I) > 1.0E-12 .and. DOS(I) > 0.0) THEN
        !  m_eff_CURR = ABS(FDD * 2.0D0 * E1 / GV2(I)) * DOS(I) * FX
        !ELSE
        !  m_eff_CURR = 0.0
        !END IF
        !
        !IF (I >= 3) THEN
        !  T  = T  + 0.5D0 * (FFOS_PREV + FFOS) * DE
        !  T1 = T1 + 0.5D0 * (FDEE_PREV + FDEE) * DE
        !  
        !  MFP  = MFP  + 0.5D0 * (MFP_PREV + MFP_CURR) * ABS(DE)
        !
        !  kappa_electron_term1 = kappa_electron_term1 + 0.5D0 * (kappa_electron_term1_PREV + kappa_electron_term1_CURR) * ABS(DE)
        !  kappa_electron_term2 = kappa_electron_term2 + 0.5D0 * (kappa_electron_term2_PREV + kappa_electron_term2_CURR) * ABS(DE)
        !
        !  IF (E1 <= 0.0) THEN
        !    T_hole     = T_hole     + 0.5D0 * (FFOS_PREV  + FFOS      ) * DE
        !    MFP_hole   = MFP_hole   + 0.5D0 * (MFP_PREV   + MFP_CURR  ) * ABS(DE)
        !    DOS_hole   = DOS_hole   + 0.5D0 * (DOS_PREV   + DOS_CURR  ) * ABS(DE)
        !    m_eff_hole = m_eff_hole + 0.5D0 * (m_eff_PREV + m_eff_CURR) * ABS(DE)
        !  ELSE
        !    T_electron     = T_electron     + 0.5D0 * (FFOS_PREV  + FFOS      ) * DE
        !    MFP_electron   = MFP_electron   + 0.5D0 * (MFP_PREV   + MFP_CURR  ) * ABS(DE)
        !    DOS_electron   = DOS_electron   + 0.5D0 * (DOS_PREV   + DOS_CURR  ) * ABS(DE)
        !    m_eff_electron = m_eff_electron + 0.5D0 * (m_eff_PREV + m_eff_CURR) * ABS(DE)
        !  END IF
        !END IF
        !
        !MFP_PREV = MFP_CURR
        !m_eff_PREV = m_eff_CURR
        !DOS_PREV = DOS_CURR
        !kappa_electron_term1_PREV = kappa_electron_term1_CURR
        !kappa_electron_term2_PREV = kappa_electron_term2_CURR
        !---------------------------------------------
        
        !---------------------------------------------
        ! Simpson's Rule
        
        FFOS_NEXT = FDD * GV2D(I) * get_tau(E1, TEM, I)
        FDEE_NEXT = FDE * GV2D(I) * get_tau(E1, TEM, I)
        
        MFP_NEXT = GV2(I)**(0.5) * get_tau(E1, TEM, I) * DOS(I) * FX
        DOS_NEXT = DOS(I) * FX
        
        kappa_electron_term1_NEXT = FDEE * E1
        kappa_electron_term2_NEXT = FDEE**2.0D0
        
        IF (GV2(I) > 1.0E-12 .and. DOS(I) > 0.0) THEN
          m_eff_NEXT = ABS(FDD * 2.0D0 * E1 / GV2(I)) * DOS(I) * FX
        ELSE
          m_eff_NEXT = 0.0
        END IF
        
        IF (I >= 4 .and. I <= LOOP_ODD) THEN
          ! Simpson's Rule
          DE = (EE(I) - EE(I-2)) / 2.0D0  ! center is EE(I-1)
          
          T  = T  + DE / 3.0D0 * (FFOS_PREV + 4.0D0 * FFOS_CURR + FFOS_NEXT)
          T1 = T1 + DE / 3.0D0 * (FDEE_PREV + 4.0D0 * FDEE_CURR + FDEE_NEXT)
          
          MFP  = MFP  + ABS(DE) / 3.0D0 * (MFP_PREV + 4.0D0 * MFP_CURR + MFP_NEXT)
          
          kappa_electron_term1  = kappa_electron_term1 + ABS(DE) / 3.0D0 * &
            & (kappa_electron_term1_PREV + 4.0D0 * kappa_electron_term1_CURR + kappa_electron_term1_NEXT)
            
          kappa_electron_term2  = kappa_electron_term2 + ABS(DE) / 3.0D0 * &
            & (kappa_electron_term2_PREV + 4.0D0 * kappa_electron_term2_CURR + kappa_electron_term2_NEXT)
          
          IF (E1 <= 0.0) THEN
            T_hole     = T_hole     +      DE / 3.0D0 * (FFOS_PREV  + 4.0D0 * FFOS_CURR  + FFOS_NEXT)
            MFP_hole   = MFP_hole   + ABS(DE) / 3.0D0 * (MFP_PREV   + 4.0D0 * MFP_CURR   + MFP_NEXT)
            DOS_hole   = DOS_hole   + ABS(DE) / 3.0D0 * (DOS_PREV   + 4.0D0 * DOS_CURR   + DOS_NEXT)
            m_eff_hole = m_eff_hole + ABS(DE) / 3.0D0 * (m_eff_PREV + 4.0D0 * m_eff_CURR + m_eff_NEXT)
          ELSE
            T_electron     = T_electron     +      DE / 3.0D0 * (FFOS_PREV  + 4.0D0 * FFOS_CURR  + FFOS_NEXT)
            MFP_electron   = MFP_electron   + ABS(DE) / 3.0D0 * (MFP_PREV   + 4.0D0 * MFP_CURR   + MFP_NEXT)
            DOS_electron   = DOS_electron   + ABS(DE) / 3.0D0 * (DOS_PREV   + 4.0D0 * DOS_CURR   + DOS_NEXT)
            m_eff_electron = m_eff_electron + ABS(DE) / 3.0D0 * (m_eff_PREV + 4.0D0 * m_eff_CURR + m_eff_NEXT)
          END IF
        ELSE IF (I > LOOP_ODD) THEN
          ! Riemann Sum
          DE = EE(I) - EE(I-1)  ! dE for integration
          
          T  = T  + 0.5D0 * (FFOS_CURR + FFOS_NEXT) * DE
          T1 = T1 + 0.5D0 * (FDEE_CURR + FDEE_NEXT) * DE
          
          MFP  = MFP  + 0.5D0 * (MFP_CURR + MFP_NEXT) * ABS(DE)
          
          IF (E1 <= 0.0) THEN
            T_hole     = T_hole     + 0.5D0 * (FFOS_CURR  + FFOS_NEXT ) * DE
            MFP_hole   = MFP_hole   + 0.5D0 * (MFP_CURR   + MFP_NEXT  ) * ABS(DE)
            DOS_hole   = DOS_hole   + 0.5D0 * (DOS_CURR   + DOS_NEXT  ) * ABS(DE)
            m_eff_hole = m_eff_hole + 0.5D0 * (m_eff_CURR + m_eff_NEXT) * ABS(DE)
          ELSE
            T_electron     = T_electron     + 0.5D0 * (FFOS_CURR  + FFOS_NEXT ) * DE
            MFP_electron   = MFP_electron   + 0.5D0 * (MFP_CURR   + MFP_NEXT  ) * ABS(DE)
            DOS_electron   = DOS_electron   + 0.5D0 * (DOS_CURR   + DOS_NEXT  ) * ABS(DE)
            m_eff_electron = m_eff_electron + 0.5D0 * (m_eff_CURR + m_eff_NEXT) * ABS(DE)
          END IF
        END IF
        
        FFOS_PREV = FFOS_CURR
        FDEE_PREV = FDEE_CURR
        MFP_PREV  = MFP_CURR
        m_eff_PREV = m_eff_CURR
        DOS_PREV = DOS_CURR
        kappa_electron_term1_PREV = kappa_electron_term1_CURR
        kappa_electron_term2_PREV = kappa_electron_term2_CURR
        
        FFOS_CURR = FFOS_NEXT
        FDEE_CURR = FDEE_NEXT
        MFP_CURR  = MFP_NEXT
        m_eff_CURR = m_eff_NEXT
        DOS_CURR = DOS_NEXT
        kappa_electron_term1_CURR = kappa_electron_term1_NEXT
        kappa_electron_term2_CURR = kappa_electron_term2_NEXT
        !---------------------------------------------
     END DO
     
     ! Calculation of sum to average
     MFP = MFP / (EE(MM - 2) - EE(2)) ! [m/s]*[s] = [m]
     IF (DOS_hole > 0.0) THEN
       MFP_hole = MFP_hole / (EE(MM - 2) - EE(2)) / DOS_hole ! [m/s]*[s] = [m]
       m_eff_hole = m_eff_hole / ABS(EE(MM - 2) - EE(2)) / DOS_hole
     END IF
     IF (DOS_electron < 0.0) THEN
       MFP_electron = MFP_electron / (EE(MM - 2) - EE(2)) / DOS_electron ! [m/s]*[s] = [m]
     END IF
     
     IF ((VEC0 - VEC) >= 0.0 .and. Nd_hole > 0.0) THEN
       Nd_hole = Nd_hole / volume * 1.0e24
       Nd_electron = -Nd_hole
       Nd_hole     =  Nd_hole + (VEC0 - VEC) / volume * 1.0e24
     ELSE IF ((VEC0 - VEC) <= 0.0 .and. Nd_electron < 0.0) THEN
       Nd_electron = Nd_electron / volume * 1.0e24
       Nd_hole     = -Nd_electron
       Nd_electron =  Nd_electron + (VEC0 - VEC) / volume * 1.0e24
     END IF
     
     Nd = Nd_hole + ABS(Nd_electron)
     CALL ComputeCarrierConcentration(m_eff, TEM, T0, alpha, beta, gamma, C, E_C, P, E_ph, N_sat, Delta_E_gap, pinningShift)
     
     temperature = TEM
     chemical_potential = CP
     mean_energy = (T1 / T)
     seebeck_coefficent = (-T1 / T / TEM * CO)
     conductivity = T
     conductivity_hole     = T_hole
     conductivity_electron = T_electron
     resistance = (1.0D0 / conductivity)
     power_factor = seebeck_coefficent**2 * conductivity
     ! kappa_electron = (Ln * T * temperature)
     kappa_electron = kappa_electron_term1 / temperature - kappa_electron_term2 / conductivity / temperature
     mean_free_path          = MFP         *1.0e10  ! [m] -> [Angstrom]
     mean_free_path_hole     = MFP_hole    *1.0e10  ! [m] -> [Angstrom]
     mean_free_path_electron = MFP_electron*1.0e10  ! [m] -> [Angstrom]
     carrier_concentration = Nc
     carrier_concentration_hole     = Nd_hole
     carrier_concentration_electron = Nd_electron
     mobility_hole     = (conductivity_hole     * (1.0D0 / (ech * Nd_hole          * 1.0e6)))
     mobility_electron = (conductivity_electron * (1.0D0 / (ech * ABS(Nd_electron) * 1.0e6)))
     ! hall_coefficient = (1.0D0 / (ech * (Nd_hole - ABS(Nd_electron)) * 1.0e6))
     hall_coefficient = (1.0D0 / ech) * &
       & (mobility_hole**2.0D0 * Nd_hole * 1.0e6 - mobility_electron**2.0D0 * ABS(Nd_electron) * 1.0e6 ) / &
       & (mobility_hole * Nd_hole*1.0e6 + mobility_electron * ABS(Nd_electron) * 1.0e6)**2.0D0
     effective_mass_hole     = ABS(m_eff_hole     * ech/ems)
     effective_mass_electron = ABS(m_eff_electron * ech/ems)
     
     IF (use_phononDOS) THEN
       Cv_DOS = compute_Cv_DOS(TEM)
       specific_heat = Cv_DOS         ! Specific heat at constant volume
       IF (vl >= 0.0 .and. vt >= 0.0 .and. tau0_phonon > 0.0) THEN
         ! [J/(mol*K)] * [mol/m^3] * ([m/s])^2 * [s] = [(J*s)/(m*K)] = [W/(m*K)]
         !kappa_phonon = (1.0D0/3.0D0) * Cv_DOS * ((N_atom/volume)*(1.0e30/6.022e23)) * ((vl + 2.0D0*vt)/3.0D0)**2.0D0 * tau0_phonon
         !
         ! phononDOS.dat + Cezairliyan
         kappa_phonon = km * ( (1.0D0/3.0D0)*(TEM/Theta_D)**2.0D0 + 2.0D0/(3.0D0*(TEM/Theta_D)) )**(-1.0D0)
         ZT = power_factor * temperature / (kappa_electron + kappa_phonon)
       END IF
     ELSE
       Cv_DOS = 0.0
       specific_heat = Cv_DOS         ! Specific heat at constant volume
       IF (Bulk_modulus > 0.0 .and. density > 0.0) THEN
         ! Cezairliyan
         kappa_phonon = km * ( (1.0D0/3.0D0)*(TEM/Theta_D_Cezairliyan_equ)**2.0D0 + &
           & 2.0D0/(3.0D0*(TEM/Theta_D_Cezairliyan_equ)) )**(-1.0D0)
         !
         ZT = power_factor * temperature / (kappa_electron + kappa_phonon)
       ELSE
         kappa_phonon = 0.0
         ZT = 0.0
       END IF
     END IF
     
     ! A(T) = T1 = mean_energy * conductivity = -1.0D0 * seebeck_coefficent * conductivity * TEM / CO
     A_T = T1  ! The numerator A(T) in Fig. 12 
     B_T = T   ! The denominator B(T) = T = conductivity
     
     ! --- Output results (skip if denominator too small) ---
     !=======================================================================
     ! Note: Precision considerations for REAL(8) (double precision)
     !-----------------------------------------------------------------------
     ! REAL(8) corresponds to IEEE 754 double precision, which provides
     ! approximately 15 - 17 significant decimal digits (typically 16 digits).
     !
     ! The machine epsilon (smallest distinguishable difference from 1.0D0)
     ! is approximately 2.220D-16. This means that differences smaller than
     ! this may be lost due to rounding errors.
     !
     ! Therefore, values like 1.0D-16 are near the precision limit, and
     ! comparisons or convergence checks should use thresholds such as:
     !
     !   IF (ABS(x - y) < 1.0D-14) THEN
     !       ! x and y are considered approximately equal
     !   END IF
     !
     ! For higher precision requirements, consider using REAL(16) or
     ! arbitrary precision libraries (e.g., MPFR).
     !=======================================================================
     IF (ABS(T) > 2.220D-16*10.0) THEN
        ! Calculate and output average energy offset <E - mu> and Seebeck coefficient
        ! T1/T is the averaged energy deviation <E - mu>
        ! -T1/T/TEM * CO gives Seebeck coefficient in muV/K
        WRITE(6,'(F8.1,1X,F10.6, 23(1X,E12.4))') &
          &   temperature &
          & , chemical_potential &
          & , mean_energy &
          & , seebeck_coefficent &
          & , conductivity &
          & , conductivity_hole &
          & , conductivity_electron &
          & , resistance &
          & , power_factor &
          & , kappa_electron &
          & , mean_free_path_hole &
          & , mean_free_path_electron &
          & , carrier_concentration &
          & , carrier_concentration_hole &
          & , carrier_concentration_electron &
          & , mobility_hole &
          & , mobility_electron &
          & , hall_coefficient &
          & , effective_mass_hole &
          & , effective_mass_electron &
          & , specific_heat &
          & , kappa_phonon &
          & , ZT &
          & , A_T &
          & , B_T
        WRITE(20,'(F8.1,1X,F10.6, 23(1X,E12.4))') &
          &   temperature &
          & , chemical_potential &
          & , mean_energy &
          & , seebeck_coefficent &
          & , conductivity &
          & , conductivity_hole &
          & , conductivity_electron &
          & , resistance &
          & , power_factor &
          & , kappa_electron &
          & , mean_free_path_hole &
          & , mean_free_path_electron &
          & , carrier_concentration &
          & , carrier_concentration_hole &
          & , carrier_concentration_electron &
          & , mobility_hole &
          & , mobility_electron &
          & , hall_coefficient &
          & , effective_mass_hole &
          & , effective_mass_electron &
          & , specific_heat &
          & , kappa_phonon &
          & , ZT &
          & , A_T &
          & , B_T
     ELSE
        ! If denominator T is too small (numerical instability), output placeholders
        WRITE(6,'(F8.1,1X,F10.6,1X,A)') TEM, CP, "-- -- -- -- -- -- -- -- -- -- -- --"
        WRITE(20,'(F8.1,1X,F10.6,1X,A)') TEM, CP, "-- -- -- -- -- -- -- -- -- -- -- --"
     END IF
  END DO
  
  WRITE(*,*)
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "-----         Calculation Finished               -----"
  WRITE(*,*) "------------------------------------------------------"
  
  WRITE(*,*)
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "-----         Unit Conversion Table              -----"
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "1 [eV] = 11604.5 [K]"
  WRITE(*,*) "1 [eV] = 241.8 [THz]"
  WRITE(*,*) "1 [eV] = 8065.544 [cm^-1]"
  WRITE(*,*) "1 [meV] = 8.0655 [cm^-1]"
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "1 [Rydberg] = 1 [Ry] = 3289.84 [THz]"
  WRITE(*,*) "1 [Rydberg] = 1 [Ry] = 13.6057 [eV]"
  WRITE(*,*) "1 [Rydberg] = 1 [Ry] = 13605.7 [meV]"
  WRITE(*,*) "------------------------------------------------------"
  !WRITE(*,*) "1 [Hartree] = 1 [Ha] = 27.2114 [eV]"
  !WRITE(*,*) "1 [Hartree] = 1 [Ha] = 2 [Ry]"
  !WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "1 [THz] = 4.135667696 [meV]"
  WRITE(*,*) "1 [THz] = 33.356 [cm^-1]"
  WRITE(*,*) "------------------------------------------------------"
  WRITE(*,*) "1 [cm^-1] = 0.123984 [eV]"
  WRITE(*,*) "------------------------------------------------------"
  !WRITE(*,*) "1 [Bohr] = 0.52918 [Angstrom]"
  !WRITE(*,*) "1 [Angstrom] = 1.8897259886 [Bohr]"
  !WRITE(*,*) "------------------------------------------------------"
  !WRITE(*,*) "1 [eV/A^3] = 160.21766208 [GPa]"
  !WRITE(*,*) "------------------------------------------------------"
  !WRITE(*,*) "1 [amu] = 1.66053906660 * 10^-27 [kg]"
  !WRITE(*,*) "1 [A^-1] = 2*PI [nm^-1]"
  !WRITE(*,*) "------------------------------------------------------"
  !WRITE(*,*) "Planck constant, h = 4.135667696 * 10^-15 [eV s]"
  !WRITE(*,*) "Reduced Planck, hbar = 6.582119569 * 10^-16 [eV s]"
  !WRITE(*,*) "Boltzmann constant, kB = 8.617333262 * 10^-5 [eV/K]"
  !WRITE(*,*) "Speed of light, c = 2.99792458 * 10^10 [cm/s]"
  !WRITE(*,*) "------------------------------------------------------"
  
  ! --- Finalize and close output file ---
  CLOSE(20)
  
  IF (use_phonon) THEN
    DEALLOCATE(broadening, lambdaArray, dosEf, omega_ln)  ! FUNCTION GetLambda(E1, T, I) RESULT(lambda)
  END IF
  IF (use_phononDOS) THEN
    DEALLOCATE(WPH, DOSPH, g_eff_omega, w_mode)           ! SUBROUTINE read_phonon_dos()
  END IF
  IF (use_a2Fdos) THEN
    DEALLOCATE(freq, a2F_total, a2F_mode)                 ! SUBROUTINE read_a2F_dos()
  END IF

END PROGRAM seebeck_analysis



SUBROUTINE ComputeCarrierConcentration(m_eff, TEM, T0, alpha, beta, gamma, C, E_C, P, E_ph, N_sat, Delta_E_gap, pinningShift)
  USE seebeck_data   ! Use external data module if necessary
  IMPLICIT NONE

  ! Input parameters
  !REAL(8), INTENT(IN) :: Nd          ! Doping concentration [cm^-3]
  !REAL(8), INTENT(IN) :: DEF         ! Energy shift relative to Fermi level [eV]
  REAL(8), INTENT(IN) :: m_eff        ! Effective mass (relative to electron mass)
  REAL(8), INTENT(IN) :: TEM          ! Temperature [K]
  REAL(8), INTENT(IN) :: T0           ! Reference temperature [K]
  REAL(8), INTENT(IN) :: alpha        ! Scattering correction factor
  REAL(8), INTENT(IN) :: beta         ! Temperature correction factor
  REAL(8), INTENT(IN) :: gamma        ! Temperature dependence coefficient
  REAL(8), INTENT(IN) :: C            ! Density of states constant
  REAL(8), INTENT(IN) :: E_C          ! Band edge energy [eV]
  REAL(8), INTENT(IN) :: P            ! Light intensity
  REAL(8), INTENT(IN) :: E_ph         ! Photon energy [eV]
  REAL(8), INTENT(IN) :: N_sat        ! Saturation doping concentration [cm^-3]
  REAL(8), INTENT(IN) :: Delta_E_gap  ! Bandgap narrowing coefficient [eV/K]
  REAL(8), INTENT(IN) :: pinningShift ! Fermi level pinning shift [eV]

  ! Output parameter
  !REAL(8), INTENT(OUT) :: Nc       ! Carrier concentration [cm^-3]

  ! Constants
  !REAL(8) :: kB                    ! Boltzmann constant
  !kB = 8.617D-5                    ! Boltzmann constant in eV/K

  ! Temporary variables
  REAL(8) :: effectiveMassTerm     ! Effective mass contribution
  REAL(8) :: boltzmannFactor       ! Boltzmann factor for thermal excitation
  REAL(8) :: scatteringCorrection  ! Correction for impurity scattering
  REAL(8) :: temperatureCorrection ! Temperature-dependent correction
  REAL(8) :: dosCorrection         ! Density of states correction
  REAL(8) :: opticalContribution   ! Optical excitation contribution
  REAL(8) :: saturationCorrection  ! Saturation effect correction
  REAL(8) :: gapCorrection         ! Bandgap narrowing adjustment
  REAL(8) :: pinningCorrection     ! Fermi level pinning adjustment

  ! Initialize carrier concentration
  Nc = 0.0D0  ! Start with zero carrier concentration

  ! ---------------------------------------------------
  ! Step 1: Base contribution from doping
  Nc = Nd
  ! Nd: Represents the base doping concentration in the material [cm^-3].

  ! ---------------------------------------------------
  ! Step 2: Effective mass contribution
  effectiveMassTerm = (m_eff / 1.0D0)**1.5
  Nc = Nc * effectiveMassTerm
  ! Adjusts carrier density based on effective mass of the material.

  ! ---------------------------------------------------
  ! Step 3: Thermal excitation (Boltzmann factor)
  !boltzmannFactor = EXP(-ABS(E_gap) / (kB * TEM))
  !Nc = Nc * boltzmannFactor
  ! Reflects the probability of carriers being thermally excited.

  ! ---------------------------------------------------
  ! Step 4: Scattering correction
  scatteringCorrection = (1.0D0 - alpha * Nd)
  IF (scatteringCorrection > 0.0D0) THEN
    Nc = Nc * scatteringCorrection
  ELSE
    scatteringCorrection = 0.0D0
  END IF
  ! Reduces carrier concentration due to impurity scattering at high doping levels.

  ! ---------------------------------------------------
  ! Step 5: Temperature correction
  temperatureCorrection = (1.0D0 + beta * (TEM / T0))
  Nc = Nc * temperatureCorrection
  ! Adjusts carrier density based on temperature-induced effects.

  ! ---------------------------------------------------
  ! Step 6: Saturation effect correction
  saturationCorrection = (1.0D0 - Nd / N_sat)
  IF (saturationCorrection > 0.0D0) THEN
    Nc = Nc * saturationCorrection
  ELSE
    saturationCorrection = 0.0D0
  END IF
  ! High doping concentrations may saturate carrier generation.

  ! ---------------------------------------------------
  ! Step 7: Bandgap narrowing adjustment
  gapCorrection = E_C - Delta_E_gap * TEM
  IF (gapCorrection > DEF) THEN
    dosCorrection = C * SQRT(gapCorrection - DEF)
    Nc = Nc + dosCorrection
  ELSE
    dosCorrection = 0.0D0
  END IF
  ! Bandgap narrowing at high temperature modifies the density of states.

  ! ---------------------------------------------------
  ! Step 8: Fermi level pinning correction
  pinningCorrection = pinningShift
  Nc = Nc + pinningCorrection
  ! Fermi level pinning may shift carrier generation under high doping or defect conditions.

  ! ---------------------------------------------------
  ! Step 9: Optical excitation contribution
  opticalContribution = alpha * P * EXP(-E_ph / (kB * TEM))
  Nc = Nc + opticalContribution
  ! Adds carrier generation due to light absorption.

  RETURN
END SUBROUTINE ComputeCarrierConcentration
