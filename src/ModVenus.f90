! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

module ModPlanet

  use ModConstants
  use ModSizeGITM
  use ModOrbital

  implicit none
! Modified (01/18/07) : SWB :   Aij, s-exponents for mutual diffusion
! Modified (06/12/08) : SWB :   ordering to species revised
! Modified (06/12/08) : SWB :   nSpecies = 6; nSpeciesTotal = 11
! Modified (04/21/17) : SWB :   nEmissions = 10; iENOUV = 1
! Majors (6):  COntrol the Pressures Gradients and winds
  integer, parameter :: nSpecies = 8
  integer, parameter :: iCO2_ = 1
  integer, parameter :: iCO_ = 2
  integer, parameter :: iO_ = 3
  integer, parameter :: iN2_ = 4
  integer, parameter :: iO2_ = 5
  integer, parameter :: iAr_ = 6
  integer, parameter :: iHe_ = 7
  integer, parameter :: iN4S_ = 8

! Minors (6) : Ride on background of mean winds derived from majors
  integer, parameter :: iH_ = 9
  integer, parameter :: iN2D_ = 10
  integer, parameter :: iNO_ = 11
  integer, parameter :: nSpeciesTotal = iNO_

! Major Ions (5):  Most Important to MWACM code
! Modified (05/21/08) : SWB :   Add N2+ to major ions list
  integer, parameter  :: iOP_ = 1
  integer, parameter  :: iO2P_ = 2
  integer, parameter  :: iCO2P_ = 3
  integer, parameter  :: iN2P_ = 4
  integer, parameter  :: iNOP_ = 5
  integer, parameter  :: ie_ = 6
  integer, parameter  :: nIons = ie_
  integer, parameter  :: nIonsAdvect = 0
  integer, parameter  :: nSpeciesAll = 16 !Ions plus neutrals

  character(len=20) :: cSpecies(nSpeciesTotal)
  character(len=20) :: cIons(nIons)

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real :: Vibration(nSpeciesTotal)

  integer, parameter :: nEmissionWavelengths = 1
  integer, parameter :: nPhotoBins = 1

  !   CP      :  HEAT CAPACITY (OR SPECIFIC HEAT) OF CO2 GAS.
  real, parameter :: HeatCapacityCO2 = 735.94              ! J/(Kg*K)

  ! When you want to program in emissions, you can use this...
  integer, parameter :: iENOUV_ = 1
! integer, parameter :: iE7320_ = 2
! integer, parameter :: iE3726_ = 3
! integer, parameter :: iE5200_ = 4
! integer, parameter :: iE10400_ = 5
! integer, parameter :: iE6300_ = 6
! integer, parameter :: iE6364_ = 7

  integer, parameter :: nEmissions = 10

  real, parameter :: GC_Venus = 8.85                    ! m/s^2
  real, parameter :: RP_Venus = 10087200                ! seconds
  real, parameter :: R_Venus = 6051.8*1000.0           ! meters
  real, parameter :: DP_Venus = 0.0

  real, parameter :: Gravitational_Constant = GC_Venus
  real, parameter :: Rotation_Period = RP_Venus
  real, parameter :: RBody = R_Venus
  real, parameter :: DipoleStrength = DP_Venus

  real, parameter :: OMEGABody = 2.00*pi/Rotation_Period  ! rad/s

  real, parameter :: HoursPerDay = Rotation_Period/3600.0
  real, parameter :: Tilt = 0

  ! This is the Vernal Equinox at Midnight (Ls = 0!!!)
  ! Revised by D. Pawlwoski:  18-FEB-2013
  ! Earth-Mars clocks are set from this epoch at vernal equinox
  integer, parameter :: iVernalYear = 1998
  integer, parameter :: iVernalMonth = 7
  integer, parameter :: iVernalDay = 14
  integer, parameter :: iVernalHour = 13
  integer, parameter :: iVernalMinute = 40
  integer, parameter :: iVernalSecond = 0

! real, parameter :: SunOrbit_A = 1.52
! real, parameter :: SunOrbit_B = 0.04
! real, parameter :: SunOrbit_C = 0.15
! real, parameter :: SunOrbit_D = 0.00
! real, parameter :: SunOrbit_E = 0.00

  real, parameter :: SunOrbit_A = semimajor_Venus
  real, parameter :: SunOrbit_B = eccentricity_Venus
  real, parameter :: SunOrbit_C = longitudePerihelion_Venus
  real, parameter :: SunOrbit_D = meanLongitude_Venus
  real, parameter :: SunOrbit_E = 210664136.06

  real :: semimajoraxis_0 = semimajor_Venus
  real :: eccentricity_0 = eccentricity_Venus
  real :: inclination_0 = inclination_Venus
  real :: longitudePerihelion_0 = longitudePerihelion_Venus
  real :: longitudeNode_0 = longitudeNode_Venus
  real :: meanLongitude_0 = meanLongitude_Venus
  real :: semimajordot = semimajordot_Venus
  real :: eccentricitydot = eccentricitydot_Venus
  real :: inclinationdot = inclinationdot_Venus
  real :: longitudePeriheliondot = longitudePeriheliondot_Venus
  real :: longitudeNodedot = longitudeNodedot_Venus
  real :: meanLongitudedot = meanLongitudedot_Venus

  real, parameter :: DaysPerYear = 1.0
  real, parameter :: SecondsPerYear = DaysPerYear*Rotation_Period

  !Used as a damping term in Vertical solver.
  real, dimension(nAlts) :: VertTau = 1.0e9

  logical :: IsEarth = .false.
  logical :: IsMars = .false.
  logical :: IsTitan = .false.
  logical :: IsVenus = .true.
  logical :: NonMagnetic = .true.
  real, parameter :: PlanetNum = 0.02
  character(len=10) :: cPlanet = "Venus"

  integer, parameter :: i3371_ = 1
  integer, parameter :: i4278_ = 2
  integer, parameter :: i5200_ = 3
  integer, parameter :: i5577_ = 4
  integer, parameter :: i6300_ = 5
  integer, parameter :: i7320_ = 6
  integer, parameter :: i10400_ = 7
  integer, parameter :: i3466_ = 8
  integer, parameter :: i7774_ = 9
  integer, parameter :: i8446_ = 10
  integer, parameter :: i3726_ = 11

!  real :: KappaTemp0 = 2.22e-4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are Modified for Mars by SWB: 1/18/07
  ! -- Most source state Dij = Dji (check)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! real, parameter, dimension(6, 6) :: Diff0 = 1.0e17 * reshape( (/ &
!  ! These are Aij coefficients from B&K (1973) formulation: Aij*1.0E+17
!  ! Use Jared Bell's Titan GITM formulation for Mars GITM
!! integer, parameter :: iCO2_    = 1
!! integer, parameter :: iCO_     = 2
!! integer, parameter :: iO_      = 3
!! integer, parameter :: iN2_     = 4
!! integer, parameter :: iO2_     = 5
!! integer, parameter :: iAr_     = 6
!     !------------------------------------------------+
!     ! i=C02      CO      O        N2      O2      Ar
!     !------------------------------------------------+
!       0.0000, 0.7762, 0.2219,  0.6580,  0.5770, 1.1920,         &  ! CO2
!       0.7762, 0.0000, 0.9466,  0.9280,  0.8300, 0.6625,         &  ! CO
!       0.2219, 0.9466, 0.0000,  0.9690,  0.9690, 0.5510,         &  ! O
!       0.6580, 0.9280, 0.9690,  0.0000,  0.7150, 0.6640,         &  ! N2
!       0.5770, 0.8300, 0.9690,  0.7150,  0.0000, 0.7170,         &  ! O2
!       1.1920, 0.6625, 0.5510,  0.6640,  0.7170, 0.000 /), (/6,6/) )! Ar
!
!  ! These are s-exponents from B&K (1973) formulation: T**s
!
!   real, parameter, dimension(6, 6) :: DiffExp = reshape( (/ &
!     !------------------------------------------------+
!     ! i=C02      CO      O     N2     O2     Ar
!     !------------------------------------------------+
!       0.000,  0.750,  0.750, 0.752, 0.749, 0.750,  &           ! CO2
!       0.750,  0.000,  0.750, 0.710, 0.724, 0.750,  &           ! CO
!       0.750,  0.750,  0.000, 0.774, 0.774, 0.841,  &           ! O
!       0.752,  0.710,  0.774, 0.000, 0.750, 0.752,  &           ! N2
!       0.749,  0.724,  0.774, 0.750, 0.000, 0.736,  &           ! O2
!       0.750,  0.750,  0.841, 0.752, 0.736, 0.000 /), (/6,6/) ) ! AR
!
!!     Arrays filled in init_radcool in data statements (np = 68)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are Modified for Mars by SWB: 1/18/07
  ! -- Most source state Dij = Dji (check)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!real, parameter, dimension(nSpecies, nSpecies) :: Diff0 = 1.0e17 * reshape( (/ &
!     !----------------------------------------------------------+
!     ! i=C02      CO      O        N2      O2      Ar    He
!     !----------------------------------------------------------+
!       0.0000, 0.7762, 0.2219,  0.6580,  0.5770, 1.1920, 2.4292, &  ! CO2
!       0.7762, 0.0000, 0.9466,  0.9280,  0.8300, 0.6625,1159.55, &  ! CO
!       0.2219, 0.9466, 0.0000,  0.9690,  0.9690, 0.5510, 3.4346, &  ! O
!       0.6580, 0.9280, 0.9690,  0.0000,  0.7150, 0.6640,1159.55, &  ! N2
!       0.5770, 0.8300, 0.9690,  0.7150,  0.0000, 0.7170, 3.2070, &  ! O2
!       1.1920, 0.6625, 0.5510,  0.6640,  0.7170, 0.0000,1000.00, &  ! Ar
!       2.4292,1159.55, 3.4346,1159.550,  3.2070,1000.00, 0.0000 /),  &  ! He
!       (/nSpecies,nSpecies/) )
!
!
!  ! These are s-exponents from B&K (1973) formulation: T**s
! real, parameter, dimension(nSpecies, nSpecies) :: DiffExp = reshape( (/ &
!     !----------------------------------------------------------+
!     ! i=C02      CO      O     N2     O2     Ar     He
!     !----------------------------------------------------------+
!       0.000,  0.750,  0.750, 0.752, 0.749, 0.750, 0.720, &     ! CO2
!       0.750,  0.000,  0.750, 0.710, 0.724, 0.750, 0.524, &     ! CO
!       0.750,  0.750,  0.000, 0.774, 0.774, 0.841, 0.749, &     ! O
!       0.752,  0.710,  0.774, 0.000, 0.750, 0.752, 0.524, &     ! N2
!       0.749,  0.724,  0.774, 0.750, 0.000, 0.736, 0.710, &     ! O2
!       0.750,  0.750,  0.841, 0.752, 0.736, 0.000, 0.524, &     ! Ar
!       0.720,  0.524,  0.749, 0.524, 0.710, 0.524, 0.000 /), &  ! He
!       (/nSpecies,nSpecies/) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real, parameter, dimension(nSpecies, nSpecies) :: &
    Diff0 = 1.0e17*reshape((/ &
                           !----------------------------------------------------------------+
                           ! i=C02      CO      O        N2      O2      Ar    He       N4S
                           !----------------------------------------------------------------+
                           0.0000, 0.7762, 0.2219, 0.6580, 0.5770, 1.1920, 2.4292, 0.2219, &  ! CO2
                           0.7762, 0.0000, 0.9466, 0.9280, 0.8300, 0.6625, 1159.55, 0.9466, &  ! CO
                           0.2219, 0.9466, 0.0000, 0.9690, 0.9690, 0.5510, 3.4346, 0.9690, &  ! O
                           0.6580, 0.9280, 0.9690, 0.0000, 0.7150, 0.6640, 1159.55, 0.9690, &  ! N2
                           0.5770, 0.8300, 0.9690, 0.7150, 0.0000, 0.7170, 3.2070, 0.9690, &  ! O2
                           1.1920, 0.6625, 0.5510, 0.6640, 0.7170, 0.0000, 1000.00, 0.5510, &  ! Ar
                           2.4292, 1159.55, 3.4346, 1159.550, 3.2070, 1000.00, 0.0000, 3.34346, &  ! He
                           0.2219, 0.9466, 0.9690, 0.9690, 0.9690, 0.5510, 3.4346, 0.0000/), &  ! N4S
                           (/nSpecies, nSpecies/))

  ! These are s-exponents from B&K (1973) formulation: T**s
  real, parameter, dimension(nSpecies, nSpecies) :: &
    DiffExp = reshape((/ &
                      !----------------------------------------------------------+
                      ! i=C02      CO      O     N2     O2     Ar     He
                      !----------------------------------------------------------+
                      0.000, 0.750, 0.750, 0.752, 0.749, 0.750, 0.720, 0.750, &  ! CO2
                      0.750, 0.000, 0.750, 0.710, 0.724, 0.750, 0.524, 0.750, &  ! CO
                      0.750, 0.750, 0.000, 0.774, 0.774, 0.841, 0.749, 0.750, &  ! O
                      0.752, 0.710, 0.774, 0.000, 0.750, 0.752, 0.524, 0.774, &  ! N2
                      0.749, 0.724, 0.774, 0.750, 0.000, 0.736, 0.710, 0.774, &  ! O2
                      0.750, 0.750, 0.841, 0.752, 0.736, 0.000, 0.524, 0.841, &  ! Ar
                      0.720, 0.524, 0.749, 0.524, 0.710, 0.524, 0.000, 0.749, &  ! He
                      0.750, 0.750, 0.750, 0.774, 0.774, 0.841, 0.749, 0.000/), &     ! O
                      (/nSpecies, nSpecies/))

  integer, parameter :: np = 68, nInAlts = 54

  !! Stuff for initial conditions

  real, Dimension(nInAlts) :: newalt
  real, Dimension(nInAlts) :: InTemp
  real, Dimension(nInAlts) :: IneTemp
  real, Dimension(nInAlts) :: InITemp
  real, Dimension(nInAlts, nSpeciesTotal) :: InNDensityS
  real, Dimension(nInAlts, nIons) :: InIDensityS

  real, allocatable :: DustLatitude(:), DustLongitude(:)
  integer, parameter:: nDustLinesMax = 400, nDustPointsMax = 200, nDustLatsMax = 50, nDustTimesMax = 20
  integer :: nDustTimes, nConrathTimes, nDustLats, nDustLons, nDustAlts = 1
  real, dimension(nDustLinesMax) :: TimeDust, TimeConrath
  real, dimension(nDustLinesMax, nLats, nLons, nBlocksMax) :: HorizontalDustProfile
  real, dimension(nDustLinesMax, nLats, nLons, nBlocksMax) :: HorizontalConrathProfile

  real :: DustPressureLevel(nDustPointsMax)
  real, dimension(nDustLinesMax, nLats, nDustPointsMax, nBlocksMax) :: &
    CumulativeTauProfile, DustMixingRatioProfile

  real, dimension(nLons, nLats, nBlocksMax) :: &
    fir, fvis, Tbot, TopL, Psurf, P125, iAltMinIono, DustDistribution, ConrathDistribution

  real, dimension(1:nLons, 1:nLats, 1:nAlts) :: MarsOrbitalDistance

!################ Nelli, April 07 ##########################
!Setting up parameters needed by the correlated k lower
!atmosphere radiation code

!     Number of atmospheric layers
  integer, PARAMETER :: LL_LAYERS = nAlts

!     Number of atmospheric levels:   2 * LL_LAYERS + 3
  integer, PARAMETER :: LL_LEVELS = 2*LL_LAYERS + 3

!C     LL_NLAYRAD is the number of radiation code layers
!C     LL_NLEVRAD is the number of radiation code levels.  Level N is the
!C
  integer, PARAMETER :: LL_NLAYRAD = LL_LAYERS + 1
  integer, PARAMETER :: LL_NLEVRAD = LL_LAYERS + 2

! Bottom layer subsurface temperature
  real, parameter :: CoreTemp = 175.0

! Surface and subsurface temperature constants
  real, parameter :: Pa = 5.927E+7
  real, parameter :: Pd = 88775.0

! Stefan-Boltzmann constant in SI
  real, parameter :: SBconstant = 5.67E-8

!      Ls variable
  real :: ell_s
!C======================================================================C
!C
!C     RADINC.H    RADiation INCludes
!C
!C     Includes for the radiation code; RADIATION LAYERS, LEVELS,
!C     number of spectral intervals. . .
!C
!C     GCM2.0  Feb 2003
!C
!C======================================================================C

!C     RADIATION parameters

!C     In radiation code, layer 1 corresponds to the stratosphere.  Level
!C     1 is the top of the stratosphere.  The dummy layer is at the same
!C     temperature as the (vertically isothermal) stratosphere, and
!C     any time it is explicitly needed, the appropriate quantities will
!C     be dealt with (aka "top". . .)

!C
!C     L_NSPECTI is the number of IR spectral intervals
!C     L_NSPECTV is the number of Visual(or Solar) spectral intervals
!C     L_NGAUSS  is the number of Gauss points for K-coefficients
!C               GAUSS POINT 9 (aka the last one) is the special case
!C     L_NWNGI   is L_NSPECTI*L_NGAUSS;  the total number of "intervals"
!C               in the IR
!C     L_NWNGV   is L_NSPECTV*L_NGAUSS;  the total number of "intervals"
!C               in the VISUAL
!C
!C     L_NPREF   is the number of reference pressures that the
!C               k-coefficients are calculated on
!C     L_PINT    is the number of Lagrange interpolated reference
!C               pressures for the CO2 k-coefficients.
!C     L_NTREF   is the number of refernce temperatures for the
!C               k-coefficients
!C     L_TAUMAX  is the largest optical depth - larger ones are set
!C               to this value.
!C
!C     L_REFH2O  The number of different water-mixing ratio values for
!C               the k-coefficients that are now CO2+H2O.
!C
!C     L_NREFV   The spectral interval number of the visible reference
!C               wavelength (i.e. the 0.67 micron band)
!C
!C
!C     L_NLTE    The number of different LTE/NLTE heating rate conversion
!C               factors
!C
!C----------------------------------------------------------------------C

  integer, parameter :: L_NSPECTI = 5
  integer, parameter :: L_NSPECTV = 7
  integer, parameter :: L_NGAUSS = 17

  integer, parameter :: L_NPREF = 11
  integer, parameter :: L_NTREF = 7
  integer, parameter :: L_TAUMAX = 35

  integer, parameter :: L_PINT = 51

  integer, parameter :: L_REFH2O = 10

  integer, parameter :: L_NREFV = 6

  integer, parameter :: L_NLTE = 37

!C----------------------------------------------------------------------C
!C
!C                             radcommon.h
!C                         FORTRAN PARAMETERS
!C                          GCM2.0  Feb 2003
!C
!C----------------------------------------------------------------------C
!C
!C  "Include" grid.h and radinc.h before this file in code that uses
!C  some or all of this common data set
!C
!C     WNOI       - Array of wavenumbers at the spectral interval
!C                  centers for the infrared.  Array is NSPECTI
!C                  elements long.
!C     DWNI       - Array of "delta wavenumber", i.e., the width,
!C                  in wavenumbers (cm^-1) of each IR spectral
!C                  interval.  NSPECTI elements long.
!C     WAVEI      - Array (NSPECTI elements long) of the wavelenght
!C                  (in microns) at the center of each IR spectral
!C                  interval.
!C     WNOV       - Array of wavenumbers at the spectral interval
!C                  center for the VISUAL.  Array is NSPECTV
!C                  elements long.
!C     DWNV       - Array of "delta wavenumber", i.e., the width,
!C                  in wavenumbers (cm^-1) of each VISUAL spectral
!C                  interval.  NSPECTV elements long.
!C     WAVEV      - Array (NSPECTV elements long) of the wavelenght
!C                  (in microns) at the center of each VISUAL spectral
!C                  interval.
!C     SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
!C                  each spectral interval.  Values are for 1 AU, and
!C                  are scaled to the Mars distance elsewhere.
!C     SOL        - Solar flux at Mars (W/M^2)
!C     TAURAY     - Array (NSPECTV elements) of the pressure-independent
!C                  part of Rayleigh scattering optical depth.
!C     PTOP       - Pressure at the top of the radiation code coordinate;
!C                  = 0.5*Ptrop
!C     FZEROI     - Fraction of zeros in the IR CO2 k-coefficients, for
!C                  each temperature, pressure, and spectral interval
!C     FZEROV     - Fraction of zeros in the VISUAL CO2 k-coefficients, for
!C                  each temperature, pressure, and spectral interval
!C
!C     AEROSOL RADIATIVE OPTICAL CONSTANTS
!C     Values are at the wavelenght interval center
!C
!C     MIE SCATTERING - Size distribution weighted
!C     Qextv    - Extinction efficiency - in the visual.
!C     QextREF  - Reference visual wavelength (.67 micron band)
!C     Qscatv   - Scattering efficiency - in the visual.
!C     WV       - Single scattering albedo - in the visual.
!C     GV       - Asymmetry parameter - in the visual.
!C
!C     Qexti    - Extinction efficiency - in the infrared.
!C     Qscati   - Scattering efficiency - in the infrared.
!C     WI       - Single scattering albedo - in the infrared.
!C     GI       - Asymmetry parameter - in the infrared.
!C
!C     VIS2IR   - VISIBLE (0.67 micron band) to IR (9 micron band) ratio.
!C
!C     XLTEFACTOR - correction factor for over-prediction of heating rates
!C                  by the LTE code.  Occurs at pressures where NLTE processes
!C                  become important
!C
!C     XLTEPRESSURE - pressure regime at which each factor is to be applied
!C
!C  "Include" grid.h and radinc.h before this file in code that uses
!C  some or all of this common data set

  REAL :: WNOI(L_NSPECTI), DWNI(L_NSPECTI), WAVEI(L_NSPECTI)
  REAL :: WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV)
  REAL :: SOLARF(L_NSPECTV), TAURAY(L_NSPECTV), SOL(L_NSPECTV)

  real :: CO2I(L_NTREF, L_PINT, L_REFH2O, L_NSPECTI, L_NGAUSS)
  real :: CO2V(L_NTREF, L_PINT, L_REFH2O, L_NSPECTV, L_NGAUSS)
  real :: FZEROI(L_NSPECTI)
  real :: FZEROV(L_NSPECTV)
  real :: PGASREF(L_NPREF), TGASREF(L_NTREF)

  real :: qextv(L_NSPECTV), qscatv(L_NSPECTV), wv(L_NSPECTV)
  real :: gv(L_NSPECTV)
  real :: QextREF, VIS2IR, Cmk, tlimit

  real :: qexti(L_NSPECTI), qscati(L_NSPECTI), wi(L_NSPECTI)
  real :: gi(L_NSPECTI)

  real :: planckir(L_NSPECTI, 8501)

  real :: PTOP, UBARI, GWEIGHT(L_NGAUSS)
  real :: PFGASREF(L_PINT)

!C  H2O and CO2 k-coefficients mixed

  real :: WREFCO2(L_REFH2O), WREFH2O(L_REFH2O)

!C  LTE/NLTE heating rate conversion factors and associated pressures

  real :: XLTEFACTOR(L_NLTE), XLTEPRESSURE(L_NLTE)

!BP: For Newtonian Cooling
  real, dimension(56) :: referenceTemperature, referenceAltitude1
  real, dimension(31) :: newtonianCoolingRate, referenceAltitude2

  real :: ModCoolingRateVenus
!  These are for the Gauss-split 0.95 case

  DATA GWEIGHT/4.8083554740D-02, 1.0563099137D-01, &
    1.4901065679D-01, 1.7227479710D-01, &
    1.7227479710D-01, 1.4901065679D-01, &
    1.0563099137D-01, 4.8083554740D-02, &
    2.5307134073D-03, 5.5595258613D-03, &
    7.8426661469D-03, 9.0670945845D-03, &
    9.0670945845D-03, 7.8426661469D-03, &
    5.5595258613D-03, 2.5307134073D-03, 0.0D0/

  DATA UBARI/0.5/

!C  These are for the CO2+H2O k-coefficients

  DATA WREFCO2/9.999999D-1, 9.99999D-1, 9.9999D-1, 9.999D-1, &
    9.99D-1, 9.9D-1, 9.0D-1, 8.0D-1, 7.0D-1, 6.0D-1/

  DATA WREFH2O/1.0D-7, 1.0D-6, 1.0D-5, 1.0D-4, 1.0D-3, 1.0D-2, &
    1.0D-1, 2.0D-1, 3.0D-1, 4.0D-1/

  DATA Cmk/3.51E+22/

  DATA TLIMIT/0.5/

! ---------------------------------------------------------------------------
!  Lopez-Valverde et al. (1998) Table 1, pg. 16809.
!
!      DATA  XLTEFACTOR / 1.007,  1.007,  1.007,  1.007,  1.008,&
!                        1.009,  1.013,  1.017,  1.021,  1.027,&
!                        1.033,  1.040,  1.047,  1.054,  1.061,&
!                        1.069,  1.078,  1.087,  1.099,  1.112,&
!                        1.129,  1.149,  1.174,  1.209,  1.263,&
!                        1.344,  1.459,  1.608,  1.796,  2.033,&
!                        2.343,  2.777,  3.404,  4.366,  5.856,&
!                        8.161, 11.908, 17.791, 27.030, 40.897,&
!                       60.798, 82.660,111.855,146.365,185.922,&
!                      237.435,284.487,346.883,422.071,513.978,&
!                      635.594 /
!
!      DATA XLTEPRESSURE / 1.122E-1, 8.822E-2, 6.912E-2, 5.397E-2, 4.202E-2,&
!                         3.261E-2, 2.523E-2, 1.946E-2, 1.497E-2, 1.149E-2,&
!                         8.793E-3, 6.713E-3, 5.115E-3, 3.890E-3, 2.953E-3,&
!                         2.239E-3, 1.696E-3, 1.284E-3, 9.719E-4, 7.355E-4,&
!                         5.566E-4, 4.212E-4, 3.187E-4, 2.412E-4, 1.825E-4,&
!                         1.381E-4, 1.045E-4, 7.908E-5, 5.984E-5, 4.530E-5,&
!                         3.430E-5, 2.600E-5, 1.973E-5, 1.500E-5, 1.145E-5,&
!                         8.761E-6, 6.731E-6, 5.192E-6, 4.024E-6, 3.141E-6,&
!                         2.478E-6, 1.981E-6, 1.601E-6, 1.306E-6, 1.076E-6,&
!                         8.939E-7, 7.462E-7, 6.242E-7, 5.234E-7, 4.397E-7,&
!                         3.702E-7 /
!
! ---------------------------------------------------------------------------
!  Bougher (2011) Calculated wrt QNLTE from MTGCM case for SZA = 0.0
!  -- Ls = 90;  F10.7 = 70; dust = weak

  DATA XLTEFACTOR/1.00, 1.000, 1.000, 1.000, 1.000, &
    1.000, 1.000, 1.000, 1.000, 1.000, &
    1.000, 1.000, 1.000, 1.000, 1.000, &
    1.000, 1.000, 1.020, 1.040, 1.081, &
    1.13225, 1.13704, 1.2334, 1.32968, 2.47989, &
    5.17858, 9.21052, 13.9471, 19.1293, 25.9336, &
    45.9666, 113.783, 146.843, 212.959, 500.000, &
    1000.00, 1000.00/
!  mbar scale
  DATA XLTEPRESSURE/1.122E-1, 8.822E-2, 6.912E-2, 5.397E-2, 4.202E-2, &
    3.261E-2, 2.523E-2, 1.946E-2, 1.497E-2, 1.149E-2, &
    8.793E-3, 6.713E-3, 5.115E-3, 3.890E-3, 2.953E-3, &
    2.239E-3, 1.696E-3, 1.284E-3, 9.719E-4, 7.982E-4, &
    4.841E-4, 2.936E-4, 1.781E-4, 1.080E-4, 6.550E-5, &
    3.970E-5, 2.410E-5, 1.460E-5, 8.900E-6, 5.400E-6, &
    3.300E-6, 2.000E-6, 1.200E-6, 7.000E-7, 5.234E-7, &
    4.397E-7, 3.702E-7/
! ---------------------------------------------------------------------------

!C  PLANCK defined variable

!      planckir

!C  SETSPI and SETSPV defined variables

!      WNOI, DWNI, WAVEI, WNOV, DWNV, WAVEV
!      SOLARF, TAURAY, PTOP, TAUREF, GWEIGHT, UBARI
!      PFGASREF

!C  SETRAD defined variables

!      CO2I, CO2V, PGASREF, TGASREF
!      QEXTV, QSCATV, WV, GV
!      QEXTI, QSCATI, WI, GI, QextREF, VIS2IR

!      fzeroi, fzerov

!      WREFCO2, WREFH2O

  real*4 :: dummyalbedo(24, 36), dummyti(24, 36)

!##########################################################################
!C
!C     Includes for the nlte_tcool code: nlte_paramdef.h, nlte_commons.h
!C
!C     November 2017; S. W. Bougher
!C
!##########################################################################
!
!       Merging of different parameter definitions for new NLTE 15um param
!       jul 2012    fgg+malv
!       Modifications for MGITM
!       nov 2017    swb
!****************************************************************************
! *** Old mz1d.par ***
! Grids parameters :

  integer, parameter :: nztabul = 79  ! # points in tabulation of Tesc & VC (ISO)

! NLTE parameters:

  integer, parameter :: nltot = 20        ! incluye el actual # alt in NLTE module
  ! y el # alturas del Tstar110
  integer, parameter ::  nl = 12        ! actual # alt in NLTE module & C.Matrix
  integer, parameter ::  nl2 = nl - 2 ! = nl-2, needed for matrix inversion (mmh2)
  integer, parameter ::  nzy = (nl - 1)*4 + 1  ! Fine grid for C.Matrix

  integer, parameter :: nl_cts = 2 + nltot - nl ! actual # alt para Tstar110
  integer, parameter ::  nzy_cts = (nl_cts - 1)*4 + 1 ! fine grid for transmit calc

!  Other NLTE parameters:
  integer, parameter ::        nisot = 4        ! number of isotopes considered
  integer, parameter ::         nb = 41                ! number of bands included
  integer, parameter ::         nhist = 36        ! # of temps in histogr.
  ! (get it from histograms!)
  integer, parameter ::   nbox_max = 4     ! max.# boxes in histogram

! *** Old tcr_15um.h ***

  integer, parameter ::  itt_cza = 13     ! Selection of NLTE scheme
  ! Upper and lower limits of
  ! NLTE model

!       Ptop_atm and Pbot_atm to Mars:  Domain of calculation applicability?
  !BP: Updates for Venus (units: bars)
  real, parameter  ::  Ptop_atm = 1.e-10, Pbottom_atm = 3.2e-2

  real*8, parameter ::  rf19 = 1.d0, rf20 = 1.d0, rf21a = 1.d0, &
                       rf21b = 1.d0, rf21c = 1.d0, rf33bc = 1.d0

! *** Old bloque_dlvr11.f ***

! data
  real :: nu(nisot, 8)
  data nu(1, 1), nu(1, 2)/667.3801, 1335.1317/
  data nu(2, 1)/662.3734/
  data nu(3, 1)/648.4784/
  data nu(4, 1)/664.7289/

  real, parameter ::  nu12_0200 = 1285.4087, nu12_1000 = 1388.1847

  integer :: indexisot(nisot)
  data indexisot/26, 28, 36, 27/

  ! ctes en el sistema cgs
!       real*8, parameter ::  vlight = 2.9979245e10, ee = 1.43876866,&
!                  hplanck = 6.6260755e-27, gamma = 1.191043934e-5
  real*8, parameter ::  vlight = 2.9979245e10, ee = 1.43876866, &
                       hplanck = 6.6260755e-27, gamma1 = 1.191043934e-5

  ! datos de marte
  real ::  imr(nisot)
  data imr/0.987, 0.00408, 0.0112, 0.000742/

!****************************************************************************
!
!       Merging of different common blocks used in the new NLTE 15um param
!       jan 2012    fgg+malv
!       Modification for MGITM usage in ModMars.f90
!       nov 2017     swb
!****************************************************************************
! *** Old datitos.cmn ***
!       common /spectralv11/ elow, deltanu
  real :: elow(nisot, nb), deltanu(nisot, nb)

!       common/nu_levs_bands_v11/ nu11, nu12, nu121,
!         &          nu21, nu31, nu41
  real*8 :: nu11, nu12, nu121, nu21, nu31, nu41

!       common /aeinstein1v11/ a1_010_000, a1_020_010
!       common /aeinstein2v11/ a2_010_000
!       common /aeinstein3v11/ a3_010_000
!       common /aeinstein4v11/ a4_010_000
  real*8 :: a1_010_000, a1_020_010, a2_010_000, &
            a3_010_000, a4_010_000

! *** Old tabulation.cmn ***
!common/input_tab_v11/ lnpnbtab,
!          &          tstar11tab, tstar21tab, tstar31tab, tstar41tab,
!          &          vc210tab, vc310tab, vc410tab
  real*8 :: lnpnbtab(nztabul), vc210tab(nztabul), vc310tab(nztabul), &
            vc410tab(nztabul)
  real*8 :: tstar11tab(nztabul), tstar21tab(nztabul), &
            tstar31tab(nztabul), tstar41tab(nztabul)

! *** Old nlte_results.cmn ***
!       common/input_avilable_from/ input_cza
  integer :: input_cza

! temperatura vibracional de entrada:
!       common/temp626/ v626t1
!       common/temp628/ v628t1
!       common/temp636/ v636t1
!       common/temp627/ v627t1
  real*8 :: v626t1(nl), v628t1(nl), v636t1(nl), v627t1(nl)

! output de cza.for
!       common /tv15um/        vt11, vt12, vt21, vt31, vt41
  real*8 ::  vt11(nl), vt12(nl), vt21(nl), vt31(nl), vt41(nl)
!       common /hr15um/        hr110,hr210,hr310,hr410,hr121
  real*8 :: hr110(nl), hr121(nl), hr210(nl), hr310(nl), hr410(nl)
!       common/sf15um/ el11,el12, el21, el31, el41
  real*8 :: el11(nl), el12(nl), el21(nl), el31(nl), el41(nl)
!       common/sl15um/ sl110,sl121, sl210,sl310,sl410
  real*8 :: sl110(nl), sl121(nl), sl210(nl), sl310(nl), sl410(nl)

! *** Old matrices.cmn ***

! curtis matrix de cza:
!       common/curtis_matrixes_15um/ c110,c121, c210,
!        &          c310,c410,
!        &          vc110,vc121,vc210,vc310,vc410
  real*8 :: c110(nl, nl), c121(nl, nl), c210(nl, nl), c310(nl, nl), &
            c410(nl, nl), vc110(nl), vc121(nl), vc210(nl), &
            vc310(nl), vc410(nl)

! for the cool-to-space formulation:
!
!      common/taustar_15um/ taustar11, taustar21, taustar31,
!        &         taustar41, taustar12, taustar11_cts
  real*8 :: taustar11(nl), taustar21(nl), taustar31(nl), &
            taustar41(nl), taustar12(nl), taustar11_cts(nl_cts)

! *** Old atmref.cmn ***

! *******************************************************************
! NLTE Subgrid
! Revision for conflicts with the rest of the code via ModMars.f90
!
!       common /atm_nl/ zl, t, pl, nt, co2, n2, co, o3p,
!          &     co2vmr, n2vmr, covmr, o3pvmr,
!          &     hrkday_factor
!       real, dimension (nl) ::  zl,t,pl,nt,co2, n2, co, o3p, &
!         co2vmr, n2vmr, covmr, o3pvmr,hrkday_factor
!       real, dimension (nl) ::  zl,tl,pl,nt,co2, n2, co, o3p, &
!         co2vmr, n2vmr, covmr, o3pvmr,hrkday_factor
  real, dimension(nl) ::  zl, tl, pl, ntl, co2, n2, co, o3p, &
                         co2vmr, n2vmr, covmr, o3pvmr, hrkday_factor
! *******************************************************************

! Subgrid Transmittances
!
!       common /atm_ny/ zy, ty, py, nty, co2y
  real :: zy(nzy), ty(nzy), py(nzy), nty(nzy), co2y(nzy)

! Grids and indexes
!       common/deltazetas/ deltaz, deltazy, deltaz_cts, deltazy_cts,
!         &    jlowerboundary, jtopboundary, jtopCTS
  real  ::    deltaz, deltazy, deltaz_cts, deltazy_cts
  integer :: jlowerboundary, jtopboundary, jtopCTS

! NLTE-CTS Subgrid
!       common /atm_nl_cts/ zl_cts, t_cts, pl_cts, nt_cts,
!        &    co2_cts, n2_cts, co_cts, o3p_cts,
!        &    co2vmr_cts, n2vmr_cts, covmr_cts, o3pvmr_cts,
!        &    hrkday_factor_cts,mmean_cts,cpnew_cts

  real :: zl_cts(nl_cts), t_cts(nl_cts), pl_cts(nl_cts), &
          nt_cts(nl_cts), co2_cts(nl_cts), &
          n2_cts(nl_cts), co_cts(nl_cts), &
          o3p_cts(nl_cts), co2vmr_cts(nl_cts), n2vmr_cts(nl_cts), &
          covmr_cts(nl_cts), o3pvmr_cts(nl_cts), &
          hrkday_factor_cts(nl_cts), mmean_cts(nl_cts), &
          cpnew_cts(nl_cts)

! CTS Subgrid Transmittances
!       common /atm_ny_cts/ zy_cts, ty_cts, py_cts, nty_cts, co2y_cts
  real :: zy_cts(nzy_cts), ty_cts(nzy_cts), py_cts(nzy_cts), &
          nty_cts(nzy_cts), co2y_cts(nzy_cts)

! *** Old rates.cmn ***
!common/rates_vt/
!       &      k19ba(4),k19bb(4),k19bc(4), k19bap(4),k19bbp(4),k19bcp(4),
!       &      k19ca(4),k19cb(4),k19cc(4), k19cap(4),k19cbp(4),k19ccp(4),
!       &      k20b(4),k20c(4), k20bp(4),k20cp(4)
  real*8 :: k19ba(4), k19bb(4), k19bc(4), k19bap(4), k19bbp(4), k19bcp(4)
  real*8 :: k19ca(4), k19cb(4), k19cc(4), k19cap(4), k19cbp(4), k19ccp(4)
  real*8 :: k20b(4), k20c(4), k20bp(4), k20cp(4)

!      common/rates_vv/
!       &              k21b(4),k21c(4), k21bp(4),k21cp(4),
!       &              k33c, k33cp(2:4)
  real*8 :: k21b(4), k21c(4), k21bp(4), k21cp(4)
  real*8 :: k33c, k33cp(2:4)

!      common/rates_last/ k23k21c, k24k21c, k34k21c, &
!                      k23k21cp, k24k21cp, k34k21cp
  real*8 :: k23k21c, k24k21c, k34k21c, k23k21cp, k24k21cp, k34k21cp

! *** Old curtis.cmn ***
!       common /ini_file/ ibcode1
  character :: ibcode1*1
!       common/block1/ alsa,alda,ka,kr
  real*8 :: ka(nbox_max), alsa(nbox_max), alda(nbox_max)
  integer :: kr
!       common/block2/ hisfile
  character :: hisfile*75
!       common/block3/ pp,ta,w
  real*8 :: pp, ta(nbox_max), w
!       common/block4/ no,sk1,xls1,xld1,thist,nbox
  integer :: nbox
  real*8        :: sk1(nhist, nbox_max), xls1(nhist, nbox_max), &
                   xld1(nhist, nbox_max), thist(nhist), &
                   no(nbox_max)
!       common/block5/eqw, aa,  cc, dd, ddbox, ccbox, mr, mr_cts
  real*8 :: eqw, aa, cc, dd
  real*8 :: ddbox(nbox_max), ccbox(nbox_max)
  real*8  :: mr(nzy), mr_cts(nzy_cts)
!       common/blockstore/no_stored, sk1_stored, xls1_stored,
!        &          xld1_stored, thist_stored, nbox_stored,
!        &          mm_stored
  real*8 :: sk1_stored(nb, nhist, nbox_max)
  real*8 :: xls1_stored(nb, nhist, nbox_max)
  real*8 :: xld1_stored(nb, nhist, nbox_max)
  real*8 :: thist_stored(nb, nhist)
  real*8 :: no_stored(nb, nbox_max)
  integer :: nbox_stored(nb), mm_stored(nb)

!****************************************************************************

contains

  subroutine init_planet

    use ModTime
    use ModIoUnit, only: UnitTmp_

    integer :: iTime(7), iiAlt, ialt

    !   Mass = AMU * mean molecular weight  (unlike TGCM codes)

    Mass(iO_) = 15.9994*AMU
    Mass(iCO_) = 12.011*AMU + Mass(iO_)
    Mass(iCO2_) = Mass(iCO_) + Mass(iO_)
    Mass(iN4S_) = 14.00674*AMU
    Mass(iN2D_) = Mass(iN4S_)
    Mass(iN2_) = Mass(iN4S_)*2
    Mass(iNO_) = Mass(iN4S_) + Mass(iO_)

    Mass(iO2_) = 2*Mass(iO_)
    Mass(iAr_) = 39.948*AMU
    Mass(iHe_) = 4.0026*AMU
    Mass(iH_) = 1.0079*AMU

    cSpecies(iO_) = "O"
    cSpecies(iO2_) = "O!D2!N"
    cSpecies(iN4S_) = "N"
    cSpecies(iN2_) = "N!D2!N"
    cSpecies(iCO_) = "CO"
    cSpecies(iCO2_) = "CO!D2!N"
    cSpecies(iNO_) = "NO"
    cSpecies(iAr_) = "Ar"
    cSpecies(iH_) = "H"
    cSpecies(iHe_) = "He"

    cIons(iO2P_) = "O!D2!U+!N"
    cIons(iCO2P_) = "CO!D2!U+!N"
    cIons(iNOP_) = "NO!U+!N"
    cIons(iOP_) = "O!U+!N"
    cIons(iN2P_) = "N!D2!U+!N"
    cIons(ie_) = "e-"

    Vibration(iCO2_) = 8.66667  ! This gives Gamma = ~1.3 (experimental value)
    !Vibration(iCO2_)  = 7.0  ! Corrected by Bougher (01/18/07)!!!!
    Vibration(iCO_) = 7.0
    Vibration(iO_) = 5.0
    Vibration(iN2_) = 7.0
    Vibration(iHe_) = 5.0
    Vibration(iN4S_) = 5.0

    MassI(iOP_) = Mass(iO_)
    MassI(iNOP_) = Mass(iO_) + Mass(iN2_)/2.0
    MassI(iCO2P_) = Mass(iCO2_)
    MassI(iO2P_) = Mass(iO2_)
    MassI(ie_) = Mass_Electron

    itime = 0
    itime(1) = iVernalYear
    itime(2) = iVernalMonth
    itime(3) = iVernalDay
    itime(4) = iVernalHour
    itime(5) = iVernalMinute
    itime(6) = iVernalSecond
    call time_int_to_real(itime, VernalTime)

    !write(*,*) 'Reading in the Mars_input.txt'

    open(UNIT=67, FILE='UA/DataIn/ALBEDO_ASCII', &
         STATUS='OLD', ACTION='READ')
    read(67, *) dummyalbedo
    close(UNIT=67)

    open(UNIT=68, FILE='UA/DataIn/THERMAL_ASCII', &
         STATUS='OLD', ACTION='READ')
    read(68, *) dummyti
    close(UNIT=68)

    open(UNIT=UnitTmp_, FILE='UA/DataIn/NewVenusAtm_2p5km.txt', &
         STATUS='OLD', ACTION='READ')

111 FORMAT(F6.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, &
           ES10.3, 1X, ES10.3, 1X, ES10.3, 1X, ES10.3, 1X, &
           ES10.3, 1X, ES10.3, 1X, ES10.3)

    InNDensityS(:, :) = 1.0e+3
    InIDensityS(:, :) = 1.0e+3

    do iiAlt = 1, 54

      read(UnitTmp_, *) &
        newalt(iiAlt), &
        InTemp(iiAlt), &
        InITemp(iiAlt), &
        IneTemp(iiAlt), &
        !
        InNDensityS(iiAlt, iCO2_), &
        InNDensityS(iiAlt, iO2_), &
        InNDensityS(iiAlt, iCO_), &
        InNDensityS(iiAlt, iN2_), &
        !
        InNDensityS(iiAlt, iO_), &
        InNDensityS(iiAlt, iAr_), &
        InIDensityS(iiAlt, ie_)

    end do

    InNDensityS = Alog(inNDensityS)
    close(Unit=UnitTmp_)

    !######## Nelli, April 07 ################################
    !
    !C             PURPOSE IS TO SET UP THE CORRELATED K LOWER ATMOPSHERE
    !C             RADIATION CODE.
    !C             INITIALIZATION CONSISTS MAINLY OF READING INPUT VALUES,
    !C             CALCULATING CONSTANTS FOR USE IN THE RADIATION CODE, AND
    !C             INITIALIZING THE SET OF VALUES MAINTAINED FOR EACH BAND PASS

    !C  Set up the radiation code stuff

    CALL RADSETUP

    !###########################################################

    call readVenusIRTable

  end subroutine init_planet

!################ Nelli, April 07 ##########################
!Filling arrays needed by the correlated k lower
!atmosphere radiation code

  subroutine radsetup

!C  GCM2.0  Feb 2003
!C
!C     PURPOSE:
!C        Bundle the new radiation code setup subroutines and call
!C     this one subroutine from main, where the three included files
!C     are also listed.  Quantities are passed between this driver
!C     and the radiation code via common (in radcommon.h).
!C
!C----------------------------------------------------------------------C

    implicit none

    REAL :: FACTOR
    INTEGER :: NW

!C======================================================================C

    call setspv
    call setspi
    call setrad

!C  Scale IR opacities (Qexti and Qscati) such that
!C  TAU(0.67 micron)/TAU(9 micron) = VIS2IR, which nominally is 2.

    QextREF = Qextv(L_NREFV)
    VIS2IR = 2.75D0

    factor = Qextv(6)/(VIS2IR*Qexti(4))

    DO NW = 1, L_NSPECTI
      Qexti(NW) = Qexti(NW)*factor
      Qscati(NW) = Qscati(NW)*factor
    END DO

    PTOP = 10.0**PFGASREF(1)

    !print*,'vis2ir=',vis2ir

    !BP: Initialize the newtonian cooling
    call readNewtonianCoolingFiles
  end subroutine radsetup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE SETSPV

!C  GCM2.0  Feb 2003
!C
!C     PURPOSE:
!C        Set up the spectral intervals in the Visual (solar).  Based
!C     on Chris McKay's SETSPV code.
!C
!C     AUTHOR
!C        Jim Schaeffer
!C
!C     UPDATES FOR
!C        Bob Haberle
!C
!C     ENVIRONMENT
!C        Sun ULTRA-2       SOLARIS 2.5.1      FORTRAN
!C
!C     REVISION HISTORY
!C        Original 9/17/1998
!C     VERSION 2.0  OCT 2001
!C
!C     INPUT PARAMETERS
!C     L_NSPECTV  - Number of spectral intervals in the Visual
!C
!C     OUTPUT PARAMETERS
!C     WNOV       - Array of wavenumbers at the spectral interval
!C                  center for the VISUAL.  Array is NSPECTV
!C                  elements long.
!C     DWNV       - Array of "delta wavenumber", i.e., the width,
!C                  in wavenumbers (cm^-1) of each VISUAL spectral
!C                  interval.  NSPECTV elements long.
!C     WAVEV      - Array (NSPECTV elements long) of the wavelenght
!C                  (in microns) at the center of each VISUAL spectral
!C                  interval.
!C     SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
!C                  each spectral interval.  Values are for 1 AU, and
!C                  are scaled to the Mars distance elsewhere.
!C     TAURAY     - Array (NSPECTV elements) of the wavelength dependent
!C                  part of Rayleigh Scattering.  The pressure dependent
!C                  part is computed elsewhere (OPTCV).
!C     CALLED BY
!C        RADIATION
!C
!C     SUBROUTINES CALLED
!C        NONE
!C
!C**********************************************************************C

    implicit none

!C     BWNV - Bin wavenumber of the edges of the VISUAL spectral bins
!C     units are inverse centimeters.  Dimension needs to be changed
!C     if the number of VISUAL bins changes.

    REAL :: BWNV(L_NSPECTV + 1)
    REAL :: SOLAR(L_NSPECTV)

    REAL ::  P0, GRAV, SCALEP, SUM, WL
    INTEGER :: N, M

!C     P0      - Rayleigh scattering reference pressure in pascals.
!C     GRAV    - Acceleration due to gravity (g) - MKS
!C     SCALEP  - multiply by 100 to convert pressure from millibars
!C               to pascals.

    DATA P0/9.423D+6/
    DATA GRAV/3.72/
    DATA SCALEP/100.0/

!C     Bin wavenumber - wavenumber [cm^(-1)] at the edges of the VISUAL
!C     spectral bins.  Go from smaller to larger wavenumbers, the same as
!C     in the IR.

    DATA BWNV/2222.22D0, 3087.37D0, 4030.63D0, 5370.57D0, &
      7651.11D0, 12500.00D0, 25000.00D0, 41666.67D0/

!C     Solar flux within each spectral interval, at 1AU (W/M^2)
!C     Sum equals 1356 W/m^2 (values from Allen, 4th edition)

    DATA SOLAR/17.0, 29.0, 52.0, 148.0, 348.0, 643.0, 118.0/

!C======================================================================C

!C     Set up mean wavenumbers and wavenumber deltas.  Units of
!C     wavenumbers is cm^(-1); units of wavelengths is microns.

    do M = 1, L_NSPECTV
      WNOV(M) = 0.5*(BWNV(M + 1) + BWNV(M))
      DWNV(M) = BWNV(M + 1) - BWNV(M)
      WAVEV(M) = 1.0E+4/WNOV(M)
    end do

!C     Sum the solar flux, and write out the result.

    sum = 0.0
    do N = 1, L_NSPECTV
      SOLARF(N) = SOLAR(N)
      sum = sum + SOLARF(N)
    end do
    !write(6,'("Solar flux at 1AU = ",f7.2," W/M^2")') sum

!C     Set up the wavelength dependent part of Rayleigh Scattering.
!C     The pressure dependent part will be computed elsewhere (OPTCV).
!C     WAVEV is in microns.  There is no Rayleigh scattering in the IR.

    do N = 1, L_NSPECTV
      WL = WAVEV(N)
      TAURAY(N) = (8.7/grav)*(1.527*(1.0 + 0.013/wl**2)/wl**4)* &
                  scalep/P0
    end do

  END SUBROUTINE SETSPV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine setspi

!C  GCM2.0  Feb 2003
!C
!C     PURPOSE:
!C        Set up the spectral intervals in the infrared.  Based on
!C     Chris McKay's SETSPI code.
!C
!C     AUTHOR
!C        Jim Schaeffer
!C
!C     UPDATES FOR
!C        Bob Haberle
!C
!C     ENVIRONMENT
!C        Sun ULTRA-2       SOLARIS 2.5.1      FORTRAN
!C
!C     REVISION HISTORY
!C        Original 9/17/1998
!C        VERSION 2  OCT 2001
!C
!C     INPUT PARAMETERS
!C     L_NSPECTI  - Number of spectral intervals in the INFRARED
!C
!C     OUTPUT PARAMETERS
!C     WNOI       - Array of wavenumbers at the spectral interval
!C                  centers for the infrared.  Array is NSPECTI
!C                  elements long.
!C     DWNI       - Array of "delta wavenumber", i.e., the width,
!C                  in wavenumbers (cm^-1) of each IR spectral
!C                  interval.  NSPECTI elements long.
!C     WAVEI      - Array (NSPECTI elements long) of the wavelenght
!C                  (in microns) at the center of each IR spectral
!C                  interval.
!C
!C     CALLED BY
!C        RADIATION
!C
!C     SUBROUTINES CALLED
!C        NONE
!C
!C**********************************************************************C

    implicit none

!C     BWNI - Bin wavenumber of the edges of the IR spectral bins
!C     units are inverse centimeters.  Dimension needs to be changed
!C     if the number of IR bins changes.

    REAL :: BWNI(L_NSPECTI + 1)

    real :: a, b, x(12), w(12), ans, y, bpa, bma, T
    real :: c1, c2, wn1, wn2, PI
    integer :: n, nw, nt, m

!C  C1 and C2 values from Goody and Yung (2nd edition)  MKS units
!C  These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4

    data c1/3.741832D-16/     ! W m^-2
    data c2/1.438786D-2/     ! m K
    data PI/3.14159265358979D0/

    data x/-0.981560634246719D0, -0.904117256370475D0, &
      -0.769902674194305D0, -0.587317954286617D0, &
      -0.367831498998180D0, -0.125233408511469D0, &
      0.125233408511469D0, 0.367831498998180D0, &
      0.587317954286617D0, 0.769902674194305D0, &
      0.904117256370475D0, 0.981560634246719D0/

    data w/0.047175336386512D0, 0.106939325995318D0, &
      0.160078328543346D0, 0.203167426723066D0, &
      0.233492536538355D0, 0.249147045813403D0, &
      0.249147045813403D0, 0.233492536538355D0, &
      0.203167426723066D0, 0.160078328543346D0, &
      0.106939325995318D0, 0.047175336386512D0/

!C======================================================================C

!C     Bin wavenumber - wavenumber [cm^(-1)] at the edges of the IR
!C     spectral bins.

    BWNI(1) = 10.000D0
    BWNI(2) = 166.667D0
    BWNI(3) = 416.667D0
    BWNI(4) = 833.333D0
    BWNI(5) = 1250.000D0
    BWNI(6) = 2500.000D0

!C     Set up mean wavenumbers and wavenumber deltas.  Units of
!C     wavenumbers is cm^(-1); units of wavelengths is microns.

    do M = 1, L_NSPECTI
      WNOI(M) = 0.5*(BWNI(M + 1) + BWNI(M))
      DWNI(M) = BWNI(M + 1) - BWNI(M)
      WAVEI(M) = 1.0E+4/WNOI(M)
    end do

!C  For each IR wavelength interval, compute the integral of B(T), the
!C  Planck function, divided by the wavelength interval, in cm-1.  The
!C  integration is in MKS units, the final answer is the same as the
!C  original planck.f; W m^-2 wavenumber^-1, where wavenumber is in CM^-1.

    DO NW = 1, L_NSPECTI
      a = 1.0D-2/BWNI(NW + 1)
      b = 1.0D-2/BWNI(NW)
      bpa = (b + a)/2.0
      bma = (b - a)/2.0
      do nt = 500, 9000
        T = dble(NT)/1.0D+1
        ans = 0.0D0
        do m = 1, 12
          y = bma*x(m) + bpa
          ans = ans + w(m)*c1/(y**5*(exp(c2/(y*T)) - 1.0D0))
        end do
        planckir(NW, nt - 499) = ans*bma/(PI*DWNI(NW))
      end do
    END DO

  end subroutine setspi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine setrad

!C  GCM2.0  Feb 2003
!C
!C     PURPOSE:
!C        Set up values used by the radiation code, such as the CO2 gas
!C     absorption coefficients.  True constants are defined, and the
!C     time-independent quantities used by the radiation code are
!C     calculated.
!C
!C     AUTHOR
!C
!C
!C     UPDATES FOR
!C        Jim Pollack
!C
!C     ENVIRONMENT
!C        NAS Cray-YMP      UNICOS       FORTRAN
!C
!C     REVISION HISTORY
!C        See notes from J. Pollack dated 2/20/90 and 3/12/90.
!C        Modified Cmk  See note from Jim Pollack dated 1/11/91.
!C        Several - See note dated 1/23/91 from Jim Pollack.
!C     VERSION 2.0  OCT 2001
!C
!C     INPUT PARAMETERS
!C     DTAU(L,M)      - Dust optical depth of layer L, and for aerosol
!C                      species M.
!C     ptrop          - Pressure of the tropopause (mb)
!C     SCALEP         - Factor to convert pressures from millibars to
!C                      Pascals
!C
!C     OUTPUT PARAMETERS
!C     PTOP       - Pressure at the TOP of the stratosphere ( or
!C                  equivalently, at the bottom of the dummy layer.
!C                  PTOP = PTROP*SDUMMY
!C
!C     AEROSOL RADIATIVE OPTICAL CONSTANTS
!C     Values are at the wavelenght interval center
!C
!C     MIE SCATTERING - Size distribution weighted
!C     Qextv    - Extinction efficiency - in the visual.
!C     Qscatv   - Scattering efficiency - in the visual.
!C     WV       - Single scattering albedo - in the visual.
!C     GV       - Asymmetry parameter - in the visual.
!C
!C     Qexti    - Extinction efficiency - in the infrared.
!C     Qscati   - Scattering efficiency - in the infrared.
!C     WI       - Single scattering albedo - in the infrared.
!C     GI       - Asymmetry parameter - in the infrared.
!C
!C     CALLED BY
!C        RAD
!C
!C     SUBROUTINES CALLED
!C        DMIESS, PLNK
!C
!C----------------------------------------------------------------------C

    implicit none

    integer :: N, NS

    real :: qev1(L_NSPECTV)
    real :: qsv1(L_NSPECTV)
    real :: gv1(L_NSPECTV)

    real :: qei1(L_NSPECTI)
    real :: qsi1(L_NSPECTI)
    real :: gi1(L_NSPECTI)

    integer :: nt, np, nw, ng

!C----------------------------------------------------------------------C

!C  Visible dust properties:  Ockert-Bell Planck-weighted values (T=6000K)

!C     Qext - Ockert-Bell values (order is increasing waveNUMBER)
!C     VISULAL WAVELENGTHS.

    data qev1/2.529D0, 2.949D0, 3.209D0, 3.337D0, 3.207D0, &
      2.938D0, 2.622D0/

!C     Qscat - Ockert-Bell values
!C     VISUAL wavelengths

    data qsv1/2.374D0, 2.637D0, 3.049D0, 3.201D0, 3.045D0, &
      2.513D0, 1.623D0/

!C     G - Ockert-Bell values
!C     VISUAL wavelengths

    data gv1/0.635D0, 0.646D0, 0.630D0, 0.630D0, 0.634D0, &
      0.700D0, 0.856D0/

!C     And now the INFRARED

!C     Qext for a modified-gamma distribution, ALPHA=2, GAMMA=0.5,
!C     Rm = 0.4 microns, using the Forget optical constants (Nr and Ni).
!C     Planck-weighted values (T=215K)
!C     INFRARED wavelengths.  (The order is increasing waveNUMBER.)

    data qei1/0.193D0, 0.867D0, 1.209D0, 2.173D0, 0.638D0/

!C     Qsca for a modified gamma-distribution, using the Forget
!C     optical constants (Nr and Ni).       INFRARED wavelengths

    data qsi1/0.027D0, 0.319D0, 0.558D0, 1.136D0, 0.237D0/

!C     g for a modified gamma-distribution, using the Forget
!C     optical constants (Nr and Ni).       INFRARED wavelengths

    data gi1/0.024D0, 0.127D0, 0.288D0, 0.423D0, 0.548D0/

!C=======================================================================

!C     Set the reference pressure and temperature arrays.  These are
!C     the pressures and temperatures at which we have k-coefficients.

    pgasref(1) = 1.0E-6
    pgasref(2) = 1.0E-5
    pgasref(3) = 1.0E-4
    pgasref(4) = 1.0E-3
    pgasref(5) = 1.0E-2
    pgasref(6) = 1.0E-1
    pgasref(7) = 1.0
    pgasref(8) = 1.0E+1
    pgasref(9) = 1.0E+2
    pgasref(10) = 1.0E+3
    pgasref(11) = 1.0E+4

    tgasref(1) = 50.0
    tgasref(2) = 100.0
    tgasref(3) = 150.0
    tgasref(4) = 200.0
    tgasref(5) = 250.0
    tgasref(6) = 300.0
    tgasref(7) = 350.0

!C     Fill the (VISUAL) arrays Qextv, Qscatv, WV, GV

    DO N = 1, L_NSPECTV
      Qextv(n) = qev1(n)
      Qscatv(n) = qsv1(n)
      IF (Qscatv(n) .GE. Qextv(n)) then
        Qscatv(n) = 0.99999*Qextv(n)
      END IF
      WV(n) = Qscatv(n)/Qextv(n)
      GV(n) = gv1(n)
    END DO

!C     Fill the (INFRARED) arrays Qexti, Qscati, WI, GI

    DO N = 1, L_NSPECTI
      Qexti(n) = qei1(n)
      Qscati(n) = qsi1(n)
      IF (Qscati(n) .GE. Qexti(n)) then
        Qscati(n) = 0.99999*Qexti(n)
      END IF
      WI(n) = Qscati(n)/Qexti(n)
      GI(n) = gi1(n)
    END DO

!C     Get CO2 k coefficients, and interpolate them to the finer
!C     pressure grid.

    call laginterp(pgasref, pfgasref)

  end subroutine setrad

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine laginterp(pgref, pint)

!C  GCM2.0  Feb 2003

!C  Lagrange interpolation (linear in log pressure) of the CO2
!C  k-coefficients in the pressure domain.  Subsequent use of these
!C  values will use a simple linear interpolation in pressure.

    implicit none

    integer :: n, nt, np, nh, ng, nw, m, i, j, k, l
    real :: co2i8(L_NTREF, L_NPREF, L_REFH2O, L_NSPECTI, L_NGAUSS)
    real :: co2v8(L_NTREF, L_NPREF, L_REFH2O, L_NSPECTV, L_NGAUSS)
    real :: pgref(L_NPREF)

    real :: x, xi(4), yi(4), ans
    real :: pint(L_PINT), pin(L_PINT), pref(L_NPREF), p

    data pin/-6.0D0, -5.8D0, -5.6D0, -5.4D0, -5.2D0, &
      -5.0D0, -4.8D0, -4.6D0, -4.4D0, -4.2D0, &
      -4.0D0, -3.8D0, -3.6D0, -3.4D0, -3.2D0, &
      -3.0D0, -2.8D0, -2.6D0, -2.4D0, -2.2D0, &
      -2.0D0, -1.8D0, -1.6D0, -1.4D0, -1.2D0, &
      -1.0D0, -0.8D0, -0.6D0, -0.4D0, -0.2D0, &
      0.0D0, 0.2D0, 0.4D0, 0.6D0, 0.8D0, &
      1.0D0, 1.2D0, 1.4D0, 1.6D0, 1.8D0, &
      2.0D0, 2.2D0, 2.4D0, 2.6D0, 2.8D0, &
      3.0D0, 3.2D0, 3.4D0, 3.6D0, 3.8D0, &
      4.0D0/

!C======================================================================!

!C  Fill pint for output from this subroutine

    do n = 1, L_PINT
      PINT(n) = PIN(n)
    end do

!C  Take log of the reference pressures

    do n = 1, L_NPREF
      pref(n) = LOG10(PGREF(n))
    end do

!C     Get CO2 k coefficients

    OPEN(Unit=25, file='UA/DataIn/CO2H2O_V_12_95_ASCII', ACTION='READ')

    DO i = 1, L_NTREF
      DO j = 1, L_NPREF
        DO k = 1, L_REFH2O
          DO l = 1, L_NSPECTV
            DO m = 1, L_NGAUSS
              read(25, *) co2v8(i, j, k, l, m)
            END DO
          END DO
        END DO
      END DO
    END DO

    DO i = 1, L_NSPECTV
      read(25, *) fzerov(i)
    END DO
    CLOSE(25)

    OPEN(unit=35, file='UA/DataIn/CO2H2O_IR_12_95_ASCII', ACTION='READ')

    DO i = 1, L_NTREF
      DO j = 1, L_NPREF
        DO k = 1, L_REFH2O
          DO l = 1, L_NSPECTI
            DO m = 1, L_NGAUSS
              read(35, *) co2i8(i, j, k, l, m)
            END DO
          END DO
        END DO
      END DO
    END DO

    DO i = 1, L_NSPECTI
      read(35, *) fzeroi(i)
    END DO
    CLOSE(35)

!!$      print*,'co2v8(3,4,1,3,2) = ',co2v8(3,4,1,3,2)
!!$      print*,'co2v8(1,2,3,4,5) = ',co2v8(1,2,3,4,5)
!!$      print*,'co2v8(3,4,1,3,3) = ',co2v8(3,4,1,3,3)
!!$
!!$      print*,'co2i8(3,4,1,3,2) = ',co2i8(3,4,1,3,2)
!!$      print*,'co2i8(1,2,3,4,5) = ',co2i8(1,2,3,4,5)
!!$      print*,'co2i8(3,4,1,3,3) = ',co2i8(3,4,1,3,3)
!C  Take Log10 of the values - we interpolate the log10 of the values,
!C  not the values themselves.   Smallest value is 1.0E-200.

    do nt = 1, L_NTREF
      do np = 1, L_NPREF
        do nh = 1, L_REFH2O
          do ng = 1, L_NGAUSS

            do nw = 1, L_NSPECTV
              if (co2v8(nt, np, nh, nw, ng) .gt. 1.0e-200) then
                co2v8(nt, np, nh, nw, ng) = log10(co2v8(nt, np, nh, nw, ng))
              else
                co2v8(nt, np, nh, nw, ng) = -200.0
              end if
            end do

            do nw = 1, L_NSPECTI
              if (co2i8(nt, np, nh, nw, ng) .gt. 1.0e-200) then
                co2i8(nt, np, nh, nw, ng) = log10(co2i8(nt, np, nh, nw, ng))
              else
                co2i8(nt, np, nh, nw, ng) = -200.0
              end if
            end do

          end do
        end do
      end do
    end do

!C  Interpolate the values:  first the IR

    do nt = 1, L_NTREF
      do nh = 1, L_REFH2O
      do nw = 1, L_NSPECTI
        do ng = 1, L_NGAUSS

!C  First, the initial interval (P=1e-6 to 1e-5)

          n = 1
          do m = 1, 5
            x = pint(m)
            xi(1) = pref(n)
            xi(2) = pref(n + 1)
            xi(3) = pref(n + 2)
            xi(4) = pref(n + 3)
            yi(1) = co2i8(nt, n, nh, nw, ng)
            yi(2) = co2i8(nt, n + 1, nh, nw, ng)
            yi(3) = co2i8(nt, n + 2, nh, nw, ng)
            yi(4) = co2i8(nt, n + 3, nh, nw, ng)
            call lagrange(x, xi, yi, ans)
            co2i(nt, m, nh, nw, ng) = 10.0**ans
          end do

          do n = 2, L_NPREF - 2
            do m = 1, 5
              i = (n - 1)*5 + m
              x = pint(i)
              xi(1) = pref(n - 1)
              xi(2) = pref(n)
              xi(3) = pref(n + 1)
              xi(4) = pref(n + 2)
              yi(1) = co2i8(nt, n - 1, nh, nw, ng)
              yi(2) = co2i8(nt, n, nh, nw, ng)
              yi(3) = co2i8(nt, n + 1, nh, nw, ng)
              yi(4) = co2i8(nt, n + 2, nh, nw, ng)
              call lagrange(x, xi, yi, ans)
              co2i(nt, i, nh, nw, ng) = 10.0**ans
            end do
          end do

!C  Now, get the last interval (P=1e+3 to 1e+4)

          n = L_NPREF - 1

          do m = 1, 5
            i = (n - 1)*5 + m
            x = pint(i)
            xi(1) = pref(n - 2)
            xi(2) = pref(n - 1)
            xi(3) = pref(n)
            xi(4) = pref(n + 1)
            yi(1) = co2i8(nt, n - 2, nh, nw, ng)
            yi(2) = co2i8(nt, n - 1, nh, nw, ng)
            yi(3) = co2i8(nt, n, nh, nw, ng)
            yi(4) = co2i8(nt, n + 1, nh, nw, ng)
            call lagrange(x, xi, yi, ans)
            co2i(nt, i, nh, nw, ng) = 10.0**ans
          end do

!C  Fill the last pressure point

          co2i(nt, L_PINT, nh, nw, ng) = 10.0**co2i8(nt, L_NPREF, nh, nw, ng)

        end do
      end do
      end do
    end do

!C  Interpolate the values:  now the Visual

    do nt = 1, L_NTREF
      do nh = 1, L_REFH2O
      do nw = 1, L_NSPECTV
        do ng = 1, L_NGAUSS

!C  First, the initial interval (P=1e-6 to 1e-5)

          n = 1
          do m = 1, 5
            x = pint(m)
            xi(1) = pref(n)
            xi(2) = pref(n + 1)
            xi(3) = pref(n + 2)
            xi(4) = pref(n + 3)
            yi(1) = co2v8(nt, n, nh, nw, ng)
            yi(2) = co2v8(nt, n + 1, nh, nw, ng)
            yi(3) = co2v8(nt, n + 2, nh, nw, ng)
            yi(4) = co2v8(nt, n + 3, nh, nw, ng)
            call lagrange(x, xi, yi, ans)
            co2v(nt, m, nh, nw, ng) = 10.0**ans
          end do

          do n = 2, L_NPREF - 2
            do m = 1, 5
              i = (n - 1)*5 + m
              x = pint(i)
              xi(1) = pref(n - 1)
              xi(2) = pref(n)
              xi(3) = pref(n + 1)
              xi(4) = pref(n + 2)
              yi(1) = co2v8(nt, n - 1, nh, nw, ng)
              yi(2) = co2v8(nt, n, nh, nw, ng)
              yi(3) = co2v8(nt, n + 1, nh, nw, ng)
              yi(4) = co2v8(nt, n + 2, nh, nw, ng)
              call lagrange(x, xi, yi, ans)
              co2v(nt, i, nh, nw, ng) = 10.0**ans
            end do
          end do

!C  Now, get the last interval (P=1e+3 to 1e+4)

          n = L_NPREF - 1

          do m = 1, 5
            i = (n - 1)*5 + m
            x = pint(i)
            xi(1) = pref(n - 2)
            xi(2) = pref(n - 1)
            xi(3) = pref(n)
            xi(4) = pref(n + 1)
            yi(1) = co2v8(nt, n - 2, nh, nw, ng)
            yi(2) = co2v8(nt, n - 1, nh, nw, ng)
            yi(3) = co2v8(nt, n, nh, nw, ng)
            yi(4) = co2v8(nt, n + 1, nh, nw, ng)
            call lagrange(x, xi, yi, ans)
            co2v(nt, i, nh, nw, ng) = 10.0**ans
          end do

!C  Fill the last pressure point

          co2v(nt, L_PINT, nh, nw, ng) = 10.0**co2v8(nt, L_NPREF, nh, nw, ng)

        end do
      end do
      end do
    end do

  end subroutine laginterp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine lagrange(x, xi, yi, ans)

!C  GCM2.0  Feb 2003
!C
!C  Lagrange interpolation - Polynomial interpolation at point x
!C  xi(1) <= x <= xi(4).  Yi(n) is the functional value at XI(n).

    implicit none

    real :: x, xi(4), yi(4), ans
    real :: fm1, fm2, fm3, fm4

!C======================================================================!

    fm1 = x - XI(1)
    fm2 = x - XI(2)
    fm3 = x - XI(3)
    fm4 = x - XI(4)

!C  Get the "answer" at the requested X

    ans = fm2*fm3*fm4*YI(1)/ &
          ((XI(1) - XI(2))*(XI(1) - XI(3))*(XI(1) - XI(4))) + &
          fm1*fm3*fm4*YI(2)/ &
          ((XI(2) - XI(1))*(XI(2) - XI(3))*(XI(2) - XI(4))) + &
          fm1*fm2*fm4*YI(3)/ &
          ((XI(3) - XI(1))*(XI(3) - XI(2))*(XI(3) - XI(4))) + &
          fm1*fm2*fm3*YI(4)/ &
          ((XI(4) - XI(1))*(XI(4) - XI(2))*(XI(4) - XI(3)))

  end subroutine lagrange

!##########################################################

  !BP (7/25/2020)
  subroutine readNewtonianCoolingFiles
    use ModIoUnit, only: UnitTmp_
    !use ModSources, only : RadCoolingRate

    integer :: iiAlt

    open(UNIT=UnitTmp_, FILE='DataIn/HausReferenceTemperature.csv', &
         STATUS='OLD', ACTION='READ')

    !Reference temperature .csv file is temperature vs. Altitude
    !Figure 4 in Haus et al. 2015
    do iiAlt = 1, 56
      read(UnitTmp_, *) referenceTemperature(iiAlt), referenceAltitude1(iiAlt)
    end do

    close(Unit=UnitTmp_)

    !Figure 24 in Haus et al. 2015
    !Contains the globally averaged cooling rates as a function of altitude
    open(UNIT=UnitTmp_, FILE='DataIn/Sub100kmGlobalAverageCooling.csv', &
         STATUS='OLD', ACTION='READ')

    do iiAlt = 1, 31
      read(UnitTmp_, *) newtonianCoolingRate(iiAlt), referenceAltitude2(iiAlt)
    end do

    close(Unit=UnitTmp_)

    newtonianCoolingRate = newtonianCoolingRate/86400
    !write(*,*) "Newtonian Cooling:", newtonianCoolingRate
  end subroutine readNewtonianCoolingFiles

  subroutine init_magheat
    return
  end subroutine init_magheat

  subroutine init_aerosol
    return
  end subroutine init_aerosol

  subroutine init_radcooling
    return
  end subroutine init_radcooling
end module ModPlanet

