module ModIEGITM

  use ModKind
  use ModCharSizeGITM
  use ModUtilities, ONLY: CON_stop, CON_set_do_test

  implicit none
  private

  integer, parameter, public :: iZero_ = 0
  ! Electric Potential Types:
  integer, parameter, public :: iSwmfPot_ = 1
  ! Auroral Types:
  integer, parameter, public :: iSwmfAur_ = 1

  integer, parameter, public :: iFTA_ = 1
  integer, parameter, public :: iFRE_ = 2
  integer, parameter, public :: iPEM_ = 3
  integer, parameter, public :: iOvationPrime_ = 4
  integer, parameter, public :: iOvationSme_ = 5
  integer, parameter, public :: iAmieAur_ = 6

  real, parameter, public :: rBadValue = -1.0e32


!  integer, external :: efield_interpret_name
!  integer, external :: aurora_interpret_name

  ! This is a IEModel_ that's used when coupling to the SWMF
  type, public :: ieModel

    logical :: isOk = .true.
    logical :: isCoupleInitialized = .false.

    integer :: iDebugLevel = 0
    integer :: iEfield_ = -1
    integer :: iAurora_ = -1

    character(len=iCharLen) :: modelDir = "none"
    character(len=iCharLen) :: northFile = "none"
    character(len=iCharLen) :: southFile = "none"
    character(len=iCharLen), allocatable, dimension(:) :: northFiles
    character(len=iCharLen), allocatable, dimension(:) :: southFiles
    integer :: nFilesNorth = 0
    integer :: nFilesSouth = 0

    integer, allocatable, dimension(:, :, :) :: IeUaInterpolationIndices
    real, allocatable, dimension(:, :, :) :: IeUaInterpolationRatios
!    integer, allocatable, dimension(:, :, :) :: IeUaInterpolationIndices
!    real, allocatable, dimension(:, :, :) :: IeUaInterpolationRatios

    ! ----------------------------------------------------------------
    ! These are the states that the models has, if we either read in
    ! a file or we are coupling to some other model.
    ! They are 3D, because there is latitude, local time, and block,
    ! where block can be north/south or whatever.
    ! ----------------------------------------------------------------
    integer :: havenLats = 0
    integer :: havenMLTs = 0
    integer :: havenBLKs = 0
    real, allocatable, dimension(:, :) :: haveLats
    real, allocatable, dimension(:, :) :: haveMLTs
    ! Field-aligned Currents:
    real, allocatable, dimension(:, :) :: haveFac
    ! Potentials:
    real, allocatable, dimension(:, :) :: havePotential
    ! Electron diffuse:
    real, allocatable, dimension(:, :) :: haveDiffuseEeFlux
    real, allocatable, dimension(:, :) :: haveDiffuseEAveE
    ! Ion diffuse:
    real, allocatable, dimension(:, :) :: haveDiffuseIeFlux
    real, allocatable, dimension(:, :) :: haveDiffuseIAveE
    ! Discrete or Monoenergetic:
    real, allocatable, dimension(:, :) :: haveMonoEeFlux
    real, allocatable, dimension(:, :) :: haveMonoEAveE
    ! Broadband or Wave-drive:
    real, allocatable, dimension(:, :) :: haveWaveEeFlux
    real, allocatable, dimension(:, :) :: haveWaveEAveE
    ! Is Polar Cap (1 if is polar cap, 0 otherwise):
    real, allocatable, dimension(:, :) :: havePolarCap

    ! ----------------------------------------------------------------
    ! These are what the code that is calling this library needs
    ! ----------------------------------------------------------------

    integer :: neednLats = 0
    integer :: neednMLTs = 0
    real, allocatable, dimension(:, :) :: needLats
    real, allocatable, dimension(:, :) :: needMLTs

    real :: needImfBz = rBadValue
    real :: needImfBy = rBadValue
    logical :: useAeForHp = .false.

    ! --------------------------------------------------------------------------
    ! Keep track of whether model has been run after time & indices are updated
    ! --------------------------------------------------------------------------

    logical :: doReadMHD = .false.
    logical :: doReadSME = .false.
    logical :: doReadKp = .false.
    logical :: doReadHPI = .false. ! not actually used yet. remove noaa hpi altogether??

    logical :: isAuroraUpdated = .false.
    logical :: isPotentialUpdated = .false.

  contains
    ! Set verbose level:
    procedure :: verbose => set_verbose

    ! Set model types:
    procedure :: efield_model => set_efield_model
    procedure :: aurora_model => set_aurora_model
    procedure :: filename_north => set_filename_north
    procedure :: filename_south => set_filename_south

    procedure :: nfiles_north => set_nfiles_north
    procedure :: nfiles_south => set_nfiles_south
    procedure :: filename_list_north => set_filename_list_north
    procedure :: filename_list_south => set_filename_list_south

    ! ! Where to find the data files for empirical models:
    procedure :: model_dir => set_model_dir

    ! ! Initialize the library:
    procedure :: init => initialize
    ! ! Initalize coupling storage
    ! ! SWMF only, Only call from wrapper and srcCoupleIE
    procedure :: initCouple => initializeCouple

    ! ! set indices to run empirical models:
    procedure :: imfBz => set_bz
    procedure :: imfBy => set_by
    procedure :: swV => set_swv
    procedure :: swN => set_swn
    procedure :: hp => set_hp
    procedure :: hpN => set_hpn
    procedure :: hpS => set_hps
    procedure :: au => set_au
    procedure :: al => set_al
    procedure :: ae => set_ae
    procedure :: kp => set_kp
    procedure :: useAeHp => set_useAeForHp
    procedure :: dontAeHp => unset_useAeForHp
    procedure :: aehp => set_hp_from_ae
    procedure :: check_indices => run_check_indices

    ! ! Grid information for the calling code:
    procedure :: nMlts => set_nMlts
    procedure :: nLats => set_nLats
    procedure :: grid => set_grid
    procedure :: mlts => set_mlts
    procedure :: lats => set_lats
    ! ! SWMF only, Only call from wrapper and srcCoupleIE
    procedure :: nIeMlts => set_ie_nMlts
    procedure :: nIeLats => set_ie_nLats
    procedure :: nIeBlks => set_ie_nBlks
    procedure :: ielats => set_ie_lats
    procedure :: iemlts => set_ie_mlts

    procedure :: time_real => set_time_real
    ! procedure :: time_ymdhms => set_time_ymdhms
    ! procedure :: check_time => run_check_time

    ! ! interpolate between IE and UA grids
    ! ! SWMF only, Only call from wrapper and srcCoupleIE
    procedure :: find_ua_point => find_ua_point
!    procedure :: find_ie_point => find_ie_point
    procedure :: ie_ua_interp_indices => set_ie_ua_interpolation_indices
!    procedure :: ua_ie_interp => set_ua_ie_interpolation_indices
    procedure :: get_ie_to_ua => get_ie_values_for_ua
!    procedure :: get_ua_to_ie => get_ua_values_for_ie

    ! ! Get model results:
    procedure :: get_potential => run_potential_model
    procedure :: get_aurora => run_aurora_model_electron_diffuse
    procedure :: get_electron_diffuse_aurora => run_aurora_model_electron_diffuse
    procedure :: get_electron_mono_aurora => run_aurora_model_electron_mono
    procedure :: get_electron_wave_aurora => run_aurora_model_electron_wave
    procedure :: get_ion_diffuse_aurora => run_aurora_model_ion_diffuse
    procedure :: get_polarcap => get_polarcap_results

  end type ieModel

  contains
      ! ------------------------------------------------------------
      ! Set the verbose level for the library:
      subroutine set_verbose(this, level)
        class(ieModel) :: this
        integer, intent(in) :: level
        this%iDebugLevel = level
      end subroutine set_verbose

  INCLUDE "grid_routines.f90"

  INCLUDE "grid_interpolation.f90"

  INCLUDE "ie_coupling_routines.f90"

  INCLUDE "ie_dummy_routines.f90"

end module ModIEGITM
