module ModIEGITM

  use ModKind
  use ModCharSizeGITM

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

    ! ----------------------------------------------------------------
    ! These are the states that the models has, if we either read in
    ! a file or we are coupling to some other model.
    ! They are 3D, because there is latitude, local time, and block,
    ! where block can be north/south or whatever.
    ! ----------------------------------------------------------------
    integer :: havenLats = 0
    integer :: havenMLTs = 0
    integer :: havenBLKs = 0
    real, allocatable, dimension(:, :, :) :: haveLats
    real, allocatable, dimension(:, :, :) :: haveMLTs
    ! Field-aligned Currents:
    real, allocatable, dimension(:, :, :) :: haveFac
    ! Potentials:
    real, allocatable, dimension(:, :, :) :: havePotential
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
    ! procedure :: mlts => set_mlts
    ! procedure :: lats => set_lats

     procedure :: time_real => set_time_real
    ! procedure :: time_ymdhms => set_time_ymdhms
    ! procedure :: check_time => run_check_time

    ! ! Get model results:
     procedure :: get_potential => run_potential_model
     procedure :: get_aurora => run_aurora_model_electron_diffuse
     procedure :: get_electron_diffuse_aurora => run_aurora_model_electron_diffuse
     procedure :: get_electron_mono_aurora => run_aurora_model_electron_mono
     procedure :: get_electron_wave_aurora => run_aurora_model_electron_wave
     procedure :: get_ion_diffuse_aurora => run_aurora_model_ion_diffuse
     procedure :: get_polarcap => get_polarcap_results

  end type ieModel

  interface
    module subroutine run_potential_model(ie, potential)
      class(ieModel) :: ie
      real, dimension(ie%neednMLTs, &
              ie%neednLats), intent(out) :: potential
    end subroutine run_potential_model
    ! ------------------------------------------------------------
    module subroutine set_efield_model(this, efield_model)
      class(ieModel) :: this
      character(len=*), intent(in) :: efield_model
    end subroutine set_efield_model
    ! ------------------------------------------------------------
    module subroutine set_aurora_model(this, aurora_model)
      class(ieModel) :: this
      character(len=*), intent(in) :: aurora_model
    end subroutine set_aurora_model
    ! ------------------------------------------------------------
    integer function efield_interpret_name(efieldString)
      use ModErrors
      implicit none
      character(len=*), intent(in) :: efieldString
      character(len=len(efieldString)) :: efieldLower
    end function efield_interpret_name
    ! ------------------------------------------------------------
    integer function aurora_interpret_name(auroraString)
      use ModErrors
      implicit none
      character(len=*), intent(in) :: auroraString
      character(len=len(auroraString)) :: auroraLower
    end function aurora_interpret_name
    ! ------------------------------------------------------------
    module subroutine set_nMlts(this, iValue)
      class(ieModel) :: this
      integer, intent(in) :: iValue
    end subroutine set_nMlts
    ! ------------------------------------------------------------
    module subroutine set_nLats(this, iValue)
      class(ieModel) :: this
      integer, intent(in) :: iValue
    end subroutine set_nLats
    ! ------------------------------------------------------------
    module subroutine set_grid(this, MltsIn, LatsIn)
!      use ModAMIE_Interface, only: set_interpolation_indices
      class(ieModel) :: this
      real, dimension(this%neednMlts, this%neednLats), intent(in) :: MltsIn
      real, dimension(this%neednMlts, this%neednLats), intent(in) :: LatsIn
    end subroutine set_grid
    ! ------------------------------------------------------------
    module subroutine set_mlts(this, MltsIn)
      class(ieModel), intent(inout) :: this
      real, dimension(this%neednMlts, this%neednLats), intent(in) :: MltsIn
      integer :: iError = 0, iMlt, iLat
    end subroutine set_mlts
    ! ------------------------------------------------------------
    module subroutine set_lats(this, LatsIn)
      class(ieModel) :: this
      real, dimension(this%neednMlts, this%neednLats), intent(in) :: LatsIn
      integer :: iError = 0, iMlt, iLat
    end subroutine set_lats
  end interface

  contains

  ! ------------------------------------------------------------
  ! Set the verbose level for the library:
  subroutine set_verbose(this, level)
    class(ieModel) :: this
    integer, intent(in) :: level
    this%iDebugLevel = level
  end subroutine set_verbose
  ! ==========================================================================
  ! filename for north
  subroutine set_filename_north(this, filename)
    class(ieModel) :: this
    character(len=*), intent(in) :: filename
  end subroutine set_filename_north
  ! South
  subroutine set_filename_south(this, filename)
    class(ieModel) :: this
    character(len=*), intent(in) :: filename
  end subroutine set_filename_south

  ! set nfiles for north
  subroutine set_nfiles_north(this, nFiles)
    class(ieModel) :: this
    integer, intent(in) :: nFiles
  end subroutine set_nfiles_north

  ! ------------------------------------------------------------
  ! filename for north
  subroutine set_filename_list_north(this, iFile, filename)
    class(ieModel) :: this
    integer, intent(in) :: iFile
    character(len=*), intent(in) :: filename
  end subroutine set_filename_list_north
  ! ------------------------------------------------------------
  ! set nfiles for south
  subroutine set_nfiles_south(this, nFiles)
    class(ieModel) :: this
    integer, intent(in) :: nFiles
  end subroutine set_nfiles_south
  ! ------------------------------------------------------------
  ! filename for south
  subroutine set_filename_list_south(this, iFile, filename)
    class(ieModel) :: this
    integer, intent(in) :: iFile
    character(len=*), intent(in) :: filename
  end subroutine set_filename_list_south
  ! ------------------------------------------------------------
  subroutine set_model_dir(this, dir)
    class(ieModel) :: this
    character(len=*), intent(in) :: dir
  end subroutine set_model_dir
  ! ------------------------------------------------------------
  subroutine run_aurora_model(ie)
    class(ieModel) :: ie
    integer :: iError = 0
    return
  end subroutine run_aurora_model
  ! ------------------------------------------------------------
  subroutine run_aurora_model_electron_diffuse(ie, eflux, avee)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: eFlux
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: AveE
    integer :: iError = 0
    return
  end subroutine run_aurora_model_electron_diffuse
  ! ------------------------------------------------------------
  subroutine run_aurora_model_electron_mono(ie, eflux, avee)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: eFlux
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: AveE
    integer :: iError = 0
    return
  end subroutine run_aurora_model_electron_mono
  ! ------------------------------------------------------------
  subroutine run_aurora_model_electron_wave(ie, eflux, avee)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: eFlux
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: AveE
    integer :: iError = 0
    return
  end subroutine run_aurora_model_electron_wave
  ! ------------------------------------------------------------
  subroutine run_aurora_model_ion_diffuse(ie, eflux, avee)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: eFlux
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: AveE
    integer :: iError = 0
    return
  end subroutine run_aurora_model_ion_diffuse
  ! ------------------------------------------------------------
  subroutine get_polarcap_results(ie, polarcap)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: polarcap
    return
  end subroutine get_polarcap_results
  ! ------------------------------------------------------------
  subroutine set_bz(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_bz
  ! ------------------------------------------------------------
  subroutine set_by(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_by
  ! ------------------------------------------------------------
  subroutine set_swv(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_swv
  ! ------------------------------------------------------------
  subroutine set_swn(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_swn
  ! ------------------------------------------------------------
  subroutine set_hp(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_hp
  ! ------------------------------------------------------------
  subroutine set_hpn(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_hpn
  ! ------------------------------------------------------------
  subroutine set_hps(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_hps
  ! ------------------------------------------------------------
  subroutine set_hp_from_ae(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_hp_from_ae
  ! ------------------------------------------------------------
  subroutine set_au(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_au
  ! ------------------------------------------------------------
  subroutine set_al(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_al
  ! ------------------------------------------------------------
  subroutine set_ae(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_ae
  ! ------------------------------------------------------------
  subroutine set_useAeForHp(this)
    class(ieModel), intent(inout) :: this
  end subroutine set_useAeForHp
  ! ------------------------------------------------------------
  subroutine unset_useAeForHp(this)
    class(ieModel), intent(inout) :: this
  end subroutine unset_useAeForHp
  ! ------------------------------------------------------------
  subroutine set_kp(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_kp
  ! ------------------------------------------------------------
  ! Initialize all of the different models
  subroutine initialize(this)
    use ModErrors
    class(ieModel) :: this
    character(len=iCharLenIE_) :: modelDirTotal
    character(len=iCharLenIE_) :: inFileNameTotal
    character(len=iCharLenIE_) :: name
    integer :: UnitTmp_ = 76
    integer :: iError = 0
    write(*,*) 'We should really figure out what this needs to do'
  end subroutine initialize
  ! ------------------------------------------------------------
  subroutine run_check_indices(this)
    class(ieModel) :: this
  end subroutine run_check_indices
  ! ------------------------------------------------------------
  subroutine set_time_real(this, ut)
    class(ieModel) :: this
    real(kind=Real8_), intent(in) :: ut
    integer :: iYear
    integer :: iMonth
    integer :: iDay
    integer, dimension(7) :: itime
    real :: rHour
  end subroutine set_time_real

end module ModIEGITM
