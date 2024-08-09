!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModAMIE_Interface

  use ModCharSize

  implicit none
  
  character (len=iCharLenIE_) :: AMIE_FileName
  integer :: AMIE_nLats, AMIE_nMlts, AMIE_nTimes, nCellsPad, AMIE_nBlks

  real :: rDummyeFlux = 1.0e-6 ! in W/m2
  real :: rDummyAveE = 2.0 ! in keV
  real :: rDummyIonAveE = 20.0 ! in keV
  
  ! For a single file
  real*4, allocatable,dimension(:) :: AMIE_Lats, AMIE_MLTs
  real*8, allocatable,dimension(:,:) :: AMIE_Time

  real*4, allocatable,dimension(:,:,:,:) :: AMIE_Value
  real, allocatable, dimension(:, :, :) :: AMIE_Interpolated

  integer, parameter :: AMIE_Closest_     = 1
  integer, parameter :: AMIE_After_       = 2
  integer, parameter :: AMIE_Interpolate_ = 3

  integer :: AMIE_iDebugLevel = 0

  integer :: AMIE_South_ = 1
  integer :: AMIE_North_ = 2
  
  integer, parameter :: iPotential_ = 1
  integer, parameter :: iPotentialy_ = 2
  integer, parameter :: iEle_diff_eflux_ = 3
  integer, parameter :: iEle_diff_avee_ = 4
  integer, parameter :: iIon_diff_eflux_ = 5
  integer, parameter :: iIon_diff_avee_ = 6
  integer, parameter :: iEle_mono_eflux_ = 7
  integer, parameter :: iEle_mono_avee_ = 8
  integer, parameter :: iEle_wave_eflux_ = 9
  integer, parameter :: iEle_wave_avee_ = 10
  logical :: useIons
  logical :: useMono
  logical :: useWave

  integer, parameter :: nValues = 10
  character (len=iCharLenIE_) :: AMIE_Names(nValues)
  integer, dimension(nValues) :: iMap_ = -1
  real, dimension(nValues) :: unitConvert
  
  ! mlts, lats, blocks, times, and variable
  real*4, allocatable,dimension(:, :, :, :, :) :: AMIE_Storage

contains

  subroutine AMIE_link_variable_names()

    AMIE_Names(iPotential_) = "Potential"
    AMIE_Names(iPotentialy_) = "PotentialY"
    AMIE_Names(iEle_diff_eflux_) = "Electron Energy Flux"
    AMIE_Names(iEle_diff_avee_) = "Electron Mean Energy"
    AMIE_Names(iIon_diff_eflux_) = "Ion Energy Flux"
    AMIE_Names(iIon_diff_avee_) = "Ion Mean Energy"
    AMIE_Names(iEle_mono_eflux_) = "ME Energy Flux"
    AMIE_Names(iEle_mono_avee_) = "ME Mean Energy"
    AMIE_Names(iEle_wave_eflux_) = "BB Energy Flux"
    AMIE_Names(iEle_wave_avee_) = "BB Mean Energy"

    ! By default, we don't need any unit conversions:
    unitConvert = 1.0
    useIons = .false.
    useMono = .false.
    useWave = .false.
    
  end subroutine AMIE_link_variable_names

  ! --------------------------------------------------
  ! Allocate all of the variables
  ! --------------------------------------------------
  
  subroutine AMIE_allocate_variables()

    implicit none

    integer :: iError = 0
    
    if (allocated(AMIE_Lats)) deallocate(AMIE_Lats)
    allocate(AMIE_Lats(AMIE_nLats + nCellsPad), stat=iError)
    if (iError /= 0) then
       write(*,*) "Error in allocating array AMIE_Lats in "
       stop
    endif

    if (allocated(AMIE_Mlts)) deallocate(AMIE_Mlts)
    allocate(AMIE_Mlts(AMIE_nMlts), stat=iError)
    if (iError /= 0) then
       write(*,*) "Error in allocating array Mlts in "
       stop
    endif
    
    ! ------------------------------------------------------------
    ! This is the variable that holds all of the variables to allo
    ! for a more generalized get variable function

    if (AMIE_iDebugLevel > 1) &
         write(*,*) 'Allocating Variables : ', &
         AMIE_nMlts, AMIE_nLats, nCellsPad, AMIE_nTimes, AMIE_nBlks
    
    allocate(AMIE_Storage( &
         nValues, &
         AMIE_nMlts, &
         AMIE_nLats + nCellsPad, &
         AMIE_nTimes, &
         AMIE_nBlks), &
         stat = iError)
    if (iError /= 0) then
       write(*,*) "Error in allocating array AMIE_Storage in readAMIEoutput"
       stop
    endif

    if (AMIE_iDebugLevel > 1) write(*,*) '  --> allocated the big array!'
    
    ! Initialize the variables to background values
    AMIE_Storage(iPotential_, :, :, :, :) = 0.0
    AMIE_Storage(iPotentialY_, :, :, :, :) = 0.0
    AMIE_Storage(iEle_diff_eflux_, :, :, :, :) = rDummyeFlux
    AMIE_Storage(iEle_diff_avee_, :, :, :, :) = rDummyAveE
    AMIE_Storage(iIon_diff_eflux_, :, :, :, :) = rDummyeFlux
    AMIE_Storage(iIon_diff_avee_, :, :, :, :) = rDummyIonAveE
    AMIE_Storage(iEle_mono_eflux_, :, :, :, :) = rDummyeFlux
    AMIE_Storage(iEle_mono_avee_, :, :, :, :) = rDummyAveE
    AMIE_Storage(iEle_wave_eflux_, :, :, :, :) = rDummyeFlux
    AMIE_Storage(iEle_wave_avee_, :, :, :, :) = rDummyAveE

    ! ------------------------
    ! Generic AMIE Value:

    allocate(AMIE_Value( &
         AMIE_nMlts, AMIE_nLats+nCellsPad, AMIE_nTimes, AMIE_nBlks), stat=iError)
    if (iError /= 0) then
       write(*,*) "Error in allocating array AMIE_Value in "
       stop
    endif

    AMIE_Value = 0.0

    ! ---------------------------------------------------
    ! AMIE Interpolated Value (this is for one time)

    allocate(AMIE_Interpolated( &
         AMIE_nMlts, AMIE_nLats+nCellsPad, AMIE_nBlks), stat=iError)
    if (iError /= 0) then
       write(*,*) "Error in allocating array AMIE_Interpolated in "
       stop
    endif

    AMIE_Value = 0.0

    ! ------------------------
    ! Time:

    allocate(AMIE_Time(AMIE_nTimes,2), stat=iError)
    if (iError /= 0) then
       write(*,*) "Error in allocating array AMIETimes in "
       stop
    endif

    if (AMIE_iDebugLevel > 1) write(*,*) '  --> done allocating!'
     
  end subroutine AMIE_allocate_variables

  subroutine AMIE_link_vars_to_keys(nFields, Fields)

    implicit none
    
    integer, parameter :: nFieldsMax = 100
    integer, intent(in) :: nFields
    character (len=30), dimension(nFieldsMax), intent(in) :: Fields
    
    integer :: iField, iVal
    do iVal = 1, nValues
       if (AMIE_iDebugLevel > 1) &
            write(*,*) 'Searching for ', trim(AMIE_Names(iVal))
       do iField = 1, nFields
          if (index(Fields(iField), trim(AMIE_Names(iVal))) > 0) then
             iMap_(iVal) = iField
             if (AMIE_iDebugLevel > 1) &
                  write(*,*) "<--- ", trim(AMIE_Names(iVal)), " Found", iField
          endif
       enddo
    enddo
    if (AMIE_iDebugLevel > 1) then
       write(*, *) 'Summary of variables found in AMIE file : '
       do iVal = 1, nValues
          if (iMap_(iVal) > 0) then
             write(*, *) 'Expected : ', trim(AMIE_Names(iVal)), &
                  '; found : ', trim(Fields(iMap_(iVal)))
          endif
       enddo
    endif

    ! Check to see if some of these exist:
    if (iMap_(iPotential_) < 1) then
       call report_error_and_crash(1, "Could not find Potential in file!")
    endif
    if ((iMap_(iEle_diff_eflux_) < 1) .or. (iMap_(iEle_diff_avee_) < 0)) then
       call report_error_and_crash(1, "Could not find Electron Diffuse in file!")
    endif
       
    if ((iMap_(iIon_diff_eflux_) > 0) .and. (iMap_(iIon_diff_avee_) > 0)) then
       useIons = .true.
       if (AMIE_iDebugLevel > 1) write(*, *) "Input Electrodynamics is using Ions!"
    else
       useIons = .false.
    endif
    if ((iMap_(iEle_mono_eflux_) > 0) .and. (iMap_(iEle_mono_avee_) > 0)) then
       useMono = .true.
       if (AMIE_iDebugLevel > 1) write(*, *) "Input Electrodynamics is using Mono!"
    else
       useMono = .false.
    endif
    if ((iMap_(iEle_wave_eflux_) > 0) .and. (iMap_(iEle_wave_avee_) > 0)) then
       useWave = .true.
       if (AMIE_iDebugLevel > 1) write(*, *) "Input Electrodynamics is using Ions!"
    else
       useWave = .false.
    endif

  end subroutine AMIE_link_vars_to_keys

end Module ModAMIE_Interface
