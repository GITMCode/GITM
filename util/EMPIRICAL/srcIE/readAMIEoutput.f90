!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf


subroutine readAMIEvariables(AmieFile, nVars, varNames, iDebugLevel)

  use ModCharSize
  implicit none
  integer, parameter :: nFieldsMax = 100
  character (len=iCharLenIE_), intent(in) :: AmieFile
  integer, intent(out) :: nVars
  character (len=30), dimension(nFieldsMax), intent(out) :: varNames
  integer, intent(in) :: iDebugLevel
  
  integer :: iUnitTmp_ = 77
  integer :: iError = 0, iVar, i
  real*4 :: dummy(1000)
  
  integer :: nlats, nmlts, ntimes
  open(iUnitTmp_, &
       file = AmieFile, &
       status = 'old', &
       form = 'UNFORMATTED', &
       iostat = iError)
  read(iUnitTmp_) nlats, nmlts, ntimes
  read(iUnitTmp_) (dummy(i), i = 1, nlats)
  read(iUnitTmp_) (dummy(i), i = 1, nmlts)
  read(iUnitTmp_) nVars
  do iVar = 1, nVars
     read(iUnitTmp_) varNames(iVar)
  enddo
  close(iUnitTmp_)

  if (iDebugLevel > 1) then
     write(*, *) "AMIE Variables : " 
     do iVar = 1, nVars
        write(*,*) iVar, '. ', trim(varNames(iVar))
     enddo
  endif
  return
end subroutine readAMIEvariables

subroutine readAMIEoutput(iBLK, IsMirror, iDebugGitm, iError)

  use ModAMIE_Interface
  use ModEIEFiles
  use ModTimeConvert, ONLY: time_int_to_real
  use ModIoUnit, ONLY : UNITTMP_
  implicit none

  integer, intent(out) :: iError
  logical, intent(in)  :: IsMirror
  integer, intent(in)  :: iBLK
  integer, intent(in)  :: iDebugGitm

  integer :: iTime, nTimesBig, nTimesTotal
  integer :: nfields
  integer :: ntemp, iyr, imo, ida, ihr, imi
  integer :: i,j, iField, iVal, iPot_, iAveE_, iEFlux_, iPotY_
  integer :: ieleMonoEFlux_ = -1
  integer :: ieleMonoAveE_ = -1
  integer :: ieleWaveEFlux_ = -1
  integer :: ieleWaveAveE_ = -1
  integer :: iIonEFlux_ = -1
  integer :: iIonAveE_ = -1
  real*4  :: swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot
  real*8  :: rtime
  integer, dimension(7) :: itime_i

  real*4, allocatable, dimension(:,:,:) :: AllData
  real*4, allocatable, dimension(:) :: TempLats
  integer, parameter :: nFieldsMax = 100
  character (len=30), dimension(nFieldsMax) :: Fields

  logical :: IsBinary, energyfluxconvert, ReverseLats = .false.

  real :: dPotential, dPotentialY
  integer :: n

  real :: factor
  nTimesBig = 0
  nTimesTotal = 0
  
  AMIE_iDebugLevel = iDebugGitm

  iError = 0
  if (AMIE_iDebugLevel >= 0) write(*,*) '> reading AMIE file : ', trim(AMIE_FileName)
  open(UnitTmp_, &
       file=AMIE_FileName, &
       status='old', &
       form='UNFORMATTED', &
       iostat=iError)
  if (iError.ne.0) then
     write(*,*) "Error opening file:", trim(AMIE_FileName)
     stop
  endif
  AMIE_nLats = 0
  IsBinary = .true.

  read(UnitTmp_, iostat = iError) AMIE_nlats, AMIE_nmlts, AMIE_ntimes
  if ((iError.ne.0).or.(AMIE_nlats.gt.500)) then
     write(*,*) "Error reading variables AMIE_nlats, AMIE_nmlts, AMIE_ntimes"
     IsBinary = .false.
  endif
  close(UnitTmp_)

  if (IsBinary) then
     call AMIE_link_variable_names()
     call readAMIEvariables(AMIE_FileName, nFields, Fields, AMIE_iDebugLevel)
     call AMIE_link_vars_to_keys(nFields, Fields)
     
     open(UnitTmp_, file = AMIE_FileName, status='old', form='UNFORMATTED')
     read(UnitTmp_) AMIE_nlats, AMIE_nmlts, AMIE_ntimes
  else
     open(UnitTmp_, file = AMIE_FileName, status='old')
     read(UnitTmp_,*) AMIE_nlats, AMIE_nmlts, AMIE_ntimes
  endif

  !\
  ! We have run into a problem with AMIE during storms.
  ! The potential is not zero at the boundary.  It is sometimes quite
  ! high.  This means that the gradient in the potential will be
  ! large - meaning that very strong flows can exist in this last
  ! cell.  This is not good.
  ! To rectify this, we will pad the AMIE results by 15 grid cells, and
  ! force the potential to go to zero linearly from the last cell to the
  ! new last cell.
  ! Since it is assumed that all of the AMIE quantities are on the same
  ! grid, we have to fill in the eflux and avee also.  We will use the
  ! last cell to fill in those value.
  ! We also have to extend the grid.
  !/

  nCellsPad = 15
  AMIE_nBlks = 2
  
  if (.not.allocated(AMIE_Storage)) call AMIE_allocate_variables()

  ! these are local variables:

  if (allocated(AllData)) deallocate(AllData)
  allocate(AllData(AMIE_nMlts, AMIE_nLats, nFields), stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array AllData in "
     stop
  endif

  ! ---------------------------------------------------------------
  ! Read and extrapolate AMIE grid

  if (AMIE_iDebugLevel > 1) write(*,*) 'Reading AMIE grid'
  
  if (IsBinary) then
     read(UnitTmp_) (AMIE_Lats(i), i = 1, AMIE_nLats)
     read(UnitTmp_) (AMIE_Mlts(i), i = 1, AMIE_nMlts)
     read(UnitTmp_) nFields
  else
     read(UnitTmp_,*) (AMIE_Lats(i), i = 1, AMIE_nLats)
     read(UnitTmp_,*) (AMIE_Mlts(i), i = 1, AMIE_nMlts)
     read(UnitTmp_,*) nFields
  endif

  if (nFields > nFieldsMax) then
     write(*,*) "Maximum number of fields in AMIE is ",nFieldsMax
     stop
  endif
  
  AMIE_Lats = 90.0 - AMIE_Lats

  !\
  ! Extrapolate the latitude grid
  !/

  if (AMIE_iDebugLevel > 1) write(*,*) 'Extrapolating AMIE grid'
  if (AMIE_Lats(AMIE_nLats) > AMIE_Lats(1)) then
     ! The AMIE data is in reverse order than what we want, so let's reverse it
     if (allocated(TempLats)) deallocate(TempLats)
     allocate(TempLats(AMIE_nLats+nCellsPad), stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array TempLats in "
        stop
     endif
     TempLats = AMIE_Lats
     do i=1,AMIE_nLats
        AMIE_Lats(i) = TempLats(AMIE_nLats+1-i)
     enddo
     ReverseLats = .true.
     deallocate(TempLats)
  endif

  do i=AMIE_nLats+1,AMIE_nLats+nCellsPad
     AMIE_Lats(i) = AMIE_Lats(i-1) + &
          (AMIE_Lats(AMIE_nLats) - AMIE_Lats(AMIE_nLats-1))
  enddo

  energyfluxconvert = .false.

  do iField=1,nfields
     if (IsBinary) then
        read(UnitTmp_) Fields(iField)
     else
        read(UnitTmp_,'(a)') Fields(iField)
     endif
  enddo

  ! Need to convert from W/m^2 to erg/cm2/s
  factor = 1.0
  if (energyfluxconvert) then
     factor = 1.0 / (1.0e-7 * 100.0 * 100.0)
     unitConvert(iEle_diff_eflux_) = factor
     unitConvert(iIon_diff_eflux_) = factor
     unitConvert(iEle_mono_eflux_) = factor
     unitConvert(iEle_wave_eflux_) = factor
  endif

  if (AMIE_iDebugLevel > 1) write(*,*) 'Reading AMIE data now...'
    
  do iTime=1,AMIE_ntimes

     if (AMIE_iDebugLevel > 1) write(*,*) ' --> iTime : ', iTime
     
     if (IsBinary) then

        read(UnitTmp_) ntemp,iyr,imo,ida,ihr,imi
        read(UnitTmp_) swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot

        do iField=1,nfields
           if (ReverseLats) then 
              read(UnitTmp_) &
                   ((AllData(j,i,iField),j=1,AMIE_nMlts),i=AMIE_nLats,1,-1)
           else
              read(UnitTmp_) &
                   ((AllData(j,i,iField),j=1,AMIE_nMlts),i=1,AMIE_nLats)
           endif
        enddo

     else

        read(UnitTmp_,*) ntemp,iyr,imo,ida,ihr,imi
        read(UnitTmp_,*) swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot

        do iField=1,nfields
           if (ReverseLats) then 
              read(UnitTmp_,*) &
                   ((AllData(j,i,iField),j=1,AMIE_nMlts),i=AMIE_nLats,1,-1)
           else
              read(UnitTmp_,*) &
                   ((AllData(j,i,iField),j=1,AMIE_nMlts),i=1,AMIE_nLats)
           endif
        enddo

     endif

     itime_i(1) = iyr
     itime_i(2) = imo
     itime_i(3) = ida
     itime_i(4) = ihr
     itime_i(5) = imi
     itime_i(6) = 0
     itime_i(7) = 0
     call time_int_to_real(itime_i,rtime)     
     AMIE_Time(iTime,iBLK) = rtime

     ! We need Potential to be in Volts
     !         AveE to be in keV
     !         EFlux to be in W/m2

     if (AMIE_iDebugLevel > 1) write(*,*) '  --> Pushing in to storage '
     do iVal = 1, nValues
        if (iMap_(iVal) > 0) then
           ! if we are mirroring, we need to do something special
           !  but only if we are talking about the potential:
           if ((IsMirror) .and. &
                ((iVal == iPotential_) .or. (iVal == iPotentialY_))) then
              do i=1,AMIE_nMlts
                 AMIE_Storage(iVal, i, 1:AMIE_nLats, iTime, iBLK) = &
                      - AllData(AMIE_nMlts+1 - i, 1:AMIE_nLats, iMap_(iVal))
              enddo
           else
              ! This is the default :
              !    (need loop for optimization issues in gfortran...)
              do i=1,AMIE_nMlts
                 AMIE_Storage(iVal, i, 1:AMIE_nLats, iTime, iBLK) = &
                      AllData(i, 1:AMIE_nLats, iMap_(iVal))
              enddo
           endif

           ! If the variable is the potential, linearly interpolate it to lower
           ! latitudes, and then smooth it in mlt:
           if ((iVal == iPotential_) .or. (iVal == iPotentialY_)) then
              ! First extend the potential:
              do i=1,AMIE_nMlts
                 dPotential = AMIE_Storage(iVal, i, AMIE_nLats, iTime, iBLK) / nCellsPad
                 do j=AMIE_nLats+1, AMIE_nLats+nCellsPad
                    AMIE_Storage(iVal, i, j, iTime, iBLK) = &
                         AMIE_Storage(iVal, i, j-1, iTime, iBLK) - dPotential
                 enddo
              enddo
              ! Then smooth the extension in MLT:
              do j=AMIE_nLats+1, AMIE_nLats+nCellsPad
                 ! We are going to smooth more and more as we go down in latitude
                 do n = 1, j-AMIE_nLats
                    i = 1
                    AMIE_Storage(iVal, i, j, iTime, iBLK) = &
                         (AMIE_Storage(iVal, AMIE_nMlts-1, j, iTime, iBLK) + &
                         2*AMIE_Storage(iVal, i, j, iTime, iBLK) + &
                         AMIE_Storage(iVal, i+1, j, iTime, iBLK))/4.0

                    do i=2,AMIE_nMlts-1
                       AMIE_Storage(iVal, i,j,iTime,iBLK) = &
                            (AMIE_Storage(iVal, i-1, j, iTime, iBLK) + &
                            2*AMIE_Storage(iVal, i, j, iTime, iBLK) + &
                            AMIE_Storage(iVal, i+1, j, iTime, iBLK))/4.0
                    enddo

                    i = AMIE_nMlts
                    AMIE_Storage(iVal, i, j, iTime, iBLK) = &
                         (AMIE_Storage(iVal, i-1, j, iTime, iBLK) + &
                         2*AMIE_Storage(iVal, i, j, iTime, iBLK) + &
                         AMIE_Storage(iVal, 2, j, iTime, iBLK))/4.0

                 enddo
              enddo
           endif
        endif
     enddo
     if (AMIE_iDebugLevel > 1) write(*,*) '  --> Done with time '
  enddo
  
  if (AMIE_iDebugLevel > 1) write(*,*) 'Done Reading AMIE file '
  close(UnitTmp_)

  ! update the number of latitudes to include padded cells
  AMIE_nLats = AMIE_nLats + nCellsPad

  deallocate(AllData, stat=iError)

end subroutine readAMIEoutput
