! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!----------------------------------------------------------------------------
! $Id: output_common.f90,v 1.63 2017/10/24 14:23:20 ridley Exp $
!
! Author: Aaron Ridley, UMichigan
!
! Comments: Routines to output binary files
!
! AGB 3/31/13: Added 1D routine to output data at a specific altitude
! AJR 8/28/13: The code was outputting data on all processors for satellite
!              files.  I corrected this to make it so that if the linear
!              interpolation routine returns -1, the processor returns.
! AGB 10/18/13: Added gravity, collision frequency, and pressure gradient
!               to 3DION output
! AGB 12/20/13: Removed gravity from 3DION output
! AJR 06/30/17: Added in 3DLST - A list of variables
!----------------------------------------------------------------------------

integer function bad_outputtype()

  use ModInputs, only: OutputType, nOutputTypes
  use ModOutputProducers, only: is_known_output_type

  implicit none

  integer :: iOutputType

  do iOutputType = 1, nOutputTypes
    if (.not. is_known_output_type(OutputType(iOutputType))) then
      bad_outputtype = iOutputType
      return
    endif
  enddo

  bad_outputtype = 0
  return

end function bad_outputtype

!----------------------------------------------------------------
! Comments: Asad added data to allow output from RCAC
!----------------------------------------------------------------

subroutine output(dir, iBlock, iOutputType)

  use ModSphereInterface, only: iStartBlk
  use ModSatellites, only: CurrentSatellitePosition, CurrentSatelliteName, &
                           CurrSat, SatAltDat
  use ModGITM
  use ModEUV
  use ModTime
  use ModInputs
  use ModSources
  ! ModUserGITM: USR type var counts now managed via init_usr_output_registry
  use ModRCMR, only: RCMRFlag
  use ModConstants, only: pi
  use ModElectrodynamics, only: nMagLats, nMagLons
  use ModOutputBackend, only: ActiveBackend
  use ModOutputContainer, only: OutputContainer, GRID_HIME
  use ModOutputProducers
  use ModGITMVersion, only: GitmVersion

  implicit none

  character(len=*), intent(in) :: dir
  integer, intent(in) :: iBlock
  integer, intent(in) :: iOutputType

  character(len=5) :: proc_str, cBlock, cType
  character(len=24) :: cTime = '', cTimeSave = ''
  integer :: iiLat, iiLon, iiAlt, nGCs, cL = 0
  integer :: iLon, iLat, iAlt, nlines, iBLK, iSpecies, i
  logical :: done, IsFirstTime = .true., IsThere, DoSaveHIMEPlot

  real :: LatFind, LonFind, AltFind
  real :: rLon, rLat, rAlt

  character(len=2) :: cYear, cMonth, cDay, cHour, cMinute, cSecond
  character(len=4) :: cYearL

  integer :: iContIdx

  !! construct naming strings

  if (iOutputType == -1) then
    cType = "1DALL"
  else if (iOutputType == -2) then
    cType = "0DALL"
  else if (iOutputType == -3) then
    cType = "1DUSR"
  else if (iOutputType == -4) then
    cType = "0DUSR"
  else
    cType = OutputType(iOutputType)
    if (cType(1:2) == "3D" .and. Is1D) then
      cType(1:2) = "1D"
    endif
  endif

  if (Is1D) then
    iiLat = 1
    iiLon = 1
    rLon = 1.0
    rLat = 1.0
  endif

  ! If there are satellites, initialize the current satellite so that
  ! the maximum from all processors will contain the real value.  This is
  ! done by setting the current value to something rediculously small for
  ! all currently known satellite input data types.

  if (CurrSat > 0 .and. RCMRFlag) then
    SatAltDat(CurrSat) = -1.0e32
  endif

  if (iOutputType <= -1) then
    LatFind = CurrentSatellitePosition(iNorth_)
    LonFind = CurrentSatellitePosition(iEast_)
    call BlockLocationIndex(LonFind, LatFind, iBlock, iiLon, iiLat, rLon, rLat)

    if (iOutputType == -2 .or. iOutputType == -4) then
      AltFind = CurrentSatellitePosition(iUp_)
      call BlockAltIndex(AltFind, iBlock, iiLon, iiLat, iiAlt, rAlt)

      if (iiAlt < 0) return
    endif

    if (iDebugLevel > 2) then
      write(*, *) 'For BlockLocationIndex:'
      write(*, *) 'LonFind, LatFind = ', LonFind, LatFind
      write(*, *) 'Found iBlock, iiLon, iiLat, rLon, rLat =', &
        iBlock, iiLon, iiLat, rLon, rLat
    endif

    if (iiLon < 0 .or. iiLat < 0) return
  endif

  call init_output_containers()
  iContIdx = get_container_idx(cType)
  if (iContIdx < 0) then
    write(*, *) "WARNING: unknown output type: ", cType
    return
  endif

  ! For regional outputs (GRID_HIME), check if this block falls within the
  ! user-defined region.  Non-zero processors return early if outside region;
  ! processor 0 continues to write the header file even outside region.
  if (containers(iContIdx)%gridKind == GRID_HIME) then
    DoSaveHIMEPlot = .false.
    ! Check if any cell of this block falls into the user-defined region
    do iLat = -1, nLats + 2
      do iLon = -1, nLons + 2
        if (Longitude(iLon, iBlock) >= HIMEPlotLonStart*pi/180. &
            .and. Longitude(iLon, iBlock) <= HIMEPlotLonEnd*pi/180. &
            .and. Latitude(iLat, iBlock) >= HIMEPlotLatStart*pi/180. &
            .and. Latitude(iLat, iBlock) <= HIMEPlotLatEnd*pi/180.) then
          DoSaveHIMEPlot = .true.
          EXIT
        endif
      enddo
    enddo
    ! Non-zero processors skip if outside region
    if (iProc /= 0 .and. .not. DoSaveHIMEPlot) return
  endif

  if ((iProc == 0 .and. iBlock == 1) .and. (iOutputType > -1)) &
    write(*, '(a,i7,i5,5i3)') &
    "Writing Output files ("//cType//") at iStep : ", &
    iStep, iTimeArray(1:6)

  if (iOutputType <= -1) &
    write(*, '(a,i7,i5,5i3)') &
    "Writing satellite file ("//trim(CurrentSatelliteName)//") at iStep : ", &
    iStep, iTimeArray(1:6)

  call calc_physics(iBlock)
  call calc_rates(iBlock)
  call calc_collisions(iBlock)
  call chapman_integrals(iBlock)
  call set_horizontal_bcs(iBlock)
  if (.not. Is1D) call calc_efield(iBlock)

  iBLK = iStartBLK + iBlock

  write(cBlock, '(a1,i4.4)') "b", iBLK

  call i2s(mod(iTimeArray(1), 100), cYear, 2)
  call i2s(iTimeArray(1), cYearL, 4)
  call i2s(iTimeArray(2), cMonth, 2)
  call i2s(iTimeArray(3), cDay, 2)
  call i2s(iTimeArray(4), cHour, 2)
  call i2s(iTimeArray(5), cMinute, 2)
  call i2s(iTimeArray(6), cSecond, 2)

  if (.not. UseSecondsInFilename) cSecond = '00'  !xianjing

  !-----------
  ! New feature - we want to be able to write to the same file over and
  ! over and over again, so we don't get 100,000,000 satellite files.
  ! So, here we are going to name the file the first time, then open
  ! that same file over and over again with an append.
  !-----------

  if (IsFirstTime .or. &
      .not. DoAppendFiles .or. &
      (iOutputType > -1 .and. .not. Is1D)) then
    if (UseCCMCFileName) then
      cTime = "GITM_"//cYearL//"-"//cMonth//"-"//cDay//"T" &
              //cHour//"-"//cMinute//"-"//cSecond
      cL = 24
    else
      cTime = "t"//cYear//cMonth//cDay//"_"//cHour//cMinute//cSecond
      cL = 14
    endif

    if (IsFirstTime) cTimeSave = cTime
  else
    cTime = cTimeSave
  endif

  !! ---------------------------------------------
  !! Write the binary data files
  !! ---------------------------------------------

  if (iOutputType <= -1) then
    ! Satellite output: always use sequential per-file I/O regardless of backend
    inquire(file=trim(dir)//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".sat", &
            EXIST=IsThere)
    if (.not. DoAppendFiles .or. tSimulation < 0.1 .or. .not. IsThere) then
      open(unit=iOutputUnit_, form="unformatted", &
           file=trim(dir)//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".sat", &
           status="unknown")
    else
      open(unit=iOutputUnit_, form="unformatted", &
           file=trim(dir)//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".sat", &
           status="unknown", position='append')
    endif
  elseif (ActiveBackend%uses_external_file) then
    ! Legacy backend: open per-block file here.
    ! For HME type output, open file only if DoSaveHIMEPlot=T. This is for iProc=0,
    ! because other iProcs with DoSaveHIMEPlot=F exit the subroutine earlier.
    if ((cType(3:5) /= 'HME' .and. (cType /= '2DMEL' .or. iBLK == 1)) .or. &
        (cType(3:5) == 'HME' .and. DoSaveHIMEPlot)) then
      inquire(file=trim(dir)//"/"//cType//"_"//cTime(1:cL)//"."//cBlock, &
              EXIST=IsThere)
      if (.not. DoAppendFiles .or. tSimulation < 0.1 .or. .not. IsThere) then
        open(unit=iOutputUnit_, form="unformatted", &
             file=trim(dir)//"/"//cType//"_"//cTime(1:cL)//"."//cBlock, &
             status="unknown")
      else
        open(unit=iOutputUnit_, form="unformatted", &
             file=trim(dir)//"/"//cType//"_"//cTime(1:cL)//"."//cBlock, &
             status="unknown", position='append')
      endif
    endif
  else
    ! Collective backends (e.g., mpiio, netcdf): the shared output file is opened
    ! collectively in write_output() around the block loop; nothing to do here.
  endif

  nGCs = containers(iContIdx)%nGhostCells
  if (containers(iContIdx)%nVars > 0) then
    if (cType == '2DMEL') then
        call fill_2dmel(containers(iCont_2DMEL), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_2DMEL), iBLK, iTimeArray)
        call containers(iCont_2DMEL)%reset()
      elseif (cType == '2DTEC') then
        call fill_2dtec(containers(iCont_2DTEC), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_2DTEC), iBLK, iTimeArray)
        call containers(iCont_2DTEC)%reset()
      elseif (cType == '2DGEL') then
        call fill_2dgel(containers(iCont_2DGEL), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_2DGEL), iBLK, iTimeArray)
        call containers(iCont_2DGEL)%reset()
      elseif (cType == '3DNEU') then
        call fill_3dneu(containers(iCont_3DNEU), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DNEU), iBLK, iTimeArray)
        call containers(iCont_3DNEU)%reset()
      elseif (cType == '3DALL') then
        call fill_3dall(containers(iCont_3DALL), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DALL), iBLK, iTimeArray)
        call containers(iCont_3DALL)%reset()
      elseif (cType == '3DION') then
        call fill_3dion(containers(iCont_3DION), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DION), iBLK, iTimeArray)
        call containers(iCont_3DION)%reset()
      elseif (cType == '3DTHM') then
        call fill_3dthm(containers(iCont_3DTHM), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DTHM), iBLK, iTimeArray)
        call containers(iCont_3DTHM)%reset()
      elseif (cType == '3DCHM') then
        call fill_3dchm(containers(iCont_3DCHM), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DCHM), iBLK, iTimeArray)
        call containers(iCont_3DCHM)%reset()
      elseif (cType == '3DLST') then
        call fill_3dlst(containers(iCont_3DLST), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DLST), iBLK, iTimeArray)
        call containers(iCont_3DLST)%reset()
      elseif (cType == '3DGLO') then
        call fill_3dglo(containers(iCont_3DGLO), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DGLO), iBLK, iTimeArray)
        call containers(iCont_3DGLO)%reset()
      elseif (cType == '3DMAG') then
        call fill_3dmag(containers(iCont_3DMAG), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DMAG), iBLK, iTimeArray)
        call containers(iCont_3DMAG)%reset()
      elseif (cType == '3DHME') then
        call fill_3dhme(containers(iCont_3DHME), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DHME), iBLK, iTimeArray)
        call containers(iCont_3DHME)%reset()
      elseif (cType == '3DMOH') then
        call fill_3dmoh(containers(iCont_3DMOH), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DMOH), iBLK, iTimeArray)
        call containers(iCont_3DMOH)%reset()
      elseif (cType == '3DMOV') then
        call fill_3dmov(containers(iCont_3DMOV), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DMOV), iBLK, iTimeArray)
        call containers(iCont_3DMOV)%reset()
      elseif (cType == '2DANC') then
        call fill_2danc(containers(iCont_2DANC), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_2DANC), iBLK, iTimeArray)
        call containers(iCont_2DANC)%reset()
      elseif (cType == '2DHME') then
        call fill_2dhme(containers(iCont_2DHME), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_2DHME), iBLK, iTimeArray)
        call containers(iCont_2DHME)%reset()
      elseif (cType == '1DALL') then
        call fill_1dall(containers(iCont_1DALL), iBlock, iiLon, iiLat, rLon, rLat)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_1DALL), iBLK, iTimeArray)
        call containers(iCont_1DALL)%reset()
      elseif (cType == '0DALL') then
        call fill_0dall(containers(iCont_0DALL), iBlock, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_0DALL), iBLK, iTimeArray)
        call containers(iCont_0DALL)%reset()
      elseif (cType == '1DGLO') then
        call fill_1dglo(containers(iCont_1DGLO), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_1DGLO), iBLK, iTimeArray)
        call containers(iCont_1DGLO)%reset()
      elseif (cType == '1DTHM') then
        call fill_1dthm(containers(iCont_1DTHM), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_1DTHM), iBLK, iTimeArray)
        call containers(iCont_1DTHM)%reset()
      elseif (cType == '1DCHM') then
        call fill_1dchm(containers(iCont_1DCHM), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_1DCHM), iBLK, iTimeArray)
        call containers(iCont_1DCHM)%reset()
      elseif (cType == '1DNEW') then
        call fill_1dnew(containers(iCont_1DNEW), iBlock, iiLon, iiLat, rLon, rLat)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_1DNEW), iBLK, iTimeArray)
        call containers(iCont_1DNEW)%reset()
      elseif (cType == '3DUSR') then
        call fill_3dusr(containers(iCont_3DUSR), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_3DUSR), iBLK, iTimeArray)
        call containers(iCont_3DUSR)%reset()
      elseif (cType == '2DUSR') then
        call fill_2dusr(containers(iCont_2DUSR), iBlock)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_2DUSR), iBLK, iTimeArray)
        call containers(iCont_2DUSR)%reset()
      elseif (cType == '1DUSR') then
        call fill_1dusr(containers(iCont_1DUSR), iBlock, iiLon, iiLat, rLon, rLat)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_1DUSR), iBLK, iTimeArray)
        call containers(iCont_1DUSR)%reset()
      elseif (cType == '0DUSR') then
        call fill_0dusr(containers(iCont_0DUSR), iBlock, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
        if (associated(ActiveBackend%write_container)) &
          call ActiveBackend%write_container(containers(iCont_0DUSR), iBLK, iTimeArray)
        call containers(iCont_0DUSR)%reset()
      endif
  endif

  if (iOutputType <= -1 .or. ActiveBackend%uses_external_file) &
    close(unit=iOutputUnit_)

  !! Now write the header file (skipped for backends where metadata lives in the data file)


  if (ActiveBackend%writes_header .and. &
      ((iProc == 0 .and. iBlock == nBlocks) .or. iOutputType <= -1)) then

    if (iOutputType <= -1) then
      inquire(file=trim(dir)//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".header", EXIST=IsThere)
      if (.not. DoAppendFiles .or. tSimulation < 0.1 .or. .not. IsThere) then
        open(unit=iOutputUnit_, &
             file=trim(dir)//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".header", &
             status="unknown")
      else
        open(unit=iOutputUnit_, &
             file=trim(dir)//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".header", &
             status="unknown", position='append')
      endif
    else
      inquire(file=trim(dir)//"/"//cType//"_"//cTime(1:cL)//".header", &
              EXIST=IsThere)
      if (.not. DoAppendFiles .or. tSimulation < 0.1 .or. .not. IsThere) then
        open(unit=iOutputUnit_, &
             file=trim(dir)//"/"//cType//"_"//cTime(1:cL)//".header", &
             status="unknown")
      else
        open(unit=iOutputUnit_, &
             file=trim(dir)//"/"//cType//"_"//cTime(1:cL)//".header", &
             status="unknown", position='append')
      endif
    endif

    call ActiveBackend%write_header(iOutputUnit_, containers(iContIdx), cType, &
                                    nAlts, nLats, nLons, nGCs, nBlocksLat, nBlocksLon, iTimeArray, GitmVersion)

    close(unit=iOutputUnit_)

  endif

  IsFirstTime = .false.

contains

  !----------------------------------------------------------------

  subroutine write_head_blocks

    ! for 1D and 0D do not write blocks
    if (cType(1:2) == "1D" .or. cType(1:2) == "0D") return

    write(iOutputUnit_, *) "BLOCKS"
    write(iOutputUnit_, "(I7,A)") 1, " nBlocksAlt"
    if (cType /= "2DMEL") then
      write(iOutputUnit_, "(I7,A)") nBlocksLat, " nBlocksLat"
      write(iOutputUnit_, "(I7,A)") nBlocksLon, " nBlocksLon"
    else
      write(iOutputUnit_, "(I7,A)") 1, " nBlocksLat"
      write(iOutputUnit_, "(I7,A)") 1, " nBlocksLon"
    endif
    write(iOutputUnit_, *) ""

  end subroutine write_head_blocks

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine write_head_time

    write(iOutputUnit_, *) "TIME"
    write(iOutputUnit_, "(I7,A)") iTimeArray(1), " Year"
    write(iOutputUnit_, "(I7,A)") iTimeArray(2), " Month"
    write(iOutputUnit_, "(I7,A)") iTimeArray(3), " Day"
    write(iOutputUnit_, "(I7,A)") iTimeArray(4), " Hour"
    write(iOutputUnit_, "(I7,A)") iTimeArray(5), " Minute"
    write(iOutputUnit_, "(I7,A)") iTimeArray(6), " Second"
    write(iOutputUnit_, "(I7,A)") iTimeArray(7), " Millisecond"
    write(iOutputUnit_, *) ""

  end subroutine write_head_time

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine write_head_version

    use ModGITMVersion

    write(iOutputUnit_, *) "VERSION"
    write(iOutputUnit_, *) trim(GitmVersion)
    write(iOutputUnit_, *) ""

  end subroutine write_head_version

end subroutine output
