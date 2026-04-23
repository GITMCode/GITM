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
  use ModOutputRegistry, only: find_output_type, nRegisteredTypes, &
                               init_output_registry

  implicit none

  integer :: iOutputType

  ! Ensure the registry is populated (safe to call multiple times;
  ! init_output_registry resets and re-registers on each call).
  if (nRegisteredTypes == 0) call init_output_registry()

  ! Validate each requested output type against the registry.
  do iOutputType = 1, nOutputTypes
    if (find_output_type(OutputType(iOutputType)) < 1) then
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
  use ModOutputRegistry, only: OutputTypeInfo, find_output_type, RegisteredTypes
  use ModOutputGather, only: gather_output
  use ModOutputBackend, only: ActiveBackend
  use ModGITMVersion, only: GitmVersion

  implicit none

  character(len=*), intent(in) :: dir
  integer, intent(in) :: iBlock
  integer, intent(in) :: iOutputType

  character(len=5) :: proc_str, cBlock, cType
  character(len=24) :: cTime = '', cTimeSave = ''
  integer :: iiLat, iiLon, iiAlt, nGCs, cL = 0
  integer :: iLon, iLat, iAlt, nVars_to_Write, nlines, iBLK, iSpecies, i
  logical :: done, IsFirstTime = .true., IsThere, DoSaveHIMEPlot

  real :: LatFind, LonFind, AltFind
  real :: rLon, rLat, rAlt

  character(len=2) :: cYear, cMonth, cDay, cHour, cMinute, cSecond
  character(len=4) :: cYearL

  ! New: registry-based output infrastructure
  type(OutputTypeInfo) :: typeInfo
  integer :: iTypeIdx
  integer :: nV, nX, nY, nZ
  real, allocatable :: outBuf(:, :, :, :)

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

  ! Xing Meng 2020-03-16: HIME type output, save output for user-defined region only
  ! It only works under the assumption of one block per processor
  if (cType(3:5) == "HME") then
    DoSaveHIMEPlot = .false.
    ! Save the current block if any cell of this block falls into the
    ! user-defined region, set DoSaveHIMEPlot to true
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
    ! Do not write output at all if the current block is outside of the
    ! user-defined region. iProc=0 writes the header file
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

  ! (Registry and backend initialized in write_output before first block loop)

  nGCs = 2

  ! Registry-based output: lookup -> gather -> backend write (all types including USR)
  iTypeIdx = find_output_type(cType)
  if (iTypeIdx < 1) then
    write(*, *) "WARNING: output type not found in registry: ", cType
  else
    typeInfo = RegisteredTypes(iTypeIdx)
    nGCs = typeInfo%nGhostCells
    nV = typeInfo%nVars
    ! Compute buffer dimensions from type metadata
    select case (typeInfo%nDims)
    case (3)
      nX = nLons + 4; nY = nLats + 4; nZ = nAlts + 4
    case (2)
      if (typeInfo%usesMagGrid) then
        nX = nMagLons + 1; nY = nMagLats
      else
        nX = nLons; nY = nLats
      endif
      nZ = 1
    case (1)
      nX = 1; nY = 1; nZ = nAlts + 4
    case default  ! 0D
      nX = 1; nY = 1; nZ = 1
    end select
    ! Gather data and write -- respecting per-type skip conditions (HME region, 2DMEL block)
    if (nV > 0) then
      if (cType(3:5) == 'HME') then
        if (DoSaveHIMEPlot) then
          allocate(outBuf(nV, nX, nY, nZ))
          call gather_output(cType, iBlock, outBuf, nV, nX, nY, nZ, &
                             iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
          call ActiveBackend%write_block(iOutputUnit_, outBuf, nV, nX, nY, nZ, iBLK)
          deallocate(outBuf)
        endif
      elseif (cType == '2DMEL') then
        if (iBLK == 1) then
          allocate(outBuf(nV, nX, nY, nZ))
          call gather_output(cType, iBlock, outBuf, nV, nX, nY, nZ, &
                             iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
          call ActiveBackend%write_block(iOutputUnit_, outBuf, nV, nX, nY, nZ, iBLK)
          deallocate(outBuf)
        endif
      else
        allocate(outBuf(nV, nX, nY, nZ))
        call gather_output(cType, iBlock, outBuf, nV, nX, nY, nZ, &
                           iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
        call ActiveBackend%write_block(iOutputUnit_, outBuf, nV, nX, nY, nZ, iBLK)
        deallocate(outBuf)
      endif
    endif
    nvars_to_write = nV
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

      ! Auto-generated header from registry metadata (all types including USR)
      call ActiveBackend%write_header(iOutputUnit_, typeInfo, cType, &
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
