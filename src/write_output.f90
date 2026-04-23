! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!-----------------------------------------------------------------------------
! $Id: write_output.f90,v 1.16 2013/10/12 04:01:00 kopmanis Exp $
!
! Author: Aaron Ridley, UMichigan
!
! Comments: Routines to write data output (both the desired model output files
!           and satellite output).  Also calls the logfile output routine.
!-----------------------------------------------------------------------------

subroutine write_output

  use ModTime
  use ModInputs
  use ModGITM
  use ModOutputRegistry, only: init_output_registry
  use ModOutputBackend, only: ActiveBackend, init_output_backend, requestedBackendName

  implicit none

  real, external :: get_timing
  real :: ProjectedTime, CompletedTime, RealTime
  integer :: i, nMLTsTmp, nLatsTmp, iBlock

  ! For cTime construction and collective open/close
  character(len=5)  :: cType
  character(len=24) :: cTime
  character(len=2)  :: cYear, cMonth, cDay, cHour, cMinute, cSecond
  character(len=4)  :: cYearL
  integer :: cL

  logical, save :: IsFirstOutput = .true.

  if (floor((tSimulation - dt)/DtReport) /= &
      floor((tsimulation)/DtReport) .and. iDebugLevel >= 0) then
    if (IsFramework) then
      if (iProc == 0) write(*, "(a,i6,a,3i2.2)") "UA:GITM2 iStep ", iStep, &
        ", Time : ", iTimeArray(4:6)
    else
      RealTime = get_timing("GITM")
      CompletedTime = (EndTime - CurrentTime)/(CurrentTime - RestartTime)
      ProjectedTime = RealTime*CompletedTime
      write(*, "(a,i8,a,i4,2i2.2,a,3i2.2,a,f9.2,a,f9.2)") "iStep: ", iStep, &
        ", Time: ", iTimeArray(1:3), ' ', iTimeArray(4:6), &
        ", WallTime: ", RealTime/60.0, " min, Proj : ", &
        ProjectedTime/60.0
    endif
  endif

  ! Initialize the output registry and I/O backend exactly once.
  if (IsFirstOutput) then
    call init_output_registry()
    call init_usr_output_registry()
    call init_output_backend(requestedBackendName)
    IsFirstOutput = .false.
  endif

  DtPlot = DtPlotSave
  if (CurrentTime >= PlotTimeChangeStart .and. &
      CurrentTime <= PlotTimeChangeEnd) then
    DtPlot = PlotTimeChangeDt
  endif

  ! Construct the filename time string (same logic as in output()).
  call i2s(mod(iTimeArray(1), 100), cYear, 2)
  call i2s(iTimeArray(1), cYearL, 4)
  call i2s(iTimeArray(2), cMonth, 2)
  call i2s(iTimeArray(3), cDay, 2)
  call i2s(iTimeArray(4), cHour, 2)
  call i2s(iTimeArray(5), cMinute, 2)
  call i2s(iTimeArray(6), cSecond, 2)
  if (.not. UseSecondsInFilename) cSecond = '00'
  if (UseCCMCFileName) then
    cTime = "GITM_"//cYearL//"-"//cMonth//"-"//cDay//"T" &
            //cHour//"-"//cMinute//"-"//cSecond
    cL = 24
  else
    cTime = "t"//cYear//cMonth//cDay//"_"//cHour//cMinute//cSecond
    cL = 14
  endif

  do i = 1, nOutputTypes
      ! Compute cType with same Is1D adjustment used inside output().
      cType = OutputType(i)
      if (cType(1:2) == "3D" .and. Is1D) cType(1:2) = "1D"
      ! For collective backends (mpiio, netcdf): collectively open the shared output
      ! file before the block loop.
      ! For legacy: open_file is a no-op; per-block files are opened inside output().
      call ActiveBackend%open_file(outputDir, cType, cTime, cL)
      do iBlock = 1, nBlocks
        call output(outputDir, iBlock, i)
      enddo
      ! For collective backends: collectively close the file after all blocks are done.
      ! For legacy: close_file is a no-op.
      call ActiveBackend%close_file()
    endif
  enddo

  call move_satellites

  if (floor((tSimulation - dt)/DtRestart) /= &
      floor((tsimulation)/DtRestart)) then
    call write_restart(restartOutDir)
  endif

  if (floor((tSimulation - dt)/DtLogfile) /= &
      floor((tsimulation)/DtLogfile)) then
    call logfile(logDir)
  endif

end subroutine write_output
