! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine finalize_gitm

  use ModInputs
  use ModSphereInterface
  use ModErrors

  implicit none

  logical :: didOk
  integer :: iError, iBlock, iOutputType
  integer :: nMLTsTmp, nLatsTmp

  call start_timing("Finalize")

  ! if (.not. Is1D) &
  !   call UA_calc_electrodynamics(nMLTsTmp, nLatsTmp)

  do iOutputType = 1, nOutputTypes
    do iBlock = 1, nBlocks
      call output("UA/data/", iBlock, iOutputType)
    enddo
  enddo

  if (IsOpenLogFile) close(iLogFileUnit_)

  if (.not. IsFrameWork) call write_restart("UA/restartOUT/")

  if (iProc == 0) then
    open(unit=iOutputUnit_, file="GITM.DONE", status="unknown")
    close(iOutputUnit_)
  endif

  if (iProc == 0) then
    call report_errors
    call report_warnings
  endif

  call end_timing("Finalize")
  call end_timing("GITM")

  if (iDebugLevel >= 0) call report_timing("all")

  if (.not. Is1D) then
    ! cleanup UAM
     !! get rid of data xfer structure
    call UAM_XFER_destroy(ok=didOk)
    if (.not. didOk) then
      call UAM_write_error()
      if (.not. IsFrameWork) call stop_gitm("problem with finalize")
    endif

  endif

  ! cleanup mpi
  if (.not. IsFrameWork) call MPI_FINALIZE(iError)

end subroutine finalize_gitm
