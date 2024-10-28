! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine finalize_gitm

  use ModInputs
  use ModSphereInterface

  implicit none

  logical :: IsOk
  integer :: iError, iBlock, iOutputType
  integer :: nMLTsTmp, nLatsTmp

  if (.not. Is1D) &
    call UA_calc_electrodynamics(nMLTsTmp, nLatsTmp)

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

  call end_timing("GITM")

  if (iDebugLevel >= 0) call report_timing("all")

  if (.not. Is1D) then
    ! cleanup UAM
     !! get rid of data xfer structure
    call UAM_XFER_destroy(ok=IsOk)
    if (.not. IsOk) then
      call UAM_write_error()
      if (.not. IsFrameWork) call stop_gitm("problem with finalize")
    endif

  endif

  ! cleanup mpi
  if (.not. IsFrameWork) call MPI_FINALIZE(iError)

end subroutine finalize_gitm
