! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine read_inputs(cFile)

  use ModGITM
  use ModInputs
  use ModMpi
  implicit none

  character(len=*), intent(in) :: cFile

  character(len=iCharLen_) :: line
  logical :: IsThere
  integer :: iError, i

  if (iProc == 0) then

    call report("Reading Inputs", 0)

    nInputLines = 1

    inquire(file=cFile, EXIST=IsThere)
    if (.not. IsThere) &
      call stop_gitm(trim(cFile)//" cannot be found by read_inputs")

    open(iInputUnit_, file=cFile, status="old")

    iError = 0
    do while (iError == 0)

      read(iInputUnit_, '(a)', iostat=iError) line

      ! A failed (EOF) read leaves line holding stale text; blank it so the
      ! slot this iteration still writes cannot be mistaken for a real input
      ! line (e.g. by parsers peeking ahead for optional trailing flags).
      if (iError /= 0) line = ""

      if (nInputLines > nInputMaxLines) &
        call stop_gitm("Too many lines of input in read_inputs")

      cInputText(nInputLines) = line
      nInputLines = nInputLines + 1

    enddo

    close(iInputUnit_)

    if (nInputLines == 0) &
      call stop_gitm("No lines of input read by read_inputs")

  endif

  ! Broadcast the number of lines and the text itself to all processors
  call MPI_Bcast(nInputLines, 1, MPI_Integer, 0, iCommGITM, ierror)

  if (iError > 0) &
    call stop_gitm("nInputLines could not be broadcast by read_inputs")

  call MPI_Bcast(cInputText, len(cInputText(1))*nInputLines, &
                 MPI_Character, 0, iCommGITM, ierror)

  if (iError > 0) &
    call stop_gitm("cInputText could not be broadcast by read_mpi")

end subroutine read_inputs
