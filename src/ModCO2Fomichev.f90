! Copyright 2026, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

module ModCO2Fomichev

  use ModGITM
  use ModInputs
  use ModMpi

  implicit None

  integer :: nFomichevLines
  character(len=iCharLen_) :: cFomichevText(nInputMaxLines) = ''

  integer, parameter :: nX0s = 43
  integer, parameter :: nColsL = 10
  integer, parameter :: nLowerTables = 8
  real :: x0(nX0s)
  real :: lowerTables(nLowerTables, nColsL, nX0s)
  
contains

  ! ----------------------------------------------------------------
  ! Read in the whole file and store it in a string
  ! ----------------------------------------------------------------
  
  subroutine read_and_store_file(cFile)

    character(len=iCharLen_) :: line  
    character(len=*), intent(in) :: cFile
    integer :: iError
    logical :: IsThere

    if (iProc == 0) then

       call report("Reading Fomichev File", 0)

       nFomichevLines = 1

       inquire(file=cFile, EXIST=IsThere)
       if (.not. IsThere) &
            call stop_gitm(trim(cFile)//" cannot be found by ModCO2Fomichev")

       open(iInputUnit_, file=cFile, status="old")

       iError = 0
       do while (iError == 0)

          read(iInputUnit_, '(a)', iostat=iError) line

          if (nFomichevLines > nInputMaxLines) &
               call stop_gitm("Too many lines of input in read_inputs")

          cFomichevText(nFomichevLines) = line
          nFomichevLines = nFomichevLines + 1

       enddo

       close(iInputUnit_)

       if (nFomichevLines == 0) &
            call stop_gitm("No lines of input read by ModCO2Fomichev")

    endif

    ! Broadcast the number of lines and the text itself to all processors
    call MPI_Bcast(nFomichevLines, 1, MPI_Integer, 0, iCommGITM, ierror)

    if (iError > 0) &
         call stop_gitm("nFomichevLines could not be broadcast ModCO2Fomichev")

    call MPI_Bcast(cFomichevText, len(cFomichevText(1))*nFomichevLines, &
         MPI_Character, 0, iCommGITM, ierror)
    
    if (iError > 0) &
         call stop_gitm("cFomichevText could not be broadcast by read_mpi")

    
  end subroutine read_and_store_file

  ! ----------------------------------------------------------------
  ! Read in a generic table
  ! ----------------------------------------------------------------
  
  subroutine parse_general_table(nRows, nCols, iStartLine, table)

    integer, intent(in) :: nRows, nCols, iStartLine
    real, intent(out) :: table(nRows, nCols)
    real :: oneRow(nCols)
    integer :: iRow, iCol, iLine, iError

    do iRow = 1, nRows
       iLine = iStartLine + iRow - 1
       read(cFomichevText(iLine), *, iostat=iError) (oneRow(iCol), iCol=1,nCols)
       table(iRow, :) = oneRow
    enddo
    
  end subroutine parse_general_table

  ! ----------------------------------------------------------------
  ! find table
  ! ----------------------------------------------------------------

  subroutine find_table(cTableName, iRow)

    character(len=*), intent(in) :: cTableName
    integer, intent(out) :: iRow
    logical :: isDone
    integer :: iLine
    
    iRow = -1
    isDone = .false.
    iLine = 1
    do while (.not. isDone)
       if (index(cFomichevText(iLine), cTableName) > 0) then
          isDone = .true.
          iRow = iLine
       endif
       iLine = iLine + 1
       if (iLine >= nFomichevLines) isDone = .true.
    enddo
    
  end subroutine find_table
  
  ! ----------------------------------------------------------------
  ! Read tables
  ! ----------------------------------------------------------------

  subroutine read_all_tables()

    integer :: iRow
    real :: lowerTable(nColsL, nX0s)

    call find_table('table 2', iRow)
    call parse_general_table(nX0s, nColsL, iRow + 2, lowerTable)

    write(*,*) lowerTable
    
  end subroutine read_all_tables
  
  ! ----------------------------------------------------------------
  ! Initialize
  ! ----------------------------------------------------------------

  subroutine initialize_fomichev_cooling

    call read_and_store_file(cFomichevFile)
    call read_all_tables
    
  end subroutine initialize_fomichev_cooling
  
end module ModCO2Fomichev
  
