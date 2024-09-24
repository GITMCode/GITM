!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine indices_set_inputs(StringInputLines)

  use ModIndices

  implicit none

  character(len=*), dimension(*), intent(in) :: StringInputLines
  character(len=iCharLenIndices_) :: StringLine
  logical :: IsDone
  integer :: iLine

  iLine = 1

  IsDone = .false.

  do while (.not. IsDone)

    StringLine = StringInputLines(iLine)

    if (StringLine(1:1) == "#") then

      if (index(StringLine, "#NGDC_INDICES") > 0) then
        call read_in_string(NameOfIndexFile)
      end if

      if (index(StringLine, "#NOAAHPI_INDICES") > 0) then
        call read_in_string(NameOfIndexFile)
      end if

      if (index(StringLine, "#MHD_INDICES") > 0) then
        call read_in_string(NameOfIndexFile)
      end if

      if (index(StringLine, "#END") > 0) then
        IsDone = .true.
      end if

      if (iLine >= MaxInputLines) then
        IsDone = .true.
      end if

    else

      iLine = iLine + 1

    end if

  end do

contains

  subroutine read_in_int(variable)
    integer, intent(out) :: variable
    iline = iline + 1
    read(StringInputLines(iline), *) variable
  end subroutine read_in_int

  subroutine read_in_logical(variable)
    logical, intent(out) :: variable
    iline = iline + 1
    read(StringInputLines(iline), *) variable
  end subroutine read_in_logical

  subroutine read_in_string(variable)
    character(len=iCharLenIndices_), intent(out) :: variable
    iline = iline + 1
    variable = StringInputLines(iline)
  end subroutine read_in_string

  subroutine read_in_real(variable)
    real :: variable
    iline = iline + 1
    read(StringInputLines(iline), *) variable
  end subroutine read_in_real

  subroutine read_in_time(variable)
    real*8 :: variable
    integer, dimension(7) :: itime
    iline = iline + 1
    read(StringInputLines(iline), *) itime(1:6)
    itime(7) = 0
!    call time_int_to_real(itime, variable)
  end subroutine read_in_time

end subroutine Indices_set_inputs
