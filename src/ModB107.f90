! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!------------------------------------------------------------------------------
! B10.7: an F10.7 -> B10.7 lookup table (#B107). When it is on, the empirical
! EUV models, MSIS and the NO photodissociation rate see B10.7 instead of
! F10.7. The conversion is always made into local variables; the f107/f107a
! globals (ModInputs) stay raw F10.7, so IRI, RCMR, the log and the restart
! files are unchanged.
!------------------------------------------------------------------------------

module ModB107

  use ModInputs, only: UseB107

  implicit none

  private
  public :: read_b107, b107_convert

  integer, parameter :: nB107Max = 1000
  integer :: nB107 = 0
  real :: F107Table(nB107Max), B107Table(nB107Max)

contains

  !----------------------------------------------------------------------------
  ! Read the table on iProc 0 and broadcast it. The file has '#' comment lines,
  ! then #START, before the rows of "F107 B107"
  !----------------------------------------------------------------------------

  subroutine read_b107(cFile)

    use ModGITM, only: iProc, iCommGITM
    use ModInputs, only: iInputUnit_
    use ModMpi

    character(len=*), intent(in) :: cFile

    character(len=200) :: cLine
    integer :: iError, i
    logical :: IsStarted
    real :: f, b

    if (iProc == 0) then

      open(iInputUnit_, file=cFile, status="old", iostat=iError)
      if (iError /= 0) call stop_gitm("Cannot open the #B107 file "//trim(cFile))

      IsStarted = .false.
      nB107 = 0
      do
        read(iInputUnit_, '(a)', iostat=iError) cLine
        if (iError /= 0) exit
        if (.not. IsStarted) then
          if (index(cLine, '#START') == 1) IsStarted = .true.
          cycle
        endif
        if (len_trim(cLine) == 0 .or. cLine(1:1) == '#') cycle
        read(cLine, *, iostat=iError) f, b
        if (iError /= 0) call stop_gitm("Bad line in the #B107 file: "//trim(cLine))
        if (nB107 == nB107Max) call stop_gitm("Too many rows in the #B107 file")
        nB107 = nB107 + 1
        F107Table(nB107) = f
        B107Table(nB107) = b
      enddo
      close(iInputUnit_)

      if (.not. IsStarted) call stop_gitm("No #START line in the #B107 file")
      if (nB107 < 2) call stop_gitm("The #B107 file needs at least 2 rows")
      do i = 2, nB107
        if (F107Table(i) <= F107Table(i - 1)) &
          call stop_gitm("F10.7 must increase down the #B107 file")
      enddo

      write(*, *) '=> B10.7 table read from ', trim(cFile), ': ', nB107, &
        ' rows, F10.7 ', F107Table(1), ' - ', F107Table(nB107)

    endif

    call MPI_Bcast(nB107, 1, MPI_INTEGER, 0, iCommGITM, iError)
    call MPI_Bcast(F107Table, nB107Max, MPI_REAL, 0, iCommGITM, iError)
    call MPI_Bcast(B107Table, nB107Max, MPI_REAL, 0, iCommGITM, iError)

    ! Show a few conversions, including across the gap at F10.7 = 130
    if (iProc == 0) write(*, '(a,7f8.2)') &
      ' => B10.7 at F10.7 = 100,129,130,131,150,250,400 : ', &
      b107_value(100.0), b107_value(129.0), b107_value(130.0), &
      b107_value(131.0), b107_value(150.0), b107_value(250.0), &
      b107_value(400.0)

  end subroutine read_b107

  !----------------------------------------------------------------------------
  ! B10.7 for one F10.7 value: linear interpolation in the table, B = F - 1
  ! below the table (continuous with its first row, 103 -> 102), and held at
  ! the last value above it (B10.7 saturates).
  !----------------------------------------------------------------------------

  real function b107_value(f)

    real, intent(in) :: f
    integer :: i

    if (f < F107Table(1)) then
      b107_value = f - 1.0
    else if (f >= F107Table(nB107)) then
      b107_value = B107Table(nB107)
    else
      i = 1
      do while (F107Table(i + 1) <= f)
        i = i + 1
      enddo
      b107_value = B107Table(i) + (B107Table(i + 1) - B107Table(i))* &
                   (f - F107Table(i))/(F107Table(i + 1) - F107Table(i))
    endif

  end function b107_value

  !----------------------------------------------------------------------------
  ! The solar index the EUV / MSIS / NO-photolysis code should use. With
  ! #B107 off this returns the inputs unchanged (bit-identical).
  ! First cut: the 81-day mean is converted directly, B(mean(F107)), not
  ! mean(B(F107))
  !----------------------------------------------------------------------------

  subroutine b107_convert(f107In, f107aIn, f107Out, f107aOut)

    real, intent(in) :: f107In, f107aIn
    real, intent(out) :: f107Out, f107aOut

    if (UseB107) then
      f107Out = b107_value(f107In)
      f107aOut = b107_value(f107aIn)
    else
      f107Out = f107In
      f107aOut = f107aIn
    endif

  end subroutine b107_convert

end module ModB107
