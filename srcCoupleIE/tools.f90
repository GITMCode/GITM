! Error/warning bookkeeping for the SWMF-coupled build.
!
! Standalone GITM gets nErrors/nWarnings from the Electrodynamics library

module ModCoupleErrorCount
  implicit none
  integer :: nErrors = 0
  integer :: nWarnings = 0
end module ModCoupleErrorCount

subroutine lower_case(String)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  character(len=*), intent(inout) :: String

  !DESCRIPTION:
  ! Change characters to lower case in String
  !EOP

  integer, parameter :: iA = ichar('A'), iZ = ichar('Z'), Di = ichar('a') - iA
  integer :: i, iC
  !--------------------------------------------------------------------------
  do i = 1, len_trim(String)
    iC = ichar(String(i:i))
    if (iC >= iA .and. iC <= iZ) String(i:i) = char(iC + Di)
  enddo

end subroutine lower_case

subroutine set_error(cError)
  use ModErrors
  use ModCoupleErrorCount
  implicit none
  character(len=*), intent(in) :: cError
  if (nErrors < nErrorsMax) then
    nErrors = nErrors + 1
    cErrorCodes(nErrors) = cError
  endif
  isOk = .false.
end subroutine set_error

subroutine report_errors()
  use ModErrors
  use ModCoupleErrorCount
  implicit none
  integer :: iError
  if (nErrors == 0) write(*, *) "No errors to report!"
  do iError = 1, min(nErrors, nErrorsMax)
    write(*, *) "--> Error : ", trim(cErrorCodes(iError))
  enddo
end subroutine report_errors

! -- This is for things that should not stop GITM, but notify user now & later. -- !
subroutine raise_warning(cWarning)
  use ModErrors
  use ModCoupleErrorCount
  implicit none
  character(len=*), intent(in) :: cWarning
  if (nWarnings < nWarningsMax) then
    nWarnings = nWarnings + 1
    cWarningCodes(nWarnings) = cWarning
  endif
  write(*, *) " -> Warning: ", trim(cWarning)
end subroutine raise_warning

subroutine report_warnings()
  use ModErrors
  use ModCoupleErrorCount
  implicit none
  integer :: iWarning
  if (nWarnings == 0) write(*, *) "No warnings to report!"
  do iWarning = 1, min(nWarnings, nWarningsMax)
    write(*, *) "--> Warning : ", trim(cWarningCodes(iWarning))
  enddo
end subroutine report_warnings
