subroutine lower_case(String)

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
    character(len=*), intent(in) :: cError
    nErrors = nErrors + 1
    cErrorCodes(nErrors) = cError
    isOk = .false.
end subroutine set_error

subroutine report_errors()
    use ModErrors
    integer :: iError
    if (nErrors == 0) write(*, *) "No errors to report!"
    do iError = 1, nErrors
        write(*, *) "--> Error : ", trim(cErrorCodes(iError))
    enddo
end subroutine report_errors

! -- This is for things that should not stop GITM, but notify user now & later. -- !
subroutine raise_warning(cWarning)
    use ModErrors
    character(len=*), intent(in) :: cWarning
    nWarnings = nWarnings + 1
    cWarningCodes(nWarnings) = cWarning
    write(*, *) " -> Warning: ", trim(cWarningCodes(nWarnings))
end subroutine raise_warning

subroutine report_warnings()
    use ModErrors
    integer :: iWarning
    if (nWarnings == 0) write(*, *) "No errors to report!"
    do iWarning = 1, nWarnings
        write(*, *) "--> Error : ", trim(cWarningCodes(iWarning))
    enddo
end subroutine report_warnings