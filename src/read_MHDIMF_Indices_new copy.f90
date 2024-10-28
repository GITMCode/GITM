! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!==============================================================================

subroutine read_MHDIMF_Indices_new(iOutputError, StartTime, EndTime)

  use ModKind
  use ModIndices
  use ModTime, only: CurrentTime
  implicit none

  integer, intent(out)     :: iOutputError
  real(Real8_), intent(in) :: EndTime, StartTime
  real(Real8_) :: EndTimeFake

  integer :: ierror, iIMF, iSW, j, npts
  logical :: done, done_inner, IsFirstLine = .true.

  ! One line of input
  character(len=iCharLenIndices_) :: line

  real(Real8_) :: TimeDelay, BufferTime = 1800.0, FirstTime, DeltaT = -1.0e32

  integer, dimension(7) :: itime

  logical :: isIgnoreEndTime = .false.
  !------------------------------------------------------------------------
  iOutputError = 0

  done = .false.

  npts = 0
  TimeDelay = 0.0

!   print*,"START<END<DETLA", StartTime, EndTime, EndTime-StartTime
  if (EndTime < StartTime) then
    isIgnoreEndTime = .true.
    EndTimeFake = 1e32
  else
    EndTimeFake = EndTime
  endif

  ! Problem - we have to figure out whether we are being called with no IMF
  ! file at all or whether we actually have an index file to read (or reread)
  ! Can't change comparability, so we can only use certain information.

  ! We know that nIndices_V(imf_by_) = 1 if there is a single index.
  ! We know that NameOfIMFFile=none if (a) we have never read the file or (b)
  !    we have no index file.

  if (NameOfIMFFile == "none" .and. nIndices_V(imf_by_) == 1) return

  call init_mod_indices

  ! If we have been here before and we read the entire file, leave
  if (.not. ReReadIMFFile .and. nIndices_V(imf_by_) > 0) return

  ! This if statement makes sure that we have been here before and
  ! want to be here again
  if (ReReadIMFFile) then
    ! If we still have a lot of data in memory, then don't bother
    ! reading more.
    if (StartTime + BufferTime < IndexTimes_TV(nIndices_V(imf_bx_), imf_by_)) return
  endif

  ! Assume that we can read the entire file
  ReReadIMFFile = .false.

  if (nIndices_V(imf_bx_) == 0) &
    NameOfIMFFile = NameOfIndexFile(:index(NameOfIndexFile, ' ') - 1)

  open(LunIndices_, file=trim(NameOfIMFFile), status="old", iostat=ierror)

  if (ierror .ne. 0) then
    iOutputError = 1
    return
  endif

  do while (.not. done)

    read(LunIndices_, '(a)', iostat=ierror) line
    !  print*, line, "READ WITH ERROR _____ ", ierror
    if (ierror /= 0) done = .true.

    if (index(line, '#DELAY') > 0) then
      read(LunIndices_, *, iostat=iError) TimeDelay
      if (iError /= 0) done = .true.
    endif

    if (index(line, '#START') > 0) then

      done_inner = .false.

      iIMF = 1
      iSW = 1

      do while (.not. done_inner)

        read(LunIndices_, *, iostat=iError) &
          (iTime(j), j=1, 7), &
          Indices_TV(iIMF, imf_bx_), &
          Indices_TV(iIMF, imf_by_), &
          Indices_TV(iIMF, imf_bz_), &
          Indices_TV(iSW, sw_vx_), &
          Indices_TV(iSW, sw_vy_), &
          Indices_TV(iSW, sw_vz_), &
          Indices_TV(iSW, sw_n_), &
          Indices_TV(iSW, sw_t_)

        if (ierror /= 0) then
          done_inner = .true.

          ! This means that the GITM time is all AFTER the first
          ! line in the file!
          if (StartTime > IndexTimes_TV(iIMF, imf_bx_)) then
            iIMF = iIMF + 1
            iSW = iSW + 1
          endif

        else

          call time_int_to_real(iTime, IndexTimes_TV(iIMF, imf_bx_))

          if (IsFirstLine) then
            FirstTime = IndexTimes_TV(iIMF, imf_bx_)
            IsFirstLine = .false.
          else
            if (DeltaT == -1.0e32) then
              DeltaT = IndexTimes_TV(iIMF, imf_bx_) - FirstTime
              if (DeltaT > BufferTime) BufferTime = 10.0*DeltaT
            endif
          endif

          IndexTimes_TV(iIMF, imf_bx_) = IndexTimes_TV(iIMF, imf_bx_) &
                                         + TimeDelay
          IndexTimes_TV(iIMF, imf_by_) = IndexTimes_TV(iIMF, imf_bx_) &
                                         + TimeDelay
          IndexTimes_TV(iIMF, imf_bz_) = IndexTimes_TV(iIMF, imf_bx_) &
                                         + TimeDelay
          IndexTimes_TV(iSW, sw_vx_) = IndexTimes_TV(iSW, imf_bx_) &
                                       + TimeDelay
          IndexTimes_TV(iSW, sw_vy_) = IndexTimes_TV(iSW, imf_bx_) &
                                       + TimeDelay
          IndexTimes_TV(iSW, sw_vz_) = IndexTimes_TV(iSW, imf_bx_) &
                                       + TimeDelay
          IndexTimes_TV(iSW, sw_n_) = IndexTimes_TV(iSW, imf_bx_) &
                                      + TimeDelay
          IndexTimes_TV(iSW, sw_t_) = IndexTimes_TV(iSW, imf_bx_) &
                                      + TimeDelay

          Indices_TV(iSW, sw_v_) = sqrt( &
                                   Indices_TV(iSW, sw_vx_)**2 + &
                                   Indices_TV(iSW, sw_vy_)**2 + &
                                   Indices_TV(iSW, sw_vz_)**2)
          IndexTimes_TV(iSW, sw_v_) = IndexTimes_TV(iSW, imf_bx_)

          ! This makes sure that we only store the values that we
          ! are really interested in

          if (IndexTimes_TV(iIMF, imf_bx_) >= StartTime - BufferTime .and. &
              IndexTimes_TV(iIMF, imf_bx_) <= EndTimeFake + BufferTime .and. &
              iIMF < MaxIndicesEntries) then

            if (abs(Indices_TV(iSW, imf_bz_)) < 200.0) iIMF = iIMF + 1
            if ((abs(Indices_TV(iSW, sw_n_)) < 900.0) .and. &
                (abs(Indices_TV(iSW, sw_v_)) < 3000.0)) iSW = iSW + 1

          else

            ! This means that the GITM time is all BEFORE the first
            ! line in the file!
            if (EndTimeFake < IndexTimes_TV(iIMF, imf_bx_) .and. iIMF == 1) then
              iIMF = iIMF + 1
              iSW = iSW + 1
              print *, "COUNTERS ADDED TO"
            endif
          endif

        endif

      enddo

      done = done_inner

    endif

  enddo

  close(LunIndices_)

  if (iIMF >= MaxIndicesEntries) ReReadIMFFile = .true.

  nIndices_V(imf_bx_) = iIMF - 2
  nIndices_V(imf_by_) = iIMF - 2
  nIndices_V(imf_bz_) = iIMF - 2
  nIndices_V(sw_vx_) = iSW - 2
  nIndices_V(sw_vy_) = iSW - 2
  nIndices_V(sw_vz_) = iSW - 2
  nIndices_V(sw_v_) = iSW - 2
  nIndices_V(sw_n_) = iSW - 2
  nIndices_V(sw_t_) = iSW - 2

  if (isIgnoreEndTime) then
    EndTimeFake = IndexTimes_TV(iIMF - 2, imf_bx_)
  endif
  ! If we have gotten to this point and we have no data,
  ! there is something wrong!
  if (nIndices_V(imf_bz_) < 2) iOutputError = 1
  if (nIndices_V(sw_v_) < 2) iOutputError = 1
  print *, CurrentTime, EndTimeFake, EndTime, iError

  call check_all_indices(CurrentTime, iError)
  print *, CurrentTime, EndTimeFake, EndTime, iError
  if (iError == 0) then
    call check_all_indices(EndTimeFake, iError)
    print *, CurrentTime, EndTimeFake, EndTime, iError
  endif
  if (iError /= 0) then
    print *, CurrentTime, EndTimeFake, EndTime, iError
    call stop_gitm("Issue with Indices! Check the file(s) times!")
  endif

end subroutine read_MHDIMF_Indices_new

