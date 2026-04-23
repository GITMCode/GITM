! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputMPIIO.f90
!
! MPI-IO backend for GITM output.
!
! All MPI ranks collectively open a single .bin file per output type per
! timestep. Each rank writes its block(s) at a byte offset computed from
! the global block index iBLK:
!
!   offset = (iBLK - 1) * nV * nX * nY * nZ * bytes_per_element
!
! Data layout: raw floats (real precision matches compilation), no Fortran
! record markers. Buffer order: buffer(iV, ix, iy, iz), iV fastest (Fortran
! column-major).
!
! A companion ASCII .header file is written by rank 0. Its content mirrors
! the legacy header format, with a FORMAT/mpiio section prepended so readers
! know this is a combined raw-binary file rather than per-block unformatted.
!
! Limitations in Phase 2:
!   - DoAppendFiles is not supported; each call creates/overwrites the file.
!   - HME output types fall back to legacy (variable block count not suitable
!     for fixed-offset addressing).

module ModOutputMPIIO

  use ModOutputRegistry, only: OutputTypeInfo
  use ModMpi

  implicit none

  integer :: mpiio_fh = -1  ! MPI file handle; -1 means no file is open

contains

  ! ------------------------------------------------------------------
  ! Open the shared output file collectively (all ranks must call).
  ! Filename: dir/cType_cTime(1:cL).bin
  ! ------------------------------------------------------------------
  subroutine mpiio_open_file(dir, cType, cTime, cL)
    use ModGITM, only: iCommGITM
    character(len=*), intent(in) :: dir, cType, cTime
    integer, intent(in) :: cL
    character(len=300) :: filename
    integer :: ierr

    write(filename, '(a,"/",a,"_",a,".bin")') &
      trim(dir), trim(cType), cTime(1:cL)

    call MPI_File_open(iCommGITM, trim(filename), &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, mpiio_fh, ierr)
    if (ierr /= MPI_SUCCESS) &
      write(*, *) 'ModOutputMPIIO: MPI_File_open failed for ', trim(filename)
  end subroutine mpiio_open_file

  ! ------------------------------------------------------------------
  ! Close the shared output file collectively.
  ! ------------------------------------------------------------------
  subroutine mpiio_close_file()
    integer :: ierr
    if (mpiio_fh /= -1) then
      call MPI_File_close(mpiio_fh, ierr)
      mpiio_fh = -1
    endif
  end subroutine mpiio_close_file

  ! ------------------------------------------------------------------
  ! Write block data at the byte offset determined by iBLK.
  !
  ! offset = (iBLK - 1) * nV * nX * nY * nZ * bytes_per_element
  !
  ! iUnit is ignored (the module holds the open file handle).
  ! ------------------------------------------------------------------
  subroutine mpiio_write_block(iUnit, buffer, nV, nX, nY, nZ, iBLK)
    integer, intent(in) :: iUnit, nV, nX, nY, nZ, iBLK
    real, intent(in) :: buffer(nV, nX, nY, nZ)
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer :: nWords, bytes_per_element, ierr
    integer :: status(MPI_STATUS_SIZE)

    call MPI_Type_size(MPI_REAL, bytes_per_element, ierr)

    nWords = nV*nX*nY*nZ
    offset = int(iBLK - 1, MPI_OFFSET_KIND)* &
             int(nWords, MPI_OFFSET_KIND)* &
             int(bytes_per_element, MPI_OFFSET_KIND)

    call MPI_File_write_at(mpiio_fh, offset, buffer, nWords, MPI_REAL, &
                           status, ierr)
    if (ierr /= MPI_SUCCESS) &
      write(*, *) 'ModOutputMPIIO: MPI_File_write_at failed, iBLK=', iBLK
  end subroutine mpiio_write_block

  ! ------------------------------------------------------------------
  ! Write the ASCII header file to an already-opened unit (iUnit).
  ! Called from output() on rank 0 after all blocks are written.
  ! Content mirrors legacy_write_header with a FORMAT/mpiio section
  ! prepended so post-processors can identify the file type.
  ! ------------------------------------------------------------------
  subroutine mpiio_write_header(iUnit, info, cType, &
                                nAlts, nLats, nLons, nGCs, &
                                nBlocksLat, nBlocksLon, &
                                iTimeArray, cVersion)
    use ModElectrodynamics, only: nMagLats, nMagLons
    integer, intent(in) :: iUnit
    type(OutputTypeInfo), intent(in) :: info
    character(len=5), intent(in) :: cType
    integer, intent(in) :: nAlts, nLats, nLons, nGCs
    integer, intent(in) :: nBlocksLat, nBlocksLon
    integer, intent(in) :: iTimeArray(7)
    character(len=*), intent(in) :: cVersion
    integer :: i

    ! FORMAT marker: post-processors use this to detect the raw-binary layout
    write(iUnit, *) "FORMAT"
    write(iUnit, *) "mpiio"
    write(iUnit, *) ""

    ! PRECISION: bytes per real element (4 = single, 8 = double)
    ! Needed by post-processor to read the raw binary correctly.
    write(iUnit, *) "PRECISION"
    write(iUnit, "(I7,A)") kind(1.0), " bytes_per_real"
    write(iUnit, *) ""

    ! --- BLOCKS section ---
    if (cType(1:2) /= "1D" .and. cType(1:2) /= "0D") then
      write(iUnit, *) "BLOCKS"
      write(iUnit, "(I7,A)") 1, " nBlocksAlt"
      if (cType /= "2DMEL") then
        write(iUnit, "(I7,A)") nBlocksLat, " nBlocksLat"
        write(iUnit, "(I7,A)") nBlocksLon, " nBlocksLon"
      else
        write(iUnit, "(I7,A)") 1, " nBlocksLat"
        write(iUnit, "(I7,A)") 1, " nBlocksLon"
      endif
      write(iUnit, *) ""
    endif

    ! --- TIME section ---
    write(iUnit, *) "TIME"
    write(iUnit, "(I7,A)") iTimeArray(1), " Year"
    write(iUnit, "(I7,A)") iTimeArray(2), " Month"
    write(iUnit, "(I7,A)") iTimeArray(3), " Day"
    write(iUnit, "(I7,A)") iTimeArray(4), " Hour"
    write(iUnit, "(I7,A)") iTimeArray(5), " Minute"
    write(iUnit, "(I7,A)") iTimeArray(6), " Second"
    write(iUnit, "(I7,A)") iTimeArray(7), " Millisecond"
    write(iUnit, *) ""

    ! --- VERSION section ---
    write(iUnit, *) "VERSION"
    write(iUnit, *) trim(cVersion)
    write(iUnit, *) ""

    ! --- NUMERICAL VALUES section ---
    write(iUnit, *) "NUMERICAL VALUES"
    write(iUnit, "(I7,6A)") info%nVars, " nvars"

    if (cType == '1DNEW') then
      write(iUnit, "(I7,7A)") nAlts, " nAltitudes"
    elseif (cType(1:2) /= "2D" .and. cType(1:2) /= "0D") then
      write(iUnit, "(I7,7A)") nAlts + 4, " nAltitudes"
    else
      write(iUnit, "(I7,7A)") 1, " nAltitudes"
    endif

    if (cType(1:2) == "1D" .or. cType(1:2) == "0D") then
      write(iUnit, "(I7,7A)") 1, " nLatitudes"
      write(iUnit, "(I7,7A)") 1, " nLongitudes"
    else
      if (info%usesMagGrid) then
        write(iUnit, "(I7,A)") nMagLats, " nLatitude"
        write(iUnit, "(I7,A)") nMagLons + 1, " nLongitudes"
        write(iUnit, *) " "
        write(iUnit, *) "NO GHOSTCELLS"
      elseif (info%nGhostCells == 0) then
        write(iUnit, "(I7,A)") nLats, " nLatitude"
        write(iUnit, "(I7,A)") nLons, " nLongitudes"
        write(iUnit, *) " "
        write(iUnit, *) "NO GHOSTCELLS"
      else
        write(iUnit, "(I7,7A)") nLats + info%nGhostCells*2, " nLatitudes"
        write(iUnit, "(I7,7A)") nLons + info%nGhostCells*2, " nLongitudes"
      endif
    endif
    write(iUnit, *) ""

    ! --- VARIABLE LIST section ---
    write(iUnit, *) "VARIABLE LIST"
    do i = 1, info%nVars
      write(iUnit, "(I7,A1,a)") i, " ", trim(info%vars(i)%name)
    enddo

    write(iUnit, *) ""
    write(iUnit, *) "END"
    write(iUnit, *) ""

  end subroutine mpiio_write_header

end module ModOutputMPIIO
