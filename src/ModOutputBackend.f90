! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
!
! I/O backend abstraction for GITM output. Uses Fortran 2003 procedure
! pointers for runtime backend selection.
!
! Contains backends for:
!   legacy  -- per-block unformatted sequential files (original format)
!   mpiio   -- single combined .bin file written collectively via MPI-IO, single header

module ModOutputBackend

  use ModOutputRegistry, only: OutputTypeInfo, OutputVar
  use ModOutputContainer, only: OutputContainer
  use ModInputs, only: iCharLen_, iOutputUnit_

  implicit none
  ! Default to legacy so runs without #OUTPUTBACKEND in UAM.in behave identically
  ! to older GITM versions.  set_inputs overrides this when #OUTPUTBACKEND is present.
  character(len=iCharLen_) :: requestedBackendName = 'legacy'

  ! ------------------------------------------------------------------
  ! Abstract interfaces for backend operations
  ! ------------------------------------------------------------------

  abstract interface
    ! Open the output file for one output type at one timestep.
    ! For legacy: no-op (output() opens per-block files itself).
    ! For mpiio:  collective MPI_File_open across all ranks.
    subroutine backend_open_file_iface(dir, cType, cTime, cL)
      character(len=*), intent(in) :: dir, cType, cTime
      integer, intent(in) :: cL
    end subroutine

    ! Close the output file.
    ! For legacy: no-op.  For mpiio: collective MPI_File_close.
    subroutine backend_close_file_iface()
    end subroutine

    ! Write one block's data buffer to the output file.
    ! iUnit -- Fortran unit (legacy only; mpiio ignores it)
    ! iBLK  -- global block index (mpiio uses for byte offset)
    subroutine backend_write_block_iface(iUnit, buffer, nV, nX, nY, nZ, iBLK)
      integer, intent(in) :: iUnit, nV, nX, nY, nZ, iBLK
      real, intent(in) :: buffer(nV, nX, nY, nZ)
    end subroutine

    subroutine backend_write_header_iface(iUnit, info, cType, &
                                          nAlts, nLats, nLons, nGCs, &
                                          nBlocksLat, nBlocksLon, &
                                          iTimeArray, cVersion)
      import :: OutputTypeInfo
      integer, intent(in) :: iUnit
      type(OutputTypeInfo), intent(in) :: info
      character(len=5), intent(in) :: cType
      integer, intent(in) :: nAlts, nLats, nLons, nGCs
      integer, intent(in) :: nBlocksLat, nBlocksLon
      integer, intent(in) :: iTimeArray(7)
      character(len=*), intent(in) :: cVersion
    end subroutine

    ! Write all variables in a container to the output file for one block.
    ! Called per block per timestep on participating ranks.
    subroutine backend_write_container_iface(c, iBlock, iTimeArray)
      import :: OutputContainer
      type(OutputContainer), intent(inout) :: c
      integer, intent(in) :: iBlock
      integer, intent(in) :: iTimeArray(7)
    end subroutine
  end interface

  ! ------------------------------------------------------------------
  ! Backend descriptor with procedure pointers
  ! ------------------------------------------------------------------

  type :: OutputBackend
    character(len=10) :: name = 'legacy'
    ! uses_external_file:
    !   .true.  => output() opens/closes the per-block data file itself (legacy).
    !   .false. => write_output() opens/closes a shared file collectively
    !              before/after the block loop; output() skips open/close.
    logical :: uses_external_file = .true.
    ! writes_header:
    !   .true.  => output() writes an ASCII .header file (legacy, mpiio).
    !   .false. => no .header is written; metadata lives in the data file itself (netcdf).
    logical :: writes_header = .true.
    procedure(backend_open_file_iface), pointer, nopass :: open_file => null()
    procedure(backend_close_file_iface), pointer, nopass :: close_file => null()
    procedure(backend_write_block_iface), pointer, nopass :: write_block => null()
    procedure(backend_write_header_iface), pointer, nopass :: write_header => null()
    procedure(backend_write_container_iface), pointer, nopass :: write_container => null()
  end type OutputBackend

  type(OutputBackend) :: ActiveBackend

contains

  ! ------------------------------------------------------------------
  ! Initialize the active backend. Call once before the first output.
  ! ------------------------------------------------------------------
  subroutine init_output_backend(backend_name)
    use ModOutputMPIIO, only: mpiio_open_file, mpiio_close_file, &
                              mpiio_write_block, mpiio_write_header, &
                              mpiio_write_container
    use ModOutputNetCDF, only: netcdf_open_file, netcdf_close_file, &
                               netcdf_write_block, netcdf_write_container
    character(len=*), intent(in) :: backend_name

    select case (trim(backend_name))
    case ('mpiio')
      ActiveBackend%name = 'mpiio'
      ActiveBackend%uses_external_file = .false.
      ActiveBackend%open_file => mpiio_open_file
      ActiveBackend%close_file => mpiio_close_file
      ActiveBackend%write_block => mpiio_write_block
      ActiveBackend%write_header => mpiio_write_header
      ActiveBackend%write_container => mpiio_write_container
    case ('netcdf')
      ActiveBackend%name = 'netcdf'
      ActiveBackend%uses_external_file = .false.
      ActiveBackend%writes_header = .false.
      ActiveBackend%open_file => netcdf_open_file
      ActiveBackend%close_file => netcdf_close_file
      ActiveBackend%write_block => netcdf_write_block
      ActiveBackend%write_container => netcdf_write_container
      ! write_header left null: output_common.f90 checks writes_header before calling it
    case default
      if (trim(backend_name) /= 'legacy') &
        write(*, *) "WARNING: unknown output backend '", trim(backend_name), &
        "', using legacy"
      ActiveBackend%name = 'legacy'
      ActiveBackend%uses_external_file = .true.
      ActiveBackend%open_file => legacy_open_file
      ActiveBackend%close_file => legacy_close_file
      ActiveBackend%write_block => legacy_write_block
      ActiveBackend%write_header => legacy_write_header
      ActiveBackend%write_container => legacy_write_container
    end select
  end subroutine init_output_backend

  ! ==================================================================
  ! Legacy backend: per-block unformatted sequential files
  ! ==================================================================

  subroutine legacy_open_file(dir, cType, cTime, cL)
    character(len=*), intent(in) :: dir, cType, cTime
    integer, intent(in) :: cL
    ! No-op: legacy opens per-block files inside output()
  end subroutine legacy_open_file

  subroutine legacy_close_file()
    ! No-op: legacy closes files inside output()
  end subroutine legacy_close_file

  ! ------------------------------------------------------------------
  ! Write a data block from buffer to an already-opened unit.
  ! Reproduces the legacy per-grid-point write pattern exactly:
  !   do iz, do iy, do ix: write(unit) buffer(:, ix, iy, iz)
  ! ------------------------------------------------------------------
  subroutine legacy_write_block(iUnit, buffer, nV, nX, nY, nZ, iBLK)
    integer, intent(in) :: iUnit, nV, nX, nY, nZ, iBLK
    real, intent(in) :: buffer(nV, nX, nY, nZ)
    integer :: ix, iy, iz

    do iz = 1, nZ
      do iy = 1, nY
        do ix = 1, nX
          write(iUnit) buffer(1:nV, ix, iy, iz)
        enddo
      enddo
    enddo
  end subroutine legacy_write_block

  ! ------------------------------------------------------------------
  ! Write a complete header file from the registry info.
  ! Reproduces the format of the legacy output_header, write_head_blocks,
  ! write_head_time, and write_head_version subroutines.
  ! ------------------------------------------------------------------
  subroutine legacy_write_header(iUnit, info, cType, &
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
      write(iUnit, "(I7,7A)") nAlts + nGCs*2, " nAltitudes"
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
      elseif (nGCs == 0) then
        write(iUnit, "(I7,A)") nLats, " nLatitude"
        write(iUnit, "(I7,A)") nLons, " nLongitudes"
        write(iUnit, *) " "
        write(iUnit, *) "NO GHOSTCELLS"
      elseif (cType(3:5) == "GEL" .or. cType(3:5) == "TEC" .or. &
              cType(1:5) == "2DANC" .or. cType(3:5) == "HME") then
        write(iUnit, "(I7,A)") nLats, " nLatitude"
        write(iUnit, "(I7,A)") nLons, " nLongitudes"
        write(iUnit, *) " "
        write(iUnit, *) "NO GHOSTCELLS"
      else
        write(iUnit, "(I7,7A)") nLats + nGCs*2, " nLatitudes"
        write(iUnit, "(I7,7A)") nLons + nGCs*2, " nLongitudes"
      endif
    endif
    write(iUnit, *) ""

    ! --- VARIABLE LIST section (auto-generated from registry) ---
    write(iUnit, *) "VARIABLE LIST"
    do i = 1, info%nVars
      write(iUnit, "(I7,A1,a)") i, " ", trim(info%vars(i)%name)
    enddo

    write(iUnit, *) ""
    write(iUnit, *) "END"
    write(iUnit, *) ""

  end subroutine legacy_write_header

  ! ---------------------------------------------------------------------------
  ! legacy_write_container — pack container vars into a single tmpbuf and write
  ! to the legacy unformatted file for the current block.
  ! Replicates the legacy per-grid-point write pattern exactly:
  !   do iz, do iy, do ix: write(unit) buffer(:, ix, iy, iz)
  ! ---------------------------------------------------------------------------
  subroutine legacy_write_container(c, iBlock, iTimeArray)
    use ModInputs, only: iOutputUnit_
    use ModOutputContainer, only: OutputContainer, output_kind
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer, intent(in) :: iTimeArray(7)

    real(output_kind), allocatable :: tmpbuf(:, :, :, :)
    integer :: nV, nX, nY, nZ
    integer :: i, ix, iy, iz, iAxis

    if (.not. c%this_rank_writes) return
    if (c%nVars == 0) return

    ! Grid dims from non-axis vars only (axis vars have degenerate dims).
    nX = 0; nY = 0; nZ = 0
    do i = 1, c%nVars
      if (.not. c%vars(i)%is_axis) then
        nX = max(nX, c%vars(i)%shape3(1))
        nY = max(nY, c%vars(i)%shape3(2))
        nZ = max(nZ, c%vars(i)%shape3(3))
      end if
    end do
    if (nX == 0 .or. nY == 0 .or. nZ == 0) return

    nV = c%nVars
    allocate(tmpbuf(nV, nX, nY, nZ))

    ! Axis vars are identified by position: first is_axis var = lon, second = lat.
    iAxis = 0
    do i = 1, nV
      if (c%vars(i)%is_axis) then
        iAxis = iAxis + 1
        if (iAxis == 1) then
          ! Longitude axis: replicate 1D data along lat and alt dimensions.
          do iz = 1, nZ
            do iy = 1, nY
              tmpbuf(i, :, iy, iz) = c%vars(i)%data(:, 1, 1)
            end do
          end do
        else
          ! Latitude axis: replicate 1D data along lon and alt dimensions.
          do iz = 1, nZ
            do ix = 1, nX
              tmpbuf(i, ix, :, iz) = c%vars(i)%data(:, 1, 1)
            end do
          end do
        end if
      else
        tmpbuf(i, :, :, :) = c%vars(i)%data
      end if
    end do

    ! Write the packed buffer to the legacy unformatted unit
    do iz = 1, nZ
      do iy = 1, nY
        do ix = 1, nX
          write(iOutputUnit_) tmpbuf(1:nV, ix, iy, iz)
        enddo
      enddo
    enddo

    deallocate(tmpbuf)
  end subroutine legacy_write_container

end module ModOutputBackend
