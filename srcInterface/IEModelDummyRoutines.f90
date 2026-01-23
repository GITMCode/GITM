

!! Set model name and interpret. It needs to be swmf though...
  subroutine set_efield_model(this, efield_model)
    class(ieModel) :: this
    character(len=*), intent(in) :: efield_model
    if (this%iDebugLevel > 0) &
      write(*, *) "=> Setting efield model to : ", trim(efield_model)
    this%iEfield_ = efield_interpret_name(efield_model)
    if (this%iDebugLevel > 0) &
      write(*, *) "=> That is model : ", this%iEfield_
  end subroutine set_efield_model

    subroutine set_aurora_model(this, aurora_model)
    class(ieModel) :: this
    character(len=*), intent(in) :: aurora_model
    if (this%iDebugLevel > 0) &
      write(*, *) "=> Setting aurora model to : ", trim(aurora_model)
    this%iEfield_ = aurora_interpret_name(aurora_model)
    if (this%iDebugLevel > 0) &
      write(*, *) "=> That is model : ", this%iAurora_
  end subroutine set_aurora_model


  integer function efield_interpret_name(efieldString)
    use ModErrors
    use ModIE
    implicit none
    character(len=*), intent(in) :: efieldString
    character(len=len(efieldString)) :: efieldLower
    
    efield_interpret_name = -1
    efieldLower = efieldString

    call lower_case(efieldLower)
    if (trim(efieldLower) == "swmf") &
      efield_interpret_name = iSwmfPot_

    if (efield_interpret_name == -1) then
      call set_error("efield model MUST be set to swmf! Received:")
      call set_error(efieldLower)
    endif
    
    return
  end function efield_interpret_name


  integer function aurora_interpret_name(auroraString)

    use ModErrors
    use ModIE
    implicit none

    character(len=*), intent(in) :: auroraString
    character(len=len(efieldString)) :: auroraLower

    aurora_interpret_name = -1

    auroraLower = auroraString
    call lower_case(auroraLower)

    if (trim(auroraLower) == "swmf") &
      aurora_interpret_name = iSwmfAur_

    if (aurora_interpret_name == -1) then
      call set_error("Aurora model MUST be set to swmf! Received:")
      call set_error(auroraLower)
    endif

    return

end function aurora_interpret_name

  ! ------------------------------------------------------------
  ! filename for north
  subroutine set_filename_north(this, filename)
    class(ieModel) :: this
    character(len=*), intent(in) :: filename
  end subroutine set_filename_north
! South
  subroutine set_filename_south(this, filename)
    class(ieModel) :: this
    character(len=*), intent(in) :: filename
  end subroutine set_filename_south

! set nfiles for north
  subroutine set_nfiles_north(this, nFiles)
    class(ieModel) :: this
    integer, intent(in) :: nFiles
  end subroutine set_nfiles_north

  ! ------------------------------------------------------------
  ! filename for north
  subroutine set_filename_list_north(this, iFile, filename)
    class(ieModel) :: this
    integer, intent(in) :: iFile
    character(len=*), intent(in) :: filename
  end subroutine set_filename_list_north

  ! ------------------------------------------------------------
  ! set nfiles for south
  subroutine set_nfiles_south(this, nFiles)
    class(ieModel) :: this
    integer, intent(in) :: nFiles
  end subroutine set_nfiles_south

  ! ------------------------------------------------------------
  ! filename for south
  subroutine set_filename_list_south(this, iFile, filename)
    class(ieModel) :: this
    integer, intent(in) :: iFile
    character(len=*), intent(in) :: filename
  end subroutine set_filename_list_south

    subroutine set_model_dir(this, dir)
    class(ieModel) :: this
    character(len=*), intent(in) :: dir
    
  end subroutine set_model_dir
