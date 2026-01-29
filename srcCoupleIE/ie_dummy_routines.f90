! ==========================================================================
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
! ------------------------------------------------------------
subroutine set_model_dir(this, dir)
  class(ieModel) :: this
  character(len=*), intent(in) :: dir
end subroutine set_model_dir
! ------------------------------------------------------------

! ------------------------------------------------------------
! These functions are for getting information that was
! derived when the auroral model was run
! ------------------------------------------------------------

! Supply the polar cap to the user:
! polarcap = 1 inside the polar cap, and 0 elsewhere.
subroutine get_polarcap_results(ie, polarcap)
  class(ieModel) :: ie
  real, dimension(ie%neednMlts, &
          ie%neednLats), intent(out) :: polarcap
  return
end subroutine get_polarcap_results
! ------------------------------------------------------------
subroutine set_bz(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_bz
! ------------------------------------------------------------
subroutine set_by(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_by
! ------------------------------------------------------------
subroutine set_swv(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_swv
! ------------------------------------------------------------
subroutine set_swn(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_swn
! ------------------------------------------------------------
subroutine set_hp(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_hp
! ------------------------------------------------------------
subroutine set_hpn(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_hpn
! ------------------------------------------------------------
subroutine set_hps(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_hps
! ------------------------------------------------------------
subroutine set_hp_from_ae(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_hp_from_ae
! ------------------------------------------------------------
subroutine set_au(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_au
! ------------------------------------------------------------
subroutine set_al(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_al
! ------------------------------------------------------------
subroutine set_ae(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_ae
! ------------------------------------------------------------
subroutine set_useAeForHp(this)
  class(ieModel), intent(inout) :: this
end subroutine set_useAeForHp
! ------------------------------------------------------------
subroutine unset_useAeForHp(this)
  class(ieModel), intent(inout) :: this
end subroutine unset_useAeForHp
! ------------------------------------------------------------
subroutine set_kp(this, value)
  class(ieModel), intent(inout) :: this
  real, intent(in) :: value
end subroutine set_kp
! ------------------------------------------------------------
! Set the current time to run the empirical model
subroutine set_time_real(this, ut)
  class(ieModel) :: this
  real(kind=Real8_), intent(in) :: ut
  integer :: iYear
  integer :: iMonth
  integer :: iDay
  integer, dimension(7) :: itime
  real :: rHour
end subroutine set_time_real
! ------------------------------------------------------------
subroutine run_check_indices(this)
  class(ieModel) :: this
end subroutine run_check_indices


