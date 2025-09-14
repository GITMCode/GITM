module ModGITMVersion

  ! Version variable is put into src/.version when compiling GITM
  ! Format is:
  ! (last commit date)_(commit hash)_(# of files changed from HEAD)
  ! - Gives info on which commit someone is working off of,
  !   and how many files were changed from that.

  INCLUDE ".version"

end module ModGITMVersion
