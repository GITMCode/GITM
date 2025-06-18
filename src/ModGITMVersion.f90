module ModGITMVersion

  ! Version information is put insto src/.version when compiling GITM
  ! Format is: 
  ! (last commit date)_(commit hash)_(# of files changed from HEAD)
  ! - Gives info on which commit someone is working off of, 
  !   and how many files were changed from that.
  
  character(25), parameter :: GitmVersion = &
  INCLUDE ".version"

end module ModGITMVersion
