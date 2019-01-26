program treecode
  use iso_fortran_env
  use constants_mod
  use treecode_system_mod
  use body_mod
  use bhtree_mod
  implicit none

  type(treecode_system) :: system
  integer :: command_count, status
  character(len=256) :: rc_file_name

  command_count = command_argument_count()

  if (command_count == 1) then
     call get_command_argument(1, rc_file_name, status = status)
  else
     print *, "There should be one command-line argument."
     print *, command_count, " found"
     stop
  end if
  

  system = treecode_system(rc_file_name)
  call system%run()

  print *, "treecode finished running"

end program treecode
