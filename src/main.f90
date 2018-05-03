program treecode
  use iso_fortran_env
  use constants_mod
  use treecode_system_mod
  use body_mod
  use bhtree_mod
  implicit none

  type(treecode_system) :: system

  call system%new_nbody_system("system_parameters.txt")
  call system%run()

end program treecode
