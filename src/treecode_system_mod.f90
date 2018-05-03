module treecode_system_mod
  use nbody_system_mod
  use bhtree_mod
  use fconfig

  type, extends(nbody_system) :: treecode_system
     type(bhtree) :: tree
   contains
     procedure :: calculate_acceleration
     procedure :: initialize
     procedure :: finalize
  end type treecode_system

contains

  subroutine initialize(system)
    class(treecode_system), intent(inout) :: system
    type(bhtree) :: tree
    type(config) :: conf

    conf = system%conf
    call tree%read_config(conf)
    
  end subroutine initialize

    subroutine finalize(system)
    class(treecode_system), intent(inout) :: system
    
  end subroutine finalize

  subroutine calculate_acceleration(system)
    class(treecode_system), intent(inout) :: system
    integer :: i, n
    !call system%tree%make_tree(system%body_array)
    !call system%tree%gravcalc()
  end subroutine calculate_acceleration

end module treecode_system_mod
