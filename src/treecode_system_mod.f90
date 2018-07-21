module treecode_system_mod
  use nbody_system_mod
  use bhtree_mod
  use fconfig

  private
  public :: treecode_system

  type, extends(nbody_system) :: treecode_system
     type(bhtree) :: tree
     character(len=:), allocatable :: config_file
   contains
     procedure :: calculate_acceleration     
     procedure :: initialize
     procedure :: step
     procedure :: finalize
  end type treecode_system

  interface treecode_system
     module procedure new_treecode_system
  end interface treecode_system      

contains

  function new_treecode_system(config_file) result(system)
    type(treecode_system) :: system
    character(*), intent(in) :: config_file    
    system%config_file = config_file
    call system%conf%read_file(config_file)
  end function new_treecode_system

  
  subroutine initialize(system)
    class(treecode_system), intent(inout) :: system
    !type(bhtree) :: tree
    type(config) :: conf

    call conf%value_from_key("dt", system%dt)
    call conf%value_from_key("t_start", system%t_start, default_value = 0d0)
    call conf%value_from_key("t_end", system%t_end)
    call conf%value_from_key("random_seed", system%random_seed, default_value = 0)

    system%conf = conf
    system%current_step = 0
    system%t = system%t_start

    system%total_steps = nint((system%t_end - system%t_start) / system%dt)

  end subroutine initialize

  
  subroutine step(system)
    class(treecode_system), intent(inout) :: system
    call system%integrate()
  end subroutine step

  
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
