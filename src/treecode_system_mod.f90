module treecode_system_mod
  use constants_mod
  use nbody_system_mod
  use bhtree_io_mod
  use bhtree_mod
  use body_mod
  use fconfig

  private
  public :: treecode_system

  type, extends(nbody_system) :: treecode_system
     type(bhtree) :: tree
     type(body_ptr), allocatable :: body_array(:)
     character(len=:), allocatable :: body_input_file
     type(config) :: conf
   contains
     procedure :: initialize
     procedure :: step
     procedure :: finalize
     procedure :: calculate_acceleration
     procedure :: checkpoint
     procedure :: write_diagnostics
     procedure :: integrate => integrate_verlet
  end type treecode_system

  type :: treecode_system_parameters
  end type treecode_system_parameters
  
  interface treecode_system
     module procedure new_treecode_system
  end interface treecode_system

contains

  function new_treecode_system(config_file) result(system)
    type(treecode_system) :: system
    character(*), intent(in) :: config_file
    call system%conf%read_file(config_file)
  end function new_treecode_system


  subroutine initialize(system)
    class(treecode_system), intent(inout) :: system
    real(prec) :: theta, epsilon
    logical :: use_quad, quick_scan    
    class(cell_opening_criteria), allocatable :: opening_criterion
    type(bh86_criterion) :: bh86

    allocate(opening_criterion, source = bh86)
    
    call system%conf%value_from_key("dt", system%dt)
    call system%conf%value_from_key("t_start", system%t_start)
    call system%conf%value_from_key("t_end", system%t_end) 
    call system%conf%value_from_key("use_quad", use_quad)
    call system%conf%value_from_key("quick_scan", quick_scan)
    call system%conf%value_from_key("theta", theta)
    call system%conf%value_from_key("body_input_file", system%body_input_file)
    call system%conf%value_from_key("dt_checkpoint", system%dt_checkpoint)
    call system%conf%value_from_key("dt_diag", system%dt_diag)
    
    system%body_array = read_bodies_from_formatted_file(system%body_input_file)

    system%current_step = 0
    system%t = system%t_start

    system%tree = bhtree(theta = theta, epsilon = epsilon, use_quad = use_quad, &
         quick_scan = quick_scan, opening_criterion = bh86)
    
    system%total_steps = nint((system%t_end - system%t_start) / system%dt)

  end subroutine initialize


  subroutine step(system)
    class(treecode_system), intent(inout) :: system
    call system%integrate()
  end subroutine step


  subroutine finalize(system)
    class(treecode_system), intent(inout) :: system
    call deallocate_body_array(system%body_array)
  end subroutine finalize

  
  subroutine calculate_acceleration(system)
    class(treecode_system), intent(inout) :: system
    integer :: i, n
    call system%tree%make_tree(system%body_array)
    call system%tree%gravcalc()
  end subroutine calculate_acceleration

  
  subroutine integrate_verlet(system)
    class(treecode_system), intent(inout) :: system
    class(body), pointer :: b
    integer :: i, n
    real(prec) :: dt

    dt = system%dt
    n = system%tree%nbodies

    do i = 1, n
       b => system%body_array(i)%ptr
       b%vel = b%vel + (b%acc * 0.5 * dt)
       b%pos = b%pos + (b%vel * dt)
    end do

    call system%calculate_acceleration()

    do i = 1, n
       b => system%body_array(i)%ptr
       b%vel = b%vel + (b%acc * 0.5 * dt)
    end do
    
  end subroutine integrate_verlet

  
  subroutine write_diagnostics(system)
    class(treecode_system), intent(inout) :: system
    call system%tree%print_diagnostics()
  end subroutine write_diagnostics

  
  subroutine checkpoint(system)
    class(treecode_system), intent(inout) :: system
  end subroutine checkpoint
  
end module treecode_system_mod
