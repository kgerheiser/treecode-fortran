module treecode_mod
  use iso_fortran_env
  use body_mod
  use fconfig
  implicit none

  private
  public :: nbody_system

  type, abstract :: nbody_system
     real(real64) :: dt, t_end, t_start, t
     integer :: nbodies, random_seed, nstep
     type(body_ptr), allocatable :: body_array(:)
     procedure(system_interface), pointer, private :: integrate => leap_frog_integrator
   contains
     procedure :: calculate_acceleration
     procedure :: set_integrator
     procedure :: should_advance
     procedure :: step_system
     procedure :: new_nbody_system
  end type nbody_system

  abstract interface
     subroutine system_interface(system)
       import
       class(nbody_system), intent(inout) :: system
     end subroutine system_interface
  end interface

contains

  subroutine new_nbody_system(system, file)
    character(*), intent(in) :: file
    class(nbody_system), intent(inout) :: system
    type(config) :: conf

    call conf%read_file(file)
    call conf%value_from_key("dt", system%dt)
    call conf%value_from_key("t_start", system%t_start)
    call conf%value_from_key("t_end", system%t_end)
    call conf%value_from_key("nstep", system%nstep, default_value = 0)
    call conf%value_from_key("random_seed", system%random_seed, default_value = 0)    
  end subroutine new_nbody_system
  
  subroutine set_integrator(system, integrator)
    class(nbody_system), intent(inout) :: system
    procedure(system_interface) :: integrator
    system%integrate => integrator
  end subroutine set_integrator

  subroutine step_system(system)
    class(nbody_system), intent(inout) :: system    
    call system%integrate()
    system%nstep = system%nstep + 1
    system%t = system%t + system%dt    
  end subroutine step_system

  subroutine run(system)
    class(nbody_system), intent(inout) :: system
    do while(system%should_advance())
       call system%step_system()
    end do
  end subroutine run

  pure logical function should_advance(system)
    class(nbody_system), intent(in) :: system
    if (system%t <= system%t_end) then
       should_advance = .false.
    else
       should_advance = .true.
    end if    
  end function should_advance

  subroutine leap_frog_integrator(system)
    class(nbody_system), intent(inout) :: system
    class(body), pointer :: b
    integer :: i, n
    real(real64) :: dt

    dt = system%dt
    n = system%nbodies

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

  end subroutine leap_frog_integrator

  subroutine calculate_acceleration(system)
    class(nbody_system), intent(inout) :: system
    integer :: i, j, n
    class(body), pointer :: b1, b2

    n = size(system%body_array)

    do i = 1, n
       b1 => system%body_array(i)%ptr
       do j = 1, n
          if (i /= j) then
             b2 => system%body_array(j)%ptr
          end if          
       end do       
    end do

  end subroutine calculate_acceleration

  
  

end module treecode_mod


module treecode_system_mod
  use treecode_mod
  use bhtree_mod
  
  type, extends(nbody_system) :: treecode_system
     type(bhtree) :: tree
   contains
     procedure :: calculate_acceleration     
  end type treecode_system

contains

  subroutine calculate_acceleration(system)
    class(treecode_system), intent(inout) :: system

    integer :: i, n

    do i = 1, system%nbodies
       system%body_array(i)%ptr%update = .true.
    end do

    call system%tree%make_tree(system%body_array)
    call system%tree%gravcalc()
    
  end subroutine calculate_acceleration
  
end module treecode_system_mod

program treecode
  use iso_fortran_env
  use constants_mod
  use treecode_mod
  use body_mod
  use bhtree_mod
  implicit none

  
  print *, "Hello, World!"

contains
  
  
end program treecode
