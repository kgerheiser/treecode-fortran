module nbody_system_mod
  use iso_fortran_env
  use body_mod
  use fconfig
  implicit none

  private
  public :: nbody_system

  type, abstract :: nbody_system
     real(real64) :: dt, t_end, t_start, t
     integer :: nbodies, random_seed, current_step, total_steps
     type(body_ptr), allocatable :: body_array(:)
     procedure(system_interface), pointer, private :: integrator
     type(config) :: conf
     logical :: verbose = .true.
   contains
     procedure(system_interface), deferred :: initialize
     procedure(system_interface), deferred :: step
     procedure(system_interface), deferred :: finalize
     procedure :: calculate_acceleration
     procedure :: set_integrator
     procedure :: integrate
     procedure :: should_advance
     procedure :: run
  end type nbody_system



  abstract interface
     subroutine system_interface(system)
       import
       class(nbody_system), intent(inout) :: system
     end subroutine system_interface
  end interface

contains


  subroutine set_integrator(system, integrator)
    class(nbody_system), intent(inout) :: system
    procedure(system_interface) :: integrator
    system%integrator => integrator
  end subroutine set_integrator


  subroutine integrate(system)
    class(nbody_system), intent(inout) :: system
    call system%integrator()
  end subroutine integrate


  subroutine run(system)
    class(nbody_system), intent(inout) :: system

    call system%initialize()

    do while(system%should_advance())
       if (system%verbose) print *, "Advancing system 1 step..."
       call system%step()
       system%current_step = system%current_step + 1
       system%t = system%t + system%dt
    end do

    call system%finalize()

  end subroutine run
  

  pure logical function should_advance(system)
    class(nbody_system), intent(in) :: system

    if (system%current_step < system%total_steps) then
       should_advance = .true.
    else
       should_advance = .false.
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

end module nbody_system_mod



