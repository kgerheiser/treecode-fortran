module nbody_system_mod
  use, intrinsic :: iso_fortran_env
  use fconfig
  implicit none

  private
  public :: nbody_system

  type, abstract :: nbody_system
     character(len=:), allocatable :: name
     real(real64) :: dt, dt_diag, dt_checkpoint, t_end, t_start, t, t_diag, t_checkpoint
     integer :: current_step, total_steps
     procedure(system_interface), pointer, private :: integrator
     logical :: verbose = .false.
   contains
     procedure(system_interface), deferred :: initialize
     procedure, private :: initialize_private
     procedure(system_interface), deferred :: step
     procedure(system_interface), deferred :: finalize
     procedure(system_interface), deferred :: checkpoint
     procedure(system_interface), deferred :: write_diagnostics
     !procedure(system_interface), deferred :: read_checkpoint
     procedure, private :: should_step, should_checkpoint
     procedure, private :: should_write_diagnostics
     procedure :: run
  end type nbody_system

  abstract interface
     subroutine system_interface(system)
       import
       class(nbody_system), intent(inout) :: system
     end subroutine system_interface
  end interface

contains


  subroutine run(system)
    class(nbody_system), intent(inout) :: system
    integer(int64) :: t_start, t_end, count_rate
    real(real64) :: elapsed_time

    if (system%verbose) print *, "initializing"
    call system%initialize()
    call system%initialize_private()

    if (system%verbose) print *, "starting"
    do while(system%should_step())
       if (system%verbose) print *, "advancing one step..."
       call system%step()
       system%current_step = system%current_step + 1
       system%t = system%t + system%dt

       if (system%should_write_diagnostics()) then
          if (system%verbose) print *, "writing diagnostics"
          call system%write_diagnostics()
          system%t_diag = system%t_diag + system%dt_diag
       end if

       if (system%should_checkpoint()) then
          if (system%verbose) print *, "checkpointing"
          call system%checkpoint()
          system%t_checkpoint = system%t_checkpoint + system%dt_checkpoint
       end if

       !print '(a, f0.2)', "time: ", system%t
    end do

    if (system%verbose) print *, "finalizing"
    call system%finalize()

  end subroutine run
  

  subroutine initialize_private(system)
    class(nbody_system), intent(inout) :: system

    system%t_end = (system%t_end + system%t_start) - (0.5d0 * system%dt)
    system%t_diag = (system%dt_diag + system%t_start) - (0.5d0 * system%dt)
    system%t_checkpoint = (system%dt_checkpoint + system%t_start) - (0.5d0 * system%dt)
    
  end subroutine initialize_private
  
  pure logical function should_step(system)
    class(nbody_system), intent(in) :: system

    if (system%t <= system%t_end) then
       should_step = .true.
    else
       should_step = .false.
    end if
  end function should_step

  
  subroutine checkpoint(system)
    class(nbody_system), intent(in) :: system
  end subroutine checkpoint

  
  pure logical function should_checkpoint(system)
    class(nbody_system), intent(in) :: system
    should_checkpoint = .false.
    if (system%t >= system%t_checkpoint) then
       should_checkpoint = .true.
    end if    
  end function should_checkpoint

  
  pure logical function should_write_diagnostics(system)
    class(nbody_system), intent(in) :: system
    should_write_diagnostics = .false.
    if (system%t >= system%t_diag) then
       should_write_diagnostics = .true.
    end if    
  end function should_write_diagnostics


end module nbody_system_mod



