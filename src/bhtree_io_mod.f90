module bhtree_io_mod
  use constants_mod
  use body_mod, only: body_ptr
  implicit none
  
contains

  function read_bodies_from_formatted_file(file) result(body_array)
    character(*), intent(in) :: file
    integer :: i, file_unit, n
    type(body_ptr), allocatable :: body_array(:)

    open(file = file, newunit = file_unit, form = 'formatted', status = 'old', action = 'read')

    read(file_unit, *) n
    allocate(body_array(n))

    do i = 1, n
       associate(b => body_array(i))
         allocate(b%ptr)
         read(file_unit, *) b%ptr%mass
         read(file_unit, *) b%ptr%pos
         read(file_unit, *) b%ptr%vel
       end associate
    end do

    close(file_unit)
           
  end function read_bodies_from_formatted_file

  
  subroutine write_bodies_to_formatted_file(file, body_array)
    character(*), intent(in) :: file
    type(body_ptr), intent(in) :: body_array(:)

    integer i, n, file_unit


    open(file = file, newunit = file_unit, form = 'formatted', status = 'replace', action = 'write')

    n = size(body_array)
    write(file_unit, *) n

    do i = 1, n
       associate(b => body_array(i))
         write(file_unit, *) b%ptr%mass
         write(file_unit, *) b%ptr%pos
         write(file_unit, *) b%ptr%vel
       end associate
    end do

    close(file_unit)

  end subroutine write_bodies_to_formatted_file

  pure function allocate_empty_body_array(n) result(body_array)
    integer, intent(in) :: n
    type(body_ptr), allocatable :: body_array(:)
    integer :: i

    allocate(body_array(n))

    do i = 1, n
       allocate(body_array(i)%ptr)
    end do
 
  end function allocate_empty_body_array

  subroutine deallocate_body_array(body_array)
    type(body_ptr), allocatable :: body_array(:)
    integer :: i, n

    n = size(body_array)

    do i = 1, n
       deallocate(body_array(i)%ptr)
    end do

    deallocate(body_array)
    
  end subroutine deallocate_body_array
  
end module bhtree_io_mod
