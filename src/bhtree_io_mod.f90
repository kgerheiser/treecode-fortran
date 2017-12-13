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
end module bhtree_io_mod
