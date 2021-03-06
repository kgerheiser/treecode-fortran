module cell_mod
  use constants_mod, only: prec, nsub, ndims
  use node_mod, only: node, node_ptr
  use body_mod, only: body
  implicit none

  private
  public :: cell, cell_ptr

  type, extends(node) :: cell
     real(prec) :: rcrit2, quad_moment(ndims,ndims), center(ndims), length
     class(node), pointer :: more => null()
     type(node_ptr) :: descendants(nsub)
   contains
     procedure :: sub_index
     procedure :: in_cell
  end type cell

  type :: cell_ptr
     class(cell), pointer :: ptr => null()
  end type cell_ptr

contains

 pure integer function sub_index(self, b) result(index)
    class(cell), intent(in) :: self
    class(node), intent(in) :: b

    integer :: k
    index = 0

    do k = 1, ndims
       if (self%pos(k) <= b%pos(k)) then
          index = index + ishft(nsub, -k)
       end if
    end do

    index = index + 1

  end function sub_index

  ! pure type(cell) function split(self, index) 
  !   class(cell), intent(in) :: self
  !   integer, intent(in) :: index    
  ! end function split

  pure logical function in_cell(self, node_)
    class(cell), intent(in) :: self
    class(node), intent(in) :: node_
    in_cell = .not. any(abs(node_%pos - self%center) > self%length / 2.0)
  end function in_cell

  
  
  
end module cell_mod
