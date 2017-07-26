module cell_mod
  use constants_mod, only: prec, nsub, ndims
  use node_mod, only: node, node_ptr, cell_type
  use body_mod, only: body
  implicit none

  private
  public :: cell, cell_ptr

  type, extends(node) :: cell
     real(prec) :: rcrit2, quad_moment(ndims,ndims)
     type(node), pointer :: more
     type(node_ptr) :: descendants(nsub)
   contains
     procedure :: sub_index => sub_index_cell
  end type cell

  type :: cell_ptr
     class(cell), pointer :: ptr => null()
  end type cell_ptr

  type :: cell_list_node
     type(cell_ptr) :: cellptr, next
  end type cell_list_node

contains

  function get_desc_ptr(self, index) result(desc_ptr)
    class(cell_ptr), intent(in) :: self
    integer, intent(in) :: index
    class(node), pointer :: desc_ptr
    desc_ptr => self%ptr%descendants(index)%ptr
  end function get_desc_ptr

  subroutine set_ptr(self, cellptr)
    class(cell_ptr), intent(inout) :: self
    class(cell_ptr), intent(in) :: cellptr

    self%ptr => cellptr%ptr
  end subroutine set_ptr

  subroutine set_desc_ptr(self, nodeptr, index)
    class(cell_ptr), intent(inout) :: self
    class(node), target, intent(in) :: nodeptr
    integer, intent(in) :: index

    self%ptr%descendants(index)%ptr => nodeptr
  end subroutine set_desc_ptr

  pure logical function has_child(self, index) result(is_associated)
    class(cell_ptr), intent(in) :: self
    integer, intent(in) :: index
    is_associated = associated(self%ptr)
  end function has_child

  pure integer function sub_index(self, b) result(index)
    class(cell_ptr), intent(in) :: self
    class(node), intent(in) :: b

    integer :: k
    index = 0
    
    do k = 1, ndims
       if (self%ptr%pos(k) <= b%pos(k)) then
          index = index + ishft(nsub, -k)
       end if
    end do
    
  end function sub_index

   pure integer function sub_index_cell(self, b) result(index)
    class(cell), intent(in) :: self
    class(node), intent(in) :: b

    integer :: i
    index = 0
    
    do i = 1, ndims
       if (self%pos(i) <= b%pos(i)) then
          index = index + ishft(nsub, -i)
       end if
    end do
    
  end function sub_index_cell

  subroutine split(c)
    class(cell_ptr), intent(inout) :: c
    
  end subroutine split

   pure logical function is_null(cellptr)
    class(cell_ptr), intent(in) :: cellptr
    is_null = .not. associated(cellptr%ptr)
  end function is_null
  
end module cell_mod
