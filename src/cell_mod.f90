module cell_mod
  use constants_mod, only: prec, nsub, ndims
  use node_mod, only: node, node_ptr, cell_type
  use body_mod, only: body
  implicit none

  private
  public :: cell, cell_ptr, cell_list

  type, extends(node) :: cell
     real(prec) :: rcrit2, quad_moment(ndims,ndims)
     type(node), pointer :: more
     type(node_ptr) :: descendants(nsub)
   contains
     procedure :: sub_index => sub_index_cell
  end type cell

  type :: cell_ptr
     class(cell), pointer :: ptr => null()
   contains
     !procedure :: set_ptr
     procedure :: get_desc_ptr
     procedure :: set_desc_ptr
     procedure :: sub_index
     procedure :: has_child
  end type cell_ptr

  type :: cell_list_node
     type(cell_ptr) :: cellptr, next
  end type cell_list_node

  type :: cell_list
     type(cell_list_node) :: head 
     type(cell_list_node) :: current
     integer :: ncells
   contains
     procedure :: initialize
     procedure :: reset
     procedure :: get_free_cell
  end type cell_list

contains

!  subroutine set_descendant_ptr(self, ptr, i)

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

    integer :: i
    index = 0
    
    do i = 1, ndims
       if (self%ptr%pos(i) <= b%pos(i)) then
          index = index + ishft(nsub, -i)
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

  subroutine initialize(list)
    class(cell_list), intent(inout) :: list
    allocate(list%head%cellptr%ptr)
    list%current%cellptr%ptr => list%head%cellptr%ptr
    list%ncells = 0
  end subroutine initialize

  subroutine reset(list)
    class(cell_list), intent(inout) :: list
    list%current%cellptr%ptr => list%head%cellptr%ptr
    list%ncells = 0
  end subroutine reset
  
  function get_free_cell(list) result(freecell)
    class(cell_list), intent(inout) :: list
    type(cell), pointer :: freecell
 
    if (associated(list%current%cellptr%ptr)) then
       freecell => list%current%cellptr%ptr
    else
       allocate(list%current%next%ptr)
       list%current%next%ptr%type = cell_type
    end if
    
    list%current%cellptr%ptr => list%current%next%ptr
  end function get_free_cell
  
end module cell_mod
