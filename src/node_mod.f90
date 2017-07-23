module node_mod
  use constants_mod, only: prec, ndims
  implicit none

  integer, parameter :: body_type = 0
  integer, parameter :: cell_type = 1

  type :: node
     logical :: update
     real(prec) :: mass, pos(ndims)
     class(node), pointer :: next => null()
     integer :: type
   contains
  end type node

  type :: node_ptr
     class(node), pointer :: ptr => null()
   contains
     procedure :: is_null
  end type node_ptr

  type :: node_list
     class(node), pointer :: head => null()
     class(node), pointer :: current => null()
     integer :: ncells
   contains
     procedure :: initialize
     procedure :: reset
     procedure :: get_free_cell
  end type node_list
  
contains

  pure logical function is_null(nodeptr)
    class(node_ptr), intent(in) :: nodeptr
    is_null = .not. associated(nodeptr%ptr)
  end function is_null

  subroutine initialize(list, kind)
    class(node_list), intent(inout) :: list
    class(node), intent(in) :: kind
    allocate(list%head, mold = kind)
    list%current => list%head
    list%ncells = 0
  end subroutine initialize

  subroutine reset(list)
    class(node_list), intent(inout) :: list
    list%current => list%head
    list%ncells = 0
  end subroutine reset
  
  function get_free_cell(list) result(freecell)
    class(node_list), intent(inout) :: list
    class(node), pointer :: freecell
 
    if (associated(list%current)) then
       freecell => list%current
    else
       allocate(list%current%next, mold = list%head)
    end if
    
    list%current => list%current%next
  end function get_free_cell

end module node_mod
