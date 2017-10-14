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
  end type node_ptr
  
contains

  pure function distance_to_squared(self, node_) result(dist2)
    class(node), intent(in) :: self, node_
    real(prec) :: dist2(ndims), diff(ndims)

    diff = node_%pos - self%pos
    dist2 = dot_product(diff, diff)
  end function distance_to_squared

  pure function distance_to(self, node_) result(dist)
    class(node), intent(in) :: self, node_
    real(prec) :: dist(ndims), diff(ndims)

    diff = node_%pos - self%pos
    dist = norm2(diff)
  end function distance_to

end module node_mod
