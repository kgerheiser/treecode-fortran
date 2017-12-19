module node_mod
  use constants_mod, only: prec, ndims
  implicit none

  integer, parameter :: body_type = 0
  integer, parameter :: cell_type = 1

  type, abstract :: node
     logical :: update = .true.
     real(prec) :: mass = 0.0_prec, pos(ndims) = 0.0_prec
     class(node), pointer :: next => null()
   contains
  end type node

  type :: node_ptr
     class(node), pointer :: ptr => null()
  end type node_ptr
  
contains

  pure function distance_to_squared(self, node_) result(dist2)
    class(node), intent(in) :: self, node_
    real(prec) :: dist2(ndims), dr(ndims)
    dr = node_%pos - self%pos
    dist2 = dot_product(dr, dr)
  end function distance_to_squared

  pure function distance_to(self, node_) result(dist)
    class(node), intent(in) :: self, node_
    real(prec) :: dist(ndims), diff(ndims)

    diff = node_%pos - self%pos
    dist = norm2(diff)
  end function distance_to

end module node_mod
