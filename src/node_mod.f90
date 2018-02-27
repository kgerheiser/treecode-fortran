module node_mod
  use constants_mod, only: prec, ndims
  implicit none

  type, abstract :: node
     logical :: update
     real(prec) :: mass, pos(ndims)
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

  pure function acceleration_to(self, b) result(acc)
    class(node), intent(in) :: self, b
    real(prec) :: acc(ndims), dr(ndims), r2, r
    dr = b%pos - self%pos
    r2 = dot_product(dr, dr)
    r = sqrt(r2)    
    acc = b%mass / (r2 * r) * dr
  end function acceleration_to

end module node_mod
