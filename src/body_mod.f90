module body_mod
  use constants_mod, only: prec, ndims, nsub
  use node_mod, only: node
  implicit none

  private
  public :: body, body_ptr, is_same_body

  type, extends(node) :: body
     real(prec) :: vel(ndims) = 0.0_prec, acc(ndims) = 0.0_prec, phi = 0.0_prec
   contains
  end type body

  type :: body_ptr
     type(body), pointer :: ptr => null()
  end type body_ptr
  
contains

  pure logical function is_same_body(node1, node2)
    class(body), intent(in) :: node1
    class(node), intent(in) :: node2
    integer i

    is_same_body = .true.

    do i = 1, ndims
       if (node1%pos(i) /= node2%pos(i)) then
          is_same_body = .false.
          exit
       end if
    end do
    
  end function is_same_body
  
end module body_mod
