module body_mod
  use constants_mod, only: prec, ndims, nsub
  use node_mod, only: node
  implicit none

  private
  public :: body, body_ptr

  type, extends(node) :: body
     real(prec) :: vel(ndims), acc(ndims), phi
   contains
  end type body

  type :: body_ptr
     type(body), pointer :: ptr => null()
  end type body_ptr
  
contains

  
end module body_mod
