module constants_mod
  use iso_fortran_env, only: real32, real64
  implicit none

  public

  integer, parameter :: ndims = 3
  integer, parameter :: nsub = 2**ndims
  integer, parameter :: prec = real32
contains
end module constants_mod
