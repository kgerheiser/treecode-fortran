
module bhtree_mod
  use constants_mod
  use node_mod, only: node, node_ptr
  use cell_mod, only: cell, cell_ptr, cell_list
  use body_mod, only: body, body_ptr
  implicit none

  private
  public :: bhtree

  type :: bhtree
     real(prec) :: theta, eps, rsize
     logical :: use_quad, bh86, sw94
     integer :: tdepth, ncells
     type(cell_list) :: cells
     type(cell), pointer :: root
   contains
     procedure :: expand_tree
  end type bhtree

contains

  subroutine init_tree(tree)
    class(bhtree), intent(inout) :: tree
    allocate(tree%root)
    tree%root%next => null()
  end subroutine init_tree


  subroutine make_tree(tree, bodies)
    type(bhtree), intent(inout) :: tree
    type(body_ptr) :: bodies(:)


    tree%root%pos = 0.0_prec

  end subroutine make_tree

  subroutine load_body(tree, b)
    type(bhtree) :: tree
    class(body), target :: b

    type(cell), pointer  :: q, c

    !fix this, should probably not be pointer
    integer :: q_index, i, index
    real(prec) :: qsize, dmax


    q => tree%root
    q_index = q%sub_index(b)
    qsize = tree%rsize

    do while(associated(q%descendants(q_index)%ptr))
       select type(child => q%descendants(q_index)%ptr)
          class is(body)
             ! if distance is 0 they are the same
          c => tree%cells%get_free_cell()

          do i = 1, ndims
             if (b%pos(i) < q%pos(i)) then
                c%pos(i) = q%pos(i) - qsize / 4.0_prec
             else
                c%pos(i) = q%pos(i) + qsize / 4.0_prec
             end if
          end do

          c%descendants(c%sub_index(child))%ptr => child
          q%descendants(q_index)%ptr => c

          q => c
          class is(cell)
          q => child
       end select

       q_index = q%sub_index(b)
       qsize = qsize / 2.0_prec
    end do

  end subroutine load_body

  recursive subroutine thread_tree(p, n)
    class(node), pointer, intent(inout) :: p, n


    integer :: i, ndesc
    type(node_ptr) :: desc(nsub + 1)

    p%next => n

    select type(p)
    type is(cell)
       ndesc = 0

       do i = 1, nsub
          if (associated(p%descendants(i)%ptr)) then
             ndesc = ndesc + 1
             desc(ndesc)%ptr => p%descendants(i)%ptr
          end if
       end do


       p%more => desc(1)%ptr
       ndesc = ndesc + 1
       desc(ndesc)%ptr =>  n

       do i = 1, ndesc
          call thread_tree(desc(i)%ptr, desc(i+1)%ptr)
       end do
    end select

  end subroutine thread_tree


  ! recursive subroutine thread_tree(p, n)
  !   class(node_ptr), intent(inout) :: p, n

  !   integer :: i, ndesc
  !   type(node_ptr) :: descendants(nsub + 1)p

  !   p%ptr%next => n%ptr

  !   select type(a => p%ptr)
  !   class is(cell)
  !      ndesc = 0
  !      do i = 1, nsub
  !         if (associated(a%descendants(i)%ptr)) then
  !            ndesc = ndesc + 1
  !            descendants(ndesc)%ptr => a%descendants(i)%ptr
  !         end if
  !      end do

  !      a%more%ptr => descendants(1)%ptr
  !      ndesc = ndesc + 1
  !      descendants(ndesc)%ptr => n%ptr

  !      do i = 1, ndesc
  !         call thread_tree(descendants(i), descendants(i+1))
  !      end do

  !   end select

  ! end subroutine thread_tree

  !class is or type is

  recursive subroutine eval_center_of_mass(p, psize, level)
    class(cell), intent(inout)  :: p
    real(prec),  intent(in)     :: psize
    integer,     intent(in)     :: level

    class(node), pointer :: q
    integer :: tree_depth, i

    tree_depth = 1
    ! cell_hist(level)++

    do i = 1, nsub
       q => p%descendants(i)%ptr
       if (associated(q)) then
          !sub_history(level)++

          select type(a => q)
             class is(cell)
             call eval_center_of_mass(a, psize/2.0_prec, level+1)
          end select
       end if

    end do

  end subroutine eval_center_of_mass

  subroutine expand_tree(tree, body_array)
    class(bhtree), intent(inout) :: tree
    class(body), intent(in) :: body_array(:)

    real(prec) :: dmax, d
    integer i, j, nbodies

    nbodies = size(body_array)
    dmax = 0d0

    do i = 1, nbodies
       do j = 1, ndims
          d = abs(body_array(i)%pos(j) - tree%root%pos(j))
          if (d > dmax) dmax = d
       end do
    end do

    do while(tree%rsize < 2.0_prec * dmax)
       tree%rsize = tree%rsize * 2.0_prec
    end do

  end subroutine expand_tree


  ! ! Fix tree values
  ! subroutine set_critical_radius(p, cm_pos, psize)
  !   type(cell_ptr), intent(inout) :: p
  !   real(prec), intent(in) :: cm_pos(ndims), psize

  !   real(prec) :: bmax2, d
  !   integer :: k

  !   if (theta == 0.0_prec) then
  !      p%ptr%rcrit2 = (2.0_prec * tree%rsize)**2
  !   else if (tree%sw94) then
  !      bmax2 = 0.0_prec
  !      do k = 1, ndims
  !         d = cm_pos(k) - p%pos(k) + psize / 2.0_prec
  !         bmax2 = bmax2 + max(d, psize - d)**2
  !      end do
  !      p%rcrit2 = bmax2 / tree%theta**2
  !   else if(tree%bh86) then
  !      p%rcrit2 = (psize / tree%theta)**2
  !   else
  !      d = norm2(cm_pos, p%pos)
  !      p%rcrit2 = (psize / tree%theta + d)**2
  !   end if

  ! end subroutine set_critical_radius


  !perhaps I should have a descendant getter instead of accessing the array directly
  
  !*****************************************************************************
  ! Recursively descend tree and find the quadrupole moment of
  ! all cells
  !*****************************************************************************
  recursive subroutine eval_quadrupole_moment(p)
    type(cell), intent(inout) :: p

    integer :: ndesc, i, ix
    class(node), pointer :: q
    type(node_ptr) :: desc(nsub)!, q
    real(prec) :: dr(ndims), drsq
    real(prec), dimension(ndims, ndims) :: drdr, I_drsq, tmpm

    ndesc = 0

    do i = 1, nsub
       if (associated(p%descendants(i)%ptr)) then
          ndesc = ndesc + 1
          desc(ndesc)%ptr => p%descendants(i)%ptr
       end if
    end do

    p%quad_moment = 0.0_prec

    do i = 1, ndesc
       q => desc(i)%ptr

       select type (q)
       type is (cell)
          call eval_quadrupole_moment(q)
       end select

       dr = q%pos - p%pos
       drdr = outer_product(dr, dr)
       drsq = dot_product(dr, dr)
       I_drsq = identity() * drsq

       tmpm = drdr * 3.0_prec
       tmpm = (tmpm - I_drsq) * q%mass

       select type(q)
       type is (cell)
          tmpm = tmpm + q%quad_moment
       end select

       p%quad_moment = p%quad_moment + tmpm

    end do

  end subroutine eval_quadrupole_moment

  
  !*****************************************************************************
  ! Outer product of two vectors
  !*****************************************************************************
  pure function outer_product(u, v) result(res)
    real(prec), intent(in) :: u(ndims), v(ndims)
    real(prec) :: res(ndims, ndims)
    integer i, j

    do i = 1, ndims
       do j = 1, ndims
          res(j,i) = v(i) * u(j)
       end do
    end do
  end function outer_product

  
  !*****************************************************************************
  ! Creates an identity matrix
  !*****************************************************************************
  pure function identity() result(I)
    real(prec) :: I(ndims, ndims)
    integer :: k

    I = 0.0_prec

    do k = 1, ndims
       I(k,k) = 1.0_prec
    end do
  end function identity


  subroutine new_tree(tree)
    type(bhtree), intent(inout) :: tree
    logical, save :: first_call = .true.

    call tree%cells%reset()
  end subroutine new_tree


end module bhtree_mod

 
