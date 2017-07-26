
module bhtree_mod
  use constants_mod
  use node_mod, only: node, node_ptr
  use cell_mod, only: cell, cell_ptr
  use body_mod, only: body, body_ptr
  implicit none

  private
  public :: bhtree

  integer, parameter :: max_depth = 32

  type :: bhtree
     real(prec) :: theta, eps, rsize
     logical :: quick_scan, use_quad, bh86, sw94, first_call = .true.
     integer :: tdepth, ncells, num_cells
     integer :: cell_hist(max_depth), subn_hist(max_depth)
     class(cell), pointer :: root
     class(node), pointer :: freecell
   contains
     procedure :: expand_box
     procedure :: make_cell
     procedure :: new_tree
     procedure :: load_body
     procedure :: eval_center_of_mass
  end type bhtree

contains

  
  !*****************************************************************************
  ! Initialize tree for hierarchical force calculation
  !*****************************************************************************
  subroutine make_tree(tree, body_array)
    class(bhtree), intent(inout) :: tree
    type(body_ptr), intent(inout) :: body_array(:)
    class(node), pointer :: null_ptr

    integer :: i

    call tree%new_tree()
    
    tree%root => tree%make_cell()
    tree%root%pos = 0.0_prec

    call tree%expand_box(body_array)

    do i = 1, size(body_array)
       call tree%load_body(body_array(i)%ptr)
    end do

    if (tree%bh86 .and. tree%sw94) then
       stop "incompatible options bh86 and sw94"
    end if

    tree%tdepth = 0
    tree%cell_hist = 0
    tree%subn_hist = 0

    call tree%eval_center_of_mass(tree%root, tree%rsize, 0)
    null_ptr => null()
    call thread_tree(tree%root, null_ptr)

    if (tree%use_quad) then
       call eval_quadrupole_moment(tree%root)
    end if
    
  end subroutine make_tree

  
  !*****************************************************************************
  ! Reclaim cells in tree, and prepare to build a new one
  !*****************************************************************************
  subroutine new_tree(tree)
    class(bhtree), intent(inout) :: tree
    type(node_ptr) :: p

    if (.not. tree%first_call) then
       p%ptr => tree%root

       do while(associated(p%ptr))
          select type(a => p%ptr)
          class is (cell)
             p%ptr%next => tree%freecell
             tree%freecell => p%ptr
             p%ptr => a%more
          class is (body)
             p%ptr => p%ptr%next
          end select
       end do
    else
       tree%first_call = .false.
    end if

    tree%root => null()
    tree%ncells = 0
    
  end subroutine new_tree

  
  !*****************************************************************************
  ! Return a pointer to a free cell
  !*****************************************************************************
  function make_cell(tree) result(freecell)
    class(bhtree), intent(inout) :: tree
    type(cell), pointer :: freecell

    integer :: i

    if (.not. associated(tree%freecell)) then
       allocate(tree%freecell, mold = freecell)
       select type (c => tree%freecell)
       class is(cell)
          freecell => c
       class default
          stop "make_cell: mismatched type, should be cell 1"
       end select
    else
       select type(c => tree%freecell)
       class is(cell)
          freecell => c
          tree%freecell => freecell%next
       class default
          stop "make_cell: mismatched type, should be cell 2"
       end select
    end if

    freecell%update = .false.
    
    do i = 1, nsub
       freecell%descendants(i)%ptr => null()
    end do
    
    tree%num_cells = tree%num_cells + 1
    
  end function make_cell

  
  !*****************************************************************************
  ! Insert body into tree 
  !*****************************************************************************
  subroutine load_body(tree, b)
    class(bhtree) :: tree
    class(body), target :: b

    type(cell), pointer  :: q, c

    integer :: q_index, i, index
    real(prec) :: qsize, dmax, dist2, dist(ndims)


    q => tree%root
    q_index = q%sub_index(b)
    qsize = tree%rsize

    do while(associated(q%descendants(q_index)%ptr))
       select type(child => q%descendants(q_index)%ptr)
       class is(body)

          dist = b%pos - child%pos
          dist2 = dot_product(dist, dist)
          if (dist2 == 0.0_prec) stop "two bodies have same position"

          c => tree%make_cell()

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

    q%descendants(q_index)%ptr => b

  end subroutine load_body
  

  !*****************************************************************************
  ! Recursive walk of tree installing next and more links
  !*****************************************************************************
  recursive subroutine thread_tree(p, n)
    class(node), intent(inout) :: p
    class(node), target, intent(inout) ::  n

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
       desc(ndesc+1)%ptr =>  n

       do i = 1, ndesc
          call thread_tree(desc(i)%ptr, desc(i+1)%ptr)
       end do
       
    end select
    
  end subroutine thread_tree

  
  !*****************************************************************************
  ! Descend tree to find cells center-of-mass, and sets critical cell radii
  !*****************************************************************************
  recursive subroutine eval_center_of_mass(tree, p, psize, level)
    class(bhtree), intent(inout)    :: tree
    class(cell), intent(inout)  :: p
    real(prec),  intent(in)     :: psize
    integer,     intent(in)     :: level

    class(node), pointer :: q
    real(prec) :: cm_pos(ndims), tmpv(ndims)
    integer :: i, k

    tree%tdepth = max(tree%tdepth, level)
    tree%cell_hist(level) = tree%cell_hist(level) + 1
    
    p%mass = 0.0_prec
    cm_pos = 0.0_prec

    do i = 1, nsub
       q => p%descendants(i)%ptr
       if (associated(q)) then
          tree%subn_hist = tree%subn_hist + 1
          
          select type(q)
             class is(cell)
             call eval_center_of_mass(tree, q, psize/2.0_prec, level+1)
          end select

          p%update = p%update .or. q%update
          p%mass = p%mass + q%mass
          cm_pos = cm_pos + (q%pos * q%mass)
       end if
    end do

    if (p%mass > 0.0_prec) then
       cm_pos = cm_pos / p%mass
    else
       cm_pos = p%pos
    end if

    do k = 1, ndims
       if (cm_pos(k) < p%pos(k) - psize / 2.0_prec .or. &
            p%pos(k) + psize / 2.0_prec <= cm_pos(k)) then
          stop "eval_center_of_mass: tree structure error"
       end if
    end do

    if (.not. tree%quick_scan) then
       call set_critical_radius(tree, p, cm_pos, psize)
    end if
    
    p%pos = cm_pos
    
  end subroutine eval_center_of_mass

  
  !*****************************************************************************
  ! Assign critical radius of cell, which determines force evaluation accuracy
  ! TODO: make this a function?
  !*****************************************************************************
  subroutine set_critical_radius(tree, p, cm_pos, psize)
    type(bhtree), intent(in) :: tree
    type(cell), intent(inout) :: p
    real(prec), intent(in) :: cm_pos(ndims), psize

    real(prec) :: bmax2, d
    integer :: k

    if (tree%theta == 0.0_prec) then
       p%rcrit2 = (2.0_prec * tree%rsize)**2
    else if (tree%sw94) then
       bmax2 = 0.0_prec
       do k = 1, ndims
          d = cm_pos(k) - p%pos(k) + psize / 2.0_prec
          bmax2 = bmax2 + max(d, psize - d)**2
       end do
       p%rcrit2 = bmax2 / tree%theta**2
    else if (tree%bh86) then
       p%rcrit2 = (psize / tree%theta)**2
    else
       d = norm2(p%pos - cm_pos)
       p%rcrit2 = (psize / tree%theta + d)**2
    end if
    
  end subroutine set_critical_radius


  !*****************************************************************************
  ! Expand size of tree so that it can fit all bodies
  !*****************************************************************************
  subroutine expand_box(tree, body_array)
    class(bhtree), intent(inout) :: tree
    type(body_ptr), intent(in) :: body_array(:)

    real(prec) :: dmax, d
    integer i, j, nbodies

    nbodies = size(body_array)
    dmax = 0.0_prec

    do i = 1, nbodies
       do j = 1, ndims
          d = abs(body_array(i)%ptr%pos(j) - tree%root%pos(j))
          if (d > dmax) dmax = d
       end do
    end do

    do while(tree%rsize < 2.0_prec * dmax)
       tree%rsize = tree%rsize * 2.0_prec
    end do
    
  end subroutine expand_box


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

end module bhtree_mod

 
