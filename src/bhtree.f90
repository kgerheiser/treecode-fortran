
module bhtree_mod
  use constants_mod
  use node_mod, only: node, node_ptr
  use cell_mod, only: cell, cell_ptr
  use body_mod, only: body, body_ptr, is_same_body
  use fconfig
  implicit none

  private
  public :: bhtree

  integer, parameter :: max_depth = 32
 
  type :: bhtree
     real(prec) :: theta, eps2, rsize
     logical :: quick_scan, use_quad, bh86, sw94, first_call = .true., print_time = .false.
     integer :: tdepth, ncells, nbbcalc, nbccalc, actmax, nbodies, actlen
     integer :: cell_hist(max_depth), subn_hist(max_depth)
     type(cell), pointer :: root => null()
     class(node), pointer :: freecell => null()
   contains
     procedure :: expand_box
     procedure :: make_cell
     procedure :: new_tree
     procedure :: load_body
     procedure :: eval_center_of_mass
     procedure :: make_tree
     procedure :: gravcalc
  end type bhtree

  type, abstract :: cell_opening_criteria
  end type cell_opening_criteria

  type, extends(cell_opening_criteria) :: bh86_criterion
  end type bh86_criterion

  type, extends(cell_opening_criteria) :: sw94_criterion     
  end type sw94_criterion

  interface bhtree
     module procedure new_bhtree
  end interface bhtree

contains
  

  pure function new_bhtree(theta, epsilon, use_quad, quick_scan, opening_criterion) result(tree)
    type(bhtree) :: tree
    real(prec), intent(in) :: theta, epsilon
    logical, optional, intent(in) :: use_quad, quick_scan
    class(cell_opening_criteria), optional, intent(in) :: opening_criterion

    tree%theta = theta
    tree%eps2 = epsilon**2

    if (present(use_quad)) then
       tree%use_quad = use_quad       
    else
       tree%use_quad = .false.
    end if

    if (present(quick_scan)) then
       tree%quick_scan = quick_scan
    else
       tree%quick_scan = .false.
    end if

    if (present(opening_criterion)) then
       select type(opening_criterion)
       type is(bh86_criterion)
          tree%bh86 = .true.
       type is(sw94_criterion)
          tree%sw94 = .true.
        class default
           tree%bh86 = .true.
       end select
    end if

  end function new_bhtree

  !*****************************************************************************
  ! Initialize tree for hierarchical force calculation
  !*****************************************************************************
  subroutine make_tree(tree, body_array)
    class(bhtree), intent(inout) :: tree
    type(body_ptr), intent(inout) :: body_array(:)
    class(node), pointer :: null_ptr

    integer :: i

    call set_update(body_array)
    call tree%new_tree()

    tree%root => tree%make_cell()
    tree%root%pos = 0.0_prec
    tree%root%center = 0.0_prec

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
    tree%nbbcalc = 0
    tree%nbccalc = 0

    call tree%eval_center_of_mass(tree%root, tree%rsize, 1)
    null_ptr => null()
    call thread_tree(tree%root, null_ptr)

    if (tree%use_quad) then
       call eval_quadrupole_moment(tree%root)
    end if

  end subroutine make_tree

  subroutine set_update(body_array)
    type(body_ptr), intent(inout) :: body_array(:)
    integer :: i, n

    n = size(body_array)

    do i = 1, n
       body_array(i)%ptr%update = .true.
    end do

  end subroutine set_update

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
       allocate(freecell)
    else
       select type(c => tree%freecell)
       class is(cell)
          freecell => c
          tree%freecell => c%next
       end select
    end if

    freecell%update = .false.

    do i = 1, nsub
       freecell%descendants(i)%ptr => null()
    end do

    tree%ncells = tree%ncells + 1

  end function make_cell


  !*****************************************************************************
  ! Insert body into tree 
  !*****************************************************************************
  subroutine load_body(tree, p)
    class(bhtree) :: tree
    class(body), target :: p

    type(cell), pointer  :: q, c
    integer :: q_index, sub_index
    real(prec) :: qsize, dist2, dr(ndims)

    q => tree%root
    q_index = q%sub_index(p)
    qsize = tree%rsize

    do while(associated(q%descendants(q_index)%ptr))

       select type(child => q%descendants(q_index)%ptr)   
       class is(body)
          dr = p%pos - child%pos
          dist2 = dot_product(dr, dr)
          if (dist2 == 0.0_prec) stop "two bodies have same position"

          c => tree%make_cell()
          c%pos = q%pos + (merge(-qsize, qsize, p%pos < q%pos) / 4.0_prec)
          c%center = q%pos + (merge(-qsize, qsize, p%pos < q%pos) / 4.0_prec)
          c%length = qsize / 2.0_prec
          sub_index = c%sub_index(child)
          c%descendants(sub_index)%ptr => child
          q%descendants(q_index)%ptr => c
       end select

       select type(child => q%descendants(q_index)%ptr)
       class is(cell)
             q => child
       class default
             stop "not a cell"
       end select

       q_index = q%sub_index(p)
       qsize = qsize / 2.0_prec
    end do

    q%descendants(q_index)%ptr => p

  end subroutine load_body


  !*****************************************************************************
  ! Recursive walk of tree installing next and more links
  !*****************************************************************************
  recursive subroutine thread_tree(p, n)
    class(node), intent(inout) :: p
    class(node), pointer, intent(inout) ::  n

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
       desc(ndesc+1)%ptr => n

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
    real(prec) :: cm_pos(ndims)
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

    real(prec) :: dmax, d, rsize
    integer i, j, nbodies

    nbodies = size(body_array)
    dmax = 0.0_prec
    rsize = 1.0_prec

    do i = 1, nbodies
       do j = 1, ndims
          d = abs(body_array(i)%ptr%pos(j) - tree%root%center(j))
          if (d > dmax) dmax = d
       end do
    end do

    do while(rsize < 2.0_prec * dmax)
       rsize = rsize * 2.0_prec
    end do

    tree%rsize = rsize
    tree%root%length = rsize

  end subroutine expand_box


  !perhaps I should have a descendant getter instead of accessing the array directly

  !*****************************************************************************
  ! Recursively descend tree and find the quadrupole moment of
  ! all cells
  !*****************************************************************************
  recursive subroutine eval_quadrupole_moment(p)
    type(cell), intent(inout) :: p

    integer :: ndesc, i
    class(node), pointer :: q
    type(node_ptr) :: desc(nsub)
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
          res(j,i) = u(j) * v(i)
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


  !*****************************************************************************
  ! Perform force calculation on all particles
  !*****************************************************************************
  subroutine gravcalc(tree)
    class(bhtree), intent(inout) :: tree
    integer(int64) :: cpustart, cpuend, cputime, clock_count
    real(prec) :: rmid(ndims)

    type(node_ptr), allocatable :: active(:)
    type(cell_ptr), allocatable :: interact(:)
    integer :: i

    tree%actlen = 216 * tree%tdepth * 5

    allocate(active(tree%actlen))
    allocate(interact(tree%actlen))

    do i = 1, tree%actlen
       allocate(interact(i)%ptr)
    end do


    call system_clock(cpustart)
    tree%actmax = 0
    tree%nbbcalc = 0
    tree%nbccalc = 0

    active(1)%ptr => tree%root
    rmid = 0.0_prec

    call walk_tree(tree, active, 1, 2, interact, 1, tree%actlen, tree%root, tree%rsize, rmid)
    call system_clock(cpuend, count_rate = clock_count)

    cputime = cpuend - cpustart
    if (tree%print_time) print *, "time: ", cputime / real(clock_count)

    deallocate(active)
    deallocate(interact)

  end subroutine gravcalc


  !*****************************************************************************
  ! Walk through the tree to calculate forces in a single pass
  !*****************************************************************************
  recursive subroutine walk_tree(tree, active, aptr, nptr, interact, cptr, bptr, p, psize, pmid)
    class(bhtree), intent(inout) :: tree
    type(node_ptr), intent(inout) :: active(:)
    type(cell_ptr), intent(inout) :: interact(:)
    integer, value :: aptr, nptr, cptr, bptr
    class(node), intent(inout) :: p
    real(prec), intent(in) :: psize, pmid(ndims)

    integer :: i, actsafe, np 
    type(cell), pointer :: c
    class(node), pointer :: q

    if (p%update) then
       np = nptr
       actsafe = tree%actlen - nsub

       do i = aptr, nptr - 1
          select type(ap => active(i)%ptr)
          type is(cell)
             if (accept(ap, psize, pmid)) then
                c => interact(cptr)%ptr
                c%mass = ap%mass
                c%pos = ap%pos
                c%quad_moment = ap%quad_moment
                cptr = cptr + 1
             else
                q => ap%more
                do while(should_loop(q, ap%next))
                   active(np)%ptr => q
                   np = np + 1
                   q => q%next
                end do
             end if
          type is(body)
             if (.not. is_same_body(ap, p)) then
                interact(bptr)%ptr%mass = ap%mass
                interact(bptr)%ptr%pos = ap%pos
                bptr = bptr - 1
             end if
          end select
       end do

       !actmax = max(actmax, np - active)
       if (np /= nptr) then
          call walk_sub(tree, active, nptr, np, interact, cptr, bptr, p, psize, pmid)
       else
          select type(p)
          type is(body)
             call sum_gravity(tree, p, interact, cptr, bptr)
             class is(cell)
             stop "recursion ended with a cell"
          end select
       end if
    end if

  contains

    !*****************************************************************************
    ! Test next level's active list against subnodes of p
    !*****************************************************************************
    recursive subroutine walk_sub(tree, active, nptr, np, interact, cptr, bptr, p, psize, pmid)
      class(bhtree), intent(inout) :: tree
      type(node_ptr), intent(inout) :: active(:)
      type(cell_ptr), intent(inout) :: interact(:)
      integer, value :: nptr, np, cptr, bptr
      class(node), intent(inout) :: p
      real(prec), intent(in) :: psize, pmid(ndims)

      real(prec) :: poff, nmid(ndims)
      class(node), pointer :: q

      poff = psize / 4.0_prec

      select type(a => p)
      type is(cell)
         q => a%more
         do while(should_loop(q, a%next))
            nmid = pmid + merge(-poff, poff, q%pos < pmid)
            call walk_tree(tree, active, nptr, np, interact, cptr, bptr, q, psize / 2.0_prec, nmid)
            q => q%next
         end do
      type is(body)
         nmid = pmid + merge(-poff, poff, p%pos < pmid)
         call walk_tree(tree, active, nptr, np, interact, cptr, bptr, p, psize / 2.0_prec, nmid)
      end select

    end subroutine walk_sub

    !*****************************************************************************
    ! Used to control looping through nodes using threading idiom. Necessary
    ! because in Fortran null() pointers are not associated. Can possibly
    ! be fixed with having root's next point to a pointer rather than
    ! just null.
    !*****************************************************************************
    logical function should_loop(q, ap) result(loop)
      class(node), pointer :: q, ap
      if (.not. associated(q) .and. .not. associated(ap)) then
         loop = .false.
      else if (associated(q, ap)) then
         loop = .false.
      else if (.not. associated(q, ap)) then
         loop = .true.
      end if
    end function should_loop

  end subroutine walk_tree


  !*****************************************************************************
  ! Calculates the acceleration and potential on a body using the interaction list
  !*****************************************************************************
  subroutine sum_gravity(tree, p0, interact, cptr, bptr)
    class(bhtree), intent(inout) :: tree
    class(body), intent(inout) :: p0
    type(cell_ptr), intent(inout) :: interact(:)
    integer, intent(in) :: cptr, bptr

    real(prec) :: acc0_bb(ndims), acc0_bc(ndims), phi0_bb, phi0_bc

    phi0_bb = 0.0_prec
    phi0_bc = 0.0_prec

    acc0_bb = 0.0_prec
    acc0_bc = 0.0_prec

!!$omp parallel sections shared(tree, interact, p0, cptr, bptr) &
!!$omp& firstprivate(phi0_bc, acc0_bc, phi0_bb, acc0_bb) num_threads(2) &
!!$omp& lastprivate(phi0_bc, acc0_bc, phi0_bb, acc0_bb)
!!$omp section
    if (tree%use_quad) then
       call sum_cell(tree, interact(:cptr-1), p0, phi0_bc, acc0_bc)
    else
       call sum_node(tree, interact(:cptr-1), p0, phi0_bc, acc0_bc)
    end if

!!$omp section
    call sum_node(tree, interact(bptr+1:), p0, phi0_bb, acc0_bb)
!!$omp end parallel sections

    p0%phi = phi0_bc + phi0_bb
    p0%acc = acc0_bc + acc0_bb

    tree%nbbcalc = tree%nbbcalc + size(interact(bptr+1:))
    tree%nbccalc = tree%nbccalc + size(interact(:cptr-1))
  end subroutine sum_gravity


  !*****************************************************************************
  ! Sum body-cell interactions using cell's center-of-mass
  !*****************************************************************************
  subroutine sum_node(tree, interact, body0, phi0, acc0)
    class(bhtree), intent(in) :: tree
    class(body), intent(in) :: body0
    type(cell_ptr), intent(in) :: interact(1:)
    real(prec), intent(inout) :: phi0, acc0(ndims)

    real(prec) :: dr(ndims)
    real(prec) :: dr2, drab, phi_p, mr3i
    type(cell), pointer :: p
    integer :: i

    do i = 1, size(interact)
       p => interact(i)%ptr
       dr = p%pos - body0%pos
       dr2 = dot_product(dr, dr) + tree%eps2
       drab = sqrt(dr2)
       phi_p = p%mass / drab
       phi0 = phi0 - phi_p
       mr3i = phi_p / dr2
       acc0 = acc0 + (dr * mr3i)
    end do

  end subroutine sum_node


  !*****************************************************************************
  ! Sum body-cell interactions using quadrupole moment
  !*****************************************************************************
  subroutine sum_cell(tree, interact, body0 , phi0, acc0)
    class(bhtree), intent(in) :: tree
    class(body), intent(in) :: body0
    type(cell_ptr), intent(in) :: interact(:)
    real(prec), intent(inout) :: phi0, acc0(ndims)

    real(prec) :: dr(ndims), qdr(ndims)
    real(prec) :: dr2, drab, phi_p, mr3i, drqdr, dr5i, phi_q
    type(cell), pointer :: p
    integer :: i

    do i = 1, size(interact) 
       p => interact(i)%ptr
       dr = p%pos - body0%pos
       dr2 = dot_product(dr, dr) + tree%eps2
       drab = sqrt(dr2)

       phi_p = p%mass / drab
       mr3i = phi_p / dr2
       qdr = matmul(p%quad_moment, dr)
       drqdr = dot_product(dr, qdr)
       dr5i = 1.0_prec / (dr2**2 * drab)
       phi_q = 0.5_prec * dr5i * drqdr
       phi0 = phi0 - (phi_p + phi_q)
       mr3i = mr3i + 5.0_prec * phi_q / dr2
       acc0 = acc0 + (dr * mr3i) + (qdr * -dr5i)
    end do

  end subroutine sum_cell


  !*****************************************************************************
  ! 
  !*****************************************************************************
  pure logical function accept(this, psize, pmid)
    class(cell), intent(in) :: this
    real(prec), intent(in) :: psize, pmid(ndims)

    real(prec) :: dmax, dsq, dk
    integer :: k

    dmax = psize
    dsq = 0.0_prec

    do k = 1, ndims
       dk = this%pos(k) - pmid(k)

       if (dk < 0.0_prec) dk = -dk
       if (dk > dmax) dmax = dk

       dk = dk - (0.5_prec * psize)

       if (dk > 0.0_prec) dsq = dsq + dk**2
    end do

    accept = (dsq > this%rcrit2) .and. (dmax > (1.5_prec * psize))
  end function accept


end module bhtree_mod
