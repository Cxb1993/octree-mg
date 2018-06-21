!> Module to use free space boundary conditions for 3D Poisson problems
module m_free_space
  use m_data_structures

  implicit none
  private

#if NDIM == 3
  type mg_free_bc_t
     logical               :: initialized = .false.
     integer               :: fft_lvl
     real(dp), allocatable :: rhs(:, :, :)
     real(dp), pointer     :: karray(:)
     real(dp)              :: inv_dr(2, mg_num_neighbors)
     real(dp)              :: r_min(2, mg_num_neighbors)
     real(dp), allocatable :: bc_x0(:, :)
     real(dp), allocatable :: bc_x1(:, :)
     real(dp), allocatable :: bc_y0(:, :)
     real(dp), allocatable :: bc_y1(:, :)
     real(dp), allocatable :: bc_z0(:, :)
     real(dp), allocatable :: bc_z1(:, :)
  end type mg_free_bc_t

  type(mg_free_bc_t) :: free_bc

  ! Public methods
  public :: mg_poisson_free_3d

contains

  subroutine mg_poisson_free_3d(mg, max_fft_frac)
    use mpi
    use poisson_solver
    use m_multigrid
    use m_restrict
    use m_prolong
    use m_ghost_cells
    type(mg_t), intent(inout)   :: mg
    !> How much smaller the fft solve has to be than the full multigrid (0.0-1.0)
    real(dp), intent(in)        :: max_fft_frac
    integer                     :: fft_lvl, lvl, n, id, nx(3), nc
    integer                     :: ix(3), ierr, n_boxes_lvl
    real(dp)                    :: dr(3)
    real(dp)                    :: max_res
    real(dp), allocatable       :: tmp(:, :, :)
    real(dp)                    :: dummy(1)
    real(dp), parameter         :: offset    = 0.0_dp
    real(dp)                    :: ehartree, eexcu, vexcu
    integer                     :: i3sd, ncomp
    character(len=*), parameter :: geocode   = 'F'
    character(len=*), parameter :: datacode  = 'G'
    integer, parameter          :: itype_scf = 8
    integer, parameter          :: ixc       = 0

    ! Correction factor for the source term 1 / (4 * pi)
    real(dp), parameter :: rhs_fac = -1 / (4 * acos(-1.0_dp))

    if (.not. free_bc%initialized) then
       ! Determine highest fully refined grid level
       do lvl = mg_highest_uniform_lvl(mg), mg%first_normal_lvl+1, -1
          ! Determine how many boxes the level contains
          n_boxes_lvl = size(mg%lvls(lvl)%ids)

          ! If the level is 'small enough', exit
          if (n_boxes_lvl <= ceiling(max_fft_frac * mg%n_boxes)) exit
       end do

       fft_lvl         = lvl
       free_bc%fft_lvl = lvl

       ! Restrict rhs to required level
       do lvl = mg%highest_lvl, fft_lvl+1, -1
          call mg_restrict_lvl(mg, mg_irhs, lvl)
       end do

       ! Add a layer of ghost cells around the domain
       nx(:) = mg%domain_size_lvl(:, fft_lvl) + 2
       dr(:) = mg%dr(:, fft_lvl)

       ! Create kernel of Green's function
       call createKernel(geocode, nx(1), nx(3), nx(3), dr(1), dr(2), dr(3),  &
            itype_scf, mg%my_rank, mg%n_cpu, free_bc%karray)

       allocate(free_bc%rhs(nx(1), nx(2), nx(3)))

       ! For interpolation of the boundary planes
       free_bc%inv_dr(:, mg_neighb_lowx) = 1 / dr(2:3)
       free_bc%r_min(:, mg_neighb_lowx)  = mg%r_min(2:3) - 0.5_dp * dr(2:3)
       free_bc%inv_dr(:, mg_neighb_lowy) = 1 / dr([1,3])
       free_bc%r_min(:, mg_neighb_lowy)  = mg%r_min([1,3]) - 0.5_dp * dr([1,3])
       free_bc%inv_dr(:, mg_neighb_lowz) = 1 / dr(1:2)
       free_bc%r_min(:, mg_neighb_lowz)  = mg%r_min(1:2) - 0.5_dp * dr(1:2)

       do n = mg_neighb_lowx+1, mg_num_neighbors, 2
          free_bc%inv_dr(:, n) = free_bc%inv_dr(:, n-1)
          free_bc%r_min(:, n)  = free_bc%r_min(:, n-1)
       end do

       ! Set boundary conditions for multigrid solver
       do n = 1, mg_num_neighbors
          mg%bc(n, mg_iphi)%boundary_cond => ghost_cells_free_bc
       end do

       allocate(tmp(nx(1), nx(2), nx(3)))
       tmp(:, :, :) = 0.0_dp

       ! Store right-hand side
       nc = mg%box_size_lvl(fft_lvl)
       do n = 1, size(mg%lvls(fft_lvl)%my_ids)
          id = mg%lvls(fft_lvl)%my_ids(n)
          ix = (mg%boxes(id)%ix - 1) * nc + 1
          tmp(ix(1)+1:ix(1)+nc, ix(2)+1:ix(2)+nc, ix(3)+1:ix(3)+nc) = &
               rhs_fac * mg%boxes(id)%cc(1:nc, 1:nc, 1:nc, mg_irhs)
       end do

       call mpi_allreduce(tmp, free_bc%rhs, product(shape(tmp)), MPI_DOUBLE, &
            MPI_SUM, mg%comm, ierr)

       i3sd  = 1
       ncomp = nx(3)

       ! Solve free-space Poisson's equation
       call PSolver(geocode, datacode, mg%my_rank, mg%n_cpu, &
            nx(1), nx(2), nx(3), ixc, dr(1), dr(2), dr(3), &
            free_bc%rhs(1, 1, i3sd), free_bc%karray, dummy, &
            ehartree, eexcu, vexcu, offset, .false., 1)

       ! Extract boundary planes by interpolation
       associate (rhs => free_bc%rhs)
         free_bc%bc_x0 = 0.5_dp * (rhs(1, :, :) + rhs(2, :, :))
         free_bc%bc_x1 = 0.5_dp * (rhs(nx(1)-1, :, :) + rhs(nx(1), :, :))
         free_bc%bc_y0 = 0.5_dp * (rhs(:, 1, :) + rhs(:, 2, :))
         free_bc%bc_y1 = 0.5_dp * (rhs(:, nx(2)-1, :) + rhs(:, nx(2), :))
         free_bc%bc_z0 = 0.5_dp * (rhs(:, :, 1) + rhs(:, :, 2))
         free_bc%bc_z1 = 0.5_dp * (rhs(:, :, nx(3)-1) + rhs(:, :, nx(3)))
       end associate

       ! Use solution as an initial guess
       nc = mg%box_size_lvl(fft_lvl)
       do n = 1, size(mg%lvls(fft_lvl)%my_ids)
          id = mg%lvls(fft_lvl)%my_ids(n)
          ix = (mg%boxes(id)%ix - 1) * nc + 1
          mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, mg_iphi) = &
               free_bc%rhs(ix(1):ix(1)+nc+1, ix(2):ix(2)+nc+1, &
               ix(3):ix(3)+nc+1)
       end do

       ! Restrict FFT solution
       do lvl = fft_lvl, mg%lowest_lvl+1, -1
          call mg_restrict_lvl(mg, mg_iphi, lvl)
       end do

       ! Prolong FFT solution
       do lvl = fft_lvl, mg%highest_lvl-1
          call mg_prolong(mg, lvl, mg_iphi, mg_iphi, mg%box_prolong, .false.)
          ! We can already use the boundary conditions from the FFT solution
          call mg_fill_ghost_cells_lvl(mg, lvl+1, mg_iphi)
       end do

       free_bc%initialized = .true.
    end if

    ! Avoid multigrid solver if we already have the full solution
    if (free_bc%fft_lvl < mg%highest_lvl) then
       ! Solve Poisson equation with free space boundary conditions
       call mg_fas_fmg(mg, .true., max_res=max_res)
    end if

  end subroutine mg_poisson_free_3d

  !> To fill ghost cells
  subroutine ghost_cells_free_bc(box, nc, iv, nb, bc_type, bc)
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< Direction
    integer, intent(out)          :: bc_type !< Type of b.c.
    double precision, intent(out) :: bc(nc, nc)
    double precision              :: rr(nc, nc, 3)
    integer                       :: ixs(2), nb_dim

    bc_type      = mg_bc_dirichlet
    nb_dim       = mg_neighb_dim(nb)
    ixs          = [1, 2]
    ixs(nb_dim:) = ixs(nb_dim:) + 1

    call mg_get_face_coords(box, nb, nc, rr)

    bc = interp_bc(free_bc, nb, rr(:, :, ixs(1)), rr(:, :, ixs(2)))
  end subroutine ghost_cells_free_bc

  elemental function interp_bc(bc, nb, x1, x2) result(val)
    type(mg_free_bc_t), intent(in) :: bc
    integer, intent(in)            :: nb
    real(dp), intent(in)           :: x1, x2
    real(dp)                       :: val
    integer                        :: ix(2)
    real(dp)                       :: frac(2), low_frac(2)
    real(dp)                       :: w(2, 2)

    frac     = ([x1, x2] - bc%r_min(:, nb)) * bc%inv_dr(:, nb)
    ix       = ceiling(frac)
    low_frac = ix - frac

    ! Bilinear interpolation
    w(1, 1) = low_frac(1) * low_frac(2)
    w(2, 1) = (1 - low_frac(1)) * low_frac(2)
    w(1, 2) = low_frac(1) * (1 - low_frac(2))
    w(2, 2) = (1 - low_frac(1)) * (1 - low_frac(2))

    select case (nb)
    case (mg_neighb_lowx)
       val = sum(w * bc%bc_x0(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_highx)
       val = sum(w * bc%bc_x1(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_lowy)
       val = sum(w * bc%bc_y0(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_highy)
       val = sum(w * bc%bc_y1(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_lowz)
       val = sum(w * bc%bc_z0(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_highz)
       val = sum(w * bc%bc_z1(ix(1):ix(1)+1, ix(2):ix(2)+1))
    end select
  end function interp_bc
#endif
end module m_free_space
