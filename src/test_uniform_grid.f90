#include "cpp_macros.h"
program test_one_level
  use mpi
  use m_data_structures
  use m_build_tree
  use m_load_balance
  use m_ghost_cells
  use m_allocate_storage
  use m_restrict
  use m_communication
  use m_prolong
  use m_multigrid

  implicit none

  integer, parameter  :: block_size  = 16
  integer, parameter  :: domain_size = 128
  real(dp), parameter :: dr          = 1.0_dp / block_size
  real(dp), parameter :: pi          = acos(-1.0_dp)

  integer    :: lvl, n, ierr
  real(dp)   :: t0, t1
  type(mg_t) :: mg

  mg%boundary_cond => my_bc
  mg%n_cycle_up    =  2
  mg%n_cycle_down  =  2
  mg%smoother_type = smoother_gsrb
  mg%residual_coarse_abs = 1e-10_dp
  mg%residual_coarse_rel = 1e-10_dp

  call comm_init(mg)
  call build_uniform_tree(mg, block_size, domain_size, dr)

  if (mg%my_rank == 0) then
     do lvl = mg%lowest_lvl, mg%highest_lvl
        print *, lvl, ":", size(mg%lvls(lvl)%ids), mg%box_size_lvl(lvl)
     end do
  end if

  call load_balance(mg)
  call allocate_storage(mg)
  call set_initial_conditions(mg)

  do n = 1, num_neighbors
     allocate(mg%bc(n)%d(mg%box_size))
     mg%bc(n)%d = 0.0_dp
     mg%bc(n)%bc_type = bc_dirichlet
  end do

  do lvl = mg%lowest_lvl, mg%highest_lvl
     call fill_ghost_cells_lvl(mg, lvl)
  end do

  call print_error(mg)

  t0 = mpi_wtime()
  do n = 1, 10
     call mg_fas_fmg(mg, n==10, n > 1)
     ! call mg_fas_vcycle(mg, .true.)
     call print_error(mg)
  end do
  t1 = mpi_wtime()
  if (mg%my_rank == 0) then
     print *, "time, time per iteration:", t1-t0, (t1-t0) / 10
     print *, "unknowns/microsec", 1e-6_dp * 10 * domain_size**NDIM / (t1-t0)
  end if
  call timers_show(mg)

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)

contains

  subroutine set_initial_conditions(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, id, lvl, nc, IJK
    real(dp)                  :: r(NDIM), sol

    do lvl = 1, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(0, mg%box_size+1)
#if NDIM == 2
             r = mg%boxes(id)%r_min + &
                  [i-0.5_dp, j-0.5_dp] * mg%dr(lvl)
             sol = product(sin(2 * pi * r))
             mg%boxes(id)%cc(i, j, i_phi) = sol
#elif NDIM == 3

             r = mg%boxes(id)%r_min + &
                  [i-0.5_dp, j-0.5_dp, k-0.5_dp] * mg%dr(lvl)
             sol = product(sin(2 * pi * r))
             mg%boxes(id)%cc(i, j, k, i_phi) = sol
#endif
          end do; CLOSE_DO

          call box_lpl(mg, id, nc, i_rhs)
          mg%boxes(id)%cc(DTIMES(:), i_phi) = 0
       end do
    end do
  end subroutine set_initial_conditions

  subroutine print_error(mg)
    type(mg_t), intent(inout) :: mg
    integer                      :: n, id, lvl, IJK, ierr
    real(dp)                     :: r(NDIM), sol, val, err, max_err

    err = 0.0_dp

    do lvl = mg%highest_lvl, mg%highest_lvl
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(1, mg%box_size)
#if NDIM == 2
             r = mg%boxes(id)%r_min + &
                  [i-0.5_dp, j-0.5_dp] * mg%dr(lvl)
#elif NDIM == 3
             r = mg%boxes(id)%r_min + &
                  [i-0.5_dp, j-0.5_dp, k-0.5_dp] * mg%dr(lvl)
#endif
             sol = product(sin(2 * pi * r))
             val = mg%boxes(id)%cc(IJK, i_phi)
             err = max(err, abs(val-sol))
             ! err = max(err, abs(mg%boxes(id)%cc(IJK, i_res)))
             ! print *, lvl, id, i, j, r, sol, val, abs(sol-val)
          end do; CLOSE_DO
       end do
    end do

    call mpi_reduce(err, max_err, 1, MPI_DOUBLE, MPI_MAX, 0, &
         mpi_comm_world, ierr)
    if (mg%my_rank == 0) print *, "max err", max_err
  end subroutine print_error

  subroutine my_bc(mg, id, nc, nb, bc_type)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: nb      !< Direction
    integer, intent(out)      :: bc_type !< Type of b.c.

    call set_bc_dirichlet_zero(mg, id, nc, nb, bc_type)
  end subroutine my_bc

end program test_one_level
