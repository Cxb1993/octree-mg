#include "cpp_macros.h"
module m_phi_bc
  use m_data_structures

  implicit none
  private

  public :: mg_phi_bc_buffer_size
  public :: mg_phi_bc_store

contains

  !> Specify minimum buffer size (per process) for communication
  subroutine mg_phi_bc_buffer_size(mg, n_send, n_recv, dsize)
    use mpi
    type(mg_t), intent(inout) :: mg
    integer, intent(out)      :: n_send(0:mg%n_cpu-1)
    integer, intent(out)      :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)      :: dsize
    integer                   :: min_lvl, ierr
    integer                   :: i, id, lvl, nb, nb_id, p_id, p_rank

    allocate(mg%comm_phibc%n_send(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))
    allocate(mg%comm_phibc%n_recv(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))

    mg%comm_phibc%n_send(:, :) = 0
    dsize = (mg%box_size/2)**(NDIM-1)
    min_lvl = max(mg%lowest_lvl+1, mg%first_normal_lvl)

    do lvl = min_lvl, mg%highest_lvl
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          do nb = 1, mg_num_neighbors
             nb_id = mg%boxes(id)%neighbors(nb)
             if (nb_id < mg_no_box) then
                p_id = mg%boxes(id)%parent
                p_rank = mg%boxes(p_id)%rank
                if (mg%my_rank /= p_rank) then
                   mg%comm_phibc%n_send(p_rank, lvl) = &
                        mg%comm_phibc%n_send(p_rank, lvl) + 1
                end if
             end if
          end do
       end do
       call mpi_alltoall(mg%comm_phibc%n_send(:, lvl), 1, MPI_INTEGER, &
            mg%comm_phibc%n_recv(:, lvl), 1, MPI_INTEGER, mg%comm, ierr)
    end do

    n_send = maxval(mg%comm_phibc%n_send, dim=2)
    n_recv = maxval(mg%comm_phibc%n_recv, dim=2)
  end subroutine mg_phi_bc_buffer_size

  subroutine mg_phi_bc_store(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: lvl, nc

    nc = mg%box_size

    do lvl = mg%first_normal_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       call mg_phi_bc_store_lvl(mg, lvl, nc)
    end do

    call mg_phi_bc_restrict(mg)

    mg%phi_bc_data_stored = .true.

  end subroutine mg_phi_bc_store

  subroutine mg_phi_bc_store_lvl(mg, lvl, nc)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer, intent(in)       :: nc
#if NDIM == 2
    real(dp)                  :: bc(nc)
#elif NDIM == 3
    real(dp)                  :: bc(nc, nc)
#endif
    integer                   :: i, id, nb, nb_id, bc_type

    do i = 1, size(mg%lvls(lvl)%my_leaves)
       id = mg%lvls(lvl)%my_leaves(i)
       do nb = 1, mg_num_neighbors
          nb_id = mg%boxes(id)%neighbors(nb)
          if (nb_id < mg_no_box) then
             ! Physical boundary
             if (associated(mg%bc(nb, mg_iphi)%boundary_cond)) then
                call mg%bc(nb, mg_iphi)%boundary_cond(mg%boxes(id), nc, &
                     mg_iphi, nb, bc_type, bc)
             else
                bc_type = mg%bc(nb, mg_iphi)%bc_type
                bc      = mg%bc(nb, mg_iphi)%bc_value
             end if

             ! Store the boundary condition type. This is not globally set in
             ! the tree, but all negative values are treated the same in
             ! other parts of the code
             mg%boxes(id)%neighbors(nb) = bc_type

             ! Store ghost cell data in the right-hand side
             call box_set_gc(mg%boxes(id), nb, nc, mg_irhs, bc)
          end if
       end do
    end do
  end subroutine mg_phi_bc_store_lvl

  subroutine mg_phi_bc_restrict(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: lvl

    do lvl = mg%highest_lvl, mg%lowest_lvl+1, -1
       call mg_phi_bc_restrict_lvl(mg, lvl)
    end do
  end subroutine mg_phi_bc_restrict

  subroutine mg_phi_bc_restrict_lvl(mg, lvl)
    use m_communication
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id, nb
    integer                   :: nc, dsize

    if (lvl <= mg%lowest_lvl) &
         error stop "cannot restrict phi_bc for lvl <= lowest_lvl"

    nc = mg%box_size_lvl(lvl)

    if (lvl >= mg%first_normal_lvl) then
       dsize = (nc/2)**(NDIM-1)

       mg%buf(:)%i_send = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          do nb = 1, mg_num_neighbors
             if (mg%boxes(id)%neighbors(nb) < mg_no_box) then
                call phi_bc_set_buffer(mg, id, nb)
             end if
          end do
       end do

       mg%buf(:)%i_recv = mg%comm_phibc%n_recv(:, lvl) * dsize
       call sort_and_transfer_buffers(mg, dsize, .false.)
       mg%buf(:)%i_recv = 0
    end if

    do i = 1, size(mg%lvls(lvl-1)%my_parents)
       id = mg%lvls(lvl-1)%my_parents(i)
       do nb = 1, mg_num_neighbors
          if (mg%boxes(id)%neighbors(nb) < mg_no_box) then
             call phi_bc_restrict_onto(mg, id, nc, mg%box_size_lvl(lvl-1), nb)
          end if
       end do
    end do
  end subroutine mg_phi_bc_restrict_lvl

  subroutine phi_bc_set_buffer(mg, id, nb)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nb
    integer                   :: n, hnc, p_id, p_rank
    integer                   :: ix

#if NDIM == 2
    integer                   :: i
    real(dp)                  :: tmp(mg%box_size/2)
#elif NDIM == 3
    integer                   :: i, j
    real(dp)                  :: tmp(mg%box_size/2, mg%box_size/2)
#endif

    p_id   = mg%boxes(id)%parent
    p_rank = mg%boxes(p_id)%rank

    if (p_rank /= mg%my_rank) then
       hnc = mg%box_size/2
       if (mg_neighb_low(nb)) then
          ix = 0
       else
          ix = mg%box_size + 1
       end if

       select case (mg_neighb_dim(nb))
#if NDIM == 2
       case (1)
          do i = 1, hnc
             tmp(i) = 0.5_dp * &
                  sum(mg%boxes(id)%cc(ix, 2*i-1:2*i, mg_irhs))
          end do
       case (2)
          do i = 1, hnc
             tmp(i) = 0.5_dp * &
                  sum(mg%boxes(id)%cc(2*i-1:2*i, ix, mg_irhs))
          end do
#elif NDIM == 3
       case (1)
          do j = 1, hnc
             do i = 1, hnc
                tmp(i, j) = 0.25_dp * &
                     sum(mg%boxes(id)%cc(ix, 2*i-1:2*i, 2*j-1:2*j, mg_irhs))
             end do
          end do
       case (2)
          do j = 1, hnc
             do i = 1, hnc
                tmp(i, j) = 0.25_dp * &
                     sum(mg%boxes(id)%cc(2*i-1:2*i, ix, 2*j-1:2*j, mg_irhs))
             end do
          end do
       case (3)
          do j = 1, hnc
             do i = 1, hnc
                tmp(i, j) = 0.25_dp * &
                     sum(mg%boxes(id)%cc(2*i-1:2*i, 2*j-1:2*j, ix, mg_irhs))
             end do
          end do
#endif
       end select

       ! Buffer
       n = size(tmp)
       i = mg%buf(p_rank)%i_send
       mg%buf(p_rank)%send(i+1:i+n) = pack(tmp, .true.)
       mg%buf(p_rank)%i_send = mg%buf(p_rank)%i_send + n

       ! To later sort the send buffer according to parent order
       i = mg%buf(p_rank)%i_ix
       n = mg_ix_to_ichild(mg%boxes(id)%ix)

       ! Sort by parent, neighbor direction and child number
       mg%buf(p_rank)%ix(i+1) = mg_num_neighbors * mg_num_children * p_id + &
            mg_num_children * nb + n
       mg%buf(p_rank)%i_ix = mg%buf(p_rank)%i_ix + 1
    end if
  end subroutine phi_bc_set_buffer

  subroutine phi_bc_restrict_onto(mg, id, nc, nc_c, nb)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: nc_c
    integer, intent(in)       :: nb
    integer                   :: hnc, dsize, i_c, c_id, n_adj_child
    integer                   :: c_rank, dix(NDIM), ix, ix_c, lvl
#if NDIM == 2
    integer                   :: i
    real(dp)                  :: tmp(nc/2)
#elif NDIM == 3
    integer                   :: i, j
    real(dp)                  :: tmp(nc/2, nc/2)
#endif

    hnc = nc/2
    lvl = mg%boxes(id)%lvl

    if (mg_neighb_low(nb)) then
       ix   = 0
       ix_c = 0
    else
       ix   = nc + 1
       ix_c = nc_c + 1
    end if

    if (lvl < mg%first_normal_lvl) then
       n_adj_child = 1
    else
       n_adj_child = size(mg_child_adj_nb, 1)
    end if

    do i_c = 1, n_adj_child
       if (lvl < mg%first_normal_lvl) then
          ! One child is neighboring all the sides
          c_id = mg%boxes(id)%children(i_c)
       else
          c_id = mg%boxes(id)%children(mg_child_adj_nb(i_c, nb))
       end if

       c_rank = mg%boxes(c_id)%rank
       dix    = mg_get_child_offset(mg, c_id)

       if (c_rank == mg%my_rank) then
          select case (mg_neighb_dim(nb))
#if NDIM == 2
          case (1)
             do i = 1, hnc
                mg%boxes(id)%cc(ix_c, dix(2)+i, mg_irhs) = &
                     0.5_dp * sum(mg%boxes(c_id)%cc(ix, 2*i-1:2*i, mg_irhs))
             end do
          case (2)
             do i = 1, hnc
                mg%boxes(id)%cc(dix(1)+i, ix_c, mg_irhs) = &
                     0.5_dp * sum(mg%boxes(c_id)%cc(2*i-1:2*i, ix, mg_irhs))
             end do
#elif NDIM == 3
          case (1)
             do j = 1, hnc
                do i = 1, hnc
                   mg%boxes(id)%cc(ix_c, dix(2)+i, dix(3)+j, mg_irhs) = 0.25_dp * &
                        sum(mg%boxes(c_id)%cc(ix, 2*i-1:2*i, 2*j-1:2*j, mg_irhs))
                end do
             end do
          case (2)
             do j = 1, hnc
                do i = 1, hnc
                   mg%boxes(id)%cc(dix(1)+i, ix_c, dix(3)+j, mg_irhs) = 0.25_dp * &
                        sum(mg%boxes(c_id)%cc(2*i-1:2*i, ix, 2*j-1:2*j, mg_irhs))
                end do
             end do
          case (3)
             do j = 1, hnc
                do i = 1, hnc
                   mg%boxes(id)%cc(dix(1)+i, dix(2)+j, ix_c, mg_irhs) = 0.25_dp * &
                        sum(mg%boxes(c_id)%cc(2*i-1:2*i, 2*j-1:2*j, ix, mg_irhs))
                end do
             end do
#endif
          end select
       else
          dsize = hnc**(NDIM-1)
          i = mg%buf(c_rank)%i_recv
#if NDIM == 2
          tmp = mg%buf(c_rank)%recv(i+1:i+dsize)
#elif NDIM == 3
          tmp = reshape(mg%buf(c_rank)%recv(i+1:i+dsize), [hnc, hnc])
#endif
          mg%buf(c_rank)%i_recv = mg%buf(c_rank)%i_recv + dsize

          select case (mg_neighb_dim(nb))
#if NDIM == 2
          case (1)
             mg%boxes(id)%cc(ix, dix(2)+1:dix(2)+hnc, mg_irhs) = tmp
          case (2)
             mg%boxes(id)%cc(dix(1)+1:dix(1)+hnc, ix, mg_irhs) = tmp
#elif NDIM == 3
          case (1)
             mg%boxes(id)%cc(ix, dix(2)+1:dix(2)+hnc, &
                  dix(3)+1:dix(3)+hnc, mg_irhs) = tmp
          case (2)
             mg%boxes(id)%cc(dix(1)+1:dix(1)+hnc, ix, &
                  dix(3)+1:dix(3)+hnc, mg_irhs) = tmp
          case (3)
             mg%boxes(id)%cc(dix(1)+1:dix(1)+hnc, &
                  dix(2)+1:dix(2)+hnc, ix, mg_irhs) = tmp
#endif
          end select
       end if
    end do

  end subroutine phi_bc_restrict_onto

  subroutine box_set_gc(box, nb, nc, iv, gc)
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)        :: nb, nc, iv
#if NDIM == 2
    real(dp), intent(in)       :: gc(nc)
#elif NDIM == 3
    real(dp), intent(in)       :: gc(nc, nc)
#endif

    select case (nb)
#if NDIM == 2
    case (mg_neighb_lowx)
       box%cc(0, 1:nc, iv)    = gc
    case (mg_neighb_highx)
       box%cc(nc+1, 1:nc, iv) = gc
    case (mg_neighb_lowy)
       box%cc(1:nc, 0, iv)    = gc
    case (mg_neighb_highy)
       box%cc(1:nc, nc+1, iv) = gc
#elif NDIM == 3
    case (mg_neighb_lowx)
       box%cc(0, 1:nc, 1:nc, iv)    = gc
    case (mg_neighb_highx)
       box%cc(nc+1, 1:nc, 1:nc, iv) = gc
    case (mg_neighb_lowy)
       box%cc(1:nc, 0, 1:nc, iv)    = gc
    case (mg_neighb_highy)
       box%cc(1:nc, nc+1, 1:nc, iv) = gc
    case (mg_neighb_lowz)
       box%cc(1:nc, 1:nc, 0, iv)    = gc
    case (mg_neighb_highz)
       box%cc(1:nc, 1:nc, nc+1, iv) = gc
#endif
    end select
  end subroutine box_set_gc

end module m_phi_bc
