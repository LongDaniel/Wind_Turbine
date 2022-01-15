!-----------------------------------------------------------
! MODULE: decomp
!
!> 2D decomposition module.
!
!> @author Anqing Xuan
!-----------------------------------------------------------
module decomp

  use MPI

  implicit none

  private

  integer, parameter, public :: wp = kind(0.0d0)
  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION

  integer, save, public :: nx_global, ny_global, nz_global  ! global size

  integer, save, public :: myid, myid1, myid2  ! local MPI rank 
  integer, save, public :: nproc, nproc1, nproc2       ! total number of processors
  logical, save, public :: istop, isbot
  integer, save, public :: topid

  !> added by plyu:
  !>   ids: whether to consider Discontinuity Smoothing for pdfx_out3 and pdfy_x_out3
  !>   idsd: whether to consider Discontinuity Smoothing for dealiasxy3
  !>   idsp: whether to consider Discontinuity Smoothing for fft_x_r2c in poisson equation
  integer, save, public :: ids, idsd, idsp

  ! 2D Cartesian topology 
  integer, save, public :: MPI_COMM_2D_CART
  integer, save, public :: MPI_COMM_2D_ROW, MPI_COMM_2D_COL

  ! neighboring blocks (for ghost-cell) 1=lower, 2=higher
  integer, save, dimension(2) :: neighbour

  ! decomposition information
  integer, save, allocatable, dimension(:), public :: xdist, ydist, zdist, ydist2
  ! staring/ending (global) index and size of data on local processor
  integer, save, dimension(3), public :: xst, xend, xsz  ! x-pencil
  integer, save, dimension(3), public :: yst, yend, ysz  ! y-pencil
  integer, save, dimension(3), public :: zst, zend, zsz  ! z-pencil

  ! Buffers for MPI_ALLTOALL transpose
  integer, save :: blk_size1 = 0, buf_size = 0
  real(wp), allocatable, dimension(:), target :: work1, work2

  ! public subroutines
  public :: decomp_init, decomp_finalize, decomp_abort
  public :: transpose_xy, transpose_yx
  public :: alloc_x, alloc_y
  public :: update_ghost

  interface transpose_xy
     module procedure transpose_xy2, transpose_xy3
  end interface transpose_xy

  interface transpose_yx
     module procedure transpose_yx2, transpose_yx3
  end interface transpose_yx

  interface update_ghost
     module procedure update_ghost_inplace1, update_ghost_inplace3
     module procedure update_ghost_outplace1, update_ghost_outplace3
  end interface update_ghost

contains

  !---------------------------------------------------------
  !> @brief Initialize the 2D domain decomposition.
  !
  !> @param[in] nx,ny,nz size of the domain
  !> @param[in] p_row,p_col the processor grid
  !---------------------------------------------------------
  subroutine decomp_init(nx,ny,nz,p_row,p_col)

    implicit none

    integer, intent(IN) :: nx,ny,nz,p_row,p_col
    
    integer :: ierror
    integer, allocatable, dimension(:) :: st, en
    integer :: coord(2)
    
    nx_global = nx
    ny_global = ny
    nz_global = nz

    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

    if (p_row==0 .and. p_col==0) then
       ! determine the best 2D processor grid
       call best_2d_grid(nproc, nproc1, nproc2)
    else
       ! check the number of processors
       if (nproc /= p_row*p_col) then
          call decomp_abort(1, &
               'nproc /= p_row*p_col')
       end if
       if (nx<p_row .or. ny<p_row .or. nz<p_col) then
          call decomp_abort(6, &
               'Make sure min(nx,ny) < p_row or min(nz) < p_col')
       end if
       nproc1 = p_row
       nproc2 = p_col
    end if
    
    ! Create 2D Catersian topology
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,[nproc2,nproc1],[.false.,.false.], &
         .false.,MPI_COMM_2D_CART, ierror)

    call MPI_CART_COORDS(MPI_COMM_2D_CART,myid,2,coord,ierror)
    myid1 = coord(2)
    myid2 = coord(1)
    
    ! derive communicators defining sub-groups for ALLTOALL
    call MPI_CART_SUB(MPI_COMM_2D_CART,(/.false.,.true./), &
         MPI_COMM_2D_COL,ierror)
    call MPI_CART_SUB(MPI_COMM_2D_CART,(/.true.,.false./), &
         MPI_COMM_2D_ROW,ierror)

    allocate(xdist(0:nproc1-1),ydist(0:nproc1-1),ydist2(0:nproc2-1),zdist(0:nproc2-1))
    ! generate global and local distribution information
    allocate(st(0:nproc1-1),en(0:nproc1-1))
    call distribute(nx,nproc1,st,en,xdist)
    yst(1)=1; yend(1)=ny; ysz(1)=ny
    yst(2)=st(myid1); yend(2)=en(myid1); ysz(2)=xdist(myid1)

    call distribute(ny,nproc1,st,en,ydist)
    xst(1)=1; xend(1)=nx; xsz(1)=nx
    xst(2)=st(myid1); xend(2)=en(myid1); xsz(2)=ydist(myid1)

    deallocate(st,en)
    allocate(st(0:nproc2-1),en(0:nproc2-1))
    call distribute(nz,nproc2,st,en,zdist)
    xst(3)=st(myid2); xend(3)=en(myid2); xsz(3)=zdist(myid2)
    yst(3)=st(myid2); yend(3)=en(myid2); ysz(3)=zdist(myid2)
    
    call distribute(ny,nproc2,st,en,ydist2)
    zst(1)=st(myid2); zend(1)=en(myid2); zsz(1)=ydist2(myid2)
    zst(2)=yst(2); zend(2)=yend(2); zsz(2)=xdist(myid1)
    zst(3)=1; zend(3)=nz; zsz(3)=nz

    ! buffer size for transpose
    blk_size1 = maxval(xdist) * maxval(ydist)
    buf_size = blk_size1*xsz(3)*nproc1

    ! allocate memory for the MPI_ALLTOALL buffers
    allocate(work1(buf_size), STAT=ierror)
    allocate(work2(buf_size), STAT=ierror)
    if (ierror /= 0) then
       call decomp_abort(2, 'Out of memory allocating transpose buffer.')
    end if

    call MPI_CART_SHIFT(MPI_COMM_2D_ROW, 0, 1, &
         neighbour(1), neighbour(2), ierror) ! lower & upper
    istop = .false.; isbot = .false.
    if (myid2 == 0) isbot = .true.
    if (myid2 == nproc2-1) istop = .true.
    topid = nproc2-1

  end subroutine decomp_init

  !---------------------------------------------------------
  !> @brief Clean up the 2D decomposition information.
  !---------------------------------------------------------
  subroutine decomp_finalize

    implicit none
    
    deallocate(xdist,ydist,ydist2,zdist)
    deallocate(work1, work2)
    
    return
  end subroutine decomp_finalize

  !---------------------------------------------------------
  !> @brief Return distribution of grid points.
  !
  !> Distributes the grid points in one dimension and return
  !! the start/end indicies and sizes on each processor.
  !
  !> @param[in] data1 size to be partitioned
  !> @param[in] proc  number of processors in that dimension
  !> @param[out] st,en array of start/end indicies
  !> @param[out] sz    array of local sizes
  !> @param[in] algorithm (optional) 1 -- higher rank get more points,
  !! 2 (default) -- lower rank get more points
  !---------------------------------------------------------
  subroutine distribute(data1,proc,st,en,sz,opt_algorithm)
  
    implicit none
    integer :: data1,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
    integer :: i,size1,nl,nu,algorithm
    integer, optional :: opt_algorithm
  
    if (present(opt_algorithm)) then
       algorithm = opt_algorithm
    else
       algorithm = 2
    end if

    select case (algorithm)
       case (1)
       ! higher rank get more grid points
       size1 = data1/proc
       nu = data1 - size1 * proc
       nl = proc - nu
       st(0) = 1
       sz(0) = size1
       en(0) = size1
       do i=1,nl-1
         st(i) = st(i-1) + size1
         sz(i) = size1
         en(i) = en(i-1) + size1
       end do
       size1 = size1 + 1
       do i=nl,proc-1
         st(i) = en(i-1) + 1
         sz(i) = size1
         en(i) = en(i-1) + size1
       end do
       en(proc-1)= data1 
       sz(proc-1)= data1-st(proc-1)+1

       case default
       ! lower rank get more grid points
       size1 = data1/proc
       nl = data1 - size1 * proc
       nu = proc - nl
       st(proc-1) = data1-size1+1
       sz(proc-1) = size1
       en(proc-1) = data1
       do i=proc-2,proc-nu,-1
         st(i) = st(i+1) - size1
         sz(i) = size1
         en(i) = en(i+1) - size1
       end do
       size1 = size1 + 1
       do i=proc-nu-1,0,-1
         en(i) = st(i+1) - 1
         sz(i) = size1
         st(i) = st(i+1) - size1
       end do
       st(0)= 1
       sz(0)= en(0)
    end select
  
    return
  end subroutine distribute


  !---------------------------------------------------------
  !> @brief Transpose array from x-pencil to y-pencil
  !
  !> @param[in]  src input array, x-pencil
  !> @param[out] dst output array, y-pencil
  !> @param[in]  nz  size of the array in 3rd dimension
  !---------------------------------------------------------
  subroutine transpose_xy3(src, dst)

    implicit none
    
    real(wp), dimension(:,:,:), intent(IN)  :: src
    real(wp), dimension(:,:,:), intent(OUT) :: dst

    integer :: j,k,m,pos,i1,i2,kmax
    integer :: ierror
    real(wp), dimension(:,:,:), pointer, contiguous :: buf

    kmax = min(size(src,3), size(dst,3))

    ! optimize for single core
    if (nproc1 == 1) then
       do k=1, kmax
          dst(:,:,k) = transpose(src(:,:,k))
       end do
       return
    end if
!print *, "trxy3, 1, ", myid;
!print *, "myid=", myid, ", blk_size1=", blk_size1, ", kmax=",kmax, &
!  ", nproc1=",nproc1, ", xsz(2)=", xsz(2), "xdist=", xdist(0:(nproc1-1))
!call MPI_Barrier(MPI_COMM_WORLD, ierror); 
    ! rearrange source array as send buffer
    buf(1:blk_size1,1:kmax,0:nproc1-1) => work1
    i2 = 0
    do m=0,nproc1-1
       i1 = i2+1
       i2 = i1+xdist(m)-1
       do k=1,kmax
          pos = 1
          do j=1,xsz(2)
             buf(pos:pos+xdist(m)-1,k,m) = src(i1:i2,j,k)
             pos = pos+xdist(m)
          end do
       end do
    end do
!print *, "trxy3, 2, ", myid; call MPI_Barrier(MPI_COMM_WORLD, ierror); 
!print *, myid, mpi_comm_2d_cart, mpi_comm_2d_col, mpi_comm_2d_row
    ! transpose using MPI_ALLTOALL
    call MPI_Alltoall(work1, blk_size1*kmax, real_type, &
                      work2, blk_size1*kmax, real_type, &
                      MPI_COMM_2D_COL, ierror)
!    print *, myid, ierror
!print *, "trxy3, 3, ", myid; call MPI_Barrier(MPI_COMM_WORLD, ierror); 

    ! rearrange receive buffer
    buf(1:blk_size1,1:kmax,0:nproc1-1) => work2
    i2 = 0
    do m=0,nproc1-1
       i1 = i2+1
       i2 = i1+ydist(m)-1
       do k=1,kmax
          pos = 1
          do j=i1,i2
             dst(j,1:ysz(2),k) = buf(pos:pos+ysz(2)-1,k,m)
             pos = pos+ysz(2)
          end do
       end do
    end do
    
  end subroutine transpose_xy3

  subroutine transpose_xy2(src, dst)

    implicit none
    
    real(wp), dimension(:,:), intent(IN)  :: src
    real(wp), dimension(:,:), intent(OUT) :: dst

    integer :: j,m,pos,i1,i2
    integer :: ierror
    real(wp), dimension(:,:), pointer, contiguous :: buf

    if (nproc1 == 1) then
       dst = transpose(src)
       return
    end if

    ! rearrange source array as send buffer
    buf(1:blk_size1,0:nproc1-1) => work1
    i2 = 0
    do m=0,nproc1-1
       i1 = i2+1
       i2 = i1+xdist(m)-1
       pos = 1
       do j=1,xsz(2)
          buf(pos:pos+xdist(m)-1,m) = src(i1:i2,j)
          pos = pos+xdist(m)
       end do
    end do

    ! transpose using MPI_ALLTOALL
    call MPI_Alltoall(work1, blk_size1, real_type, &
                      work2, blk_size1, real_type, &
                      MPI_COMM_2D_COL, ierror)

    ! rearrange receive buffer
    buf(1:blk_size1,0:nproc1-1) => work2
    i2 = 0
    do m=0,nproc1-1
       i1 = i2+1
       i2 = i1+ydist(m)-1
       pos = 1
       do j=i1,i2
          dst(j,1:ysz(2)) = buf(pos:pos+ysz(2)-1,m)
          pos = pos+ysz(2)
       end do
    end do
    
  end subroutine transpose_xy2

  !---------------------------------------------------------
  !> @brief Transpose array from y-pencil to x-pencil
  !
  !> @param[in]  src input array, y-pencil
  !> @param[out] dst output array, x-pencil
  !> @param[in]  nz  size of the array in 3rd dimension
  !---------------------------------------------------------
  subroutine transpose_yx3(src, dst)

    implicit none
    
    real(wp), dimension(:,:,:), intent(IN)  :: src
    real(wp), dimension(:,:,:), intent(OUT) :: dst

    integer :: j,k,m,pos,i1,i2,kmax
    integer :: ierror
    real(wp), dimension(:,:,:), pointer, contiguous :: buf

    kmax = min(size(src,3), size(dst,3))
    if (nproc1 == 1) then
       do k=1, kmax
          dst(:,:,k) = transpose(src(:,:,k))
       end do
       return
    end if

    ! rearrange source array as send buffer
    buf(1:blk_size1,1:kmax,0:nproc1-1) => work1
    i2 = 0
    do m=0,nproc1-1
       i1 = i2+1
       i2 = i1+ydist(m)-1
       do k=1,kmax
          pos = 1
          do j=i1,i2
             buf(pos:pos+ysz(2)-1,k,m) = src(j,1:ysz(2),k) 
             pos = pos+ysz(2)
          end do
       end do
    end do

    ! transpose using MPI_ALLTOALL
    call MPI_Alltoall(work1, blk_size1*kmax, real_type, &
                      work2, blk_size1*kmax, real_type, &
                      MPI_COMM_2D_COL, ierror)

    ! rearrange receive buffer
    buf(1:blk_size1,1:kmax,0:nproc1-1) => work2
    i2 = 0
    do m=0,nproc1-1
       i1 = i2+1
       i2 = i1+xdist(m)-1
       do k=1,kmax
          pos = 1
          do j=1,xsz(2)
             dst(i1:i2,j,k) = buf(pos:pos+xdist(m)-1,k,m)
             pos = pos+xdist(m)
          end do
       end do
    end do
    
  end subroutine transpose_yx3

  subroutine transpose_yx2(src, dst)

    implicit none
    
    real(wp), dimension(:,:), intent(IN)  :: src
    real(wp), dimension(:,:), intent(OUT) :: dst

    integer :: j,m,pos,i1,i2
    integer :: ierror
    real(wp), dimension(:,:), pointer, contiguous :: buf

    if (nproc1 == 1) then
       dst = transpose(src)
       return
    end if

    ! rearrange source array as send buffer
    buf(1:blk_size1,0:nproc1-1) => work1
    i2 = 0
    do m=0,nproc1-1
       i1 = i2+1
       i2 = i1+ydist(m)-1
       pos = 1
       do j=i1,i2
          buf(pos:pos+ysz(2)-1,m) = src(j,1:ysz(2)) 
          pos = pos+ysz(2)
       end do
    end do

    ! transpose using MPI_ALLTOALL
    call MPI_Alltoall(work1, blk_size1, real_type, &
                      work2, blk_size1, real_type, &
                      MPI_COMM_2D_COL, ierror)

    ! rearrange receive buffer
    buf(1:blk_size1,0:nproc1-1) => work2
    i2 = 0
    do m=0,nproc1-1
       i1 = i2+1
       i2 = i1+xdist(m)-1
       pos = 1
       do j=1,xsz(2)
          dst(i1:i2,j) = buf(pos:pos+xdist(m)-1,m)
          pos = pos+xdist(m)
       end do
    end do
    
  end subroutine transpose_yx2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Auto-tuning algorithm to select the best 2D processor grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine best_2d_grid(iproc, best_p_row, best_p_col)

    implicit none

    integer, intent(IN) :: iproc
    integer, intent(OUT) :: best_p_row, best_p_col

    integer, allocatable, dimension(:) :: factors
    double precision :: t1, t2, best_time
    integer :: nfact, i, row, col, ierror
    integer :: best_row, best_col

    real(wp), allocatable, dimension(:,:,:) :: u1, u2

    if (myid==0) write(*,*) 'In auto-tuning mode......'

    best_time = huge(t1)
    best_row = -1
    best_col = -1
    
    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if (myid==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

    do i=1, nfact

       row = factors(i)
       col = iproc / row

       ! enforce the limitation of 2D decomposition
       if (min(nx_global,ny_global)>=row .and. &
            min(ny_global,nz_global)>=col) then

          ! 2D Catersian topology
          call decomp_init(nx_global,ny_global,nz_global,row,col)

          ! arrays for X,Y and Z-pencils
          allocate(u1(xsz(1),xsz(2),xsz(3)))
          allocate(u2(ysz(1),ysz(2),ysz(3)))

          ! timing the transposition routines
          t1 = MPI_Wtime()
          call transpose_xy(u1,u2)
          call transpose_yx(u2,u1)
          t2 = MPI_Wtime() - t1

          deallocate(u1,u2)

          call decomp_finalize()

          call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                   MPI_COMM_WORLD,ierror)
          t1 = t1 / dble(nproc)

          if (myid==0) then
             write(*,*) 'processor grid', row, ' by ', col, ' time=', t1
          end if

          if (best_time > t1) then
             best_time = t1
             best_row = row
             best_col = col
          end if

       end if
       
    end do ! loop through processer grid

    deallocate(factors)

    if (best_row/=-1) then
       best_p_row = best_row
       best_p_col = best_col
       if (myid==0) then
          write(*,*) 'the best processor grid is probably ', &
               best_p_row, ' by ', best_p_col
       end if
    else
       call decomp_abort(9, &
            'The processor-grid auto-tuning failed. ')
    end if

    return
  end subroutine best_2d_grid

#include "factor.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Ghost cell subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "ghost.f90"


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Error handling
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_abort(errorcode, msg)

    implicit none

    integer, intent(IN) :: errorcode
    character(len=*), intent(IN), optional :: msg

    integer :: ierror
    
    if (present(msg)) then
       if (myid==0) then
          write(*,'(A6,I2.2,A,A)') 'ERROR(',errorcode,'): ',msg
       end if
    end if
    call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)

    return
  end subroutine decomp_abort


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routines to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "alloc.f90"
    
  
end module decomp

