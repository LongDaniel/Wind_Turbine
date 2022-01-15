!-------------------------------------------------------------------------
! MODULE: fft
!
!> FFT transform library based on FFTW.
!
!> @author Anqing Xuan modified by William Xuanting Hao for HOS applications.
!-------------------------------------------------------------------------
module fft_hos

  use hos_param, only : wp, nxhos, nyhos, ncpu_hos, myid, numprocs, ierr
  use, intrinsic :: iso_c_binding

  implicit none

  public

  integer(C_INT), parameter :: FFTW_MEASURE = 0
  integer(C_INT), parameter :: FFTW_DESTROY_INPUT = 1
  integer(C_INT), parameter :: FFTW_UNALIGNED = 2
  integer(C_INT), parameter :: FFTW_EXHAUSTIVE = 8
  integer(C_INT), parameter :: FFTW_PRESERVE_INPUT = 16
  integer(C_INT), parameter :: FFTW_PATIENT = 32
  integer(C_INT), parameter :: FFTW_ESTIMATE = 64
  integer(C_INT), parameter :: FFTW_ESTIMATE_PATIENT = 128

  interface
     type(C_PTR) function fftw_plan_many_dft_r2c(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,flags) &
                   bind(C, name='fftw_plan_many_dft_r2c')
       import
       integer(C_INT), value :: rank
       integer(C_INT), dimension(*), intent(in) :: n
       integer(C_INT), value :: howmany
       real(C_DOUBLE), dimension(*), intent(out) :: in
       integer(C_INT), dimension(*), intent(in) :: inembed
       integer(C_INT), value :: istride
       integer(C_INT), value :: idist
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
       integer(C_INT), dimension(*), intent(in) :: onembed
       integer(C_INT), value :: ostride
       integer(C_INT), value :: odist
       integer(C_INT), value :: flags
     end function fftw_plan_many_dft_r2c

     type(C_PTR) function fftw_plan_many_dft_c2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,flags) &
                   bind(C, name='fftw_plan_many_dft_c2r')
       import
       integer(C_INT), value :: rank
       integer(C_INT), dimension(*), intent(in) :: n
       integer(C_INT), value :: howmany
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
       integer(C_INT), dimension(*), intent(in) :: inembed
       integer(C_INT), value :: istride
       integer(C_INT), value :: idist
       real(C_DOUBLE), dimension(*), intent(out) :: out
       integer(C_INT), dimension(*), intent(in) :: onembed
       integer(C_INT), value :: ostride
       integer(C_INT), value :: odist
       integer(C_INT), value :: flags
     end function fftw_plan_many_dft_c2r

     subroutine fftw_destroy_plan(p) bind(C, name='fftw_destroy_plan')
       import
       type(C_PTR), value :: p
     end subroutine fftw_destroy_plan

     subroutine my_fftw_execute_dft_r2c(p,in,out) bind(C, name='fftw_execute_dft_r2c')
       import
       type(C_PTR), value :: p
       real(C_DOUBLE), dimension(*), intent(in) :: in
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
     end subroutine my_fftw_execute_dft_r2c

     subroutine my_fftw_execute_dft_c2r(p,in,out) bind(C, name='fftw_execute_dft_c2r')
       import
       type(C_PTR), value :: p
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
       real(C_DOUBLE), dimension(*), intent(out) :: out
     end subroutine my_fftw_execute_dft_c2r
  end interface

  ! FFTW plan flag: FFTW_MEASURE or FFTW_ESTIMATE, see FFTW3 documentation
  integer(C_INT), parameter :: flag_r2c = ior(FFTW_ESTIMATE, FFTW_PRESERVE_INPUT)
  integer(C_INT), parameter :: flag_r2c_d = ior(FFTW_ESTIMATE, FFTW_DESTROY_INPUT)
  integer(C_INT), parameter :: flag_c2r = ior(FFTW_ESTIMATE, FFTW_DESTROY_INPUT)

  ! FFTW plans
  type(C_PTR), save :: plan_r2c_x_hos, plan_c2r_x_hos, plan_r2c_y_hos, plan_c2r_y_hos

  ! FFT temporary array
  real(wp), save, allocatable, dimension(:,:), target :: buf_x_hos, buf_y_hos
  complex(wp), save, contiguous, dimension(:,:), pointer :: p_buf_x_hos, p_buf_y_hos
  real(wp), save, allocatable, dimension(:,:), public :: bufb_hos


  public :: fft_init_hos, fft_finalize_hos
  public :: fft_for_hos, fft_bac_hos
  public :: fft_for_x_hos, fft_bac_x_hos
  public :: fft_for_y_hos, fft_bac_y_hos
  public :: fft_for_xy_hos, fft_bac_xy_hos
  public :: transpose_2d
  
  interface fft_for_x_hos
     module procedure fft_for_x_1_hos, fft_for_x_2_hos
     module procedure fft_for_x_3_hos, fft_for_x_4_hos
  end interface

  interface fft_bac_x_hos
     module procedure fft_bac_x_1_hos, fft_bac_x_2_hos
     module procedure fft_bac_x_3_hos, fft_bac_x_4_hos
  end interface

  interface fft_for_y_hos
     module procedure fft_for_y_1_hos, fft_for_y_2_hos
     module procedure fft_for_y_3_hos, fft_for_y_4_hos
  end interface

  interface fft_bac_y_hos
     module procedure fft_bac_y_1_hos, fft_bac_y_2_hos
     module procedure fft_bac_y_3_hos, fft_bac_y_4_hos
  end interface

  interface fft_for_xy_hos
     module procedure fft_for_xy_1_hos, fft_for_xy_2_hos
     module procedure fft_for_xy_3_hos, fft_for_xy_4_hos
  end interface

  interface fft_bac_xy_hos
     module procedure fft_bac_xy_1_hos, fft_bac_xy_2_hos
     module procedure fft_bac_xy_3_hos, fft_bac_xy_4_hos
  end interface

  interface transpose_2d
     module procedure old_2dtrans
  end interface

contains

  !========================================================
  !> @create plan
  !========================================================
  subroutine fft_init_hos
    
    implicit none
    
    real(wp), allocatable, dimension(:,:) :: phy_x_hos, phy_y_hos

    integer xsz(2), ysz(2)

    xsz(1) = nxhos
    xsz(2) = nyhos / ncpu_hos
    ysz(1) = nyhos
    ysz(2) = nxhos / ncpu_hos

    allocate(phy_x_hos(nxhos,nyhos/ncpu_hos))
    allocate(phy_y_hos(nyhos,nxhos/ncpu_hos))
    allocate(buf_x_hos(nxhos+2,nyhos/ncpu_hos))
    allocate(buf_y_hos(nyhos+2,nxhos/ncpu_hos))
    allocate(bufb_hos(nyhos,nxhos/ncpu_hos))

    call c_f_pointer(c_loc(buf_x_hos), p_buf_x_hos, [size(buf_x_hos,1)/2, size(buf_x_hos,2)])
    call c_f_pointer(c_loc(buf_y_hos), p_buf_y_hos, [size(buf_y_hos,1)/2, size(buf_y_hos,2)])

    plan_r2c_x_hos = fftw_plan_many_dft_r2c(1, xsz(1), xsz(2), &
         phy_x_hos, [xsz(1)], 1, xsz(1), &
         p_buf_x_hos, [xsz(1)/2+1], 1, xsz(1)/2+1, &
         flag_r2c)
    
    plan_c2r_x_hos = fftw_plan_many_dft_c2r(1, xsz(1), xsz(2), &
         p_buf_x_hos, [xsz(1)/2+1], 1, xsz(1)/2+1, &
         phy_x_hos, [xsz(1)],  1, xsz(1), &
         flag_c2r)

    plan_r2c_y_hos = fftw_plan_many_dft_r2c(1, ysz(1), ysz(2), &
         phy_y_hos, [ysz(1)], 1, ysz(1), &
         p_buf_y_hos, [ysz(1)/2+1], 1, ysz(1)/2+1, &
         flag_r2c)

    plan_c2r_y_hos = fftw_plan_many_dft_c2r(1, ysz(1), ysz(2), &
         p_buf_y_hos, [ysz(1)/2+1], 1, ysz(1)/2+1, &
         phy_y_hos, [ysz(1)],  1, ysz(1), &
         flag_c2r)

  end subroutine fft_init_hos
  
  !========================================================
  !> @destroy plan
  !========================================================
  subroutine fft_finalize_hos
    
    implicit none

    call dfftw_destroy_plan(plan_r2c_x_hos)
    call dfftw_destroy_plan(plan_c2r_x_hos)
    call dfftw_destroy_plan(plan_r2c_y_hos)
    call dfftw_destroy_plan(plan_c2r_y_hos)       

    deallocate(buf_x_hos)
    deallocate(buf_y_hos)
    deallocate(bufb_hos)
  end subroutine fft_finalize_hos
  
  !========================================================
  !> @fft 
  !========================================================
   
  !---------------------------------------------------------
  
  subroutine fft_for_x_1_hos(inout)
    
    implicit none
    
    real(wp), intent(inout), contiguous, dimension(:,:)::inout  

    if (size(inout,1) == nxhos) then    
       call my_fftw_execute_dft_r2c(plan_r2c_x_hos,inout,p_buf_x_hos)
       inout(:,:) = buf_x_hos(1:nxhos,1:nyhos/ncpu_hos) / nxhos
    else if (size(inout,1) == nyhos) then
       call my_fftw_execute_dft_r2c(plan_r2c_y_hos,inout,p_buf_y_hos)
       inout(:,:) = buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) / nyhos
    else
       print *,"This subroutine is only for FFT of array types (nxhos,:) or (nyhos,:)"
       print *,"To correct this bug, use 'call fft_for_x_hos(a,size(a,1),size(a,2))' instead!"
       stop
    end if


  end subroutine fft_for_x_1_hos
 
  !---------------------------------------------------------

  !---------------------------------------------------------

  subroutine fft_for_x_2_hos(in,out)

    implicit none

    real(wp), intent(in), contiguous, dimension(:,:):: in
    real(wp), intent(out), contiguous, dimension(:,:) :: out

    if (size(in,1) == nxhos) then
       call my_fftw_execute_dft_r2c(plan_r2c_x_hos,in,p_buf_x_hos)
       out(:,:) = buf_x_hos(1:nxhos,1:nyhos/ncpu_hos) / nxhos
    else if (size(in,1) == nyhos) then
       call my_fftw_execute_dft_r2c(plan_r2c_y_hos,in,p_buf_y_hos)
       out(:,:) = buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) / nyhos
    else
       print *,"This subroutine is only for FFT of array types (nxhos,:) or (nyhos,:)"
       print *,"To correct this bug, use 'call fft_for_x_hos(a,size(a,1),size(a,2))' instead!"
       stop
    end if

  end subroutine fft_for_x_2_hos

  !---------------------------------------------------------

  subroutine fft_for_x_3_hos(inout,nx,ny)

    implicit none

    integer, intent(in) :: nx, ny
    real(wp), intent(inout), dimension(:,:) :: inout

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    integer*8 :: tmp_plan

    integer :: i, k

    allocate(tmp_buf(nx),tmp_p_buf(nx/2+1))
    call dfftw_plan_dft_r2c_1d(tmp_plan,nx,tmp_buf,tmp_p_buf,FFTW_ESTIMATE)

    do k = 1, ny
       tmp_buf = inout(1:nx,k)
       inout(:,k) = 0.0
       call dfftw_execute_dft_r2c(tmp_plan, tmp_buf, tmp_p_buf)
       do i = 1, nx / 2
          inout(i*2-1,k) = real(tmp_p_buf(i)) / nx
          inout(i*2,k) = aimag(tmp_p_buf(i)) / nx
       end do
    end do

    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf)

  end subroutine fft_for_x_3_hos

  !---------------------------------------------------------

  subroutine fft_for_x_4_hos(in,out,nx,ny)

    implicit none

    integer, intent(in) :: nx, ny
    real(wp), intent(in), dimension(:,:)::in
    real(wp), intent(out), dimension(:,:)::out

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    integer*8 :: tmp_plan

    integer :: i, k

    allocate(tmp_buf(nx),tmp_p_buf(nx/2+1))
    call dfftw_plan_dft_r2c_1d(tmp_plan,nx,tmp_buf,tmp_p_buf,FFTW_ESTIMATE)

    do k = 1, ny
       tmp_buf = in(1:nx,k)
       call dfftw_execute_dft_r2c(tmp_plan, tmp_buf, tmp_p_buf)
       do i = 1, nx / 2
          out(i*2-1,k) = real(tmp_p_buf(i)) / nx
          out(i*2,k) = aimag(tmp_p_buf(i)) / nx
       end do
    end do

    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf)

  end subroutine fft_for_x_4_hos

  !---------------------------------------------------------
  !---------------------------------------------------------

  subroutine fft_for_hos(in,out,n)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in), dimension(:)::in
    real(wp), intent(out), dimension(:)::out

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    integer*8 :: tmp_plan

    integer :: i

    allocate(tmp_buf(n),tmp_p_buf(n/2+1))
    call dfftw_plan_dft_r2c_1d(tmp_plan,n,tmp_buf,tmp_p_buf,FFTW_ESTIMATE)

    tmp_buf = in(1:n)
    call dfftw_execute_dft_r2c(tmp_plan, tmp_buf, tmp_p_buf)
    do i = 1, n / 2
       out(i*2-1) = real(tmp_p_buf(i)) / n
       out(i*2) = aimag(tmp_p_buf(i)) / n
    end do

    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf)

  end subroutine fft_for_hos
  
  !---------------------------------------------------------

  subroutine fft_bac_x_1_hos(inout)

    implicit none

    real(wp), intent(inout), dimension(:,:)::inout

    if (size(inout,1) == nxhos) then
       buf_x_hos = 0
       buf_x_hos(1:nxhos,1:nyhos/ncpu_hos) = inout
       call my_fftw_execute_dft_c2r(plan_c2r_x_hos,p_buf_x_hos,inout)
    else if (size(inout,1) == nyhos) then
       buf_y_hos = 0
       buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) = inout
       call my_fftw_execute_dft_c2r(plan_c2r_y_hos,p_buf_y_hos,inout)
    else
       print *,"This subroutine is only for inverse FFT of array types (nxhos,:) or (nyhos,:)"
       print *,"To correct this bug, use 'call fft_bac_x_hos(a,size(a,1),size(a,2))' instead!"
       stop
    end if

  end subroutine fft_bac_x_1_hos
  !---------------------------------------------------------

  subroutine fft_bac_x_2_hos(in, out)

    implicit none

    real(wp), intent(in), dimension(:,:) :: in
    real(wp), intent(out), dimension(:,:) :: out

    if (size(in,1) == nxhos) then
       buf_x_hos = 0
       buf_x_hos(1:nxhos,1:nyhos/ncpu_hos) = in
       call my_fftw_execute_dft_c2r(plan_c2r_x_hos,p_buf_x_hos,out)
    else if (size(in,1) == nyhos) then
       buf_y_hos = 0
       buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) = in
       call my_fftw_execute_dft_c2r(plan_c2r_y_hos,p_buf_y_hos,out)
    else
       print *,"This subroutine is only for inverse FFT of array types (nxhos,:) or (nyhos,:)"
       print *,"To correct this bug, use 'call fft_bac_x_hos(a,size(a,1),size(a,2))' instead!"
       stop
    end if

  end subroutine fft_bac_x_2_hos


  !---------------------------------------------------------

  subroutine fft_bac_x_3_hos(inout,nx,ny)

    implicit none

    integer, intent(in) :: nx, ny
    real(wp), intent(inout), dimension(:,:)::inout

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    integer*8 :: tmp_plan

    integer :: i, k

    allocate(tmp_buf(nx),tmp_p_buf(nx/2+1))
    call dfftw_plan_dft_c2r_1d(tmp_plan,nx,tmp_p_buf,tmp_buf,FFTW_ESTIMATE)

    do k = 1, ny
       tmp_p_buf = 0.0
       do i = 1, nx / 2
          tmp_p_buf(i) = cmplx(inout(i*2-1,k),inout(i*2,k),wp)
       end do
       tmp_p_buf(nx / 2+1) = 0.0
       call dfftw_execute_dft_c2r(tmp_plan, tmp_p_buf, tmp_buf)
       inout(:,k) = tmp_buf
    end do

    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf,tmp_p_buf)

  end subroutine fft_bac_x_3_hos
  !---------------------------------------------------------

  subroutine fft_bac_x_4_hos(in, out, nx, ny)

    implicit none

    integer, intent(in) :: nx, ny
    real(wp), intent(in), dimension(:,:) :: in
    real(wp), intent(inout), dimension(:,:) :: out

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    integer*8 :: tmp_plan

    integer :: i, k

    allocate(tmp_buf(nx),tmp_p_buf(nx/2+1))
    call dfftw_plan_dft_c2r_1d(tmp_plan,nx,tmp_p_buf,tmp_buf,FFTW_ESTIMATE)

    do k = 1, ny
       tmp_p_buf = 0.0
       do i = 1, nx / 2
          tmp_p_buf(i) = cmplx(in(i*2-1,k),in(i*2,k),wp)
       end do
       tmp_p_buf(nx/2+1) = 0.0
       call dfftw_execute_dft_c2r(tmp_plan, tmp_p_buf, tmp_buf)
       out(:,k) = tmp_buf
    end do

    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf)

  end subroutine fft_bac_x_4_hos

  !---------------------------------------------------------

  !---------------------------------------------------------

  subroutine fft_bac_hos(in, out, n)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in), dimension(:) :: in
    real(wp), intent(inout), dimension(:) :: out

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    integer*8 :: tmp_plan

    integer :: i

    allocate(tmp_buf(n),tmp_p_buf(n/2+1))
    call dfftw_plan_dft_c2r_1d(tmp_plan,n,tmp_p_buf,tmp_buf,FFTW_ESTIMATE)

    tmp_p_buf = 0.0
    do i = 1, n / 2
       tmp_p_buf(i) = cmplx(in(i*2-1),in(i*2),wp)
    end do
    tmp_p_buf(n/2+1) = 0.0
    call dfftw_execute_dft_c2r(tmp_plan, tmp_p_buf, tmp_buf)
    out(:) = tmp_buf

    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf)

  end subroutine fft_bac_hos

  !---------------------------------------------------------

  subroutine fft_for_y_1_hos(inout)

    implicit none

    real(wp), intent(inout), dimension(:,:)::inout

    if (size(inout,2) == nyhos/ncpu_hos) then
       call transpose_2d(inout,bufb_hos,nxhos,nyhos/ncpu_hos)
       call my_fftw_execute_dft_r2c(plan_r2c_y_hos,bufb_hos,p_buf_y_hos)
       bufb_hos(:,:) = buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) / nyhos
       call transpose_2d(bufb_hos,inout,nyhos,nxhos/ncpu_hos)
    else
       print *,"This subroutine fft_for_y_hos(inout) is only for FFT of array types (nxhos,:) in y direction!"
       print *,"To correct this bug, use 'call fft_for_y_hos(a,size(a,1),size(a,2))' instead!"
       stop
    end if

  end subroutine fft_for_y_1_hos

  !---------------------------------------------------------

  subroutine fft_for_y_2_hos(in, out)

    implicit none

    real(wp), intent(in), dimension(:,:) :: in
    real(wp), intent(out), dimension(:,:) :: out

    if (size(in,2) == nyhos/ncpu_hos) then
       call transpose_2d(in,bufb_hos,nxhos,nyhos/ncpu_hos)
       call my_fftw_execute_dft_r2c(plan_r2c_y_hos,bufb_hos,p_buf_y_hos)
       bufb_hos(:,:) = buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) / nyhos
       call transpose_2d(bufb_hos,out,nyhos,nxhos/ncpu_hos)
    else
       print *,"This subroutine fft_for_y_hos(in, out) is only for FFT of array types (nxhos,:) in y direction!"
       print *,"To correct this bug, use 'call fft_for_y_hos(a,size(a,1),size(a,2))' instead!"
       stop
    end if

  end subroutine fft_for_y_2_hos


  !---------------------------------------------------------

  subroutine fft_for_y_3_hos(inout,nx,ny)

    implicit none

    integer, intent(in) :: nx, ny
    real(wp), intent(inout), dimension(:,:) :: inout

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    real(wp), allocatable, dimension(:,:) :: at,b,bt
    integer*8 :: tmp_plan

    integer :: i, k

    allocate(tmp_buf(ny*ncpu_hos),tmp_p_buf(ny*ncpu_hos/2+1))
    allocate(at(nx,ny),b(ny*ncpu_hos,nx/ncpu_hos),bt(ny*ncpu_hos,nx/ncpu_hos))
    call dfftw_plan_dft_r2c_1d(tmp_plan,ny*ncpu_hos,tmp_buf,tmp_p_buf,FFTW_ESTIMATE)

    call transpose_2d(inout,b,nx,ny)

    do k = 1, nx / ncpu_hos
       tmp_buf = b(:,k)
       call dfftw_execute_dft_r2c(tmp_plan, tmp_buf, tmp_p_buf)
       do i = 1, ny * ncpu_hos / 2
          b(i*2-1,k) = real(tmp_p_buf(i)) / ny / ncpu_hos
          b(i*2,k) = aimag(tmp_p_buf(i)) / ny / ncpu_hos
       end do
    end do

    call transpose_2d(b,inout,ny*ncpu_hos,nx/ncpu_hos)
    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf,tmp_p_buf,at,b,bt)

  end subroutine fft_for_y_3_hos

  !---------------------------------------------------------
  subroutine fft_for_y_4_hos(in, out,nx,ny)

    implicit none

    integer, intent(in) :: nx, ny
    real(wp), intent(in), dimension(:,:) :: in
    real(wp), intent(out), dimension(:,:) :: out

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    real(wp), allocatable, dimension(:,:) :: at,b,bt
    integer*8 :: tmp_plan

    integer :: i, k

    allocate(tmp_buf(ny*ncpu_hos),tmp_p_buf(ny*ncpu_hos/2+1))
    allocate(at(nx,ny),b(ny*ncpu_hos,nx/ncpu_hos),bt(ny*ncpu_hos,nx/ncpu_hos))
    call dfftw_plan_dft_r2c_1d(tmp_plan,ny*ncpu_hos,tmp_buf,tmp_p_buf,FFTW_ESTIMATE)

    call transpose_2d(in,b,nx,ny)

    do k = 1, nx / ncpu_hos
       tmp_buf = b(:,k)
       call dfftw_execute_dft_r2c(tmp_plan, tmp_buf, tmp_p_buf)
       do i = 1, ny * ncpu_hos / 2
          b(i*2-1,k) = real(tmp_p_buf(i)) / ny / ncpu_hos
          b(i*2,k) = aimag(tmp_p_buf(i)) / ny / ncpu_hos
       end do
    end do

    call transpose_2d(b,out,ny*ncpu_hos,nx/ncpu_hos)
    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf,tmp_p_buf,at,b,bt)

  end subroutine fft_for_y_4_hos

  !---------------------------------------------------------

  subroutine fft_bac_y_1_hos(inout)

    implicit none

    real(wp), intent(inout), dimension(:,:)::inout

    if (size(inout,2) == nyhos/ncpu_hos) then
       call transpose_2d(inout,bufb_hos,nxhos,nyhos/ncpu_hos)
       buf_y_hos = 0
       buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) = bufb_hos      
       call my_fftw_execute_dft_c2r(plan_c2r_y_hos,p_buf_y_hos,bufb_hos)
       call transpose_2d(bufb_hos,inout,nyhos,nxhos/ncpu_hos)
    else
       print *,"This subroutine fft_bac_y_hos(inout) is only for FFT of array types (nxhos,:) in y direction!"
       print *,"To correct this bug, use 'call fft_bac_y_hos(a,size(a,1),size(a,2))' instead!"
       stop
    end if

  end subroutine fft_bac_y_1_hos

  !---------------------------------------------------------

  subroutine fft_bac_y_2_hos(in, out)

    implicit none

    real(wp), intent(in), dimension(:,:) :: in
    real(wp), intent(out), dimension(:,:) :: out

    if (size(in,2) == nyhos/ncpu_hos) then
       call transpose_2d(in,bufb_hos,nxhos,nyhos/ncpu_hos)
       buf_y_hos = 0
       buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) = bufb_hos
       call my_fftw_execute_dft_c2r(plan_c2r_y_hos,p_buf_y_hos,bufb_hos)
       call transpose_2d(bufb_hos,out,nyhos,nxhos/ncpu_hos)
    else
       print *,"This subroutine fft_bac_y_hos(in,out) is only for FFT of array types (nxhos,:) in y direction!"
       print *,"To correct this bug, use 'call fft_bac_y_hos(a,size(a,1),size(a,2))' instead!"
       stop
    end if

  end subroutine fft_bac_y_2_hos


  !---------------------------------------------------------

  subroutine fft_bac_y_3_hos(inout,nx,ny)

    implicit none

    integer, intent(in) :: nx, ny
    real(wp), intent(inout), dimension(:,:) :: inout

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    real(wp), allocatable, dimension(:,:) :: at,b,bt
    integer*8 :: tmp_plan

    integer :: i, k

    allocate(tmp_buf(ny*ncpu_hos),tmp_p_buf(ny*ncpu_hos/2+1))
    allocate(at(nx,ny),b(ny*ncpu_hos,nx/ncpu_hos),bt(ny*ncpu_hos,nx/ncpu_hos))
    call dfftw_plan_dft_c2r_1d(tmp_plan,ny*ncpu_hos,tmp_p_buf,tmp_buf,FFTW_ESTIMATE)

    call transpose_2d(inout,b,nx,ny)

    do k = 1, nx / ncpu_hos
       tmp_p_buf = 0.0
       do i = 1, ny * ncpu_hos/2
          tmp_p_buf(i) = cmplx(b(i*2-1,k),b(i*2,k),wp)
       end do
       call dfftw_execute_dft_c2r(tmp_plan, tmp_p_buf, tmp_buf)
       b(:,k) = tmp_buf
    end do

    call transpose_2d(b,inout,ny*ncpu_hos,nx/ncpu_hos)
    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf,tmp_p_buf,at,b,bt)

  end subroutine fft_bac_y_3_hos

  !---------------------------------------------------------

  subroutine fft_bac_y_4_hos(in,out,nx,ny)

    implicit none

    integer, intent(in) :: nx, ny
    real(wp), intent(inout), dimension(:,:) :: in
    real(wp), intent(out), dimension(:,:) :: out

    real(wp), allocatable, dimension(:) :: tmp_buf
    complex(wp), allocatable, dimension(:) :: tmp_p_buf
    real(wp), allocatable, dimension(:,:) :: b,bt
    integer*8 :: tmp_plan

    integer :: i, k

    allocate(tmp_buf(ny*ncpu_hos),tmp_p_buf(ny*ncpu_hos/2+1))
    allocate(b(ny*ncpu_hos,nx/ncpu_hos),bt(ny*ncpu_hos,nx/ncpu_hos))
    call dfftw_plan_dft_c2r_1d(tmp_plan,ny*ncpu_hos,tmp_p_buf,tmp_buf,FFTW_ESTIMATE)

    call transpose_2d(in,b,nx,ny)

    do k = 1, nx / ncpu_hos
       tmp_p_buf = 0.0
       do i = 1, ny * ncpu_hos/2
          tmp_p_buf(i) = cmplx(b(i*2-1,k),b(i*2,k),wp)
       end do
       call dfftw_execute_dft_c2r(tmp_plan, tmp_p_buf, tmp_buf)
       b(:,k) = tmp_buf
    end do

    call transpose_2d(b,out,ny*ncpu_hos,nx/ncpu_hos)
    call dfftw_destroy_plan(tmp_plan)
    deallocate(tmp_buf,tmp_p_buf,b,bt)

  end subroutine fft_bac_y_4_hos

  !---------------------------------------------------------

  subroutine fft_for_xy_1_hos(inout)
    
    implicit none

    real(wp), intent(inout), dimension(:,:) :: inout
    
    call my_fftw_execute_dft_r2c(plan_r2c_x_hos,inout,p_buf_x_hos)
    inout(:,:) = buf_x_hos(1:nxhos,1:nyhos/ncpu_hos) / nxhos
    call transpose_2d(inout,bufb_hos,nxhos,nyhos/ncpu_hos)
    call my_fftw_execute_dft_r2c(plan_r2c_y_hos,bufb_hos,p_buf_y_hos)
    bufb_hos(:,:) = buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) / nyhos
    call transpose_2d(bufb_hos,inout,nyhos,nxhos/ncpu_hos)

  end subroutine fft_for_xy_1_hos

  !---------------------------------------------------------

  subroutine fft_for_xy_2_hos(in, out)

    implicit none

    real(wp), intent(in), dimension(:,:) :: in
    real(wp), intent(out), dimension(:,:)::out

    call my_fftw_execute_dft_r2c(plan_r2c_x_hos,in,p_buf_x_hos)
    out(:,:) = buf_x_hos(1:nxhos,1:nyhos/ncpu_hos) / nxhos
    call transpose_2d(out,bufb_hos,nxhos,nyhos/ncpu_hos)
    call my_fftw_execute_dft_r2c(plan_r2c_y_hos,bufb_hos,p_buf_y_hos)
    bufb_hos(:,:) = buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) / nyhos
    call transpose_2d(bufb_hos,out,nyhos,nxhos/ncpu_hos)

  end subroutine fft_for_xy_2_hos

  !---------------------------------------------------------

  subroutine fft_for_xy_3_hos(inout, nx, ny)

    implicit none

    real(wp), intent(inout), dimension(:,:) :: inout
    integer, intent(in) :: nx,ny

    real(wp), allocatable, dimension(:,:) :: tmp

    allocate(tmp(ny,nx))

    call fft_for_x_hos(inout,nx,ny)
    tmp = transpose(inout)
    call fft_for_x_hos(tmp,ny,nx)
    inout = transpose(tmp)

    deallocate(tmp)

  end subroutine fft_for_xy_3_hos

  !---------------------------------------------------------

  subroutine fft_for_xy_4_hos(in, out, nx, ny)

    implicit none

    real(wp), intent(in), dimension(:,:) :: in
    real(wp), intent(out), dimension(:,:)::out
    integer, intent(in) :: nx,ny

    real(wp), allocatable, dimension(:,:) :: tmpin, tmpout

    allocate(tmpin(nx,ny))
    allocate(tmpout(ny,nx))

    tmpin = in
    call fft_for_x_hos(tmpin,nx,ny)
    tmpout = transpose(tmpin)
    call fft_for_x_hos(tmpout,ny,nx)
    out = transpose(tmpout)
    
    deallocate(tmpin,tmpout)

  end subroutine fft_for_xy_4_hos


  !---------------------------------------------------------

  subroutine fft_bac_xy_1_hos(inout)
    
    implicit none
    
    real(wp), intent(inout), dimension(:,:)::inout
    
    call transpose_2d(inout,bufb_hos,nxhos,nyhos/ncpu_hos)
    buf_y_hos = 0
    buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) = bufb_hos
    call my_fftw_execute_dft_c2r(plan_c2r_y_hos,p_buf_y_hos,bufb_hos)
    call transpose_2d(bufb_hos,inout,nyhos,nxhos/ncpu_hos)
    buf_x_hos = 0
    buf_x_hos(1:nxhos,1:nyhos/ncpu_hos) = inout

    call my_fftw_execute_dft_c2r(plan_c2r_x_hos,p_buf_x_hos,inout)


  end subroutine fft_bac_xy_1_hos

  !---------------------------------------------------------

  subroutine fft_bac_xy_2_hos(in, out)

    implicit none

    real(wp), intent(in), dimension(:,:) :: in
    real(wp), intent(out), dimension(:,:)::out

    call transpose_2d(in,bufb_hos,nxhos,nyhos/ncpu_hos)
    buf_y_hos = 0
    buf_y_hos(1:nyhos,1:nxhos/ncpu_hos) = bufb_hos
    call my_fftw_execute_dft_c2r(plan_c2r_y_hos,p_buf_y_hos,bufb_hos)
    call transpose_2d(bufb_hos,out,nyhos,nxhos/ncpu_hos)
    buf_x_hos = 0
    buf_x_hos(1:nxhos,1:nyhos/ncpu_hos) = out

    call my_fftw_execute_dft_c2r(plan_c2r_x_hos,p_buf_x_hos,out)

  end subroutine fft_bac_xy_2_hos

  !---------------------------------------------------------

  subroutine fft_bac_xy_3_hos(inout, nx, ny)

    implicit none

    real(wp), intent(inout), dimension(:,:) :: inout
    integer, intent(in) :: nx,ny

    real(wp), allocatable, dimension(:,:) :: tmp

    allocate(tmp(ny,nx))

    tmp = transpose(inout)
    call fft_bac_x_hos(tmp,ny,nx)
    inout = transpose(tmp)
    call fft_bac_x_hos(inout,nx,ny)

    deallocate(tmp)

  end subroutine fft_bac_xy_3_hos

  !---------------------------------------------------------

  subroutine fft_bac_xy_4_hos(in, out, nx, ny)

    implicit none

    real(wp), intent(in), dimension(:,:) :: in
    real(wp), intent(out), dimension(:,:)::out
    integer, intent(in) :: nx,ny

    real(wp), allocatable, dimension(:,:) :: tmpin, tmpout

    allocate(tmpin(ny,nx))
    allocate(tmpout(nx,ny))
   
    tmpin = transpose(in)
    call fft_bac_x_hos(tmpin,ny,nx)
    tmpout = transpose(tmpin)
    call fft_bac_x_hos(tmpout,nx,ny)
    out = transpose(tmpout)

    deallocate(tmpin,tmpout)

  end subroutine fft_bac_xy_4_hos


  !---------------------------------------------------------


  subroutine old_2dtrans(a,b,ndx,ndy)

!    use hos, only: wp,myid, numprocs, ierr
    !     by di yang 2005                                                                                             
    !     input array a(ndx,ndy). ndx and ndy are the numbers of points                                               
    !     in x and y direction in one processor. ncpu_hos is the number of                                                
    !     of processors in computation. output b(ndy*ncpu_hos,ndx/ncpu_hos) which                                             
    !     is the transpose array of a(ndx,ndy) 
    use MPI
    implicit none
!    include 'mpif.h'                                                                                                 
    integer, intent(in) :: ndx,ndy
    real(wp), intent(in) :: a(ndx,ndy)
    real(wp), intent(inout) :: b(ndy*ncpu_hos,ndx/ncpu_hos)
    real(wp), allocatable, dimension(:,:) :: tmpa, at

    integer i,j,n
    integer nbb,ndyglob
    integer m1,n1,k1,m2,n2,k2,ii

    ndyglob=ndy*ncpu_hos
    nbb=ndy*ndx/ncpu_hos

    allocate(tmpa(ndx,ndy), at(ndx,ndy))
    tmpa = a
    do j=1,ndy
       do i=1,ndx
          at(i,j)=tmpa(i,j)
       enddo
    enddo

    do j=1,ndy
       do i=1,ndx
          k1=j+(i-1)*ndy
          k2=i+(j-1)*ndx
          call getij_hos(m1,n1,k1,ndx)
          call getij_hos(m2,n2,k2,ndx)
          tmpa(m1,n1)=at(m2,n2)
          at(m2,n2)=0.
       enddo
    enddo

    call mpi_alltoall(tmpa,nbb,mpi_double_precision,at,nbb,mpi_double_precision,mpi_comm_world,ierr)

    do j=1,ndx/ncpu_hos
       ii=0
       do n=1,numprocs
          do i=1,ndy
             ii=ii+1
             k1=ii+(j-1)*ndyglob
             k2=i+(n-1)*nbb+(j-1)*ndy
             call getij_hos(m1,n1,k1,ndyglob)
             call getij_hos(m2,n2,k2,ndx)
             b(m1,n1)=at(m2,n2)
          enddo
       enddo
    enddo

    deallocate(tmpa, at)
    return
  end subroutine old_2dtrans

  subroutine getij_hos(i,j,k,cols)

    !     input the offset of the element to the start location                                                       
    !     in the old array in the physical memory. input the                                                          
    !     number of columns in the new array. get the position                                                        
    !     (i,j) of the element in the new array.                                                                      

    implicit none
    integer i,j,k,cols
    i=mod(k-1,cols)+1
    j=(k-1)/cols+1
    return
  end subroutine getij_hos

  
end module fft_hos
