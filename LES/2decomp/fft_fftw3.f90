!-----------------------------------------------------------
! MODULE: fft
!
!> FFT transform library based on 2D decomposition module and FFTW.
!
!> @author Anqing Xuan
!-----------------------------------------------------------
module fft
  
  use decomp  ! 2D decomposition module
  use, intrinsic :: iso_c_binding

  implicit none

  private

  include "fftw3_interface.f90"

  ! FFTW plan flag: FFTW_MEASURE or FFTW_ESTIMATE, see FFTW3 documentation
  integer(C_INT), parameter :: flag_r2c_out = ior(FFTW_ESTIMATE, FFTW_PRESERVE_INPUT)
  integer(C_INT), parameter :: flag_r2c_in = ior(FFTW_ESTIMATE, FFTW_DESTROY_INPUT)
  integer(C_INT), parameter :: flag_c2r = ior(FFTW_ESTIMATE, FFTW_DESTROY_INPUT)

  ! FFTW plans
  type(C_PTR), save :: plan_r2c_x_out, plan_r2c_x_in, plan_r2c_y_out, plan_r2c_y_in
  type(C_PTR), save :: plan_c2r_x, plan_c2r_y

  ! FFT temporary array
  real(wp), save, allocatable, dimension(:,:), target :: buf_x, buf_y
  complex(wp), save, contiguous, dimension(:,:), pointer :: p_buf_x, p_buf_y

  public :: fft_init, fft_finalize
  public :: fft_r2c_x, fft_c2r_x
  public :: fft_r2c_y, fft_c2r_y
  public :: fft_r2c_xy, fft_c2r_xy

  !> added by plyu
  real(wp), dimension(:), allocatable, public, save :: co_s32_x, co_s32_y
  public :: fft_r2c_xy_out_s32

  interface fft_r2c_x
     module procedure fft_r2c_x_in, fft_r2c_x_in2
     module procedure fft_r2c_x_out, fft_r2c_x_out2
     module procedure fft_r2c_x_out_cplx, fft_r2c_x_out2_cplx
  end interface fft_r2c_x

  interface fft_r2c_y
     module procedure fft_r2c_y_in, fft_r2c_y_in2
     module procedure fft_r2c_y_out, fft_r2c_y_out2
     module procedure fft_r2c_y_out_cplx, fft_r2c_y_out2_cplx
  end interface fft_r2c_y

  interface fft_r2c_xy
     module procedure fft_r2c_xy_out, fft_r2c_xy_out2
  end interface fft_r2c_xy

  interface fft_c2r_x
     module procedure fft_c2r_x_in, fft_c2r_x_in2
     module procedure fft_c2r_x_out, fft_c2r_x_out2
     module procedure fft_c2r_x_out_cplx, fft_c2r_x_out2_cplx
  end interface fft_c2r_x

  interface fft_c2r_y
     module procedure fft_c2r_y_in, fft_c2r_y_in2
     module procedure fft_c2r_y_out, fft_c2r_y_out2
     module procedure fft_c2r_y_out_cplx, fft_c2r_y_out2_cplx
  end interface fft_c2r_y

  interface fft_c2r_xy
     module procedure fft_c2r_xy_out, fft_c2r_xy_out2
  end interface fft_c2r_xy

contains

  !---------------------------------------------------------
  !> @brief Initialize the 2D decomposition FFTW plan.
  !---------------------------------------------------------
  subroutine fft_init

    implicit none

    real(wp), allocatable, dimension(:,:) :: phy_x, phy_y

    !> added by plyu
    integer :: i

    ! Check FFT dimensions
    if (mod(xsz(1),2) /= 0 .or. mod(ysz(1),2) /= 0) then
       call decomp_abort(20, 'FFT size must be even.')
    end if

    allocate(phy_x(xsz(1),xsz(2)))
    allocate(phy_y(ysz(1),ysz(2)))
    allocate(buf_x(xsz(1)+2,xsz(2)))
    allocate(buf_y(ysz(1)+2,ysz(2)))
    call c_f_pointer(c_loc(buf_x), p_buf_x, [size(buf_x,1)/2, size(buf_x,2)])
    call c_f_pointer(c_loc(buf_y), p_buf_y, [size(buf_y,1)/2, size(buf_y,2)])

    plan_r2c_x_out = fftw_plan_many_dft_r2c(1, xsz(1), xsz(2),              &
                                        phy_x, [xsz(1)],     1, xsz(1),     &
                                      p_buf_x, [xsz(1)/2+1], 1, xsz(1)/2+1, &
                                        flag_r2c_out)
    plan_r2c_x_in = fftw_plan_many_dft_r2c(1, xsz(1), xsz(2),               &
                                        phy_x, [xsz(1)],     1, xsz(1),     &
                                      p_buf_x, [xsz(1)/2+1], 1, xsz(1)/2+1, &
                                        flag_r2c_in)
    plan_c2r_x = fftw_plan_many_dft_c2r(1, xsz(1), xsz(2),                  &
                                      p_buf_x, [xsz(1)/2+1], 1, xsz(1)/2+1, &
                                        phy_x, [xsz(1)],     1, xsz(1),     &
                                        flag_c2r)
    plan_r2c_y_out = fftw_plan_many_dft_r2c(1, ysz(1), ysz(2),              &
                                        phy_y, [ysz(1)],     1, ysz(1),     &
                                      p_buf_y, [ysz(1)/2+1], 1, ysz(1)/2+1, &
                                        flag_r2c_out)
    plan_r2c_y_in = fftw_plan_many_dft_r2c(1, ysz(1), ysz(2),               &
                                        phy_y, [ysz(1)],     1, ysz(1),     &
                                      p_buf_y, [ysz(1)/2+1], 1, ysz(1)/2+1, &
                                        flag_r2c_in)
    plan_c2r_y = fftw_plan_many_dft_c2r(1, ysz(1), ysz(2),                  &
                                      p_buf_y, [ysz(1)/2+1], 1, ysz(1)/2+1, &
                                        phy_y, [ysz(1)],     1, ysz(1),     &
                                        flag_c2r)

    !> added by plyu
    allocate(co_s32_x(xsz(1)/2), co_s32_y(ysz(1)/2))
    do i = 1, xsz(1)/2
      co_s32_x(i) = exp(-36.0*(real(i,wp)/real(xsz(1)/2,wp))**36)
    enddo
    do i = 1, ysz(1)/2
      co_s32_y(i) = exp(-36.0*(real(i,wp)/real(ysz(1)/2,wp))**36)
    enddo
 
  end subroutine fft_init

  !---------------------------------------------------------
  !> @brief Clean up the FFTW plan.
  !---------------------------------------------------------
  subroutine fft_finalize

    implicit none

    call fftw_destroy_plan(plan_r2c_x_in)
    call fftw_destroy_plan(plan_r2c_x_out)
    call fftw_destroy_plan(plan_c2r_x)
    call fftw_destroy_plan(plan_r2c_y_in)
    call fftw_destroy_plan(plan_r2c_y_out)
    call fftw_destroy_plan(plan_c2r_y)

    deallocate(buf_x)
    deallocate(buf_y)

  end subroutine fft_finalize


  !---------------------------------------------------------
  !> @brief Multiple r2c transform in x direction, unnormalized.
  !
  !> Multiple r2c transform in x direction. The output is 
  !! unnormalized and only keeps modes from 0 to N/2-1, ie, 
  !! the last mode N/2 is discarded. The input data is kept
  !! after transformation.
  !
  !> @param[in]  input   data to be transformed, x-pencil
  !> @param[out] output  transformed data, x-pencil, unnormalized
  !> @param[in]  nz      number of grids to be transformed in z direction
  !---------------------------------------------------------
  subroutine fft_r2c_x_out(input, output)

    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(IN)  :: input
    real(wp), dimension(:,:,:), contiguous, intent(OUT) :: output

    integer :: k

    do k=1, size(input,3)
       call my_fftw_execute_dft_r2c(plan_r2c_x_out,input(:,:,k),p_buf_x)
       output(1:xsz(1),:,k) = buf_x(1:xsz(1),:)
    end do
    output(xsz(1)+1:,:,:) = 0

  end subroutine fft_r2c_x_out

  subroutine fft_r2c_x_out_cplx(input, output)

    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(IN)  :: input
    complex(wp), dimension(:,:,:), contiguous, intent(OUT) :: output

    integer :: k

    do k=1, size(input,3)
       call my_fftw_execute_dft_r2c(plan_r2c_x_out,input(:,:,k),output(:,:,k))
    end do
    output(xsz(1)/2+1:,:,:) = 0

  end subroutine fft_r2c_x_out_cplx

  subroutine fft_r2c_x_in(inout)

    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(INOUT) :: inout

    integer :: k

    do k=1, size(inout,3)
       call my_fftw_execute_dft_r2c(plan_r2c_x_in,inout(:,:,k),p_buf_x)
       inout(:,:,k) = buf_x(1:xsz(1),1:xsz(2))
    end do

  end subroutine fft_r2c_x_in

  subroutine fft_r2c_x_out2(input, output)

    implicit none

    real(wp), dimension(:,:), contiguous, intent(IN)  :: input
    real(wp), dimension(:,:), contiguous, intent(OUT) :: output

    call my_fftw_execute_dft_r2c(plan_r2c_x_out,input,p_buf_x)
    output = buf_x(1:xsz(1),1:xsz(2))
    output(xsz(1)+1:,:) = 0

  end subroutine fft_r2c_x_out2

  subroutine fft_r2c_x_out2_cplx(input, output)

    implicit none

    real(wp), dimension(:,:), contiguous, intent(IN)  :: input
    complex(wp), dimension(:,:), contiguous, intent(OUT) :: output

    call my_fftw_execute_dft_r2c(plan_r2c_x_out,input,output)
    output(xsz(1)/2+1:,:) = 0

  end subroutine fft_r2c_x_out2_cplx

  subroutine fft_r2c_x_in2(inout)

    implicit none

    real(wp), dimension(:,:), contiguous, intent(INOUT) :: inout

    call my_fftw_execute_dft_r2c(plan_r2c_x_in,inout,p_buf_x)
    inout = buf_x(1:xsz(1),1:xsz(2))

  end subroutine fft_r2c_x_in2

  !---------------------------------------------------------
  !> @brief Multiple r2c transform in y direction, unnormalized.
  !
  !> Multiple r2c transform in y direction. The output is 
  !! unnormalized and only keeps modes from 0 to N/2-1, ie, 
  !! the last mode N/2 is discarded. The input data is kept
  !! after transformation.
  !
  !> @param[in]  input   data to be transformed, y-pencil
  !> @param[out] output  transformed data, y-pencil, unnormalized
  !> @param[in]  nz      number of grids to be transformed in z direction
  !---------------------------------------------------------
  subroutine fft_r2c_y_out(input, output)

    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(IN)  :: input
    real(wp), dimension(:,:,:), contiguous, intent(OUT) :: output

    integer :: k

    do k=1, size(input,3)
       call my_fftw_execute_dft_r2c(plan_r2c_y_out,input(:,:,k),p_buf_y)
       output(1:ysz(1),1:ysz(2),k) = buf_y(1:ysz(1),1:ysz(2))
    end do
    output(ysz(1)+1:,:,:) = 0

  end subroutine fft_r2c_y_out

  subroutine fft_r2c_y_out_cplx(input, output)

    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(IN)  :: input
    complex(wp), dimension(:,:,:), contiguous, intent(OUT) :: output

    integer :: k

    do k=1, size(input,3)
       call my_fftw_execute_dft_r2c(plan_r2c_y_out,input(:,:,k),output(:,:,k))
    end do
    output(ysz(1)/2+1:,:,:) = 0

  end subroutine fft_r2c_y_out_cplx

  subroutine fft_r2c_y_in(inout)

    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(INOUT) :: inout

    integer :: k

    do k=1, size(inout,3)
       call my_fftw_execute_dft_r2c(plan_r2c_y_in,inout(:,:,k),p_buf_y)
       inout(:,:,k) = buf_y(1:ysz(1),1:ysz(2))
    end do

  end subroutine fft_r2c_y_in

  subroutine fft_r2c_y_out2(input, output)

    implicit none

    real(wp), dimension(:,:), contiguous, intent(IN)  :: input
    real(wp), dimension(:,:), contiguous, intent(OUT) :: output

    call my_fftw_execute_dft_r2c(plan_r2c_y_out,input,p_buf_y)
    output = buf_y(1:ysz(1),1:ysz(2))
    output(ysz(1)+1:,:) = 0

  end subroutine fft_r2c_y_out2

  subroutine fft_r2c_y_out2_cplx(input, output)

    implicit none

    real(wp), dimension(:,:), contiguous, intent(IN)  :: input
    complex(wp), dimension(:,:), contiguous, intent(OUT) :: output

    call my_fftw_execute_dft_r2c(plan_r2c_y_out,input,output)
    output(ysz(1)/2+1:,:) = 0

  end subroutine fft_r2c_y_out2_cplx

  subroutine fft_r2c_y_in2(inout)

    implicit none

    real(wp), dimension(:,:), contiguous, intent(INOUT) :: inout

    call my_fftw_execute_dft_r2c(plan_r2c_y_in,inout,p_buf_y)
    inout = buf_y(1:ysz(1),1:ysz(2))

  end subroutine fft_r2c_y_in2

  !---------------------------------------------------------
  !> @brief Multiple c2r transform in x direction.
  !
  !> Multiple c2r transform in x direction. The input should have
  !! modes from 0 to N/2-1, and the last mode N/2 is assumed to be
  !! zero. The input data is kept after transformation.
  !
  !> @param[in]  input   data to be transformed, x-pencil
  !> @param[out] output  transformed data, x-pencil
  !> @param[in]  nz      number of grids to be transformed in z direction
  !---------------------------------------------------------
  subroutine fft_c2r_x_out(input, output)
    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(IN)  :: input
    real(wp), dimension(:,:,:), contiguous, intent(OUT) :: output

    integer :: k

    do k=1, size(input,3)
       buf_x(1:xsz(1),1:xsz(2)) = input(1:xsz(1),1:xsz(2),k)
       buf_x(xsz(1)+1:,:) = 0
       call my_fftw_execute_dft_c2r(plan_c2r_x,p_buf_x,output(:,:,k))
    end do

  end subroutine fft_c2r_x_out

  subroutine fft_c2r_x_out_cplx(input, output)
    implicit none

    complex(wp), dimension(:,:,:), contiguous, intent(INOUT)  :: input
    real(wp), dimension(:,:,:), contiguous, intent(OUT) :: output

    integer :: k

    do k=1, size(input,3)
       call my_fftw_execute_dft_c2r(plan_c2r_x,input(:,:,k),output(:,:,k))
    end do

  end subroutine fft_c2r_x_out_cplx

  subroutine fft_c2r_x_in(inout)
    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(INOUT) :: inout

    integer :: k

    do k=1, size(inout,3)
       buf_x(1:xsz(1),1:xsz(2)) = inout(1:xsz(1),1:xsz(2),k)
       buf_x(xsz(1)+1:,:) = 0
       call my_fftw_execute_dft_c2r(plan_c2r_x,p_buf_x,inout(:,:,k))
    end do

  end subroutine fft_c2r_x_in

  subroutine fft_c2r_x_out2(input, output)

    implicit none

    real(wp), dimension(:,:), contiguous, intent(IN)  :: input
    real(wp), dimension(:,:), contiguous, intent(OUT) :: output

    buf_x(1:xsz(1),1:xsz(2)) = input(1:xsz(1),1:xsz(2))
    buf_x(xsz(1)+1:,:) = 0
    call my_fftw_execute_dft_c2r(plan_c2r_x,p_buf_x,output)

  end subroutine fft_c2r_x_out2

  subroutine fft_c2r_x_out2_cplx(input, output)

    implicit none

    complex(wp), dimension(:,:), contiguous, intent(INOUT)  :: input
    real(wp), dimension(:,:), contiguous, intent(OUT) :: output

    call my_fftw_execute_dft_c2r(plan_c2r_x,input,output)

  end subroutine fft_c2r_x_out2_cplx

  subroutine fft_c2r_x_in2(inout)

    implicit none

    real(wp), dimension(:,:), contiguous, intent(INOUT) :: inout

    buf_x(1:xsz(1),1:xsz(2)) = inout(1:xsz(1),1:xsz(2))
    buf_x(xsz(1)+1:,:) = 0
    call my_fftw_execute_dft_c2r(plan_c2r_x,p_buf_x,inout)

  end subroutine fft_c2r_x_in2

  !---------------------------------------------------------
  !> @brief Multiple c2r transform in y direction.
  !
  !> Multiple c2r transform in y direction. The input should have
  !! modes from 0 to N/2-1, and the last mode N/2 is assumed to be
  !! zero. The input data is kept after transformation.
  !
  !> @param[in]  input   data to be transformed, y-pencil
  !> @param[out] output  transformed data, y-pencil
  !> @param[in]  nz      number of grids to be transformed in z direction
  !---------------------------------------------------------
  subroutine fft_c2r_y_out(input, output)
    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(IN)  :: input
    real(wp), dimension(:,:,:), contiguous, intent(OUT) :: output

    integer :: k

    do k=1, size(input,3)
       buf_y(1:ysz(1),1:ysz(2)) = input(1:ysz(1),1:ysz(2),k)
       buf_y(ysz(1)+1:,:) = 0
       call my_fftw_execute_dft_c2r(plan_c2r_y,p_buf_y,output(:,:,k))
    end do

  end subroutine fft_c2r_y_out

  subroutine fft_c2r_y_out_cplx(input, output)
    implicit none

    complex(wp), dimension(:,:,:), contiguous, intent(INOUT)  :: input
    real(wp), dimension(:,:,:), contiguous, intent(OUT) :: output

    integer :: k

    do k=1, size(input,3)
       call my_fftw_execute_dft_c2r(plan_c2r_y,input(:,:,k),output(:,:,k))
    end do

  end subroutine fft_c2r_y_out_cplx

  subroutine fft_c2r_y_in(inout)
    implicit none

    real(wp), dimension(:,:,:), contiguous, intent(INOUT) :: inout

    integer :: k

    do k=1, size(inout,3)
       buf_y(1:ysz(1),1:ysz(2)) = inout(1:ysz(1),1:ysz(2),k)
       buf_y(ysz(1)+1:,:) = 0
       call my_fftw_execute_dft_c2r(plan_c2r_y,p_buf_y,inout(:,:,k))
    end do

  end subroutine fft_c2r_y_in

  subroutine fft_c2r_y_out2(input, output)
    implicit none

    real(wp), dimension(:,:), contiguous, intent(IN)  :: input
    real(wp), dimension(:,:), contiguous, intent(OUT) :: output

    buf_y(1:ysz(1),1:ysz(2)) = input(1:ysz(1),1:ysz(2))
    buf_y(ysz(1)+1:,:) = 0
    call my_fftw_execute_dft_c2r(plan_c2r_y,p_buf_y,output(:,:))

  end subroutine fft_c2r_y_out2

  subroutine fft_c2r_y_out2_cplx(input, output)
    implicit none

    complex(wp), dimension(:,:), contiguous, intent(INOUT)  :: input
    real(wp), dimension(:,:), contiguous, intent(OUT) :: output

    call my_fftw_execute_dft_c2r(plan_c2r_y,input(:,:),output(:,:))

  end subroutine fft_c2r_y_out2_cplx

  subroutine fft_c2r_y_in2(inout)

    implicit none

    real(wp), dimension(:,:), contiguous, intent(INOUT) :: inout

    buf_y(1:ysz(1),1:ysz(2)) = inout(1:ysz(1),1:ysz(2))
    buf_y(ysz(1)+1:,:) = 0
    call my_fftw_execute_dft_c2r(plan_c2r_y,p_buf_y,inout)

  end subroutine fft_c2r_y_in2

  !---------------------------------------------------------
  !> @brief Multiple r2c transform in x and then y direction, unnormalized.
  !
  !> Multiple r2c transform in x and y direction. The output is 
  !! unnormalized and only keeps modes from 0 to N/2-1, ie, 
  !! the last mode N/2 is discarded. The input data is <b>overwritten</b>
  !! after transformation.
  !
  !> @param[in]  input   data to be transformed, x-pencil
  !> @param[out] output  transformed data, y-pencil, unnormalized
  !> @param[in]  nz      number of grids to be transformed in z direction
  !---------------------------------------------------------
  subroutine fft_r2c_xy_out(input, output)

     implicit none

     real(wp), dimension(:,:,:), contiguous, intent(INOUT) :: input
     real(wp), dimension(:,:,:), contiguous, intent(OUT)   :: output

     call fft_r2c_x(input)
     call transpose_xy(input, output)
     call fft_r2c_y(output)

  end subroutine fft_r2c_xy_out
  
  subroutine fft_r2c_xy_out_s32(input, output)
     implicit none

     real(wp), dimension(:,:,:), contiguous, intent(INOUT) :: input
     real(wp), dimension(:,:,:), contiguous, intent(OUT)   :: output

     integer :: i

     call fft_r2c_x(input)
     
     do i = 1, xsz(1)/2
       input((2*i-1):(2*i),:,:) = co_s32_x(i) * input((2*i-1):(2*i),:,:)
     enddo

     call transpose_xy(input, output)
     call fft_r2c_y(output)

     do i = 1, ysz(1)/2
       output((2*i-1):(2*i),:,:) = co_s32_y(i) * output((2*i-1):(2*i),:,:)
     enddo

  end subroutine fft_r2c_xy_out_s32

  subroutine fft_r2c_xy_out2(input, output)

     implicit none

     real(wp), dimension(:,:), contiguous, intent(INOUT) :: input
     real(wp), dimension(:,:), contiguous, intent(OUT)   :: output

     call fft_r2c_x(input)
     call transpose_xy(input, output)
     call fft_r2c_y(output)

  end subroutine fft_r2c_xy_out2

  !---------------------------------------------------------
  !> @brief Multiple c2r transform in y and then x direction.
  !
  !> Multiple c2r transform in y and then x direction. The input should have
  !! modes from 0 to N/2-1, and the last mode N/2 is assumed to be
  !! zero. The input data is overwritten after transformation.
  !
  !> @param[in]  input   data to be transformed, y-pencil
  !> @param[out] output  transformed data, x-pencil
  !> @param[in]  nz      number of grids to be transformed in z direction
  !---------------------------------------------------------
  subroutine fft_c2r_xy_out(input, output)

     implicit none

     real(wp), dimension(:,:,:), contiguous, intent(INOUT) :: input
     real(wp), dimension(:,:,:), contiguous, intent(OUT) :: output

     call fft_c2r_y(input)
     call transpose_yx(input, output)
     call fft_c2r_x(output)

  end subroutine fft_c2r_xy_out

  subroutine fft_c2r_xy_out2(input, output)

     implicit none

     real(wp), dimension(:,:), contiguous, intent(INOUT) :: input
     real(wp), dimension(:,:), contiguous, intent(OUT) :: output

     call fft_c2r_y(input)
     call transpose_yx(input, output)
     call fft_c2r_x(output)

  end subroutine fft_c2r_xy_out2

end module fft
