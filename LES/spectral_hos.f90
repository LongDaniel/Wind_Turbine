module spectral_hos

  ! use hos_param
  use hos_param, only: myid_hos,ierr_hos
  use decomp, only: mpi_comm_2d_col
  use fft_hos

  implicit none

  private

  public :: pdfx_hos, pdfy_hos
  public :: pdfxx_hos, pdfyy_hos
  public :: onetoall, alltoone, dealiasxy_hos
  public :: spec_x, spec_y
  public :: get_amp_phase
    
  interface dealiasxy_hos
     module procedure dealiasxy_arbitrary, dealiasxy_twothirds
  end interface

  interface pdfx_hos
     module procedure pdfx_hos_1, pdfx_hos_2
  end interface

  interface pdfy_hos
     module procedure pdfy_hos_1, pdfy_hos_2
  end interface

  interface onetoall
     module procedure new_onetoall
  end interface
  
  interface alltoone
     module procedure new_alltoone
  end interface

contains

  !---------------------------------------------------------

  subroutine dealiasxy_arbitrary(f,ratio)

    implicit none

    real(wp), intent(in) :: ratio
    real(wp), intent(inout), dimension(:,:) :: f

    integer nxpat, nypat
    
    nxpat = int(nxhos * ratio / 2.0)
    nypat = int(nyhos * ratio / 2.0)


    call fft_for_x_hos(f)       
    f(nxpat + 1 : nxhos,:) = 0.0
    call transpose_2d(f,bufb_hos, nxhos, nyhos/ncpu_hos)
    call fft_for_x_hos(bufb_hos)
    bufb_hos(nypat + 1 : nyhos,:) = 0.0
    
    call fft_bac_x_hos(bufb_hos)
    call transpose_2d(bufb_hos,f, nyhos, nxhos/ncpu_hos)
    call fft_bac_x_hos(f)       
    
    return
  end subroutine dealiasxy_arbitrary

  !---------------------------------------------------------

  subroutine dealiasxy_twothirds(f)

    implicit none

    real(wp), intent(inout), dimension(:,:) :: f

    call fft_for_x_hos(f)
    f(nxhos/3*2 + 1 : nxhos,:) = 0.0
    call transpose_2d(f,bufb_hos, nxhos, nyhos/ncpu_hos)
    call fft_for_x_hos(bufb_hos)
    bufb_hos(nyhos/3*2 + 1 : nyhos,:) = 0.0
    
    call fft_bac_x_hos(bufb_hos)
    call transpose_2d(bufb_hos,f, nyhos, nxhos/ncpu_hos)
    call fft_bac_x_hos(f)
    
    return
  end subroutine dealiasxy_twothirds

  !---------------------------------------------------------

  !---------------------------------------------------------
  subroutine pdfx_hos_1(a,ax,pex)

    implicit none

    real(wp), intent(in) :: pex
    real(wp), intent(in), dimension(:,:) :: a
    real(wp), intent(inout), dimension(:,:) :: ax

    complex(wp) :: pmodex
    integer :: i

    call my_fftw_execute_dft_r2c(plan_r2c_x_hos,a,p_buf_x_hos)

    do i = 1, nxhos / 2
       pmodex = cmplx(0.0,(i-1)*pex,wp)
       p_buf_x_hos(i,:) = p_buf_x_hos(i,:) * pmodex / nxhos
    end do
    
    p_buf_x_hos(nxhos/2+1,:) = 0.0
    call my_fftw_execute_dft_c2r(plan_c2r_x_hos,p_buf_x_hos,ax)

  end subroutine pdfx_hos_1

  !---------------------------------------------------------


  subroutine pdfx_hos_2(a,ax,pex,nx,ny)

    implicit none

    integer, intent(in) :: nx,ny
    real(wp), intent(in) :: pex
    real(wp), intent(in), dimension(:,:) :: a
    real(wp), intent(inout), dimension(:,:) :: ax

    real(wp) tmp1, tmp2
    integer :: i,k,nmod

    call fft_for_x_hos(a,ax,nx,ny)

    do i = 1, nx, 2
       nmod = (i-1)/1
       do k = 1, ny
          tmp1 = -pex*nmod*ax(i+1,k)
          tmp2 = pex*nmod*ax(i,k)
          ax(i,k) = tmp1
          ax(i+1,k) = tmp2
       end do
    end do

    call fft_bac_x_hos(ax,nx,ny)

  end subroutine pdfx_hos_2

  !---------------------------------------------------------

  subroutine pdfxx_hos(a,axx,pex)

    implicit none

    real(wp), intent(in) :: pex
    real(wp), intent(in), dimension(:,:) :: a
    real(wp), intent(inout), dimension(:,:) :: axx

    integer :: i

    call my_fftw_execute_dft_r2c(plan_r2c_x_hos,a,p_buf_x_hos)

    do i = 1, nxhos / 2 + 1
       p_buf_x_hos(i,:) = -p_buf_x_hos(i,:) * ((i-1)*pex)**2 / nxhos
    end do

    p_buf_x_hos(nxhos/2+1,:) = 0.0
    call my_fftw_execute_dft_c2r(plan_c2r_x_hos,p_buf_x_hos,axx)

  end subroutine pdfxx_hos

  !---------------------------------------------------------

  subroutine pdfy_hos_1(a,ay,pey)

    implicit none
    
    real(wp), intent(in) :: pey
    real(wp), intent(in), dimension(:,:) :: a
    real(wp), intent(inout), dimension(:,:) :: ay

    complex(wp) :: pmodey
    integer :: i

    call transpose_2d(a,bufb_hos,nxhos,nyhos/ncpu_hos)
    call my_fftw_execute_dft_r2c(plan_r2c_y_hos,bufb_hos,p_buf_y_hos)

    do i = 1, nyhos / 2
       pmodey = cmplx(0,(i - 1) * pey, wp)
       p_buf_y_hos(i,:) = p_buf_y_hos(i,:) * pmodey  / nyhos
    end do
    p_buf_y_hos(nyhos/2+1,:) = 0.0
    call my_fftw_execute_dft_c2r(plan_c2r_y_hos,p_buf_y_hos,bufb_hos)

    call transpose_2d(bufb_hos,ay,nyhos,nxhos/ncpu_hos)

  end subroutine pdfy_hos_1

  !---------------------------------------------------------

  subroutine pdfy_hos_2(a,ay,pey,nx,ny)

    implicit none

    integer, intent(in) :: nx,ny
    real(wp), intent(in) :: pey
    real(wp), intent(in), dimension(:,:) :: a
    real(wp), intent(inout), dimension(:,:) :: ay

    real(wp), allocatable, dimension(:,:) :: tmpa, tmpay
    real(wp) tmp1, tmp2
    integer :: j,k,nmod

    allocate(tmpa(ny,nx),tmpay(ny,nx))

    tmpa = transpose(a)
    call pdfx_hos_2(tmpa,tmpay,pey,ny,nx)
    ay = transpose(tmpay)

    deallocate(tmpa,tmpay)

  end subroutine pdfy_hos_2

  !---------------------------------------------------------

  subroutine pdfyy_hos(a,ayy,pey)

    implicit none

    real(wp), intent(in) :: pey
    real(wp), intent(in), dimension(:,:) :: a
    real(wp), intent(inout), dimension(:,:) :: ayy

    integer :: i

    call transpose_2d(a,bufb_hos,nxhos,nyhos/ncpu_hos)
    call my_fftw_execute_dft_r2c(plan_r2c_y_hos,bufb_hos,p_buf_y_hos)

    do i = 1, nyhos / 2
       p_buf_y_hos(i,:) = -p_buf_y_hos(i,:) * ((i - 1) * pey)**2 / nyhos
    end do
    p_buf_y_hos(nyhos/2+1,:) = 0.0
    call my_fftw_execute_dft_c2r(plan_c2r_y_hos,p_buf_y_hos,bufb_hos)

    call transpose_2d(bufb_hos,ayy,nyhos,nxhos/ncpu_hos)

  end subroutine pdfyy_hos

  !---------------------------------------------------------

  !========================================================

  !========================================================
  !> @brief Allocate and initialize variables in module navier.
  !========================================================

   subroutine new_alltoone(f,fa)

     implicit none
     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:,:) :: fa

     integer j,js
     real(wp), dimension(nxhos,nyhos) :: f0
     f0 = 0
     fa = 0
     !print *, 'plyudebug, ato 1'
     do js = 1, nyhos / ncpu_hos
        j = myid_hos * nyhos / ncpu_hos + js
        !print *, 'plyudebug_nl', j, js
	f0(:,j) = f(:,js)
     end do
     !print *, 'plyudebug, ato 2'
     call mpi_allreduce(f0,fa,nxhos*nyhos,mpi_double_precision,mpi_sum,mpi_comm_2d_col,ierr_hos)
     !print *, 'plyudebug, ato 3'
   end subroutine new_alltoone

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine new_onetoall(fa,f)

     implicit none
     real(wp), intent(inout), dimension(:,:) :: f,fa

     integer j,js
     real(wp), dimension(nxhos,nyhos) :: f0
     f = 0
     f0 = 0
     if (myid_hos /= 0) then
        fa = 0
     endif

     call mpi_allreduce(fa,f0,nxhos*nyhos,mpi_double_precision,mpi_sum,mpi_comm_2d_col,ierr_hos)

     do j = 1, nyhos / ncpu_hos
        js = myid_hos * nyhos / ncpu_hos + j
        f(:,j) = f0(:,js)
     end do

   end subroutine new_onetoall

!--------------------------------------------------------------------

   subroutine spec_x(f,skx)

     implicit none
     
     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:) :: skx

     real(wp), allocatable, dimension(:,:) :: tmp      
     real(wp) ftn
     integer l,m

     allocate(tmp(nxhos,nyhos))

     skx = 0
     call fft_for_x_hos(f,tmp,nxhos,nyhos)

     skx = 0.0
     do l = 3, nxhos, 2
        do m = 1, nyhos
           ftn = tmp(l,m)**2 + tmp(l+1,m)**2
           skx((l-1)/2) = skx((l-1)/2) + 2 * ftn / pex_hos / nyhos
        end do
     end do

     deallocate(tmp)

   end subroutine spec_x

!-------------------------------------------------------------------------

   subroutine spec_y(f,sky)

     implicit none

     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:) :: sky

     real(wp), allocatable, dimension(:,:) :: tmp
     real(wp) ftn
     integer l,m

     allocate(tmp(nyhos,nxhos))

     sky = 0
     
     tmp = transpose(f)
     call fft_for_x_hos(tmp,nyhos,nxhos)

     sky = 0.0
     do m = 3, nyhos, 2
        do l = 1, nxhos
           ftn = tmp(m,l)**2 + tmp(m+1,l)**2
           sky((m-1)/2) = sky((m-1)/2) + 2 * ftn / pey_hos / nxhos
        end do
     end do

     deallocate(tmp)

   end subroutine spec_y

   subroutine get_amp_phase(f,amp,theta,nx)
     
     implicit none

     integer, intent(in) :: nx
     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:) :: amp, theta

     real(wp), allocatable, dimension(:) :: tmp
     integer i
     real(wp) tmtmp(nx,1)

     allocate(tmp(nx))
     call fft_for_x_hos(f,tmtmp,nx,1)
     
     tmp = tmtmp(:,1)

     do i = 2, nx/2
        amp(i)=2.0*sqrt(tmp(2*i-1)**2+tmp(2*i)**2)
        if(tmp(2*i-1) > 0) then
           theta(i)=atan(tmp(2*i)/tmp(2*i-1))
        elseif(tmp(2*i-1) < 0) then
           if(tmp(2*i) > 0) then
              theta(i)=atan(tmp(2*i)/tmp(2*i-1))+twopi/2
           else
              theta(i)=atan(tmp(2*i)/tmp(2*i-1))-twopi/2
           endif
        elseif(tmp(2*i) > 0) then
           theta(i)=twopi/2/2.0
        elseif(tmp(2*i) < 0) then
           theta(i)=-twopi/2/2.0
        else
           theta(i)=0.0
        endif
        theta(i)=-theta(i)
     end do

     deallocate(tmp)
     
   end subroutine get_amp_phase

end module spectral_hos
