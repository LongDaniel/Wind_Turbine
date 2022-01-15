module spectral_hos

  use hos_param
  use fft_hos

  implicit none

  private

  public :: pdfx_hos, pdfy_hos
  public :: pdfxx_hos, pdfyy_hos
  public :: onetoall, alltoone, dealiasxy_hos
  public :: spec_x, spec_y, spec_xy, spec_1dk, spec_1df

  interface pdfx_hos
     module procedure pdfx_hos_1, pdfx_hos_2
  end interface

  interface pdfy_hos
     module procedure pdfy_hos_1, pdfy_hos_2
  end interface
    
  interface dealiasxy_hos
     module procedure dealiasxy_arbitrary, dealiasxy_twothirds
  end interface

  interface onetoall
     module procedure new_onetoall, onetoall_arb
  end interface
  
  interface alltoone
     module procedure new_alltoone, alltoone_arb
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

  subroutine smooth_1_hos(eta,vps,ratio)

    implicit none
    
    real(wp), intent(in) :: ratio
    real(wp), intent(inout), dimension(:,:) :: eta, vps

    integer nxpat, nypat

    nxpat = int(nxhos * ratio / 2.0)
    nypat = int(nyhos * ratio / 2.0)
    call fft_for_x_hos(eta)
    eta(nxpat + 1 : nxhos,:) = 0.0
    call transpose_2d(eta,bufb_hos, nxhos, nyhos/ncpu_hos)
    call fft_for_x_hos(bufb_hos,nyhos,nxhos/ncpu_hos)
    bufb_hos(nypat + 1 : nyhos,:) = 0.0

    call fft_bac_x_hos(bufb_hos,nxhos,nyhos/ncpu_hos)
    call transpose_2d(bufb_hos,eta, nyhos, nxhos/ncpu_hos)
    call fft_bac_x_hos(eta)

    call fft_for_x_hos(vps)
    vps(nxpat + 1 : nxhos,:) = 0.0
    call transpose_2d(vps,bufb_hos, nxhos, nyhos/ncpu_hos)
    call fft_for_x_hos(bufb_hos,nyhos,nxhos/ncpu_hos)
    bufb_hos(nypat + 1 : nyhos,:) = 0.0

    call fft_bac_x_hos(bufb_hos,nxhos,nyhos/ncpu_hos)
    call transpose_2d(bufb_hos,vps, nyhos, nxhos/ncpu_hos)
    call fft_bac_x_hos(vps)

  end subroutine smooth_1_hos

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
    integer :: nmod

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

     do js = 1, nyhos / ncpu_hos
        j = myid * nyhos / ncpu_hos + js
        f0(:,j) = f(:,js)
     end do

     call mpi_allreduce(f0,fa,nxhos*nyhos,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

   end subroutine new_alltoone

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine alltoone_arb(f,fa,nx,ny)

     implicit none
     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:,:) :: fa
     integer, intent(in) :: nx,ny

     integer j,js
     real(wp), allocatable, dimension(:,:) :: f0
     
     allocate(f0(nx,ny))
     f0 = 0
     fa = 0

     do js = 1, ny / ncpu_hos
        j = myid * ny / ncpu_hos + js
        f0(:,j) = f(:,js)
     end do

     call mpi_allreduce(f0,fa,nx*ny,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
     deallocate(f0)

   end subroutine alltoone_arb

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
     if (myid /= 0) then
        fa = 0
     endif

     call mpi_allreduce(fa,f0,nxhos*nyhos,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

     do j = 1, nyhos / ncpu_hos
        js = myid * nyhos / ncpu_hos + j
        f(:,j) = f0(:,js)
     end do

   end subroutine new_onetoall

!--------------------------------------------------------------------

   subroutine onetoall_arb(fa,f,nx,ny)

     implicit none
     real(wp), intent(inout), dimension(:,:) :: f,fa
     integer, intent(in) :: nx,ny

     integer j,js
     real(wp), allocatable, dimension(:,:) :: f0
     f = 0
     allocate(f0(nx,ny))
     f0 = 0
     if (myid /= 0) then
        fa = 0
     endif

     call mpi_allreduce(fa,f0,nx*ny,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

     do j = 1, ny / ncpu_hos
        js = myid * ny / ncpu_hos + j
        f(:,j) = f0(:,js)
     end do

     deallocate(f0)
   end subroutine onetoall_arb

!--------------------------------------------------------------------


   subroutine spec_xy(f,skxy)

     implicit none
     
     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:,:) :: skxy

     real(wp), allocatable, dimension(:,:) :: tmp
     real(wp) ftn, wkx, wky, con, wav
     integer i,j,l,m,lp1,mp1,root

     allocate(tmp(nxhos,nyhos))

     call fft_for_xy_hos(f,tmp,nxhos,nyhos)

     skxy = 0.0_wp
     con = 0.0_wp

     do l = 1, nxhos - 1, 2
        lp1 = l+1
        wkx = (l - 1) / 2 * pex_hos
        i = (l + 1) / 2
        do m = 1, nyhos - 1, 2
           mp1 = m+1
           j = (m + 1) / 2
           wky = (m - 1) / 2 * pey_hos
           wav = sqrt(wkx**2 + wky**2)
                      
           ftn = tmp(l,m)**2 + tmp(l,mp1)**2 + tmp(lp1,m)**2 + tmp(lp1,mp1)**2
           con = 4.0_wp
           if (l == 1 .or. m == 1) con = 2.0_wp
           if (l == 1 .and. m == 1) con = 1.0_wp
           skxy(i,j) = skxy(i,j) + con * ftn / pex_hos / pey_hos
        end do
     end do

     root = 0
     call mpi_bcast(skxy, nxhos*nyhos , mpi_double_precision, root, mpi_comm_world,ierr)

     deallocate(tmp)

   end subroutine spec_xy

!--------------------------------------------------------------------

   subroutine spec_x(f,skx)

     implicit none
     
     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:) :: skx

     real(wp), allocatable, dimension(:,:) :: tmp      
     real(wp) ftn
     integer l,m,root

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

     root = 0
     call mpi_bcast(skx, nxhos/2 , mpi_double_precision, root, mpi_comm_world,ierr)

   end subroutine spec_x

!-------------------------------------------------------------------------

   subroutine spec_y(f,sky)

     implicit none

     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:) :: sky

     real(wp), allocatable, dimension(:,:) :: tmp
     real(wp) ftn
     integer l,m,root

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

     root = 0
     call mpi_bcast(sky, nyhos/2 , mpi_double_precision, root, mpi_comm_world,ierr)

   end subroutine spec_y

!-------------------------------------------------------------------------

   subroutine spec_1dk(f, sk, nmax, dwk)
     implicit none

     integer, intent(in) :: nmax
     real(wp), intent(in) :: dwk
     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:) :: sk

     real(wp), allocatable, dimension(:,:) :: tmp
     real(wp) fkn, wkx, wky, wav, con
     integer l,m,n,lp1,mp1,root

     allocate(tmp(nxhos,nyhos))

     sk = 0
     call fft_for_xy_hos(f,tmp,nxhos,nyhos)
     do l=1,nxhos-1,2
        lp1=l+1
        wkx=(l-1)/2*pex_hos
        do m=1,nyhos-1,2
           wky=(m-1)/2*pey_hos
           wav=sqrt(wkx**2+wky**2)
           do n=1,nmax
              if(wav >= ((n-0.5)*dwk).and.wav < ((n+0.5)*dwk))then
!                 nwk(n)=nwk(n)+1
                 exit
              endif
           enddo
           mp1=m+1
           fkn=tmp(l,m)**2+tmp(l,mp1)**2+tmp(lp1,m)**2+tmp(lp1,mp1)**2
           con=4.
           if(l.eq.1.or.m.eq.1) con=2.
           if(l.eq.1.and.m.eq.1) con=1.
           if (n == nmax+1) n = nmax
           sk(n)=sk(n)+con*fkn
        enddo
     enddo
     
     do n=1,nmax
        sk(n)=sk(n)/dwk
     enddo
     deallocate(tmp)
   end subroutine spec_1dk

   subroutine spec_1df(f, sf, nmax, dwf)
     implicit none

     integer, intent(in) :: nmax
     real(wp), intent(in) :: dwf
     real(wp), intent(in), dimension(:,:) :: f
     real(wp), intent(out), dimension(:) :: sf

     real(wp), allocatable, dimension(:,:) :: tmp
     real(wp) fkn, wkx, wky, wav, con, freq
     integer l,m,n,lp1,mp1,root

     allocate(tmp(nxhos,nyhos))

     sf = 0
     call fft_for_xy_hos(f,tmp,nxhos,nyhos)
     do l=1,nxhos-1,2
        lp1=l+1
        wkx=(l-1)/2*pex_hos
        do m=1,nyhos-1,2
           wky=(m-1)/2*pey_hos
           wav=sqrt(wkx**2+wky**2)
           freq = sqrt(wav/fr2)/twopi
           do n=1,nmax
              if(freq >= ((n-0.5)*dwf) .and. freq < ((n+0.5)*dwf))then
!                 nwk(n)=nwk(n)+1
                 exit
              endif
           enddo
           mp1=m+1
           fkn=tmp(l,m)**2+tmp(l,mp1)**2+tmp(lp1,m)**2+tmp(lp1,mp1)**2
           con=4.
           if(l.eq.1.or.m.eq.1) con=2.
           if(l.eq.1.and.m.eq.1) con=1.
           if (n == nmax+1) n = nmax
           sf(n)=sf(n)+con*fkn
        enddo
     enddo
     
     do n=1,nmax
        sf(n)=sf(n)/dwf
     enddo
     deallocate(tmp)
   end subroutine spec_1df

end module spectral_hos
