!all checked
module solver_common
   use decomp
   implicit none

   private

   real(wp), allocatable, dimension(:,:,:), public :: tmp_x1, tmp_x2, tmp_x3, tmp_x4
   real(wp), allocatable, dimension(:,:,:), public :: tmp_x5, tmp_x6, tmp_x7!, tmp_x8
   real(wp), allocatable, dimension(:,:,:), public :: tmp_y1!, tmp_y2

   !real(wp), allocatable, dimension(:,:)   :: uzeta_coef, wzeta_coef

   public :: calc_uzeta, calc_wzeta
   public :: les_filter
   interface les_filter
      module procedure les_filter_2d, les_filter_3d
   end interface les_filter

contains

   !========================================================
   !> @brief Calculate \f$ \partial f/\partial \zeta \f$ for
   !!  f located at u nodes.
   !
   !> Second order scheme. No boundary conditions are imposed.
   !
   !> @param[in]  f      the function to calculate
   !> @param[in]  expand the number of ghost points of input f
   !> @param[out] fz     the derivative of f
   !========================================================

  !checked
  !input: f, expand, 
  !output: fz
  subroutine calc_uzeta(f, fz, expand)
      use grid, only : dz

      implicit none

      integer, intent(IN) :: expand
      real(wp), dimension(xsz(1),xsz(2),1-expand:xsz(3)+expand), intent(IN) :: f
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: fz

      integer :: k

      do k=1,xsz(3)
         fz(:,:,k)=(f(:,:,k+1)-f(:,:,k-1))/(dz(k)+dz(k-1))
      end do

      if (isbot) then
         fz(:,:,2) = (f(:,:,3)-f(:,:,2))/dz(2)
      end if      

      ! note: k = nz-1, finite difference schemes are not same in subroutine fun_u
      ! for terms: u_zeta and t11_zeta/t12_zeta
      
      if (istop) then
         fz(:,:,xsz(3)-1) = (f(:,:,xsz(3)-1)-f(:,:,xsz(3)-2))/(2.0_wp*dz(xsz(3)-2))
         fz(:,:,xsz(3)) = (f(:,:,xsz(3)+1)-f(:,:,xsz(3)-1))/dz(xsz(3)-2)
      end if

   end subroutine calc_uzeta

   !========================================================
   !> @brief Calculate \f$ \partial f/\partial \zeta \f$ for
   !!  f located at w nodes.
   !
   !> Second order scheme. No boundary conditions are imposed.
   !
   !> @param[in]  f      the function to calculate
   !> @param[in]  expand the number of ghost points of input f
   !> @param[out] fz     the derivative of f
   !========================================================

   !checked
   !input: f, expand, dz
   !output: fz
   subroutine calc_wzeta(f, fz, expand)
      !use grid, only : dz, dzw
     use grid, only: dz
      implicit none

      integer, intent(IN) :: expand
      real(wp), dimension(xsz(1),xsz(2),1-expand:xsz(3)+expand), intent(IN) :: f
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: fz

      integer :: k

      do k=1,xsz(3)
         fz(:,:,k)=(f(:,:,k+1)-f(:,:,k-1)) /(2.0_wp*dz(k))
      end do

   end subroutine calc_wzeta

   !same as Eillott
   !checked
   !input/output: a
   !input: mfilt, deltax, deltay, dealias(dealiasing indicator)
   
   subroutine les_filter_2d(a, mfilt, deltax, deltay, dealias)
      use param, only : pex, pey
      use fft
      use utils
      use constants, only : PI
      implicit none

      integer,  intent(IN) :: mfilt
      real(wp), intent(IN) :: deltax, deltay
      real(wp), dimension(:,:), intent(INOUT) :: a
      logical,  intent(IN), optional :: dealias

      integer  :: i, j, ckxi, ckyi
      real(wp) :: ckx, cky
      real(wp), allocatable, dimension(:,:) :: tmp
      logical  :: dealias_
      
      allocate(tmp(ysz(1), ysz(2)))
      if (present(dealias)) then
         dealias_ = dealias
      else
         dealias_ = .false.
      end if

      if (mfilt == 2) then
         ckx = PI/2/deltax
         cky = PI/2/deltay
         ckxi = ceiling(ckx/pex)
         ckyi = ceiling(cky/pey)
      end if

      call fft_r2c_x(a)

      if (mfilt == 1) then
         do j=1, xsz(2)
            do i=0, xsz(1)/2-1
               a(2*i+1,j) = a(2*i+1,j) * exp(-(pex*i*2*deltax)**2/24)
               a(2*i+2,j) = a(2*i+2,j) * exp(-(pex*i*2*deltax)**2/24)
            end do
         end do
      else if (mfilt == 2) then
         a(2*ckxi+1:,:) = 0
      else
         call error_abort(12, 'Invalid MFILT.')
      end if
      if (dealias_) then
         a(((xsz(1)/2)*2/3)*2+1:,:) = 0
      end if
      a = a/xsz(1)

      call transpose_xy(a, tmp)

      call fft_r2c_y(tmp)

      if (mfilt == 1) then
         do j=1, ysz(2)
            do i=0, ysz(1)/2-1
               tmp(2*i+1,j) = tmp(2*i+1,j) * exp(-(pey*i*2*deltay)**2/24)
               tmp(2*i+2,j) = tmp(2*i+2,j) * exp(-(pey*i*2*deltay)**2/24)
            end do
         end do
      else if (mfilt == 2) then
         tmp(2*ckyi+1:,:) = 0
      end if
      if (dealias_) then
         tmp(((ysz(1)/2)*2/3)*2+1:,:) = 0
      end if
      tmp = tmp/ysz(1)

      call fft_c2r_y(tmp)

      call transpose_yx(tmp, a)

      call fft_c2r_x(a)
      
      deallocate(tmp)

   end subroutine les_filter_2d

   ! same as Eillott
   ! checked
   ! input/output: a
   ! input: mfilt, deltax, deltay, dealias(dealiasing indicator)
   subroutine les_filter_3d(a, mfilt, deltax, deltay, dealias)
      use param, only : pex, pey
      use fft
      use utils
      use constants, only : PI
      implicit none

      integer,  intent(IN) :: mfilt
      real(wp), intent(IN) :: deltax, deltay
      real(wp), dimension(:,:,:), intent(INOUT) :: a
      logical,  intent(IN), optional :: dealias

      integer  :: i, j, k, ckxi, ckyi
      real(wp) :: ckx, cky
      real(wp), allocatable, dimension(:,:,:) :: tmp
      logical  :: dealias_

      allocate(tmp(ysz(1), ysz(2), ysz(3)))
      if (present(dealias)) then
         dealias_ = dealias
      else
         dealias_ = .false.
      end if

      if (mfilt == 2) then
         ckx = PI/2/deltax
         cky = PI/2/deltay
         ckxi = ceiling(ckx/pex)
         ckyi = ceiling(cky/pey)
      end if

      call fft_r2c_x(a)

      if (mfilt == 1) then
         do k=1, xsz(3)
            do j=1, xsz(2)
               do i=0, xsz(1)/2-1
                  a(2*i+1,j,k) = a(2*i+1,j,k) * exp(-(pex*i*2*deltax)**2/24)
                  a(2*i+2,j,k) = a(2*i+2,j,k) * exp(-(pex*i*2*deltax)**2/24)
               end do
            end do
         end do
      else if (mfilt == 2) then
         a(2*ckxi+1:,:,:) = 0
      else
         call error_abort(12, 'Invalid MFILT.')
      end if
      if (dealias_) then
         a(((xsz(1)/2)*2/3)*2+1:,:,1:xsz(3)) = 0
      end if
      a = a/xsz(1)

      call transpose_xy(a, tmp)

      call fft_r2c_y(tmp)

      if (mfilt == 1) then
         do k=1, ysz(3)
            do j=1, ysz(2)
               do i=0, ysz(1)/2-1
                  tmp(2*i+1,j,k) = tmp(2*i+1,j,k) * exp(-(pey*i*2*deltay)**2/24)
                  tmp(2*i+2,j,k) = tmp(2*i+2,j,k) * exp(-(pey*i*2*deltay)**2/24)
               end do
            end do
         end do
      else if (mfilt == 2) then
         tmp(2*ckyi+1:,:,:) = 0
      end if
      if (dealias_) then
         tmp(((ysz(1)/2)*2/3)*2+1:,:,1:ysz(3)) = 0
      end if
      tmp = tmp/ysz(1)

      call fft_c2r_y(tmp)

      call transpose_yx(tmp, a)

      call fft_c2r_x(a)

      deallocate(tmp)
   end subroutine les_filter_3d

end module solver_common
