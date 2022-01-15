!-----------------------------------------------------------
! MODULE: spectral
!
!> Base subroutines for spectral method, including dealising
!!  and differentiation.
!
!> @author Anqing Xuan
!
!-----------------------------------------------------------
module spectral
   use decomp
   use fft

   implicit none

   private
   complex(wp), dimension(:,:), allocatable, save :: bufx_cplx, bufy_cplx
   
   !> changed from private to public by plyu
   real(wp), dimension(:,:,:), allocatable, public, save :: bufy
   !> added by plyu
   real(wp), dimension(:,:,:), allocatable, save :: bufy2, bufx
   !real(wp), dimension(:), allocatable, public, save :: co_s32_x, co_s32_y    

   public :: spectral_init
   public :: pdfx, pdfxx, pdfy, pdfy_x, pdfyy
   public :: pdfx_, pdfxx_, pdfy_, pdfyy_
   public :: dealiasxy, dealiasx
   public :: dealiasxy_2, dealiasx_2
   public :: phase_shift_x_2d

   interface dealiasxy
      module procedure dealiasxy2, dealiasxy3
      module procedure dealiasxy3_s !> added by plyu
   end interface dealiasxy

   interface dealiasxy_2
      module procedure dealiasxy2_2, dealiasxy3_2
   end interface dealiasxy_2

   interface pdfx
      module procedure pdfx_out2, pdfx_out3
      module procedure pdfx_in2, pdfx_in3
      module procedure pdfx_out3_s !> added by plyu
   end interface pdfx

   interface pdfxx
      module procedure pdfxx_out2, pdfxx_out3
      module procedure pdfxx_in2, pdfxx_in3
   end interface pdfxx

   interface pdfy
      module procedure pdfy_out2, pdfy_out3
      module procedure pdfy_in2, pdfy_in3
   end interface pdfy

   interface pdfyy
      module procedure pdfyy_out2, pdfyy_out3
      module procedure pdfyy_in2, pdfyy_in3
   end interface pdfyy

   interface pdfy_x
      module procedure pdfy_x_out2, pdfy_x_out3
      module procedure pdfy_x_in2, pdfy_x_in3
      module procedure pdfy_x_out3_s !> added by plyu
   end interface pdfy_x

contains

   subroutine spectral_init
      implicit none
      integer :: i

      allocate(bufy(ysz(1),ysz(2),ysz(3)))
      allocate(bufx_cplx(0:xsz(1)/2,xsz(2)))
      allocate(bufy_cplx(0:ysz(1)/2,ysz(2)))
      
      !> added by plyu
      allocate(bufy2(ysz(1),ysz(2),ysz(3)))
      allocate(bufx(xsz(1),xsz(2),xsz(3)))
      !> following lines has been moved to module fft
      !allocate(co_s32_x(xsz(1)/2), co_s32_y(ysz(1)/2))
      !do i = 1, xsz(1)/2
      !  co_s32_x(i) = exp(-36.0*(real(i,wp)/real(xsz(1)/2,wp))**36)
      !enddo
      !do i = 1, ysz(1)/2
      !  co_s32_y(i) = exp(-36.0*(real(i,wp)/real(ysz(1)/2,wp))**36)
      !enddo

   end subroutine spectral_init

   !--------------------------------------------------------
   !> @brief Calculate x partial derivative.
   !
   !> @param[in] input Input data, physical domain, x-pencil
   !> @param[out] output Output data, physical domain, x-pencil
   !> @param[in] nz
   !> @param[in] pex   unit wavenumber
   !--------------------------------------------------------
   subroutine pdfx_out3(input, output, pex)
     !> hacked by plyu
     use discontinuity_smooth, only: tbn_lim_x
     implicit none

     real(wp), intent(IN) :: pex
     real(wp), dimension(:,:,:), intent(IN),  contiguous :: input
     real(wp), dimension(:,:,:), intent(OUT), contiguous :: output

     integer  :: i, j, k
     real(wp) :: invn

     if ((ids==1 .or. ids==2) .and. ASSOCIATED(tbn_lim_x)) then
       call pdfx_out3_s(input, output, pex, tbn_lim_x)
     elseif (ids==0) then
       !> original version:
       invn = pex/real(xsz(1),wp)
       do k = 1, size(input,3)
          call fft_r2c_x(input(:,:,k), bufx_cplx(:,:))
          do j = 1, xsz(2)
             do i = 0, xsz(1)/2-1
                ! tmp = output(i*2+1,j,k)
                ! output(i*2+1,j,k) = -output(i*2+2,j,k)*i*invn
                ! output(i*2+2,j,k) = tmp*i*invn
                bufx_cplx(i,j) = cmplx(0, i*invn, wp) * bufx_cplx(i,j)
             end do
             bufx_cplx(xsz(1)/2,j) = 0
          end do
          call fft_c2r_x(bufx_cplx, output(:,:,k))
       end do
     endif

   end subroutine pdfx_out3

   subroutine pdfx_out2(input, output, pex)
      implicit none

      real(wp), intent(IN) :: pex
      real(wp), dimension(:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:), intent(OUT), contiguous :: output

      integer  :: i, j
      real(wp) :: invn

      call fft_r2c_x(input, bufx_cplx)

      invn = pex/real(xsz(1),wp)
      do j = 1, xsz(2)
         do i = 0, xsz(1)/2-1
            ! tmp = output(i*2+1,j)
            ! output(i*2+1,j) = -output(i*2+2,j)*i*invn
            ! output(i*2+2,j) = tmp*i*invn
            bufx_cplx(i,j) = cmplx(0, i*invn, wp) * bufx_cplx(i,j)
         end do
         bufx_cplx(xsz(1)/2,j) = 0
      end do

      call fft_c2r_x(bufx_cplx, output)

   end subroutine pdfx_out2

   subroutine pdfx_in3(inout, pex)
      implicit none

      real(wp), intent(IN) :: pex
      real(wp), dimension(:,:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j, k
      real(wp) :: invn

      invn = pex/real(xsz(1),wp)
      do k = 1, size(inout,3)
         call fft_r2c_x(inout(:,:,k), bufx_cplx)
         do j = 1, xsz(2)
            do i = 0, xsz(1)/2-1
               ! tmp = inout(i*2+1,j,k)
               ! inout(i*2+1,j,k) = -inout(i*2+2,j,k)*i*invn
               ! inout(i*2+2,j,k) = tmp*i*invn
               bufx_cplx(i,j) = cmplx(0, i*invn, wp) * bufx_cplx(i,j)
            end do
            bufx_cplx(xsz(1)/2,j) = 0
         end do
      call fft_c2r_x(bufx_cplx, inout(:,:,k))
      end do

   end subroutine pdfx_in3

   subroutine pdfx_in2(inout, pex)
      implicit none

      real(wp), intent(IN) :: pex
      real(wp), dimension(:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j
      real(wp) :: invn

      call fft_r2c_x(inout, bufx_cplx)

      invn = pex/real(xsz(1),wp)
      do j = 1, xsz(2)
         do i = 0, xsz(1)/2-1
            ! tmp = inout(i*2+1,j)
            ! inout(i*2+1,j) = -inout(i*2+2,j)*i*invn
            ! inout(i*2+2,j) = tmp*i*invn
            bufx_cplx(i,j) = cmplx(0, i*invn, wp) * bufx_cplx(i,j)
         end do
         bufx_cplx(xsz(1)/2,j) = 0
      end do

      call fft_c2r_x(bufx_cplx, inout)

   end subroutine pdfx_in2

   !> added by plyu. lim3d contains info about which location in (y,z) plane needs
   !! smoothing and where the discontinuity occurs in that x-direction line at
   !! certain (y,z) location.
   subroutine pdfx_out3_s(input, output, pex, lim3d)
      use fll_mods_m
      use discontinuity_smooth, only : cubic_smooth_3d, postproc_smooth_3d, &
        i_div_ustar, i_grad_p, n_pdfx
      implicit none

      real(wp), intent(IN) :: pex
      real(wp), dimension(:,:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:,:), intent(OUT), contiguous :: output
      type(dnode), pointer, intent(in) :: lim3d
      
      integer  :: i, j, k
      real(wp) :: invn, dx, twopi_local
      type(dnode), pointer :: ptmp, p_2d_indices, p_2d_index, lim1d
      
      !> for debug
      logical :: ok
      type(func_data_set) :: fpar
      integer :: i_dbg, nt0, i0, j0, k0, i_print

      i_dbg = 0 !> set by user to enable debug mode    
      i_print = 0 !> it could be determined by the program

      if (i_dbg > 0 .and. i_grad_p>0) then
        i_print = 1
        n_pdfx = n_pdfx + 1

        nt0 = 10; i0 = 67; j0 = 48; k0 = 17
      endif
       
      if (i_print > 0 .and. n_pdfx == nt0) then 
        !> original version:
        invn = pex/real(xsz(1),wp)
        do k = 1, size(input,3)
           call fft_r2c_x(input(:,:,k), bufx_cplx(:,:))
           do j = 1, xsz(2)
              do i = 0, xsz(1)/2-1
                 ! tmp = output(i*2+1,j,k)
                 ! output(i*2+1,j,k) = -output(i*2+2,j,k)*i*invn
                 ! output(i*2+2,j,k) = tmp*i*invn
                 bufx_cplx(i,j) = cmplx(0, i*invn, wp) * bufx_cplx(i,j)
              end do
              bufx_cplx(xsz(1)/2,j) = 0
           end do
           call fft_c2r_x(bufx_cplx, output(:,:,k))
        end do
           
        !> print result of original version
        if (myid.eq.0) print *, 'Times of calling pdfx_debug:', n_pdfx, ', original dy/dx'
        if (xst(2)<=j0 .and. (xst(2)+xsz(2)-1)>=j0 .and. xst(3)<=k0 &
          .and. (xst(3)+xsz(3)-1)>=k0) then
          open(1050, file='dbg_3_original_du.dat', action='write')
          print *, 'VARIABLS = X, dU'            
          do i = i0-5, i0+5
            print *, i, output(i, j0-xst(2)+1, k0-xst(3)+1)
          enddo
          close(1050)
        endif
      endif

      !> plyudebug: write a block before smoothing
      if (i_print > 0) then
        if(myid.eq.0) print *, 'Times of calling pdfx_debug:', n_pdfx, ', before_pdf'
        if ( n_pdfx .eq. nt0) then
          if (xst(2)<=j0 .and. (xst(2)+xsz(2)-1)>=j0 .and. xst(3)<=k0 &
            .and. (xst(3)+xsz(3)-1)>=k0) then
            open(1050, file='dbg_0_before_pdf.dat', action='write')
            print *, 'VARIABLS = X, U'            
            do i = i0-5, i0+5
              print *, i, input(i, j0-xst(2)+1, k0-xst(3)+1)
            enddo
            close(1050)
          endif
        endif
      endif

      bufx(:,:,:) = input(:,:,:)
      call cubic_smooth_3d (bufx, lim3d)
      
      if (i_print > 0) then
        if(myid.eq.0) print *, 'Times of calling pdfx_debug:', n_pdfx, ', after_smooth'
        if ( n_pdfx .eq. nt0) then
          if (xst(2)<=j0 .and. (xst(2)+xsz(2)-1)>=j0 .and. xst(3)<=k0 &
            .and. (xst(3)+xsz(3)-1)>=k0) then
            open(1050, file='dbg_1_after_smooth_u.dat', action='write')
            print *, 'VARIABLS = X, U'            
            do i = i0-5, i0+5
              print *, i, bufx(i, j0-xst(2)+1, k0-xst(3)+1)
            enddo
            close(1050)
          endif
        endif
      endif

      invn = pex/real(xsz(1),wp)
      do k = 1, size(input,3)
         call fft_r2c_x(bufx(:,:,k), bufx_cplx(:,:))
         do j = 1, xsz(2)
            do i = 0, xsz(1)/2-1
               ! tmp = output(i*2+1,j,k)
               ! output(i*2+1,j,k) = -output(i*2+2,j,k)*i*invn
               ! output(i*2+2,j,k) = tmp*i*invn
               bufx_cplx(i,j) = cmplx(0, i*invn, wp) * bufx_cplx(i,j)
            end do
            bufx_cplx(xsz(1)/2,j) = 0
         end do
         call fft_c2r_x(bufx_cplx, output(:,:,k))
      end do
      
      if (i_print > 0) then
        if(myid.eq.0) print *, 'Times of calling pdfx_debug:', n_pdfx, ', after pdfx'
        if ( n_pdfx .eq. nt0) then
          if (xst(2)<=j0 .and. (xst(2)+xsz(2)-1)>=j0 .and. xst(3)<=k0 &
            .and. (xst(3)+xsz(3)-1)>=k0) then
            open(1050, file='dbg_2_after_pdfx_u.dat', action='write')
            print *, 'VARIABLS = X, dU'            
            do i = i0-5, i0+5
              print *, i, output(i, j0-xst(2)+1, k0-xst(3)+1)
            enddo
            close(1050)
          endif
        endif
      endif
      
      if (ids==1) then
        !> do 2nd order central difference for the smoothed part
        twopi_local = 3.141592653589793238_wp * 2.0_wp
        dx = twopi_local / pex / real(xsz(1), wp)
        call postproc_smooth_3d (input, output, lim3d, 1, dx)
        
        if (i_print > 0) then
          if (myid.eq.0) print *, 'Times of calling pdfx_debug:', n_pdfx, ', after fdmx'
          if ( n_pdfx .eq. nt0) then
            if (xst(2)<=j0 .and. (xst(2)+xsz(2)-1)>=j0 .and. xst(3)<=k0 &
              .and. (xst(3)+xsz(3)-1)>=k0) then
              open(1050, file='dbg_3_after_fdm_u.dat', action='write')
              print *, 'VARIABLS = X, dU'            
              do i = i0-5, i0+5
                print *, i, output(i, j0-xst(2)+1, k0-xst(3)+1)
              enddo
              close(1050)
            endif
          endif
        endif

      endif
      !> ignore the precision loss if ids==2
   end subroutine pdfx_out3_s
   
   !--------------------------------------------------------
   !> @brief Calculate x partial derivative in spectral.
   !
   !> @param[in] input Input data, spectral domain, x-pencil
   !> @param[out] output Output data, spectral domain, x-pencil
   !> @param[in] nz
   !> @param[in] pex   unit wavenumber
   !--------------------------------------------------------
   subroutine pdfx_(input, output, nz, pex)
      implicit none

      integer,  intent(IN) :: nz
      real(wp), intent(IN) :: pex
      real(wp), dimension(xsz(1),xsz(2),nz) :: input, output

      integer  :: i, j, k
      real(wp) :: tmp, invn

      invn = pex/real(xsz(1),wp)
      do k = 1, nz
         do j = 1, xsz(2)
            do i = 0, xsz(1)/2-1
               tmp = input(i*2+1,j,k)
               output(i*2+1,j,k) = -input(i*2+2,j,k)*i*invn
               output(i*2+2,j,k) = tmp*i*invn
            end do
         end do
      end do

      end subroutine pdfx_

   !--------------------------------------------------------
   !> @brief Calculate second-order x partial derivative.
   !
   !> @param[in] input Input data, physical domain, x-pencil
   !> @param[out] output Output data, physical domain, x-pencil
   !> @param[in] nz
   !> @param[in] pex   unit wavenumber
   !--------------------------------------------------------
   subroutine pdfxx_out3(input, output, pex)
      implicit none

      real(wp), intent(IN) :: pex
      real(wp), dimension(:,:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:,:), intent(OUT), contiguous :: output

      integer  :: i, j, k
      real(wp) :: invn

      invn = -pex**2/real(xsz(1),wp)
      do k = 1, size(input,3)
         call fft_r2c_x(input(:,:,k), bufx_cplx)
         do j = 1, xsz(2)
            do i = 0, xsz(1)/2-1
               ! output(i*2+1,j,k) = output(i*2+1,j,k)*(i**2)*invn
               ! output(i*2+2,j,k) = output(i*2+2,j,k)*(i**2)*invn
               bufx_cplx(i,j) = (i**2)*invn * bufx_cplx(i,j)
            end do
            bufx_cplx(xsz(1)/2,j) = 0
         end do
         call fft_c2r_x(bufx_cplx, output(:,:,k))
      end do

   end subroutine pdfxx_out3

   subroutine pdfxx_out2(input, output, pex)
      implicit none

      real(wp), intent(IN) :: pex
      real(wp), dimension(:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:), intent(OUT), contiguous :: output

      integer  :: i, j
      real(wp) :: invn

      call fft_r2c_x(input, bufx_cplx)
      
      invn = -pex**2/real(xsz(1),wp)
      do j = 1, xsz(2)
         do i = 0, xsz(1)/2-1
            ! output(i*2+1,j) = output(i*2+1,j)*(i**2)*invn
            ! output(i*2+2,j) = output(i*2+2,j)*(i**2)*invn
            bufx_cplx(i,j) = (i**2)*invn * bufx_cplx(i,j)
         end do
         bufx_cplx(xsz(1)/2,j) = 0
      end do

      call fft_c2r_x(bufx_cplx, output)

   end subroutine pdfxx_out2

   subroutine pdfxx_in3(inout, pex)
      implicit none

      real(wp), intent(IN) :: pex
      real(wp), dimension(:,:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j, k
      real(wp) :: invn

      invn = -pex**2/real(xsz(1),wp)
      do k = 1, size(inout,3)
         call fft_r2c_x(inout(:,:,k), bufx_cplx)
         do j = 1, xsz(2)
            do i = 0, xsz(1)/2-1
               ! inout(i*2+1,j,k) = inout(i*2+1,j,k)*(i**2)*invn
               ! inout(i*2+2,j,k) = inout(i*2+2,j,k)*(i**2)*invn
               bufx_cplx(i,j) = (i**2)*invn * bufx_cplx(i,j)
            end do
            bufx_cplx(xsz(1)/2,j) = 0
         end do
         call fft_c2r_x(bufx_cplx, inout(:,:,k))
      end do

      call fft_c2r_x(inout)

   end subroutine pdfxx_in3

   subroutine pdfxx_in2(inout, pex)
      implicit none

      real(wp), intent(IN) :: pex
      real(wp), dimension(:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j
      real(wp) :: invn

      call fft_r2c_x(inout, bufx_cplx)
      
      invn = -pex**2/real(xsz(1),wp)
      do j = 1, xsz(2)
         do i = 0, xsz(1)/2-1
            ! inout(i*2+1,j) = inout(i*2+1,j)*(i**2)*invn
            ! inout(i*2+2,j) = inout(i*2+2,j)*(i**2)*invn
            bufx_cplx(i,j) = (i**2)*invn * bufx_cplx(i,j)
         end do
         bufx_cplx(xsz(1)/2,j) = 0
      end do

      call fft_c2r_x(bufx_cplx, inout)

   end subroutine pdfxx_in2

   !--------------------------------------------------------
   !> @brief Calculate second order x partial derivative in spectral.
   !
   !> @param[in] input Input data, spectral domain, x-pencil
   !> @param[out] output Output data, spectral domain, x-pencil
   !> @param[in] nz
   !> @param[in] pex   unit wavenumber
   !--------------------------------------------------------
   subroutine pdfxx_(input, output, nz, pex)
      implicit none

      integer,  intent(IN) :: nz
      real(wp), intent(IN) :: pex
      real(wp), dimension(xsz(1),xsz(2),nz) :: input, output

      integer  :: i, j, k
      real(wp) :: invn

      invn = -pex**2/real(xsz(1),wp)
      do k = 1, nz
         do j = 1, xsz(2)
            do i = 0, xsz(1)/2-1
               output(i*2+1,j,k) = input(i*2+1,j,k)*(i**2)*invn
               output(i*2+2,j,k) = input(i*2+2,j,k)*(i**2)*invn
            end do
         end do
      end do

      end subroutine pdfxx_

   !--------------------------------------------------------
   !> @brief Calculate y partial derivative.
   !
   !> @param[in] input Input data, physical domain, y-pencil
   !> @param[out] output Output data, physical domain, y-pencil
   !> @param[in] nz
   !> @param[in] pey   unit wavenumber
   !--------------------------------------------------------
   subroutine pdfy_out3(input, output, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:,:), intent(OUT), contiguous :: output

      integer  :: i, j, k
      real(wp) :: invn

      invn = pey/real(ysz(1),wp)
      do k = 1, size(input,3)
         call fft_r2c_y(input(:,:,k), bufy_cplx)
         do j = 1, ysz(2)
            do i = 0, ysz(1)/2-1
               ! tmp = output(i*2+1,j,k)
               ! output(i*2+1,j,k) = -output(i*2+2,j,k)*i*invn
               ! output(i*2+2,j,k) = tmp*i*invn
               bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
            end do
            bufy_cplx(ysz(1)/2,j) = 0
         end do
         call fft_c2r_y(bufy_cplx, output(:,:,k))
      end do

   end subroutine pdfy_out3

   subroutine pdfy_out2(input, output, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:), intent(OUT), contiguous :: output

      integer  :: i, j
      real(wp) :: invn

      call fft_r2c_y(input, bufy_cplx)

      invn = pey/real(ysz(1),wp)
      do j = 1, ysz(2)
         do i = 0, ysz(1)/2-1
            ! tmp = output(i*2+1,j)
            ! output(i*2+1,j) = -output(i*2+2,j)*i*invn
            ! output(i*2+2,j) = tmp*i*invn
            bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
         end do
         bufy_cplx(ysz(1)/2,j) = 0
      end do

      call fft_c2r_y(bufy_cplx, output)

   end subroutine pdfy_out2

   subroutine pdfy_in3(inout, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j, k
      real(wp) :: invn

      invn = pey/real(ysz(1),wp)
      do k = 1, size(inout,3)
         call fft_r2c_y(inout(:,:,k), bufy_cplx)
         do j = 1, ysz(2)
            do i = 0, ysz(1)/2-1
               ! tmp = inout(i*2+1,j,k)
               ! inout(i*2+1,j,k) = -inout(i*2+2,j,k)*i*invn
               ! inout(i*2+2,j,k) = tmp*i*invn
               bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
            end do
            bufy_cplx(ysz(1)/2,j) = 0
         end do
         call fft_c2r_y(bufy_cplx, inout(:,:,k))
      end do

   end subroutine pdfy_in3

   subroutine pdfy_in2(inout, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j
      real(wp) :: invn

      call fft_r2c_y(inout, bufy_cplx)

      invn = pey/real(ysz(1),wp)
      do j = 1, ysz(2)
         do i = 0, ysz(1)/2-1
            ! tmp = inout(i*2+1,j)
            ! inout(i*2+1,j) = -inout(i*2+2,j)*i*invn
            ! inout(i*2+2,j) = tmp*i*invn
            bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
         end do
         bufy_cplx(ysz(1)/2,j) = 0
      end do

      call fft_c2r_y(bufy_cplx, inout)

   end subroutine pdfy_in2

   !--------------------------------------------------------
   !> @brief Calculate y partial derivative.
   !
   !> @param[in] input Input data, physical domain, x-pencil
   !> @param[out] output Output data, physical domain, x-pencil
   !> @param[in] nz
   !> @param[in] pey   unit wavenumber
   !--------------------------------------------------------
   subroutine pdfy_x_out3(input, output, pey)
     !> hacked by plyu
     use discontinuity_smooth, only : tbn_lim_y
     implicit none
     real(wp), intent(IN) :: pey
     real(wp), dimension(:,:,:), intent(IN),  contiguous :: input
     real(wp), dimension(:,:,:), intent(OUT), contiguous :: output
      
     integer  :: i, j, k
     real(wp) :: invn
     
     if ((ids==1 .or. ids==2) .and. ASSOCIATED(tbn_lim_y)) then     
       call pdfy_x_out3_s(input, output, pey, tbn_lim_y)
     else 
       !> original version: 
       call transpose_xy(input, bufy)
 
       invn = pey/real(ysz(1),wp)
       do k = 1, size(input,3)
          call fft_r2c_y(bufy(:,:,k), bufy_cplx)
          do j = 1, ysz(2)
             do i = 0, ysz(1)/2-1
                ! tmp = bufy(i*2+1,j,k)
                ! bufy(i*2+1,j,k) = -bufy(i*2+2,j,k)*i*invn
                ! bufy(i*2+2,j,k) = tmp*i*invn
                bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
             end do
             bufy_cplx(ysz(1)/2,j) = 0
          end do
          call fft_c2r_y(bufy_cplx, bufy(:,:,k))
       end do

       call transpose_yx(bufy(:,:,1:size(input,3)), output)
     endif

   end subroutine pdfy_x_out3

   subroutine pdfy_x_out2(input, output, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:), intent(OUT), contiguous :: output

      integer  :: i, j
      real(wp) :: invn

      call transpose_xy(input, bufy(:,:,1))

      call fft_r2c_y(bufy(:,:,1), bufy_cplx)

      invn = pey/real(ysz(1),wp)
      do j = 1, ysz(2)
         do i = 0, ysz(1)/2-1
            ! tmp = bufy(i*2+1,j,1)
            ! bufy(i*2+1,j,1) = -bufy(i*2+2,j,1)*i*invn
            ! bufy(i*2+2,j,1) = tmp*i*invn
            bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
         end do
         bufy_cplx(ysz(1)/2,j) = 0
      end do

      call fft_c2r_y(bufy_cplx, bufy(:,:,1))

      call transpose_yx(bufy(:,:,1), output)

   end subroutine pdfy_x_out2

   subroutine pdfy_x_in3(inout, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j, k
      real(wp) :: invn

      call transpose_xy(inout, bufy)

      invn = pey/real(ysz(1),wp)
      do k = 1, size(inout,3)
         call fft_r2c_y(bufy(:,:,k), bufy_cplx)
         do j = 1, ysz(2)
            do i = 0, ysz(1)/2-1
               ! tmp = bufy(i*2+1,j,k)
               ! bufy(i*2+1,j,k) = -bufy(i*2+2,j,k)*i*invn
               ! bufy(i*2+2,j,k) = tmp*i*invn
               bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
            end do
            bufy_cplx(ysz(1)/2,j) = 0
         end do
         call fft_c2r_y(bufy_cplx, bufy(:,:,k))
      end do

      call transpose_yx(bufy(:,:,1:size(inout,3)), inout)

   end subroutine pdfy_x_in3

   subroutine pdfy_x_in2(inout, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j
      real(wp) :: invn

      call transpose_xy(inout, bufy(:,:,1))

      call fft_r2c_y(bufy(:,:,1), bufy_cplx)

      invn = pey/real(ysz(1),wp)
      do j = 1, ysz(2)
         do i = 0, ysz(1)/2-1
            ! tmp = bufy(i*2+1,j,1)
            ! bufy(i*2+1,j,1) = -bufy(i*2+2,j,1)*i*invn
            ! bufy(i*2+2,j,1) = tmp*i*invn
            bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
         end do
         bufy_cplx(ysz(1)/2,j) = 0
      end do

      call fft_c2r_y(bufy_cplx, bufy(:,:,1))

      call transpose_yx(bufy(:,:,1), inout)

   end subroutine pdfy_x_in2
   
   subroutine pdfy_x_out3_s (input, output, pey, lim3d)
      use fll_mods_m
      use discontinuity_smooth, only : cubic_smooth_3d, postproc_smooth_3d, &
        i_grad_p, n_pdfy
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:,:), intent(OUT), contiguous :: output
      type(dnode), pointer, intent(in) :: lim3d

      integer  :: i, j, k
      real(wp) :: invn, dy, twopi_local
      type(dnode), pointer :: ptmp, p_2d_indices, p_2d_index, lim1d

      !> for debug
      integer :: i_dbg, nt0, i0, j0, k0, i_print

      i_dbg = 0  !> set it to 1 if you want to see the process in detail 
      i_print = 0

      if (i_dbg > 0 .and. i_grad_p>0) then
        i_print = 1
        n_pdfy = n_pdfy + 1

        nt0 = 40; i0 = 66; j0=11; k0=17 !> the desired timestep and location
      endif

      if (i_print > 0 .and. n_pdfy == nt0) then
        !> original version: 
        call transpose_xy(input, bufy)
 
        invn = pey/real(ysz(1),wp)
        do k = 1, size(input,3)
           call fft_r2c_y(bufy(:,:,k), bufy_cplx)
           do j = 1, ysz(2)
              do i = 0, ysz(1)/2-1
                 ! tmp = bufy(i*2+1,j,k)
                 ! bufy(i*2+1,j,k) = -bufy(i*2+2,j,k)*i*invn
                 ! bufy(i*2+2,j,k) = tmp*i*invn
                 bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
              end do
              bufy_cplx(ysz(1)/2,j) = 0
           end do
           call fft_c2r_y(bufy_cplx, bufy(:,:,k))
        end do

        !> print result of original version
        if (myid.eq.0) print *, 'Times of calling pdfy_x_debug:', n_pdfy, ', original dy/dx'
        if (yst(2)<=i0 .and. (yst(2)+ysz(2)-1)>=i0 .and. yst(3)<=k0 &
          .and. (yst(3)+ysz(3)-1)>=k0) then
          open(1050, file='data_du_pdf_orig.dat', action='write')
          print *, 'VARIABLS = Y, dU'            
          do j = 35, 75
            print *, j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
            write(1050, *) j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
          enddo
          close(1050)
        endif
        
        !call transpose_yx(bufy(:,:,1:size(input,3)), output)

      endif

      call transpose_xy(input, bufy)
        
      if (i_print > 0 .and. n_pdfy .eq. nt0) then
        !> write the input before smoothing 
        if (myid.eq.0) print *, 'Times of calling pdfy_x_debug:', n_pdfy, ', before pdfy'
        if (yst(2)<=i0 .and. (yst(2)+ysz(2)-1)>=i0 .and. yst(3)<=k0 &
          .and. (yst(3)+ysz(3)-1)>=k0) then
          open(1051, file='data_u_orig.dat', action='write')
          print *, 'VARIABLS = Y, dU'            
          do j = 35, 75
            print *, j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
            write(1051, *) j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
          enddo
          close(1051)
        endif
      endif

      bufy2(:,:,:) = bufy(:,:,:) !> back up the original value
      call cubic_smooth_3d (bufy, lim3d)

      if (i_print > 0 .and. n_pdfy .eq. nt0) then
        !> write the input after smoothing 
        if (myid.eq.0) print *, 'Times of calling pdfy_x_debug:', n_pdfy, ', after smooth'
        if (yst(2)<=i0 .and. (yst(2)+ysz(2)-1)>=i0 .and. yst(3)<=k0 &
          .and. (yst(3)+ysz(3)-1)>=k0) then
          open(1052, file='data_u_smooth.dat', action='write')
          print *, 'VARIABLS = Y, dU'            
          do j = 35, 75
            print *, j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
            write(1052, *) j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
          enddo
          close(1052)
        endif
      endif
      
      invn = pey/real(ysz(1),wp)
      do k = 1, size(input,3)
         call fft_r2c_y(bufy(:,:,k), bufy_cplx)
         do j = 1, ysz(2)
            do i = 0, ysz(1)/2-1
               ! tmp = bufy(i*2+1,j,k)
               ! bufy(i*2+1,j,k) = -bufy(i*2+2,j,k)*i*invn
               ! bufy(i*2+2,j,k) = tmp*i*invn
               bufy_cplx(i,j) = cmplx(0, i*invn, wp) * bufy_cplx(i,j)
            end do
            bufy_cplx(ysz(1)/2,j) = 0
         end do
         call fft_c2r_y(bufy_cplx, bufy(:,:,k))
      end do

      if (i_print > 0 .and. n_pdfy .eq. nt0) then
        !> write the input after smoothing 
        if (myid.eq.0) print *, 'Times of calling pdfy_x_debug:', n_pdfy, ', after pdfy'
        if (yst(2)<=i0 .and. (yst(2)+ysz(2)-1)>=i0 .and. yst(3)<=k0 &
          .and. (yst(3)+ysz(3)-1)>=k0) then
          open(1053, file='data_du_pdf.dat', action='write')
          print *, 'VARIABLS = Y, dU'            
          do j = 35, 75
            print *, j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
            write(1053, *) j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
          enddo
          close(1053)
        endif
      endif
      
      if (ids==1) then
        !> do 2nd order central difference for the smoothed part
        twopi_local = 3.141592653589793238_wp * 2.0_wp
        dy = twopi_local / pey / real(ysz(1), wp)
        call postproc_smooth_3d(bufy2, bufy, lim3d, 1, dy)
        
        if (i_print > 0 .and. n_pdfy .eq. nt0) then
          if (myid.eq.0) print *, 'Times of calling pdfy_x_debug:', n_pdfy, ', after fdmy'
          if (yst(2)<=i0 .and. (yst(2)+ysz(2)-1)>=i0 .and. yst(3)<=k0 &
            .and. (yst(3)+ysz(3)-1)>=k0) then
            open(1054, file='data_du_fdm.dat', action='write')
            print *, 'VARIABLS = Y, dU'            
            do j = 35, 75
              print *, j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
              write(1054, *) j, bufy(j, i0-yst(2)+1, k0-yst(3)+1)
            enddo
            close(1054)
          endif
        endif
      endif

      !> ignore the precision loss if ids==2

      call transpose_yx(bufy(:,:,1:size(input,3)), output)

   end subroutine pdfy_x_out3_s

   !--------------------------------------------------------
   !> @brief Calculate y partial derivative in spectral.
   !
   !> @param[in] input Input data, spectral domain, y-pencil
   !> @param[out] output Output data, spectral domain, y-pencil
   !> @param[in] nz
   !> @param[in] pey   unit wavenumber
   !--------------------------------------------------------
   subroutine pdfy_(input, output, nz, pey)
      implicit none

      integer,  intent(IN) :: nz
      real(wp), intent(IN) :: pey
      real(wp), dimension(ysz(1),ysz(2),nz) :: input,output

      integer  :: i, j, k
      real(wp) :: tmp, invn

      invn = pey/real(ysz(1),wp)
      do k = 1, nz
         do j = 1, ysz(2)
            do i = 0, ysz(1)/2-1
               tmp = input(i*2+1,j,k)
               output(i*2+1,j,k) = -input(i*2+2,j,k)*i*invn
               output(i*2+2,j,k) = tmp*i*invn
            end do
         end do
      end do

   end subroutine pdfy_

   !--------------------------------------------------------
   !> @brief Calculate second-order y partial derivative.
   !
   !> @param[in] input Input data, physical domain, y-pencil
   !> @param[out] output Output data, physical domain, y-pencil
   !> @param[in] nz
   !> @param[in] pey   unit wavenumber
   !--------------------------------------------------------
   subroutine pdfyy_out3(input, output, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:,:), intent(OUT), contiguous :: output

      integer  :: i, j, k
      real(wp) :: invn

      invn = -pey**2/real(ysz(1),wp)
      do k = 1, size(input,3)
         call fft_r2c_y(input(:,:,k), bufy_cplx)
         do j = 1, ysz(2)
            do i = 0, ysz(1)/2-1
               ! output(i*2+1,j,k) = output(i*2+1,j,k)*(i**2)*invn
               ! output(i*2+2,j,k) = output(i*2+2,j,k)*(i**2)*invn
               bufy_cplx(i,j) = (i**2)*invn * bufy_cplx(i,j)
            end do
            bufy_cplx(ysz(1)/2,j) = 0
         end do
         call fft_c2r_y(bufy_cplx, output(:,:,k))
      end do

   end subroutine pdfyy_out3

   subroutine pdfyy_out2(input, output, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:), intent(IN),  contiguous :: input
      real(wp), dimension(:,:), intent(OUT), contiguous :: output

      integer  :: i, j
      real(wp) :: invn

      call fft_r2c_y(input, bufy_cplx)
      
      invn = -pey**2/real(ysz(1),wp)
      do j = 1, ysz(2)
         do i = 0, ysz(1)/2-1
            ! output(i*2+1,j) = output(i*2+1,j)*(i**2)*invn
            ! output(i*2+2,j) = output(i*2+2,j)*(i**2)*invn
            bufy_cplx(i,j) = (i**2)*invn * bufy_cplx(i,j)
         end do
         bufy_cplx(ysz(1)/2,j) = 0
      end do

      call fft_c2r_y(bufy_cplx, output)

   end subroutine pdfyy_out2

   subroutine pdfyy_in3(inout, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j, k
      real(wp) :: invn

      invn = -pey**2/real(ysz(1),wp)
      do k = 1, size(inout,3)
         call fft_r2c_y(inout(:,:,k), bufy_cplx)
         do j = 1, ysz(2)
            do i = 0, ysz(1)/2-1
               ! inout(i*2+1,j,k) = inout(i*2+1,j,k)*(i**2)*invn
               ! inout(i*2+2,j,k) = inout(i*2+2,j,k)*(i**2)*invn
               bufy_cplx(i,j) = (i**2)*invn * bufy_cplx(i,j)
            end do
            bufy_cplx(ysz(1)/2,j) = 0
         end do
         call fft_c2r_y(bufy_cplx, inout(:,:,k))
      end do

   end subroutine pdfyy_in3

   subroutine pdfyy_in2(inout, pey)
      implicit none

      real(wp), intent(IN) :: pey
      real(wp), dimension(:,:), intent(INOUT), contiguous :: inout

      integer  :: i, j
      real(wp) :: invn

      call fft_r2c_y(inout, bufy_cplx)
      
      invn = -pey**2/real(ysz(1),wp)
      do j = 1, ysz(2)
         do i = 0, ysz(1)/2-1
            ! inout(i*2+1,j) = inout(i*2+1,j)*(i**2)*invn
            ! inout(i*2+2,j) = inout(i*2+2,j)*(i**2)*invn
            bufy_cplx(i,j) = (i**2)*invn * bufy_cplx(i,j)
         end do
         bufy_cplx(ysz(1)/2,j) = 0
      end do

      call fft_c2r_y(bufy_cplx, inout)

   end subroutine pdfyy_in2

   !--------------------------------------------------------
   !> @brief Calculate second-order y partial derivative in spectral.
   !
   !> @param[in] input Input data, spectral domain, y-pencil
   !> @param[out] output Output data, spectral domain, y-pencil
   !> @param[in] nz
   !> @param[in] pey   unit wavenumber
   !--------------------------------------------------------
   subroutine pdfyy_(input, output, nz, pey)
      implicit none

      integer,  intent(IN) :: nz
      real(wp), intent(IN) :: pey
      real(wp), dimension(ysz(1),ysz(2),nz) :: input, output

      integer  :: i, j, k
      real(wp) :: invn

      invn = -pey**2/real(ysz(1),wp)
      do k = 1, nz
         do j = 1, ysz(2)
            do i = 0, ysz(1)/2-1
               output(i*2+1,j,k) = input(i*2+1,j,k)*(i**2)*invn
               output(i*2+2,j,k) = input(i*2+2,j,k)*(i**2)*invn
            end do
         end do
      end do

   end subroutine pdfyy_


   !--------------------------------------------------------
   !> @brief Dealias data
   !
   !> @param[in,out] input Data to be dealiased, physical 
   !!                domain, x-pencil. The input data is 
   !!                overwritten by dealiased one.
   !> @param[in] nz
   !--------------------------------------------------------
   subroutine dealiasxy3(input)
      !use mpi
      use discontinuity_smooth, only : tbn_lim_x, tbn_lim_y
      implicit none

      real(wp), dimension(:,:,:), intent(INOUT) :: input

      integer :: k, kc_x, kc_y, ierr_wt

      if (idsd==1 .and. ASSOCIATED(tbn_lim_x) .and. ASSOCIATED(tbn_lim_y)) then
        call dealiasxy3_s (input, tbn_lim_x, tbn_lim_y)
      elseif (idsd==2) then
        call dealiasxy3_s32 (input)
      else    
        ! return
        !print *, "dexy3, 1"
        kc_x = ((xsz(1)/2)*2/3)*2+1
        kc_y = ((ysz(1)/2)*2/3)
        !print *, "dexy3, 2"

        call fft_r2c_x(input)
        input(kc_x:,:,:) = 0
        !call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
        !print *, "dexy3, 3"
        !print *, "input(1:5,1:5)=", input(1:5,1:5,2)
        call transpose_xy(input, bufy)
        !print *, "bufy(1:5, 1:5)=", bufy(1:5, 1:5, 2)
        !call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
        !print *, "dexy3, 4"

        do k=1, size(input, 3)
           call fft_r2c_y(bufy(:,:,k), bufy_cplx)

           bufy_cplx(kc_y:,:) = 0

           call fft_c2r_y(bufy_cplx, bufy(:,:,k))
        end do
        !call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
        !print *, "dexy3, 5"

        call transpose_yx(bufy(:,:,1:size(input,3)), input)
        !call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
        !print *, "dexy3, 6"

        call fft_c2r_x(input)
        !print *, "dexy3, 7"

        input = input/nx_global/ny_global
      endif

   end subroutine dealiasxy3
   
   subroutine dealiasxy3_s32(input)
      !use mpi
      implicit none

      real(wp), dimension(:,:,:), intent(INOUT) :: input

      integer :: i, k, kc_x, kc_y, ierr_wt

        kc_x = ((xsz(1)/2)*2/3)*2+1
        kc_y = ((ysz(1)/2)*2/3)

        call fft_r2c_x(input)

        do i = 1, xsz(1)/2
          input((2*i-1):(2*i),:,:) = co_s32_x(i) * input((2*i-1):(2*i),:,:)
        enddo

        call transpose_xy(input, bufy)

        do k=1, size(input, 3)
           call fft_r2c_y(bufy(:,:,k), bufy_cplx)

           do i = 1, ysz(1)/2
             bufy_cplx(i,:) = co_s32_y(i) * bufy_cplx(i,:)
           enddo

           call fft_c2r_y(bufy_cplx, bufy(:,:,k))
        end do

        call transpose_yx(bufy(:,:,1:size(input,3)), input)

        call fft_c2r_x(input)

        input = input/nx_global/ny_global

   end subroutine dealiasxy3_s32
   
   subroutine dealiasxy3_s(input, lim3d_x, lim3d_y)
      !use mpi
      use fll_mods_m
      use discontinuity_smooth, only : cubic_smooth_3d, postproc_smooth_3d
      implicit none

      real(wp), dimension(:,:,:), intent(INOUT) :: input
      type(dnode), pointer, intent(in) :: lim3d_x, lim3d_y

      integer :: k, kc_x, kc_y, ierr_wt

      !> at the start, do smoothing for x direction
      bufx(:,:,:) = input(:,:,:)
      call cubic_smooth_3d (input, lim3d_x)
      !call transpose_xy(input, bufy)
      !call cubic_smooth_3d (bufy, lim3d_y)
      !call transpose_yx(bufy, input)

      kc_x = ((xsz(1)/2)*2/3)*2+1
      kc_y = ((ysz(1)/2)*2/3)

      call fft_r2c_x(input)
      input(kc_x:,:,:) = 0
      call transpose_xy(input, bufy)

      !> do smoothing for y direction
      bufy2(:,:,:) = bufy(:,:,:)
      call cubic_smooth_3d (bufy, lim3d_y)

      do k=1, size(input, 3)
         call fft_r2c_y(bufy(:,:,k), bufy_cplx)

         bufy_cplx(kc_y:,:) = 0

         call fft_c2r_y(bufy_cplx, bufy(:,:,k))
      end do

      !> after the dealias in y, do recovering for y
      call postproc_smooth_3d (bufy2, bufy, lim3d_y, 0)

      call transpose_yx(bufy(:,:,1:size(input,3)), input)

      call fft_c2r_x(input)

      input = input/nx_global/ny_global

      !> do recovering for x
      call postproc_smooth_3d (bufx, input, lim3d_x, 0) 

   end subroutine dealiasxy3_s

   subroutine dealiasxy2(input)
      implicit none

      real(wp), dimension(:,:), intent(INOUT) :: input

      integer :: kc_x, kc_y

      ! return

      kc_x = ((xsz(1)/2)*2/3)*2+1
      kc_y = ((ysz(1)/2)*2/3)

      call fft_r2c_x(input)
      input(kc_x:,:) = 0

      call transpose_xy(input, bufy(:,:,1))

      call fft_r2c_y(bufy(:,:,1), bufy_cplx)

      bufy_cplx(kc_y:,:) = 0

      call fft_c2r_y(bufy_cplx, bufy(:,:,1))

      call transpose_yx(bufy(:,:,1), input)

      call fft_c2r_x(input)

      input = input/nx_global/ny_global

   end subroutine dealiasxy2

   subroutine dealiasxy3_2(input)
      implicit none

      real(wp), dimension(:,:,:), intent(INOUT) :: input

      integer :: k, kc_x, kc_y

      kc_x = ((xsz(1)/2)*2/3)*2+1
      kc_y = ((ysz(1)/2)*2/3)

      call fft_r2c_x(input)
      input(kc_x:,:,:) = 0

      call transpose_xy(input, bufy)

      do k=1, size(input, 3)
         call fft_r2c_y(bufy(:,:,k), bufy_cplx)

         bufy_cplx(kc_y:,:) = 0

         call fft_c2r_y(bufy_cplx, bufy(:,:,k))
      end do

      call transpose_yx(bufy(:,:,1:size(input,3)), input)

      call fft_c2r_x(input)

      input = input/nx_global/ny_global

   end subroutine dealiasxy3_2

   subroutine dealiasxy2_2(input)
      implicit none

      real(wp), dimension(:,:), intent(INOUT) :: input

      integer :: kc_x, kc_y

      kc_x = ((xsz(1)/2)*2/3)*2+1
      kc_y = ((ysz(1)/2)*2/3)

      call fft_r2c_x(input)
      input(kc_x:,:) = 0

      call transpose_xy(input, bufy(:,:,1))

      call fft_r2c_y(bufy(:,:,1), bufy_cplx)

      bufy_cplx(kc_y:,:) = 0

      call fft_c2r_y(bufy_cplx, bufy(:,:,1))

      call transpose_yx(bufy(:,:,1), input)

      call fft_c2r_x(input)

      input = input/nx_global/ny_global

   end subroutine dealiasxy2_2
   !--------------------------------------------------------
   !> @brief Dealias data only in x direction
   !
   !> @param[in,out] input Data to be dealiased, physical 
   !!                domain, x-pencil. The input data is 
   !!                overwritten by dealiased one.
   !> @param[in] nz
   !--------------------------------------------------------
   subroutine dealiasx(input, nz)
      implicit none

      integer, intent(IN) :: nz
      real(wp), dimension(xsz(1),xsz(2),nz) :: input

      integer :: kc_x, kc_y

      return

      kc_x = ((xsz(1)/2)*2/3)*2+1
      kc_y = ((ysz(1)/2)*2/3)*2+1

      call fft_r2c_x(input)

      input(kc_x:,:,:) = 0
      input = input/nx_global

      call fft_c2r_x(input)

   end subroutine dealiasx

   subroutine dealiasx_2(input, nz)
      implicit none

      integer, intent(IN) :: nz
      real(wp), dimension(xsz(1),xsz(2),nz) :: input

      integer :: kc_x, kc_y

      kc_x = ((xsz(1)/2)*2/3)*2+1
      kc_y = ((ysz(1)/2)*2/3)*2+1

      call fft_r2c_x(input)

      input(kc_x:,:,:) = 0
      input = input/nx_global

      call fft_c2r_x(input)

   end subroutine dealiasx_2

   !------------------------------------------------
   ! Phase shift a 2D data in X direction
   ! Suppose the input is f(x), the output satisfies
   !   g(x)=f(x-theta)
   !------------------------------------------------
   subroutine phase_shift_x_2d(input, output, theta)
      implicit none

      real(wp), dimension(:,:), intent(IN) :: input
      real(wp), dimension(:,:), intent(OUT) :: output
      real(wp), intent(IN) :: theta

      integer :: i, j

      call fft_r2c_x(input, bufx_cplx(:,:))

      do j=1, size(input,2)
         do i = 1, xsz(1)/2-1
            bufx_cplx(i,j) = bufx_cplx(i,j)*cmplx(cos(i*theta),-sin(i*theta),wp)
         end do
         bufx_cplx(xsz(1)/2,j) = 0
      end do
      bufx_cplx = bufx_cplx/size(input,1)

      call fft_c2r_x(bufx_cplx(:,:), output)

   end subroutine phase_shift_x_2d

end module spectral
