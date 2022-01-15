!all checked
!input&output checked
module grid
   use decomp
   use constants
   use param, only : pex, pey, hbar, xl, yl, zl, dx, dy, nz, clbeta,&
        clgama, timewavy, iwavy, tcoef
   
   implicit none

   real(wp), allocatable, dimension(:)   :: zz, dz, zw, dzw
   real(wp), allocatable, dimension(:,:) :: hh, ht, eta, eta0
   real(wp), allocatable, dimension(:,:) :: ex, exx, ey, eyy, exy
   real(wp), allocatable, dimension(:,:) :: hx, hxx, hy, hyy, hxy
   real(wp), allocatable, dimension(:,:) :: ehx, ehx2, ehy, ehy2
   real(wp), allocatable, dimension(:,:) :: reh, reh2, her, her2
   real(wp), allocatable, dimension(:,:,:) :: zetax, zetay, zetat
   real(wp), allocatable, dimension(:,:,:) :: zetax_w, zetay_w, zetat_w

   real(wp), allocatable, dimension(:,:) :: et


   real(wp), allocatable, dimension(:,:) :: exr, eyr, hxr, hyr

   integer :: level

   public :: grid_init
   public :: nl_coef, coef_et
   public :: grid_gen
   !public :: bottom2
   !public :: bottom_hos_les

contains

   !--------------------------------------------------------
   !> @brief Allocate and initialize variables in module grid.
   !--------------------------------------------------------

  !checked
  subroutine grid_init(level_)
      implicit none

      integer, intent(IN) :: level_

      level = level_
      !print *, 'plyudebug, xsz:', xsz(1), xsz(2) 
      allocate(zz(1-level:xsz(3)+level), zw(1-level:xsz(3)+level))
      allocate(dz(1-level:xsz(3)+level),dzw(1-level:xsz(3)+level))
      allocate(eta(xsz(1),xsz(2)));  allocate(eta0(xsz(1),xsz(2)))
      allocate(hh(xsz(1),xsz(2)));   allocate(ht(xsz(1),xsz(2)))
      allocate(ex(xsz(1),xsz(2)));   allocate(exx(xsz(1),xsz(2)))
      allocate(ey(xsz(1),xsz(2)));   allocate(eyy(xsz(1),xsz(2)))
      allocate(hx(xsz(1),xsz(2)));   allocate(hxx(xsz(1),xsz(2)))
      allocate(hy(xsz(1),xsz(2)));   allocate(hyy(xsz(1),xsz(2)))
      allocate(exy(xsz(1),xsz(2)));  allocate(hxy(xsz(1),xsz(2)))
      allocate(ehx(xsz(1),xsz(2)));  allocate(ehx2(xsz(1),xsz(2)))
      allocate(ehy(xsz(1),xsz(2)));  allocate(ehy2(xsz(1),xsz(2)))
      allocate(reh(xsz(1),xsz(2)));  allocate(reh2(xsz(1),xsz(2)))
      allocate(her(xsz(1),xsz(2)));  allocate(her2(xsz(1),xsz(2)))
      allocate(exr(xsz(1),xsz(2)));  allocate(eyr(xsz(1),xsz(2)))
      allocate(hxr(xsz(1),xsz(2)));  allocate(hyr(xsz(1),xsz(2)))
      allocate(zetax(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(zetay(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(zetat(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(zetax_w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(zetay_w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(zetat_w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(et(xsz(1),xsz(2)))

      hh = 0; ht=0; hx = 0; hxx = 0; hy = 0; hyy = 0; hxy = 0

   end subroutine grid_init

   !--------------------------------------------------------
   !> @brief Calculate the grid transformation coefficients.
   !
   !> The coefficients to calculate are the following: \n
   !> ex(\f$\eta'_x\f$), exx(\f$\eta'_{xx}\f$), 
   !! ey(\f$\eta'_y\f$), eyy(\f$\eta'_{yy}\f$), exy(\f$\eta'_{xy}\f$), \n
   !> ehx(\f$\eta_x\f$), ehx2(\f$\eta_x^2\f$), 
   !! ehy(\f$\eta_y\f$), ehy2(\f$\eta_y^2\f$), \n
   !> reh(\f$1/(1+\eta_x^2+\eta_y^2)\f$), reh2(\f$1/(1+\eta_x^2+\eta_y^2)^2\f$),
   !! her(\f$1/(\eta'+\bar{H})\f$), her2(\f$1/(\eta'+\bar{H})^2\f$), \n
   !> zetax(\f$\zeta_x\f$), zetay(\f$\zeta_y\f$) at u,v nodes,
   !! zetax_w(\f$\zeta_x\f$), zetay_w(\f$\zeta_y\f$) at w nodes.
   !--------------------------------------------------------

   !checked
   !------------------------------------------------------------
   !input: pex, pey, hbar, zz
   !input: eta, hx, hy, 
   !output: ex, exx, ey, eyy, exy
   !output: ehx, ehy, ehx2, ehy2, reh, reh2, her, her2
   !output: exr, eyr, hxr, hyr, zetax, zetay, zetax_w, zetay_w
   !------------------------------------------------------------
   subroutine nl_coef
      use fft
      use spectral

      implicit none

      real(wp), dimension(ysz(1),ysz(2)) :: tmpy1, tmpy2

      integer :: k

      !print *, "nl, 1"
      call dealiasxy(eta)

      call fft_r2c_x(eta,ex)
      exx = ex
      !print *, "eta=",eta(1:5,1), ", ex=", ex(1:5,1)
      
      ! compute eta_x and eta_xx
      call pdfx_(ex,ex,1,pex)
      call pdfxx_(exx,exx,1,pex)
      call fft_c2r_x(ex)
      call fft_c2r_x(exx)

      !print *, "nl, 2"
      call transpose_xy(eta,tmpy1)
      !tmpy1 saves the eta in spectral domain in y-pencil
      call fft_r2c_y(tmpy1)
      ! compute eta_y
      call pdfy_(tmpy1,tmpy2,1,pey)
      call fft_c2r_y(tmpy2)
      call transpose_yx(tmpy2,ey)
      ! compute eta_yy
      call pdfyy_(tmpy1,tmpy2,1,pey)
      call fft_c2r_y(tmpy2)
      call transpose_yx(tmpy2,eyy)

      ! compute eta_xy
      !print *, "nl, 3"
      call transpose_xy(ex,tmpy2)
      call pdfy(tmpy2,pey)
      call transpose_yx(tmpy2,exy)

      !print *, "nl, 4"
      ehx = ex-hx
      ehy = ey-hy
      ehx2 = ehx**2
      ehy2 = ehy**2
      call dealiasxy(ehx2)
      call dealiasxy(ehy2)

      !print *, "nl, 5"
      !print *, "ehx2=", ehx2(1:5,1), "ehy2=", ehy2(1:5,1)
      reh = 1/(1+ehx2+ehy2)
      call dealiasxy(reh)
      !print *, "nl, 5.1"
      reh2 = reh**2
      call dealiasxy(reh2)

      !print *, "nl, 5.2"
      her = 1/(eta+hbar)
      call dealiasxy(her)
      !print *, "nl, 5.3"
      her2 = her**2
      call dealiasxy(her2)

      exr = ex*her
      eyr = ey*her
      hxr = hx*her
      hyr = hy*her
      ! call dealiasxy(exr)
      ! call dealiasxy(eyr)
      ! call dealiasxy(hxr)
      ! call dealiasxy(hyr)

      !print *, "nl, 6"
      do k=lbound(zetax,3),ubound(zetax,3)
         zetax(:,:,k) = hxr(:,:)-zz(k)*exr(:,:)
         zetay(:,:,k) = hyr(:,:)-zz(k)*eyr(:,:)
      end do
      do k=lbound(zetax_w,3),ubound(zetax_w,3)
         zetax_w(:,:,k) = hxr(:,:)-zw(k)*exr(:,:)
         zetay_w(:,:,k) = hyr(:,:)-zw(k)*eyr(:,:)
      end do

      ! if(myid==0) print*,'hx', hx(:,1)
      ! if(myid==0) print*,'eta', eta(:,1)
      ! if(myid==0) print*,'hxr', hxr(:,1)

   end subroutine nl_coef

   !--------------------------------------------------------
   !> @brief Calculate coeffecient \f$\zeta_t\f$.
   !
   !> @param[in] utop,vtop,wtop velocity at the surface
   !--------------------------------------------------------

   !checked
   !--------------------------------------------------------
   !input: ht, her, zz, zw
   !output: et(set zero), zetat, zetat_w
   !--------------------------------------------------------
   subroutine coef_et
     use spectral
     use mpi
     
     implicit none
     
     integer :: k
     
     real(wp), allocatable, dimension(:,:) :: tmp
     
     et = 0
     allocate(tmp(xsz(1),xsz(2)))
     do k=lbound(zetat,3),ubound(zetat,3)
        tmp = (1-zz(k))*ht-zz(k)*et
        call dealiasxy(tmp)
        zetat(:,:,k) = tmp*her 
        call dealiasxy(zetat(:,:,k))
     end do

     do k=lbound(zetat_w,3),ubound(zetat_w,3)
        tmp = (1-zw(k))*ht-zw(k)*et
        call dealiasxy(tmp)
        zetat_w(:,:,k) = tmp*her 
        call dealiasxy(zetat_w(:,:,k))
     end do
     deallocate(tmp)

   end subroutine coef_et
   
   !checked
   !--------------------------------------------------------
   !input: nz
   !output: zz, dz, zw, dzw
   !--------------------------------------------------------
   subroutine grid_gen
     
     implicit none
     integer  :: k, k1, k2
     real(wp), allocatable, dimension(:) :: zz1tmp, zz2tmp, zztmp, zwtmp, dztmp, dzwtmp
     
     allocate(zz1tmp(nz+1), zz2tmp(nz+1), zztmp(nz+1), zwtmp(nz+1), dztmp(nz+1), dzwtmp(nz+1))
     
     !> original unifrom grid    
     !call cluster_nl_1(zz1tmp, zz2tmp, zztmp, zwtmp, dztmp)
     !> modified by plyu for turbine cases.
     call cluster_nl_3(zz1tmp, zz2tmp, zztmp, zwtmp, dztmp)

     
     !     call cluster_nl(nz_global,zl,clbeta,clgama,zz1tmp,zz2tmp,zztmp,zwtmp,dztmp)
     !     call cluster_nl_2(nz_global,zl,clbeta,clgama,zz1tmp,zz2tmp,zztmp,zwtmp,dztmp)     
     if (isbot) then
        k1 = 1; k2 = ubound(zz,1)
     else if (istop) then
        k1 = lbound(zz,1); k2 = xsz(3)+1
     else
        k1 = lbound(zz,1); k2 = ubound(zz,1)
     end if
     do k=k1, k2
        zz(k) = zztmp(xst(3)+k-1)
        dz(k) = dztmp(xst(3)+k-1)
        zw(k) = zwtmp(xst(3)+k-1)
     end do
     
     do k=1, nz-1
        dztmp(k) = zwtmp(k+1) - zwtmp(k)
     enddo
     if (isbot) then 
        k1 = 1; k2 = ubound(zz,1)
     else if (istop) then 
        k1 = lbound(zz,1); k2 = xsz(3)-1
     else 
        k1 = lbound(zz,1); k2 = ubound(zz,1)
     endif
     do k=k1, k2
        dzw(k) = dztmp(xst(3)+k-1)
     enddo
     if (istop) dzw(xsz(3))=dz(xsz(3)) 
     if (isbot) dzw(0) = dz(1)/2.0_wp
     
   end subroutine grid_gen
   
   !checked
   !--------------------------------------------------------
   !input: nz, clbeta, zl
   !output: zz1, zz2, zz, zw, dz
   !--------------------------------------------------------
   subroutine cluster_nl_1(ZZ1,ZZ2,ZZ,ZW,DZ)
     
     !     CLGAMA: LARGEST GRID RATION
     !     CLBETA: CLUSTERING PAPAMETER: SMALLER VALUE -> MORE CLUSTERED
     
     use param, only: nx, ny 
     implicit none
     
     real(wp), dimension(*), intent(out) :: zz1,zz2,zz,zw,dz
     
     integer  :: neven=4, k, i, j
     real(wp) :: a, etak
     
     a = 0.5
     do k=neven+1,nz-neven-1
        etak=(k-neven)/(nz-2.0_wp*neven) / 2.
        zz1(k)=0.5*(1.+sinh(clbeta*(etak-a))/sinh(clbeta*a))
     end do
     zz1(neven)=0.
     zz1(nz-neven)=0.5
     
     do k=neven,nz-neven
        dz(k)=1.
     end do
     
     do k=1,neven-1
        dz(k)=dz(neven)
     end do
     
     do k=nz-neven+1,nz
        dz(k)=dz(nz-neven)
     end do
     
     dz(1)=dz(1)/2.
     dz(nz-1)=dz(nz-1)/2.
     dz(nz)=dz(nz)/2.
     
     zz2(1)=0.
     do k=1,nz
        zz2(k+1)=zz2(k)+dz(k)
     end do
     
     do k=1,nz+1
        zz(k)=(zz2(k)-zz2(1))/(zz2(nz)-zz2(1))*zl
     end do
     
     do k=1,nz
        dz(k)=zz(k+1)-zz(k)
     end do
     
     do k=1,nz
        if(k .eq. 1)then
           zw(k)=zz(k)
        else if(k .eq. nz-1) then
           zw(k)=zz(k+1)
        else if(k .eq. nz) then
           zw(k)=zz(k+1)+dz(k)
        else
           zw(k)=(zz(k)+zz(k+1))/2.
        endif
     enddo
     zw(nz+1) = 2.0*zw(nz) - zw(nz-1)

     if (myid == 0) then
        open(40)
        write(40,499)
499     format('variables=k,zz, dz, ratio, zw') ! modified by plyu
        
        do k=nz,1,-1
           write(40,599)float(k), zz(k), dz(k), dz(k+1)/dz(k), zw(k) ! modified by plyu
        enddo
599     format(5e24.14)
        close(40)
     end if

     ! added by plyu: export xyz.dat for inflow generation
     if (myid ==0) then
       open(41, file='xyz.dat', action='write')
       write(41,*) ny, nz, nx
       do j = 1, ny
         write(41,*) yl/ny*j, zz(1), xl/nx*1
       enddo
       do k = 1, nz
         write(41,*) yl/ny*1, zz(k), xl/nx*1
       enddo
       do i = 1, nx
         write(41,*) yl/ny*1, zz(1), xl/nx*i
       enddo 
       close(41)
     endif

     
   END SUBROUTINE  cluster_nl_1 
   
   
   ! added by plyu to refine grid from bottom to turbine blade tip
   !--------------------------------------------------------
   !input: nz, nz1, zl, L1, rat
   !output: zz1, zz2, zz, zw, dz
   !--------------------------------------------------------
   subroutine cluster_nl_3(ZZ1,ZZ2,ZZ,ZW,DZ)
     
     !     nz1: grid number from bottom to blade tip
     !     L1: length from bottom to blade tip
     !     rat: grading ratio from blade tip to top
     
     use mpi
     use param, only: nx, ny 
     implicit none
     
     real(wp), dimension(*), intent(out) :: zz1,zz2,zz,zw,dz
     
     integer  :: neven=2, k, i, j
     real(wp) :: a, etak

     integer :: nz1, nz2, nz3, ierr
     real(wp) :: L1, L2, L3, rat1, rat2, rat3, q1, q2, q3, a1, a2, a3

     if (myid .eq. 0) then
       open(15010, file='z_grid.inp', action='read', status='old')
       read(15010,*) nz1, L1, rat1
       read(15010,*) nz2, L2, rat2
       read(15010,*) nz3, L3, rat3
       close(15010)
     endif

     if ((nz1+nz2+nz3).ne.nz .or. dabs(L1+L2+L3-ZL)>1.0e-6) then
       if(myid.eq.0) print *, 'grid number or block size incorrect in z_grid.inp'
     endif

     nz3 = nz - nz1 - nz2
     L3 = zl - L1 - L2

     call mpi_bcast(nz1, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr)
     call mpi_bcast(L1, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr)
     call mpi_bcast(rat1, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr)
     call mpi_bcast(nz2, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr)
     call mpi_bcast(L2, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr)
     call mpi_bcast(rat2, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr)
     call mpi_bcast(nz3, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr)
     call mpi_bcast(L3, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr)
     call mpi_bcast(rat3, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr)
     
     !if(myid.eq.0) print *,'plyudebug,', nz1, L1, rat
     
     !nz2 = nz - nz1
     !L2 = zl - L1

     !do k = 1, nz1
     !  dz(k) = L1 / nz1;! print *,'k,dz=',k,dz(k)
     !enddo

     q1 = rat1**(1.0_wp / (nz1 - 1 - neven))
     a1 = L1 / ((1.0_wp-q1**(nz1-neven))/(1.0_wp-q1)+2.0_wp)
     do k = 1, neven
       dz(k) = a1
     enddo
     do k = 1, nz1-neven
       dz(neven+k) = a1 * q1**(k-1);! print *,'k,dz=',k,dz(k)
     enddo
     
     q2 = rat2**(1.0_wp / (nz2 - 1))
     a2 = L2 /((1.0_wp-q2**nz2)/(1.0_wp-q2))
     do k = 1, nz2
       dz(nz1+k) = a2 * q2**(k-1);! print *,'k,dz=',k,dz(k)
     enddo
     
     q3 = rat3**(1.0_wp / (nz3 - 1 - neven))
     a3 = L3 / ((1.0_wp-q3**(nz3-neven))/(1.0_wp-q3)+2.0_wp*q3**(nz3-neven-1))
     do k = 1, nz3-neven
       dz(k+nz1+nz2) = a3 * q3**(k-1);! print *,'k,dz=',k,dz(k)
     enddo
     
     do k = nz1+nz2+nz3-neven+1, nz1+nz2+nz3
       dz(k) = dz(nz-neven)
     enddo
     
     if(myid.eq.0) then
       print *,'plyudebug, z_grid para: q1, a1, nz1, L1=',q1,a1,nz1,L1
       print *,'plyudebug, z_grid para: q2, a2, nz2, L2=',q2,a2,nz2,L2
       print *,'plyudebug, z_grid para: q3, a3, nz3, L3=',q3,a3,nz3,L3
     endif

     zz1(1) = 0.0_wp
     do k = 2, nz+1
       zz1(k) = zz1(k-1) + dz(k-1); !print *, 'k,zz1=', k,zz1(k)
     enddo
     
     dz(1) = dz(1) /2.0_wp
     dz(nz-1) = dz(nz-1) / 2.0_wp
     dz(nz) = dz(nz) /2.0_wp
     
     zz2(1) = 0.0_wp
     do k = 2, nz+1
       zz2(k) = zz2(k-1) + dz(k-1); !print *, 'k,zz2=',k,zz2(k)
     enddo
     
     do k = 1, nz+1
       zz(k) = (zz2(k)-zz2(1))/(zz2(nz)-zz2(1))*zl;! print*,'k,zz=',k,zz(k)
     enddo 
     
     do k=1,nz
        dz(k)=zz(k+1)-zz(k)
     end do
     
     do k=1,nz
        if(k .eq. 1)then
           zw(k)=zz(k)
        else if(k .eq. nz-1) then
           zw(k)=zz(k+1)
        else if(k .eq. nz) then
           zw(k)=zz(k+1)+dz(k)
        else
           zw(k)=(zz(k)+zz(k+1))/2.
        endif
     enddo
     zw(nz+1) = 2.0*zw(nz) - zw(nz-1)

     if (myid == 0) then
        open(40)
        write(40,499)
499     format('variables=k, zz, dz, ratio, zw')
        
        do k=nz,1,-1
           write(40,599) float(k), zz(k), dz(k), dz(k+1)/dz(k), zw(k)
        enddo
599     format(5e24.14)
        close(40)
     end if
     
     ! added by plyu: export xyz.dat for inflow generation
     if (myid ==0) then
       open(41, file='xyz.dat', action='write')
       write(41,*) ny, nz, nx
       do j = 1, ny
         write(41,*) yl/ny*j, zz(1), xl/nx*1
       enddo
       do k = 1, nz
         write(41,*) yl/ny*1, zz(k), xl/nx*1
       enddo
       do i = 1, nx
         write(41,*) yl/ny*1, zz(1), xl/nx*i
       enddo 
       close(41)
     endif
     
   END SUBROUTINE  cluster_nl_3
   !> cluster_nl_3 is added by plyu to refine grid for turbine 

   subroutine cluster_nl(nz,zl,clbeta,clgama,zz1,zz2,zz,zw,dz)
   !     clgama: largest grid ration
   !     clbeta: clustering papameter: smaller value -> more clustered
   !     zl domain size in the vertical direction
   !     zz vertical coordinate in the transformed domain
   !     zw staggered grid in vertial coordinate
      implicit none

      integer,  intent(IN) :: nz
      real(wp), intent(IN) :: zl, clbeta, clgama
      real(wp), dimension(*), intent(OUT):: zz1,zz2,zz,zw,dz
 
      integer  :: neven=4, k
      real(wp) :: a, etak 
      
      do k=neven,nz-neven
         etak=(k-neven)/(nz-2.0_wp*neven)
         zz1(k)=((1.0_wp+clbeta)*((clbeta+1.0_wp)/(clbeta-1.0_wp))**(2.0_wp*etak-1.0_wp) &
         +1.0_wp-clbeta)/(2.0_wp*(1.0_wp+((clbeta+1.0_wp)/(clbeta-1.0_wp))**(2.0_wp*etak-1.0_wp)))
      end do
      
      
      a=1.0_wp/(2.0_wp*clbeta)*log((1.0_wp+(exp(clbeta)-1.0_wp)*0.5_wp) &
      /(1.0_wp+(exp(-clbeta)-1.0_wp)*0.5_wp))
      do k=neven+1,nz-neven-1
         etak=(k-neven)/(nz-2.0_wp*neven) 
         zz1(k)=0.5_wp*(1.0_wp+sinh(clbeta*(etak-a))/sinh(clbeta*a))
      end do
      zz1(neven)=0.0_wp
      zz1(nz-neven)=1.0_wp
      
      do k=neven,nz-neven
         dz(k)=(1.0_wp-cos(zz1(k)*twopi))*(clgama-1.0_wp)/2.0_wp+1.0_wp
      end do
      
      do k=1,neven-1
         dz(k)=dz(neven)
      end do
      
      do k=nz-neven+1,nz
         dz(k)=dz(nz-neven)
      end do
      
      dz(1)=dz(1)/2.0_wp
      dz(nz-1)=dz(nz-1)/2.0_wp
      dz(nz)=dz(nz)/2.0_wp
      
      zz2(1)=0.0_wp
      do k=1,nz
         zz2(k+1)=zz2(k)+dz(k)
      end do
      
      do k=1,nz+1
         zz(k)=(zz2(k)-zz2(1))/(zz2(nz)-zz2(1))*zl
      end do
      
      do k=1,nz
         dz(k)=zz(k+1)-zz(k)
      end do
      
      do k=1,nz
         if(k == 1) then
            zw(k)=zz(k)
         else if(k == nz-1) then
            zw(k)=zz(k+1)
         else if(k == nz) then
            zw(k)=zz(k+1)+dz(k)
         else
            zw(k)=(zz(k)+zz(k+1))/2
         endif
      enddo
   
  end subroutine cluster_nl 

  subroutine cluster_nl_2(nz,zl,clbeta,clgama,zz1,zz2,zz,zw,dz)
    !     clgama: largest grid ration
    !     clbeta: clustering papameter: smaller value -> more clustered
    !     zl domain size in the vertical direction
    !     zz vertical coordinate in the transformed domain
    !     zw staggered grid in vertial coordinate
    implicit none
    
    integer,  intent(IN) :: nz
    real(wp), intent(IN) :: zl, clbeta, clgama
    real(wp), dimension(*), intent(OUT):: zz1,zz2,zz,zw,dz
    
    integer  :: neven1=4, neven2, k
    real(wp) :: a, etak 
    
    neven2=4
    
    do k=neven1,nz-neven2
       etak=(k-neven1)/(nz-1.0_wp*neven1-1.0_wp*neven2)
       zz1(k)=((1.0_wp+clbeta)*((clbeta+1.0_wp)/(clbeta-1.0_wp))**(2.0_wp*etak-1.0_wp) &
            +1.0_wp-clbeta)/(2.0_wp*(1.0_wp+((clbeta+1.0_wp)/(clbeta-1.0_wp))**(2.0_wp*etak-1.0_wp)))
    end do
    
    
    a=1.0_wp/(2.0_wp*clbeta)*log((1.0_wp+(exp(clbeta)-1.0_wp)*0.5_wp) &
         /(1.0_wp+(exp(-clbeta)-1.0_wp)*0.5_wp))
    do k=neven1+1,nz-neven2-1
       etak=(k-neven1)/(nz-1.0_wp*neven1-1.0_wp*neven2) 
       zz1(k)=0.5_wp*(1.0_wp+sinh(clbeta*(etak-a))/sinh(clbeta*a))
    end do
    zz1(neven1)=0.0_wp
    zz1(nz-neven2)=1.0_wp
    
    do k=neven1,nz-neven2
       ! do k=neven,nz/2
       dz(k)=(1.0_wp-cos(zz1(k)*twopi/2))*(clgama-1.0_wp)/2.0_wp+1.0_wp
    end do
    
    ! do k=nz/2+1, nz-neven
    !    dz(k)=dz(nz/2)
    ! end do 
    
    do k=1,neven1-1
       dz(k)=dz(neven1)
    end do
    
    do k=nz-neven2+1,nz
       dz(k)=dz(nz-neven2)
    end do
    
    dz(1)=dz(1)/2.0_wp
    dz(nz-1)=dz(nz-1)/2.0_wp
    dz(nz)=dz(nz)/2.0_wp
      
    zz2(1)=0.0_wp
    do k=1,nz
       zz2(k+1)=zz2(k)+dz(k)
    end do
    
    do k=1,nz+1
       zz(k)=(zz2(k)-zz2(1))/(zz2(nz)-zz2(1))*zl
    end do
    
    do k=1,nz
       dz(k)=zz(k+1)-zz(k)
    end do
    
    do k=1,nz
       if(k == 1) then
          zw(k)=zz(k)
       else if(k == nz-1) then
          zw(k)=zz(k+1)
       else if(k == nz) then
          zw(k)=zz(k+1)+dz(k)
       else
          zw(k)=(zz(k)+zz(k+1))/2
       endif
    enddo
    
!    write(40,499)
499 format('variables=k,zz, dz, ratio')
    
!    do k=nz,1,-1
!       write(40,599)float(k), zz(k), dz(k), dz(k+1)/dz(k)     
!    enddo
599 format(4e24.14)
    
  end subroutine cluster_nl_2
end module grid


