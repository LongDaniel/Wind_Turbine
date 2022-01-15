module hos

   use constants
   use decomp, only : wp, mpi_comm_2d_col

   use hos_param
   use fft_hos
   use smooth
   use spectral_hos
   implicit none

   public :: my_mpi_init, my_mpi_finalize, input_hos, hos_init
   public :: jonswap_2d, jonswap_3d

!   public :: mpe_2dtrans, getij

contains

  subroutine my_mpi_init()
    implicit none

    call mpi_init(ierr_hos)
    call mpi_comm_rank(mpi_comm_2d_col,myid_hos,ierr_hos)
    ! call mpi_comm_size(mpi_comm_2d_col,numprocs,ierr)

  end subroutine my_mpi_init

  subroutine my_mpi_finalize()
    
    implicit none

    call mpi_finalize(ierr_hos)

  end subroutine my_mpi_finalize

   subroutine input_hos()

     implicit none

      open(12)
      read(12,*) nxhos,nyhos,npw,ncpu_hos
      read(12,*) istart_hos
      read(12,*) iswell
      read(12,*) ntime_hos,noutd_hos,noutc_hos
      read(12,*) ntp
      read(12,*) pex_hos,pey_hos
      read(12,*) nswavex
      read(12,*) aka_hos,nswellx,nswelly
      read(12,*) fr2_hos
      read(12,*) ustar,u10,uss
      read(12,*) gamma,fetch,g,phi
      read(12,*) ist_hos,weber
      close(12)

      ! if (ncpu_hos /= numprocs) then
      !    print *,"CPU asked does not equal to that provided!"
      !    stop
      ! end if

      dx_hos = twopi / pex_hos / nxhos
      dy_hos = twopi / pey_hos / nyhos

   end subroutine input_hos

  !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================
   subroutine hos_init()
      implicit none

      allocate(eta_hos(nxhos,nyhos/ncpu_hos),ex_hos(nxhos,nyhos/ncpu_hos),ey_hos(nxhos,nyhos/ncpu_hos))
      allocate(vps_hos(nxhos,nyhos/ncpu_hos),vpsx_hos(nxhos,nyhos/ncpu_hos),vpsy_hos(nxhos,nyhos/ncpu_hos))
      allocate(feta(nxhos,nyhos/ncpu_hos,4),fvps(nxhos,nyhos/ncpu_hos,4))
      allocate(u_hos(nxhos,nyhos/ncpu_hos),v_hos(nxhos,nyhos/ncpu_hos),w_hos(nxhos,nyhos/ncpu_hos))
      allocate(pa_hos(nxhos,nyhos/ncpu_hos),pa_st(nxhos,nyhos/ncpu_hos))
      allocate(pa0_hos(nxhos,nyhos/ncpu_hos))
      allocate(wvn_hos(nxhos,nyhos/ncpu_hos,npw),r_hos(nxhos,nyhos/ncpu_hos,npw),zp_hos(nxhos,nyhos/ncpu_hos,npw))
      !print *,nxhos,nyhos,ncpu_hos
      !print *,size(eta_hos,1),size(eta_hos,2)

      eta_hos = 0.0
      vps_hos = 0.0
      time_hos = 0.0
      
   end subroutine hos_init

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine linear_cap

     implicit none

     integer i,j
     real(wp)::theta,omega,amp,x,phase

     real(wp), allocatable, dimension(:,:) :: eo,vo

     allocate(eo(nxhos,nyhos),vo(nxhos,nyhos))
     if (myid_hos == 0) then
        call random_number(theta)
        theta = theta * twopi
        omega = (nswavex * pex_hos / fr2_hos + (nswavex * pex_hos)**3 / weber)**0.5

        amp = aka_hos / nswavex / pex_hos
        do i = 1, nxhos
           x = (i - 1) * dx_hos
           phase = nswavex * pex_hos * x + theta
           do j = 1, nyhos
              eo(i,j) = amp * cos(phase)
              vo(i,j) = (amp / fr2_hos / omega) * sin(phase)
           enddo
        enddo
     end if

     call onetoall(eo,eta_hos)
     call onetoall(vo,vps_hos)

     deallocate(eo,vo)

     return
   end subroutine linear_cap


   subroutine stokes_2d()

     implicit none

     integer i,j
     real(wp)::theta,omega,amp,x,phase

     real(wp), allocatable, dimension(:,:) :: eo,vo

     allocate(eo(nxhos,nyhos),vo(nxhos,nyhos))
     if (myid_hos == 0) then
        call random_number(theta)
        theta = theta * twopi
        omega = (nswavex * pex_hos / fr2_hos)**0.5
        
        amp = aka_hos / nswavex / pex_hos
        do i = 1, nxhos
           x = (i - 1) * dx_hos
           phase = nswavex * pex_hos * x + theta
           do j = 1, nyhos 
              eo(i,j) = amp * (cos(phase) + 0.5 * aka_hos * cos(2*phase) + 0.375 * aka_hos**2 * cos(3*phase))
              vo(i,j) = (amp / fr2_hos / omega) * exp(nswavex*pex_hos * eo(i,j)) * sin(phase)
           enddo
        enddo
     end if

     call onetoall(eo,eta_hos)
     call onetoall(vo,vps_hos)
     
     deallocate(eo,vo)

     return
   end subroutine stokes_2d

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================


   subroutine jonswap_init()

     implicit none

     alpha = 0.076 * (u10**2 / fetch / g)**0.22
     omega_p = 22.0 * (g**2 / u10 / fetch)**(1.0/3.0)
     lambda_p = twopi * g / omega_p**2
     omega0 = (nswavex * pex_hos / fr2_hos)**0.5
     

   end subroutine jonswap_init

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine jonswap_2d()

     implicit none

     integer i,j,k
     real(wp)::theta,omega,sw,aa,x,sigma

     do k = 1, (nxhos + 2) / 2
        call random_number(theta)
        theta = theta * twopi
        omega = (k * pex_hos / fr2_hos)**0.5
        if (omega < omega0) then
           sigma = 0.07
        else
           sigma = 0.09
        end if
        sw = alpha / fr2_hos**2 / omega**5 * exp(-1.25*(omega/omega0)**(-4.0)) &
             * gamma**(exp(-(omega-omega0)**2/2/sigma**2/omega0**2))
        sw = sw / fr2_hos / 2.0 / omega
        aa = (2.0 * sw * pex_hos)**0.5
        do i = 1, nxhos
           x = (i - 1) * dx_hos
           do j = 1, nyhos / ncpu_hos
              eta_hos(i,j) = eta_hos(i,j) + aa * cos(k*pex_hos*x+theta)
              vps_hos(i,j) = vps_hos(i,j) + aa / fr2_hos / omega * sin(k*pex_hos*x+theta)
           enddo
        enddo
     end do
     
     return
   end subroutine jonswap_2d



   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine jonswap_3d

     implicit none
     
     real(wp), allocatable, dimension(:,:) :: eta1,vps1
     integer :: kx, ky, kyr, ky1, ky2, kyr1, kyr2, kyy1, kyy2
     real(wp) :: k, ck, theta, theta1, sw, omega, aa, sigma

     allocate(eta1(nxhos,nyhos),vps1(nxhos,nyhos))
     eta1 = 0
     vps1 = 0
     
     print*, pex_hos, pey_hos

     do kx = 0, nxhos / 2 - 1
        do ky = 1, nyhos, 2
           kyr = (ky + 1) / 2 - 1

           if (kyr < nyhos / 2) then
              if (kx .ne. 0 .or. kyr .ne. 0) then
                 k = ((kx * pex_hos)**2 + (kyr * pey_hos)**2)**0.5
                 ck = kx * pex_hos / k
                 call random_number(theta)
                 call random_number(theta1)
                 theta = theta * twopi
                 theta1 = theta1 * twopi
                 omega = (k / fr2_hos)**0.5
                 if ( omega .lt. omega0 ) then
                    sigma = 0.07_wp
                 else
                    sigma = 0.09_wp
                 end if
                 sw = alpha / fr2_hos**2 / omega**5 * exp(-1.25*(omega/omega0)**(-4.))&
                      * gamma**(exp(-(omega-omega0)**2/2/sigma**2/omega0**2)) * (4. / twopi) * ck**2
                 sw = sw / fr2_hos**2 / 2.0_wp / omega**3

                 aa = ( 2.0_wp * sw * pex_hos * pey_hos )**0.5
                 
                 omega = ( k / fr2_hos )**0.5_wp
                 ky1 = 1
                 ky2 = nyhos
                 kyr1 = 2 * kyr + 1
                 kyr2 = 2 * kyr + 2
                 kyy1 = kyr1
                 kyy2 = kyr2
                 if ( kyr1 .eq. 1 .and. kyr2 .eq. 2 ) then
                    theta1 = theta
                 end if
                 aa = aa / 4
                 if(kx.eq.0.or.kyr.eq.0) then
                    aa=aa*sqrt(2.0_wp)
                 endif
                 if ( kyr1 .ge. ky1 .and. kyr1 .le. ky2 ) then
                    eta1(kx*2+1,kyy1) = aa * cos(theta) + aa * cos(theta1)
                    eta1(kx*2+2,kyy1) = aa * sin(theta) + aa * sin(theta1)
                    vps1(kx*2+1,kyy1) = aa / fr2_hos / omega * sin(theta) + aa / fr2_hos / omega * sin(theta1)
                    vps1(kx*2+2,kyy1) = - aa / fr2_hos / omega * cos(theta) - aa / fr2_hos / omega * cos(theta1)
                 end if
                 
                 if ( kyr2 .ge. ky1 .and. kyr2 .le. ky2 ) then
                    eta1(kx*2+1,kyy2) = aa * sin(theta) - aa * sin(theta1)
                    eta1(kx*2+2,kyy2) = - aa * cos(theta) + aa * cos(theta1)
                    vps1(kx*2+1,kyy2) = - aa / fr2_hos / omega * cos(theta) + aa / fr2_hos / omega * cos(theta1)
                    vps1(kx*2+2,kyy2) = - aa / fr2_hos / omega * sin(theta) + aa / fr2_hos / omega * sin(theta1)
                 end if
              end if
           end if
        end do
     end do

     call onetoall(eta1,eta_hos)
     call onetoall(vps1,vps_hos)

     call fft_bac_xy_hos(eta_hos)
     call fft_bac_xy_hos(vps_hos)

     deallocate(eta1,vps1)

   end subroutine jonswap_3d

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================
   subroutine hos_wave_3d(eta,vps,dt,pa)

     implicit none

     real(wp), intent(in) :: dt
     real(wp), intent(in), dimension(:,:) :: pa
     real(wp), intent(inout), dimension(:,:) :: eta, vps

     real(wp) :: ratio
     integer :: irk
     real(wp), allocatable, dimension(:,:) :: eta0, vps0
     real(wp) :: beta1,beta2

     allocate(eta0(nxhos,nyhos/ncpu_hos), vps0(nxhos,nyhos/ncpu_hos))

!     ratio = 0.5_wp
     ratio = 0.8_wp
     beta1 = 8.0_wp
     beta2 = 30.0_wp
     pa_st = 0.0_wp
      !print *, 'surf_rk4,0'
      !print *, eta(1:32,1)
     call surf_rk4(eta,vps,eta0,vps0,feta,fvps,dt,1)
      !print *, 'surf_rk4,1'
      !print *, eta(1:32,1)

!     if (myid_hos == 0) then
!        print *,"check eta before fft"
!        print '(" ",E25.16," ")',eta(1,1)
!     end if

     if (ist_hos == 1) then
        call pa_surten(eta,ex_hos,ey_hos,pa_st)
     end if

     call righ(eta,ex_hos,ey_hos,vpsx_hos,vpsy_hos,w_hos,feta,fvps,pa,pa_st,1)

     do irk = 2,4
        call surf_rk4(eta,vps,eta0,vps0,feta,fvps,dt,irk)
        !print *, 'surf_rk4,',irk
        !print *, eta(1:32,1)

        ! call filter_lp(eta,vps,beta1,beta2)
!        call filter_lp(eta,vps,ratio)

        call derivh(eta,vps,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

        if (ist_hos == 1) then
           call pa_surten(eta,ex_hos,ey_hos,pa_st)
        end if

        call zeta(eta,zp_hos)

        call boundvp(vps,r_hos,zp_hos)

        call wsurf(w_hos,r_hos,zp_hos,wvn_hos)

        call righ(eta,ex_hos,ey_hos,vpsx_hos,vpsy_hos,w_hos,feta,fvps,pa,pa_st,irk)

     end do

     call surf_update(eta,vps,eta0,vps0,feta,fvps,dt)
      !print *, 'surf_update'
      !print *, eta(1:32,1)
     call filter_lp(eta,vps,beta1,beta2)
        ! print*, 'eta', eta(1:32,1)
        ! print*, 'vps', vps(1:32,1)

     !call filter_lp(eta,vps,ratio)
! lp: low pass. ratio: remained part
      !print *, 'filter_lp'
      !print *, eta(1:32,1) 
    
     call derivh(eta,vps,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

     call zeta(eta,zp_hos)
     
     call boundvp(vps,r_hos,zp_hos)
     
     call wsurf(w_hos,r_hos,zp_hos,wvn_hos)

     call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

     deallocate(eta0,vps0)
   end subroutine hos_wave_3d

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine surf_rk4(eta,vps,eta0,vps0,feta,fvps,dt,nrk4)

     implicit none
     
     real(wp), intent(in) :: dt
     integer, intent(in) :: nrk4
     real(wp), intent(in), dimension(:,:,:) :: feta, fvps
     real(wp), intent(inout), dimension(:,:) :: eta0, eta, vps0, vps

     real(wp) :: fac

     if (nrk4 == 1) then
        eta0(:,:) = eta(:,:)
        vps0(:,:) = vps(:,:)
     else if (nrk4 == 2 .or. nrk4 == 3 .or. nrk4 == 4) then
        fac = 0.5_wp * dt
        if (nrk4 == 4) fac = dt

        eta(:,:) = eta0(:,:) + fac * feta(:,:,nrk4 - 1)
        vps(:,:) = vps0(:,:) + fac * fvps(:,:,nrk4 - 1)
     else

         PRINT*, "Invalid value for IRK in SURF_RK4 !"
         PRINT*, "IRK=",NRK4
         STOP
     end if

   end subroutine surf_rk4

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine surf_update(eta,vps,eta0,vps0,feta,fvps,dt)

     implicit none

     real(wp), intent(in) :: dt
     real(wp), intent(in), dimension(:,:,:) :: feta, fvps
     real(wp), intent(inout), dimension(:,:) :: eta0, eta, vps0, vps

     eta(:,:) = eta0(:,:) + dt * (feta(:,:,1) + 2 * feta(:,:,2) + 2 * feta(:,:,3) + feta(:,:,4)) / 6
     vps(:,:) = vps0(:,:) + dt * (fvps(:,:,1) + 2 * fvps(:,:,2) + 2 * fvps(:,:,3) + fvps(:,:,4)) / 6

   end subroutine surf_update

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine zeta(eta,zp_hos)

     implicit none
     
     real(wp), intent(in), dimension(:,:) :: eta
     real(wp), intent(inout), dimension(:,:,:) :: zp_hos
     
     integer k
     real(wp) :: ratio

     zp_hos(:,:,1) = eta(:,:)
     ratio = 1.0
     do k = 2, npw - 1
        zp_hos(:,:,k) = zp_hos(:,:,k-1) * zp_hos(:,:,1) / (k * 1.0)
        call dealiasxy_hos(zp_hos(:,:,k))
     end do

   end subroutine zeta

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine boundvp(vps,r,zp_hos)

     implicit none
     
     real(wp), intent(in), dimension(:,:) :: vps
     real(wp), intent(inout), dimension(:,:,:) :: r, zp_hos

     integer k, k1, k2
     real(wp), dimension(nxhos,nyhos/ncpu_hos) :: dr

     r(:,:,1) = vps

     call fft_for_xy_hos(r(:,:,1))

     do k = 2, npw
        r(:,:,k) = 0.0

        do k1 = 1, k-1
           k2 = k - k1
           dr(:,:) = r(:,:,k2) * wvn_hos(:,:,k1)
           call fft_bac_xy_hos(dr)
           r(:,:,k) = r(:,:,k) - zp_hos(:,:,k1) * dr(:,:)
        end do

        call fft_for_xy_hos(r(:,:,k))
     end do

     do k = 1, npw
        call fft_bac_xy_hos(r(:,:,k))
     end do
     ! print*, 'r', r(1:6,3,3)

   end subroutine boundvp

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine wavenum(wvn_hos)

     implicit none
     
     real(wp), intent(inout), dimension(:,:,:) :: wvn_hos

     integer :: l,m,k,modex,modey
     real(wp) :: ax,ay,an 
     
     do k = 1, npw
        do m = 1, nyhos / ncpu_hos
           modey = (myid_hos * nyhos / ncpu_hos + m - 1) / 2
           ay = pey_hos * modey
           do l = 1, nxhos
              modex = (l - 1) / 2
              ax = pex_hos * modex
              an = (ax**2 + ay**2)**0.5
              wvn_hos(l,m,k) = an**k
           end do
        end do
     end do

   end subroutine wavenum

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine wsurf(ws,r,zp_hos,wvn_hos)

     implicit none

     real(wp), intent(in), dimension(:,:,:) :: r,zp_hos,wvn_hos
     real(wp), intent(inout), dimension(:,:) :: ws
     real(wp), dimension(nxhos,nyhos/ncpu_hos,npw) :: pk
     real(wp), dimension(nxhos,nyhos/ncpu_hos) :: tmp
     integer :: k,k1
     real(wp) :: ratio

     pk(:,:,1) = r(:,:,1)

     do k = 2, npw
        pk(:,:,k) = pk(:,:,k-1) + r(:,:,k)
     end do

     ws = 0

     do k = 1, npw - 1
        k1 = npw - k

        call fft_for_xy_hos(pk(:,:,k1))

        tmp(:,:) = pk(:,:,k1) * wvn_hos(:,:,k+1)

        call fft_bac_xy_hos(tmp)

        ws(:,:) = ws(:,:) + zp_hos(:,:,k) * tmp(:,:)
                
     end do

     call fft_for_xy_hos(pk(:,:,npw))

     tmp(:,:) = pk(:,:,npw) * wvn_hos(:,:,1)

     call fft_bac_xy_hos(tmp)

     ws = ws + tmp
     ratio = 1.0_wp
     call dealiasxy_hos(ws)
     ! print*,'ws', ws(1:32,1)
        
   end subroutine wsurf

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine derivh(eta,vps,ex,ey,vpsx,vpsy)

     implicit none

     real(wp), intent(in), dimension(:,:) :: eta,vps
     real(wp), intent(inout), dimension(:,:) :: ex,ey,vpsx,vpsy

     call pdfx_hos(eta,ex,pex_hos)
     call pdfy_hos(eta,ey,pey_hos)
     call pdfx_hos(vps,vpsx,pex_hos)
     call pdfy_hos(vps,vpsy,pey_hos)

   end subroutine derivh

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine uvsurf(us,vs,ws,ex,ey,vpsx,vpsy)

     implicit none

     real(wp), intent(in), dimension(:,:) :: ws,ex,ey,vpsx,vpsy
     real(wp), intent(inout), dimension(:,:) :: us,vs
     real(wp) :: ratio

     us = vpsx - ws * ex
     vs = vpsy - ws * ey

     ratio = 1.0_wp
     call dealiasxy_hos(us)
     call dealiasxy_hos(vs)

      ! print*, us(1:32,1)
      ! print*, vs(1:32,1)
   end subroutine uvsurf

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine righ(eta,ex,ey,vpsx,vpsy,ws,feta,fvps,pa,pa_st,nrk4)

     implicit none

     integer, intent(in) :: nrk4
     real(wp), intent(in), dimension(:,:) :: eta,ex,ey,vpsx,vpsy,ws,pa,pa_st
     real(wp), intent(inout), dimension(:,:,:) :: feta,fvps 
     
     real(wp), allocatable, dimension(:,:) :: t1, t2
     real(wp) :: ratio

     ratio = 1.0_wp
     allocate(t1(nxhos,nyhos/ncpu_hos),t2(nxhos,nyhos/ncpu_hos))
     if(nrk4.le.0.or.nrk4.gt.4) then
        print*, "invalid value for irk in righ !"
        print*, "irk=",nrk4
        stop
     endif

     t1 = 1.0_wp + ex * ex + ey * ey
     call dealiasxy_hos(t1)

     feta(:,:,nrk4) = -(vpsx * ex + vpsy * ey) + t1 * ws 
     call dealiasxy_hos(feta(:,:,nrk4))

     t2 = ws * ws
     call dealiasxy_hos(t2)
     fvps(:,:,nrk4) = -eta / fr2_hos - 0.5_wp * (vpsx**2 + vpsy**2) + 0.5_wp * t1 * t2 - pa - pa_st 
     call dealiasxy_hos(fvps(:,:,nrk4))

     deallocate(t1,t2)

   end subroutine righ

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine pa_surten(eta,ex,ey,pa_st)

     implicit none

     real(wp), intent(in), dimension(:,:) :: eta,ex,ey
     real(wp), intent(out), dimension(:,:) :: pa_st

     real(wp), allocatable, dimension(:,:) :: exx,eyy,exy

     ! integer i,j


     allocate(exx(nxhos,nyhos/ncpu_hos))
     allocate(eyy(nxhos,nyhos/ncpu_hos))
     allocate(exy(nxhos,nyhos/ncpu_hos))

     call pdfx_hos(ex,exx,pex_hos)
     call pdfy_hos(ey,eyy,pey_hos)
     call pdfx_hos(ey,exy,pex_hos)

     pa_st = -(2/weber) * ((1+ey**2)*exx - 2*ex*ey*exy + (1+ex**2) * eyy) / 2 / (1+ex**2 + ey**2)**1.5

     deallocate(exx,eyy,exy)

   end subroutine pa_surten

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

end module hos
