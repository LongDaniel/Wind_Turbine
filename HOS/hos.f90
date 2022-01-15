module hos

   use hos_param
   use fft_hos
   use smooth
   use spectral_hos
   implicit none

   public :: my_mpi_init, my_mpi_finalize, input_hos, hos_init
   public :: jonswap_2d, jonswap_3d, jonswap_dir

!   public :: mpe_2dtrans, getij

contains

  subroutine my_mpi_init()
    implicit none

    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,numprocs,ierr)

  end subroutine my_mpi_init

  subroutine my_mpi_finalize()
    
    implicit none

    call mpi_finalize(ierr)

  end subroutine my_mpi_finalize

   subroutine input_hos()

     implicit none

      open(12)
      read(12,*) nxhos,nyhos,npw,ncpu_hos
      read(12,*) istart_hos
      read(12,*) iswell
      read(12,*) ntime_hos,noutd_hos,noutc_hos
      read(12,*) ntp, dt_hos
      read(12,*) pex_hos,pey_hos
      read(12,*) nswavex
      read(12,*) aka,nswellx,nswelly
      read(12,*) fr2
      read(12,*) ustar,u10,uss
      read(12,*) gamma,fetch,g,phi
      read(12,*) ist_hos,bond,icap_hos
      close(12)

      if (ncpu_hos /= numprocs) then
         print *,"CPU asked does not equal to that provided!"
         stop
      end if

      if (mod(nyhos, ncpu_hos) /= 0) then
         print *,"Mod(nyhos, ncpu_hos) /= 0!"
         stop
      end if

      dx_hos = twopi / pex_hos / nxhos
      dy_hos = twopi / pey_hos / nyhos

   end subroutine input_hos

  !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================
   subroutine hos_init()
      implicit none

      allocate(eta_hos(nxhos,nyhos/ncpu_hos))
      allocate(ex_hos(nxhos,nyhos/ncpu_hos)) 
      allocate(ey_hos(nxhos,nyhos/ncpu_hos))
      allocate(vps_hos(nxhos,nyhos/ncpu_hos))
      allocate(vpsx_hos(nxhos,nyhos/ncpu_hos))
      allocate(vpsy_hos(nxhos,nyhos/ncpu_hos))
      allocate(feta(nxhos,nyhos/ncpu_hos,4),fvps(nxhos,nyhos/ncpu_hos,4))
      allocate(u_hos(nxhos,nyhos/ncpu_hos),v_hos(nxhos,nyhos/ncpu_hos),w_hos(nxhos,nyhos/ncpu_hos),pa_hos(nxhos,nyhos/ncpu_hos))
      allocate(wvn(nxhos,nyhos/ncpu_hos,npw),r_hos(nxhos,nyhos/ncpu_hos,npw),zp(nxhos,nyhos/ncpu_hos,npw))

      if (ist_hos == 1) then
         allocate(pa_st(nxhos,nyhos/ncpu_hos))
      end if

      eta_hos = 0.0
      vps_hos = 0.0
      time_hos = 0.0
      
   end subroutine hos_init

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine stokes_2d()

     implicit none

     integer i,j,k
     real(wp)::theta,omega,amp,x,sigma,phase

     real(wp), allocatable, dimension(:,:) :: eo,vo

     allocate(eo(nxhos,nyhos),vo(nxhos,nyhos))
     if (myid == 0) then
        !call random_number(theta)
        !theta = 0.0
        !theta = theta * twopi
        
        omega = (nswavex * pex_hos / fr2)**0.5
        theta = 71970.462834274076 * omega
        
        amp = aka / nswavex / pex_hos

        print *, "stokes_2d: omega=", omega, ", amp=", amp
        do i = 1, nxhos
           x = (i - 1) * dx_hos
           phase = nswavex * pex_hos * x + theta
           do j = 1, nyhos 
              eo(i,j) = amp * (cos(phase) + 0.5 * aka * cos(2*phase) + 0.375 * aka**2 * cos(3*phase))
              vo(i,j) = (amp / fr2 / omega) * exp(nswavex*pex_hos * eo(i,j)) * sin(phase)
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

   subroutine linear_2d

     implicit none

     integer i,j,k
     real(wp)::theta,omega,amp,x,sigma,phase

     real(wp), allocatable, dimension(:,:) :: eo,vo, eta0, vps0

     allocate(eo(nxhos,nyhos),vo(nxhos,nyhos))
     allocate(eta0(nxhos,nyhos/ncpu_hos),vps0(nxhos,nyhos/ncpu_hos))

     if (myid == 0) then        
        call random_number(theta)
        theta = theta * twopi
        if (ist_hos == 0) then
           omega = (nswellx * pex_hos / fr2)**0.5
        else
           omega = (nswellx * pex_hos * (1/fr2 + (nswellx*pex_hos)**2/bond/fr2))**0.5
        end if
        
        amp = aka / nswellx / pex_hos
        do i = 1, nxhos
           x = (i - 1) * dx_hos
           phase = nswellx * pex_hos * x + theta
           do j = 1, nyhos
              eo(i,j) = amp * cos(phase)
              vo(i,j) = (amp / fr2 / omega) * exp(nswellx*pex_hos * eo(i,j)) * sin(phase)
           enddo
        enddo
     end if

     call onetoall(eo,eta0)
     call onetoall(vo,vps0)

     eta_hos = eta_hos + eta0
     vps_hos = vps_hos + vps0

     deallocate(eo,vo)
     deallocate(eta0,vps0)

     return
   end subroutine linear_2d

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine cap_2d

     implicit none

     integer i,j,k
     real(wp)::theta,omega,amp,x,sigma,phase

     real(wp), allocatable, dimension(:,:) :: eo,vo
     real(wp), allocatable,dimension(:) :: tmpx, eta1,vps1
     real(wp) hlambda,tmpa,tmpc
     integer nt

     allocate(eo(nxhos,nyhos),vo(nxhos,nyhos))
     allocate(tmpx(nxhos/nswavex), eta1(nxhos/nswavex),vps1(nxhos/nswavex))

     if (myid == 0) then
        call random_number(theta)
        theta = theta * twopi        

        ! obtain dimensionless eta1, vps1 in one wave length
        hlambda = aka / pi
        tmpc = sqrt((nswavex * pex_hos / fr2 / bond) / sqrt(1 + 0.25 * pi**2 * hlambda**2))
        tmpa = (-2/hlambda/pi) * (1+sqrt(1+0.25 * pi**2 * hlambda**2))       
        do i = 1, nxhos / nswavex
           tmpx(i) = (i-1.0) / nxhos
        end do
        
        nt = 50
        vps1 = tmpx
        do i = 1, nt
           vps1 = -tmpx + 2*tmpa*sin(twopi*vps1)/pi/(1+2*tmpa*cos(twopi*vps1)+tmpa**2)
        end do

        eta1 = 2*tmpa*(tmpa+cos(twopi*vps1))/pi/(1+2*tmpa*cos(twopi*vps1)+tmpa**2)
        eta1 = eta1 - sum(eta1) / nxhos
        vps1 = vps1 + tmpx 

        ! extend to the entire domain
        do k = 0, nswavex - 1
           do j = 1, nyhos
              eo(k*nxhos/nswavex+1:(k+1)*nxhos/nswavex,j) = eta1 * (twopi / pex_hos / nswavex)
              vo(k*nxhos/nswavex+1:(k+1)*nxhos/nswavex,j) = vps1 * (twopi / pex_hos / nswavex) * tmpc
           end do
        end do

     end if

     call onetoall(eo,eta_hos)
     call onetoall(vo,vps_hos)

     eta_hos = eta_hos

     deallocate(eo,vo)
     deallocate(tmpx,eta1,vps1)

     return
   end subroutine cap_2d

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine jonswap_init()

     implicit none

     alpha = 0.076 * (u10**2 / fetch / g)**0.22
     omega_p = 22.0 * (g**2 / u10 / fetch)**(1.0/3.0)
     lambda_p = twopi * g / omega_p**2
     omega0 = (nswavex * pex_hos / fr2)**0.5
     
     if (myid == 0) then
        print *,"alpha=", alpha
        print *,"omega_p=", omega_p
        print *,"lambda_p=", lambda_p
        print *,"omega0=", omega0
     end if

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
        omega = (k * pex_hos / fr2)**0.5
        if (omega < omega0) then
           sigma = 0.07
        else
           sigma = 0.09
        end if
        sw = alpha / fr2**2 / omega**5 * exp(-1.25*(omega/omega0)**(-4.0)) &
             * gamma**(exp(-(omega-omega0)**2/2/sigma**2/omega0**2))
        sw = sw / fr2 / 2.0 / omega
        aa = (2.0 * sw * pex_hos)**0.5
        do i = 1, nxhos
           x = (i - 1) * dx_hos
           do j = 1, nyhos / ncpu_hos
              eta_hos(i,j) = eta_hos(i,j) + aa * cos(k*pex_hos*x+theta)
              vps_hos(i,j) = vps_hos(i,j) + aa / fr2 / omega * sin(k*pex_hos*x+theta)
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
     
     real(wp), allocatable, dimension(:,:) :: eta1,vps1, eta0, vps0
     integer :: kx, ky, kxr, kyr, ky1, ky2, kyr1, kyr2, kyy1, kyy2
     real(wp) :: k, ck, theta, theta1, sw, omega, aa, sigma

     allocate(eta1(nxhos,nyhos),vps1(nxhos,nyhos))
     allocate(eta0(nxhos,nyhos/ncpu_hos),vps0(nxhos,nyhos/ncpu_hos))

     eta1 = 0.0
     vps1 = 0.0
     
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
                 if (ist_hos == 0) then
                    omega = (k / fr2)**0.5
                 else if (ist_hos == 1) then
                    omega = (k * (1/fr2 + k**2/bond/fr2))**0.5
                 end if
                       
                 if ( omega .lt. omega0 ) then
                    sigma = 0.07
                 else
                    sigma = 0.09
                 end if
                 sw = alpha / fr2**2 / omega**5 * exp(-1.25*(omega/omega0)**(-4.))&
                      * gamma**(exp(-(omega-omega0)**2/2/sigma**2/omega0**2)) * (4. / twopi) * ck**2
                 if (ist_hos == 0) then
                    sw = sw / fr2**2 / 2. / omega**3
                 else if (ist_hos == 1) then
                    sw = sw * (1 / fr2 / 2. / omega / k) * (1 + 3*k**2/bond)
                 end if

                 aa = ( 2. * sw * pex_hos * pey_hos )**0.5
                 

                 ky1 = 1
                 ky2 = nyhos
                 kyr1 = 2 * kyr + 1
                 kyr2 = 2 * kyr + 2
                 kyy1 = kyr1
                 kyy2 = kyr2
                 if ( kyr1 .eq. 1 .and. kyr2 .eq. 2 ) then
                    theta1 = theta
                 end if
                 aa = aa / 4.
                 if(kx.eq.0.or.kyr.eq.0) then
                    aa=aa*sqrt(2.)
                 endif
                 if ( kyr1 .ge. ky1 .and. kyr1 .le. ky2 ) then
                    eta1(kx*2+1,kyy1) = aa * cos(theta) + aa * cos(theta1)
                    eta1(kx*2+2,kyy1) = aa * sin(theta) + aa * sin(theta1)
                    vps1(kx*2+1,kyy1) = aa / fr2 / omega * sin(theta) + aa / fr2 / omega * sin(theta1)
                    vps1(kx*2+2,kyy1) = - aa / fr2 / omega * cos(theta) - aa / fr2 / omega * cos(theta1)
                 end if
                 
                 if ( kyr2 .ge. ky1 .and. kyr2 .le. ky2 ) then
                    eta1(kx*2+1,kyy2) = aa * sin(theta) - aa * sin(theta1)
                    eta1(kx*2+2,kyy2) = - aa * cos(theta) + aa * cos(theta1)
                    vps1(kx*2+1,kyy2) = - aa / fr2 / omega * cos(theta) + aa / fr2 / omega * cos(theta1)
                    vps1(kx*2+2,kyy2) = - aa / fr2 / omega * sin(theta) + aa / fr2 / omega * sin(theta1)
                 end if
              end if
           end if
        end do
     end do

     call onetoall(eta1,eta0)
     call onetoall(vps1,vps0)

     call fft_bac_xy_hos(eta0)
     call fft_bac_xy_hos(vps0)

     call dealiasxy_hos(eta0)
     call dealiasxy_hos(vps0)

     eta_hos = eta_hos + eta0
     vps_hos = vps_hos + vps0

     deallocate(eta1,vps1)
     deallocate(eta0,vps0)

   end subroutine jonswap_3d

   subroutine jonswap_dir(fk,nfhos,nthos,u10_dir,fetch_dir,dfq, dth,ctime)

     implicit none

     real(wp), intent(in) :: u10_dir,fetch_dir,dfq, dth,ctime
     integer, intent(in) :: nfhos,nthos
     real(wp), intent(out), dimension(:,:) :: fk

     integer :: i,j
     real(wp) :: k, ck, theta, omega, sigma
     real(wp) :: alpha_dir, omegap_dir

     fk = 0

     alpha_dir = 0.076 * (u10_dir**2 / fetch_dir / g)**0.22
     omegap_dir = 22.0 * (g**2 / u10_dir / fetch_dir)**(1.0/3.0) * ctime
     
     do i = 1, nfhos
        omega = dfq * i * twopi
        do j = 1, nthos / ncpu_hos
           theta = (j+myid*nthos/ncpu_hos-1) * dth - pi
           if (abs(theta) >= pi/2) then
              ck = 0
           else
              ck = cos(theta)
           end if

           if ( omega .lt. omegap_dir ) then
              sigma = 0.07
           else
              sigma = 0.09
           end if
           fk(i,j) = twopi * alpha / fr2**2 / omega**5 * exp(-1.25*(omega/omegap_dir)**(-4.))&
                * gamma**(exp(-(omega-omegap_dir)**2/2/sigma**2/omegap_dir**2)) * (4. / twopi) * ck**2
        end do
     end do
     
     fk(1,:) = 0

   end subroutine jonswap_dir

   subroutine const_spec_pm_3d(eta,vps,time,phase)

     implicit none

     real(wp), intent(in) :: time
     real(wp), intent(in), dimension(:,:,:) :: phase
     real(wp), intent(out), dimension(:,:) :: eta,vps

     real(wp), allocatable, dimension(:,:) :: eta1,vps1
     real(wp) :: inv_cu
     integer :: kx, ky, kxr, kyr, ky1, ky2, kyr1, kyr2, kyy1, kyy2
     real(wp) :: k, ck, theta, theta1, sw, omega, aa, sigma, beta, capgam

     allocate(eta1(nxhos,nyhos),vps1(nxhos,nyhos))
     eta1 = 0.0
     vps1 = 0.0
     
     alpha = 8.1e-3_wp
!     alpha = 0.006 * (inv_cu)**0.55
     omega0 = (nswavex * pex_hos / fr2)**0.5

!     if (inv_cu < 1) then
!        gamma = 1.7
!     else 
!        gamma = 1.7 + 6.0 * log(inv_cu)
!     end if

!     sigma = 0.08 * (1 + 4 / inv_cu**3)

     do kx = 0, nxhos / 2 - 1
        do ky = 1, nyhos, 2
           kyr = (ky + 1) / 2 - 1

           if (kyr < nyhos / 2) then
              if (kx .ne. 0 .or. kyr .ne. 0) then
                 k = ((kx * pex_hos)**2 + (kyr * pey_hos)**2)**0.5
!                 ck = acos(kx * pex_hos / k)
                 ck = kx * pex_hos / k
                 
                 omega = (k / fr2)**0.5
                 theta = phase(kx+1,kyr+1,1) - omega * time
                 theta1 = phase(kx+1,kyr+1,2)- omega * time

!                 beta = 1.24
!                 if (omega/omega0 < 0.95 .and. omega/omega0 > 0.56) then
!                    beta = 2.61 * (omega/omega0)**1.3
!                 end if
!                 if (omega/omega0 < 1.6 .and. omega/omega0 >= 0.95) then
!                    beta = 2.28 * (omega/omega0)**(-1.3)
!                 end if

!                 sw = alpha / fr2**2 / omega**5 * (omega/omega0) * exp(-(omega/omega0)**(-4.))&
!                      * gamma**(exp(-(omega-omega0)**2/2/sigma**2/omega0**2)) &
!                      * 0.5 * beta / cosh(beta*ck)**2
                 sw = alpha / fr2**2 / omega**5 * exp(-1.25*(omega0/omega)**4)&
                      * (4. / twopi) * ck**2
                 sw = sw / fr2**2 / 2. / omega**3

                 aa = ( 2. * sw * pex_hos * pey_hos )**0.5
                 
                 ky1 = 1
                 ky2 = nyhos
                 kyr1 = 2 * kyr + 1
                 kyr2 = 2 * kyr + 2
                 kyy1 = kyr1
                 kyy2 = kyr2
                 if ( kyr1 .eq. 1 .and. kyr2 .eq. 2 ) then
                    theta1 = theta
                 end if
                 aa = aa / 4.
                 if(kx.eq.0.or.kyr.eq.0) then
                    aa=aa*sqrt(2.)
                 endif
                 if ( kyr1 .ge. ky1 .and. kyr1 .le. ky2 ) then
                    eta1(kx*2+1,kyy1) = aa * cos(theta) + aa * cos(theta1)
                    eta1(kx*2+2,kyy1) = aa * sin(theta) + aa * sin(theta1)
                    vps1(kx*2+1,kyy1) = aa / fr2 / omega * sin(theta) + aa / fr2 / omega * sin(theta1)
                    vps1(kx*2+2,kyy1) = - aa / fr2 / omega * cos(theta) - aa / fr2 / omega * cos(theta1)
                 end if
                 
                 if ( kyr2 .ge. ky1 .and. kyr2 .le. ky2 ) then
                    eta1(kx*2+1,kyy2) = aa * sin(theta) - aa * sin(theta1)
                    eta1(kx*2+2,kyy2) = - aa * cos(theta) + aa * cos(theta1)
                    vps1(kx*2+1,kyy2) = - aa / fr2 / omega * cos(theta) + aa / fr2 / omega * cos(theta1)
                    vps1(kx*2+2,kyy2) = - aa / fr2 / omega * sin(theta) + aa / fr2 / omega * sin(theta1)
                 end if
              end if
           end if
        end do
     end do

     call onetoall(eta1,eta)
     call onetoall(vps1,vps)

     call fft_bac_xy_hos(eta)
     call fft_bac_xy_hos(vps)

     call dealiasxy_hos(eta)
     call dealiasxy_hos(vps)

     deallocate(eta1,vps1)     

   end subroutine const_spec_pm_3d

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================
   subroutine hos_wave_3d(eta,vps,dt,pa,period)

     implicit none

     real(wp), intent(in) :: period, dt
     real(wp), intent(in), dimension(:,:) :: pa
     real(wp), intent(inout), dimension(:,:) :: eta, vps

     real(wp) :: ratio
     integer :: irk
     real(wp), allocatable, dimension(:,:) :: eta0, vps0
     real(wp) :: beta1,beta2

     allocate(eta0(nxhos,nyhos/ncpu_hos), vps0(nxhos,nyhos/ncpu_hos))

     ratio = 0.5
     beta1 = 8.0
     beta2 = 30.0

     call surf_rk4(eta,vps,eta0,vps0,feta,fvps,dt,1)

     if (ist_hos == 1) then
        call pa_surten(eta,ex_hos,ey_hos,pa_st)
     end if

     call righ(eta,ex_hos,ey_hos,vpsx_hos,vpsy_hos,w_hos,feta,fvps,pa,pa_st,1)

     do irk = 2,4
        call surf_rk4(eta,vps,eta0,vps0,feta,fvps,dt,irk)

!        call filter_lp(eta,vps,ratio)

!        call derivh(eta,vps,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
!        call filter_lp(eta,vps,ex_hos,ey_hos)
!        call filter_lp(eta,vps,beta1,beta2)
!        call filter_lp(eta,vps,ratio)

        call derivh(eta,vps,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

        if (ist_hos == 1) then
           call pa_surten(eta,ex_hos,ey_hos,pa_st)
        end if

        call zeta(eta,zp)

        call boundvp(vps,r_hos,zp)

        call wsurf(w_hos,r_hos,zp,wvn)

        call righ(eta,ex_hos,ey_hos,vpsx_hos,vpsy_hos,w_hos,feta,fvps,pa,pa_st,irk)

     end do

     call surf_update(eta,vps,eta0,vps0,feta,fvps,dt)

     call filter_lp(eta,vps,beta1,beta2)
!     call filter_lp(eta,vps,ratio)
     
     call derivh(eta,vps,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
     
     !> uncomment later
     !call filter_lp(eta,vps,ex_hos,ey_hos)
     !call derivh(eta,vps,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

     call zeta(eta,zp)
     
     call boundvp(vps,r_hos,zp)
     
     call wsurf(w_hos,r_hos,zp,wvn)

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
        fac = 0.5 * dt
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

     eta(:,:) = eta0(:,:) + dt * (feta(:,:,1) + 2.0 * feta(:,:,2) + 2.0 * feta(:,:,3) + feta(:,:,4)) / 6.0
     vps(:,:) = vps0(:,:) + dt * (fvps(:,:,1) + 2.0 * fvps(:,:,2) + 2.0 * fvps(:,:,3) + fvps(:,:,4)) / 6.0

   end subroutine surf_update

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine zeta(eta,zp)

     implicit none
     
     real(wp), intent(in), dimension(:,:) :: eta
     real(wp), intent(inout), dimension(:,:,:) :: zp
     
     integer k
     real(wp) :: ratio

     zp(:,:,1) = eta(:,:)
     ratio = 1.0
     do k = 2, npw - 1
        zp(:,:,k) = zp(:,:,k-1) * zp(:,:,1) / (k * 1.0)
        call dealiasxy_hos(zp(:,:,k))
     end do

   end subroutine zeta

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine boundvp(vps,r,zp)

     implicit none
     
     real(wp), intent(in), dimension(:,:) :: vps
     real(wp), intent(inout), dimension(:,:,:) :: r, zp

     integer k, k1, k2
     real(wp), dimension(nxhos,nyhos/ncpu_hos) :: dr

     r(:,:,1) = vps

     call fft_for_xy_hos(r(:,:,1))

     do k = 2, npw
        r(:,:,k) = 0.0

        do k1 = 1, k-1
           k2 = k - k1
           dr(:,:) = r(:,:,k2) * wvn(:,:,k1)
           call fft_bac_xy_hos(dr)
           r(:,:,k) = r(:,:,k) - zp(:,:,k1) * dr(:,:)
        end do

        call fft_for_xy_hos(r(:,:,k))
     end do

     do k = 1, npw
        call fft_bac_xy_hos(r(:,:,k))
     end do

   end subroutine boundvp

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine wavenum(wvn)

     implicit none
     
     real(wp), intent(inout), dimension(:,:,:) :: wvn

     integer :: l,m,k,modex,modey
     real(wp) :: ax,ay,an 
     
     do k = 1, npw
        do m = 1, nyhos / ncpu_hos
           modey = (myid * nyhos / ncpu_hos + m - 1) / 2
           ay = pey_hos * modey
           do l = 1, nxhos
              modex = (l - 1) / 2
              ax = pex_hos * modex
              an = (ax**2 + ay**2)**0.5
              wvn(l,m,k) = an**k
           end do
        end do
     end do

   end subroutine wavenum

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

   subroutine wsurf(ws,r,zp,wvn)

     implicit none

     real(wp), intent(in), dimension(:,:,:) :: r,zp,wvn
     real(wp), intent(inout), dimension(:,:) :: ws
     real(wp), allocatable, dimension(:,:,:) :: pk
     real(wp), allocatable, dimension(:,:) :: tmp
     integer :: k,k1
     real(wp) :: ratio

     allocate(pk(nxhos,nyhos/ncpu_hos,npw),tmp(nxhos,nyhos/ncpu_hos))

     pk(:,:,1) = r(:,:,1)

     do k = 2, npw
        pk(:,:,k) = pk(:,:,k-1) + r(:,:,k)
     end do

     ws = 0

     do k = 1, npw - 1
        k1 = npw - k

        call fft_for_xy_hos(pk(:,:,k1))

        tmp(:,:) = pk(:,:,k1) * wvn(:,:,k+1)

        call fft_bac_xy_hos(tmp)

        ws(:,:) = ws(:,:) + zp(:,:,k) * tmp(:,:)
                
     end do

     call fft_for_xy_hos(pk(:,:,npw))

     tmp(:,:) = pk(:,:,npw) * wvn(:,:,1)

     call fft_bac_xy_hos(tmp)

     ws = ws + tmp
     ratio = 1.0
     call dealiasxy_hos(ws)
        
     deallocate(pk,tmp)

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

     ratio = 1.0
     call dealiasxy_hos(us)
     call dealiasxy_hos(vs)

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

     ratio = 1.0
     allocate(t1(nxhos,nyhos/ncpu_hos),t2(nxhos,nyhos/ncpu_hos))
     if(nrk4.le.0.or.nrk4.gt.4) then
        print*, "invalid value for irk in righ !"
        print*, "irk=",nrk4
        stop
     endif

     t1 = 1.0 + ex * ex + ey * ey
     call dealiasxy_hos(t1)

     feta(:,:,nrk4) = -(vpsx * ex + vpsy * ey) + t1 * ws 
     call dealiasxy_hos(feta(:,:,nrk4))

     t2 = ws * ws
     call dealiasxy_hos(t2)
     fvps(:,:,nrk4) = - 0.5 * (vpsx**2 + vpsy**2) + 0.5 * t1 * t2 - pa

     if (ist_hos == 1) then
        fvps(:,:,nrk4) = fvps(:,:,nrk4) - pa_st
     end if


     if (icap_hos == 0) then
        fvps(:,:,nrk4) = fvps(:,:,nrk4) - eta / fr2 
     end if

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

     integer i,j


     allocate(exx(nxhos,nyhos/ncpu_hos))
     allocate(eyy(nxhos,nyhos/ncpu_hos))
     allocate(exy(nxhos,nyhos/ncpu_hos))

     call pdfx_hos(ex,exx,pex_hos)
     call pdfy_hos(ey,eyy,pey_hos)
     call pdfx_hos(ey,exy,pex_hos)

     pa_st = -(2/bond/fr2) * ((1+ey**2)*exx - 2*ex*ey*exy + (1+ex**2) * eyy) / 2 / (1+ex**2 + ey**2)**1.5

!     pa_st = (1/bond/fr2) * exx / (1+ex**2)**1.5

     deallocate(exx,eyy,exy)

   end subroutine pa_surten

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================

end module hos
