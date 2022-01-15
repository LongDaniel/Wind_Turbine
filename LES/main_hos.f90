! program grid_test

!   use MPI
!   use hos
!   use fft_hos
!   implicit none

!   call my_mpi_init
!   call input_hos
!   call hos_init
!   call fft_init_hos
!   call init_random_seed(0)

!   call run_hos
! !  call test
  
!   call fft_finalize_hos
!   if (myid_hos == 0) then
!      print *,"program ends successfully!"
!   end if
!   call my_mpi_finalize

! end program grid_test

!-------------------------------------------------------------

! subroutine test
  
!   use hos
!   use fft_hos
!   implicit none

!   real(wp) x,y,pa_ex(nxhos,nyhos/ncpu_hos)
!   real(wp) pao(nxhos,nyhos), pao_ex(nxhos,nyhos), eo(nxhos,nyhos)

!   integer i,j

!   do i = 1, nxhos
!      x = (i - 1) * dx_hos
!      do j = 1, nyhos/ncpu_hos
!         y = (j + myid_hos * nyhos/ncpu_hos - 1) * dy_hos
!         eta_hos(i,j) = cos(x) + sin(y)
!         pa_ex(i,j) = (2 / weber) * (-cos(x)-cos(x)*cos(y)**2-sin(y)+sin(y)*sin(x)**2) &
! / 2 / (1+sin(x)**2+cos(y)**2)**1.5
!      end do
!   end do

!   call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

!   call pa_surten(eta_hos,ex_hos,ey_hos,pa_st)

!   call alltoone(pa_st,pao)
!   call alltoone(eta_hos,eo)

!   if (myid_hos == 0) then

!      write(93,*) ' ZONE T="',0.0,'" I=', nxhos, ' J=',nyhos,' F=POINT'
!      do j = 1, nyhos
!         y = (j - 1) * dy_hos
!         do i = 1, nxhos
!            x = (i - 1) * dx_hos
!            write(93,*) x,y,pao(i,j),eo(i,j)
        ! end do
     ! end do
  ! end if
        

  ! print *,"test end!"

  ! return
! end subroutine test

!------------------------------------------------------------------------------

subroutine run_hos

  use hos
  use io_hos
  implicit none

  integer it,icon
  real(wp) :: sigma, period
  real(wp) :: tmp

  icon = 16
!  call init_old_data

  if (istart_hos .eq. 0) then
     time_hos = 0.0
     pa_hos(1:nxhos,1:nyhos/ncpu_hos) = 0.0
     call jonswap_init
     call jonswap_3d
!     call stokes_2d
!     call linear_cap
!     call save_hos(eta_hos,vps_hos,pa_hos)

     ! ioutc_hos = noutc_hos - 1
     ! call outsurf(eta_hos,ioutd_hos,ioutc_hos,time_hos/period)
     return
  end if

!  call read_hos(eta_hos,vps_hos,pa_hos)

  if (myid_hos == 0) then
     print *,"reread successful!"
  end if

  call wavenum(wvn_hos)
  call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
  call zeta(eta_hos,zp_hos)  
  tmp = 1.0

  call boundvp(vps_hos,r_hos,zp_hos)

  call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
  call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

  sigma = (nswavex * pex_hos / fr2_hos)**0.5
  period = twopi / sigma
  dt_hos = period / ntp

  ! ioutc_hos = noutc_hos - 1
  ! call outsurf(eta_hos,ioutd_hos,ioutc_hos,time_hos/period)

  
  do it = 1, ntime_hos

     time_hos = time_hos + dt_hos

     call hos_wave_3d(eta_hos,vps_hos,dt_hos,pa_hos)

     tmp = sqrt(sum(eta_hos*eta_hos))
     if (isnan(sum(eta_hos))) then
        print *,"code blows up!"
        stop
     end if
     
     ! call outsurf(eta_hos,ioutd_hos,ioutc_hos,time_hos/period)

     if (myid_hos == 0) print *,"time/Tp=",time_hos/period

  end do

  if (myid_hos == 0) then
     print *,"check new eta here:"
     print '(" ",E25.16," ")',eta_hos(1:5,1)
  end if

!  call save_hos(eta_hos,vps_hos,pa_hos)

end subroutine run_hos

!-------------------------------------------------------------------------------

SUBROUTINE INIT_RANDOM_SEED(MYSEED)

     ! BY XUANTING HAO

     ! THIS SUBROUTINE INITIALIZES THE RANDOM NUMBER SEED
     ! MYSEED IS USED TO CONTROL THE DIFFERENCE BETWEEN DIFFERENT CPUS.

     ! THE SUBROUTINE MUST BE CALLED AT THE BEGINING OF PROGRAM.

     ! 10/1/2014 FIRST EDITION


      IMPLICIT NONE
      INTEGER :: I, N, CLOCK
      INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
      INTEGER MYSEED

      CALL RANDOM_SEED(SIZE = N)
      ALLOCATE(SEED(N))

      CALL SYSTEM_CLOCK(COUNT=CLOCK)

      SEED = CLOCK + MYSEED * (/ (I - 1, I = 1, N) /)
      CALL RANDOM_SEED(PUT = SEED)

      DEALLOCATE(SEED)
END SUBROUTINE INIT_RANDOM_SEED
