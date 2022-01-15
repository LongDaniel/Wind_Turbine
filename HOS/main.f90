program grid_test

  use MPI
  use hos
  use fft_hos
  use src_stat
  use post_proc
  implicit none
  integer ipost

  call my_mpi_init
  call input_hos

  call hos_init
  call fft_init_hos
  call init_random_seed(0)

  !call src_init
  read(15,*) ipost
  if (ipost == 0) then
     call run_hos
  else
     call post_wind_wave
  end if
  
  call fft_finalize_hos
  if (myid == 0) then
     print *,"program ends successfully!"
  end if
  call my_mpi_finalize

end program grid_test

!-------------------------------------------------------------


!------------------------------------------------------------------------------

subroutine run_hos

  use hos
  use io_hos
  use src_stat
  implicit none

  integer it,icon
  real(wp) :: sigma, period
  real(wp) :: tmp
  real(wp), dimension(nxhos/2,nyhos/2,2) :: phase
  real(wp) :: eo(nxhos,nyhos), skxy(nxhos/2,-nyhos/2:nyhos/2)
  real(wp), dimension(nfhos,nthos/ncpu_hos) :: fk, snl
  real(wp) :: omega, sw, ck, x, y, dwk
  real(wp), dimension(nfhos,nthos) :: fkall, snlall
  integer i,j, nmax
  real(wp), allocatable, dimension(:) :: sk, sf
!  real(wp) :: err

  icon = 16
!print *,1
  if (istart_hos .eq. 0) then
     time_hos = 0.0
     pa_hos(1:nxhos,1:nyhos/ncpu_hos) = 0.0

!     call jonswap_init
!     call jonswap_3d
!     call linear_2d
     call stokes_2d
!     call cap_2d
!     call save_hos(eta_hos,vps_hos,pa_hos)
!print *, 2
     call save_hos(time_hos,ioutd_hos,ioutc_hos,eta_hos,vps_hos,pa_hos)
!print *, 3
     ioutc_hos = noutc_hos - 1
     if (ist_hos == 1) then
        call outsurf(eta_hos,pa_st,ioutd_hos,ioutc_hos,time_hos)
     else
        call outsurf(ioutd_hos,ioutc_hos,time_hos)
     end if
!print *,4
     !test fkk2omg
     call alltoone(eta_hos, eo)
     call fullspec_xy(eo,skxy)
!     call fkk2omg(skxy, fk)
!print *, 5
     if (nfhos >= 32 .and. nthos >= 32) then
!        call jonswap_dir(fk,nfhos,nthos,dfq,dth)
        !call snl_webb(fk,snl)
     end if
!print *, 6
     call alltoone(snl,snlall,nfhos,nthos)
     call alltoone(fk,fkall,nfhos,nthos)

!     nmax = max(nxhos,nyhos)
!     dwk = sqrt((nxhos/2*pex_hos)**2+(nyhos/2*pey_hos)**2) / nmax
     nmax = nxhos/2-1
     dwk = pex_hos
     allocate(sk(nmax), sf(nfhos))
     call spec_1dk(eo, sk, nmax, dwk)
!     call spec_1df(eo, sf, nfhos, dfq)
!     call spec_1df_2(eo, sf, nfhos, dfq, sk, nmax, dwk)
     
     !export kx-ky spectrum
     !if (myid == 0) then
     !   open(384, action='write')
     !   write(384,*) " VARIABLES = freq,theta,ftheory,snl"
     !   write(384,"(A,F16.8,A,I5,A,I5,A)") ' ZONE T="',0.0,'" I=', nfhos, ' J=',nthos,' F=POINT'
     !   do j = 1, nthos
     !      do i = 1, nfhos
     !         write(384,'(25e12.4)') i*dfq, j*dth-pi,fkall(i,j),snlall(i,j)
     !      end do
     !   end do
     !   close(384)
     !end if
     
     !print *, 7
     return
  end if

  call read_hos(time_hos,ioutd_hos,ioutc_hos,eta_hos,vps_hos,pa_hos)

  if (myid == 0) then
     print *,"reread successful!"
  end if

  call wavenum(wvn)
  call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
  call zeta(eta_hos,zp)  
  tmp = 1.0

  call boundvp(vps_hos,r_hos,zp)

  call wsurf(w_hos,r_hos,zp,wvn)
  call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

  if (ist_hos == 1) then
     sigma = (nswavex * pex_hos * (1/fr2 + (nswavex*pex_hos)**2/bond/fr2))**0.5
  else
     sigma = (nswavex * pex_hos / fr2)**0.5
  end if

  period = twopi / sigma

  if (dt_hos <= 0) then
     dt_hos = period / ntp
  end if

  if (myid == 0) then
     print *,"dt_hos=",dt_hos
  end if

  do it = 1, ntime_hos

     time_hos = time_hos + dt_hos

     call hos_wave_3d(eta_hos,vps_hos,dt_hos,pa_hos,period)

     tmp = sqrt(sum(eta_hos*eta_hos))
     if (isnan(sum(eta_hos))) then
        print *,"code blows up!"
        stop
     end if

     if (ist_hos == 1) then
        call outsurf(eta_hos,pa_st,ioutd_hos,ioutc_hos,time_hos/period)
     else
        call outsurf(ioutd_hos,ioutc_hos,time_hos/period)
     end if

     if (myid == 0) print *,"Step=", it, ", time/Tp=",time_hos/period

  end do

  call save_hos(time_hos,ioutd_hos,ioutc_hos,eta_hos,vps_hos,pa_hos)

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
