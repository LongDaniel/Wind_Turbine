program fst

  use post_proc

!#ifndef POST
!  call run_sim
!#else
  call post_les_hos_turbine
  call post_tke_budget
  call post_vortex_dynamics
!#endif

!  call post_les_hos
!  call post_les_hos_wave_coherent_analysis
!  call post_les_hos_k_omg_spectrum
end program fst

subroutine run_sim

   use iso_fortran_env, only : INT64
   use param
   use mpi
   use decomp
   use fft
   use spectral
   use grid
   use navier
   use constants
   use io
   use utils
   !HOS modules
   !need to check
   use hos
   use hos_param
   use smooth
   use spectral_hos
   use fft_hos
   use io_hos
   use wavecontrol

   use turbine_model !> plyu
   !end
   
   implicit none

   integer :: ioutc, ioutd
   real(wp) :: time

   integer :: ierror, ierr

   double precision :: t1, t2

   !HOS variables
   real(wp) :: tmp
   integer icon

   ! save HOS 
   character (len=128) :: fileid, fparam_hos, fdata_hos
   
   integer(INT64) :: id
   integer i,j,k
   
   !> added by plyu
   integer :: irlm ! Remove Lateral Mean flow.
   
   call MPI_INIT(ierror)
   call input_iht
   call input_les('LES.IN')
   call decomp_init(nx,ny,nz,np1,np2)
   call fft_init
   call spectral_init
   call init_random_seed(myid)

   call grid_init(1)
   call navier_init(1)
   call les_init(1)

   call input_hos_par

   call fft_init_hos

   icon = 16

   !> added by plyu. They might be updated in turbine_model.
   ids = 0
   idsd = 0
   idsp = 0
    
   if (istart==0) then

      call grid_gen
     
      !> add by plyu
      call read_turbine_model_param(iturbine, time)

      time = 0; !print *, myid, 0.0
      call hos_init; !print *, myid, 0.1
      call bottom_hos_les(time); !print *, myid, 0.2
      call nl_coef; !print *, myid, 0.3
      !> plyunote: the cutoff() in add_turbulenc_velocity may lead to memory
      !!           overflow, so it may need be commented for large mesh.
      call add_turbulence_velocity(pu, pv, pw, u, v, w, 1);!print *,myid,0.4
      call add_log_velocity(u,1); !print *, myid, 0.5
     call mpi_barrier(mpi_comm_world, ierr);!print *, myid, 1 
      call dealiasxy(u(:,:,1:xsz(3)))
      call dealiasxy(v(:,:,1:xsz(3)))
      call dealiasxy(w(:,:,1:xsz(3)))
      call update_ghost(u,1)
      call update_ghost(v,1)
      call update_ghost(w,1)
     call mpi_barrier(mpi_comm_world, ierr);!print *, myid, 2 
      uzfs = 0
      vzfs = 0
      pp = 0
      
      call bc_lnr
      
      ioutc = 0; ioutd = 0
      
      call output_xy_average(time)
      call output_cut(time)
      call saverestart(time, ioutd, ioutc)
     call mpi_barrier(mpi_comm_world, ierr);!print *, myid, 3 
   else
      !call les_init(1)
      call reread(time, ioutd, ioutc)
      call grid_gen
      
      !> add by plyu
      call read_turbine_model_param(iturbine, time)

      if(iwavy.eq.1 .or.iwavy.eq.3)then
         
         call bottom2(time)
         
      elseif(iwavy.eq.7)then
         !need to add hos initialization subroutines below
         if(myid==0) print *, 'HOS initialization started!'
         if(isbot)then
            call hos_init
            call read_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos)
            !print*,'read_hos'
            !print*, eta_hos(1:32,1)
            
            call wavenum(wvn_hos)
            call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
            call zeta(eta_hos,zp_hos)
            tmp=1.0_wp
            call boundvp(vps_hos,r_hos,zp_hos)
            call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
            call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
            
            call bottom_hos_les(time)
            
            if(myid==0) print *, 'HOS initialized!'
         endif
      else
         if(myid==0) print *,"invalid iwavy"
         stop
      endif
      
      ! transfer data to the upper cpus
      call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
      call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
      call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
      call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
      call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
      call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
      call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
      call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
      
      ! plyunote: commented by me
      !eta=hh+0
      
      call init_p_poisson
      call nl_coef
      
      ! kfilt = 0
      if (myid==0) open(20)
      do i=1,ntime
         time = time+dt
         
         !> added by plyu
         !> sometime we may observe undesired lateral flow in simulation.
         !! if you dislike it, use this simple function to remove it.
         irlm = 0  
         if (irlm==1) then
           call remove_plane_mean(v)
         endif

         call nlfs_rk2_nn_les(time)
         
         ioutc = mod(ioutc+1,noutc)
         ioutd = mod(ioutd+1,noutd)
         
         if (ioutc == 0) then
            if (iturbine .ne. -100 .or. (iturbine.eq.-100 .and. ioutd==0)) then
              call output_xy_average(time)
              call output_cut(time)
              !call output_surface(time)
              ! if(isbot) call outsurf(eta_hos,time_hos)
            endif

            !> plyunote: output for turbine models
            if (iturbine .ne. 0) call output_turbine

            !> plyunote: for inlet plane timehistory output
            if (iturbine .eq. -100) call output_inlet(time)
         end if
         if (ioutd == 0) then
            !> plyunote: modified
            call output_all_2(time, id)
            if (iwavy .eq. 7) then
              write(fileid, '(i0.8)') nint(time/dt)
              fileid = trim(adjustl(fileid))
              fparam_hos = trim("restart_param_hos.dat"//fileid)
              fdata_hos = trim("restart_hos.h5"//fileid)
              if (isbot) then
                call save_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos,&
                  fparam_hos, fdata_hos)
              endif
            endif
            if (myid == 0) write(20,*) id
             
             !> added by plyu
             !irlm = 0
             !if (irlm==1) then
             !  call remove_plane_mean(v)
             !endif
         end if
      end do
      if (myid==0) close(20)
      
      !save hos data
      if(iwavy.eq.7)then
         if(isbot)then
            call save_hos(time, ioutd, ioutc, eta_hos,vps_hos,pa0_hos)
         endif
      endif
      
      call saverestart(time, ioutd, ioutc)
      
   end if
   
   ! call spectral_finalize
   call fft_finalize
   call decomp_finalize
   !HOS
   call fft_finalize_hos

   !> plyunote
   if (myid .eq. 0) print *, "LES ends successfully."
   call mpi_finalize(ierror)
 end subroutine run_sim


#include "initial.f90"


subroutine nlfs_rk2_nn_les(time)
   use navier
   use grid
   use wavecontrol
   use utils
   use constants, only : TWOPI
   use param, only : pex, pey, iwcontrol,  dt, Fr2, nth, rwe
   use hos
   use spectral
   
   implicit none

   real(wp), intent(IN) :: time
   !real(wp) :: th
   integer  :: ith

   !integer, dimension(3) :: maxidx
   !real(wp) :: rtmp

   real(wp) :: as, at, si, co, sigma
   real(wp), allocatable, dimension(:,:) :: etatmp, ettmp, etat, ett, padd

   integer :: nx, ny

   real(wp) :: temp
   integer :: ierr

   allocate(etatmp(xsz(1),xsz(2)))
   allocate(ettmp(xsz(1),xsz(2)))
   allocate(etat(xsz(1),xsz(2)))
   allocate(ett(xsz(1),xsz(2)))
   allocate(padd(xsz(1),xsz(2)))
   
   if (myid==0) then
      print '(1X,A)','======================'
      print '(1X,A,F15.7)', ' Time=', time
      print '(1X,A)','======================'
   end if

   !--------------
   ! Wave control
   !--------------
   if (iwcontrol==1)then
      padd = 0
      !if (istop .and. aa>1d-3) then
      if (isbot .and. aa>1d-3) then
         nx = 4
         ny = 0
         
         call DIVTRAVSTAND(nx, ny, U(:,:,1),V(:,:,1), &
              w(:,:,1),AS,AT,SI,CO,ETATMP,ETTMP,ETAT,ETT)
         call KILLSTANDINGWAVE_ONE(nx, ny, etatmp,ettmp,padd)
         call KILLSTANDINGWAVE_ALL(nx, ny, at,W(:,:,1),padd)
         
         at1 = (at+at1*exp(-dt)) / (1+exp(-dt))
         at = at1
         
         if (myid1==0) then
            write(999,*) time, as, aa, at, si, co
         end if
         
         sigma = sqrt(ak/Fr2 + ak**3 * RWe)
         temp = (aa**2-at**2)*twopi**2/pex/pey/Fr2/2
         call MAINTAIN3(ETAt,temp,AT,SIGMA,AK,PAdd)
      end if
   endif

   !-----------------
   ! Prediction step
   !-----------------
   !checked @grid.f
   call coef_et
   
   !get u,v,w
   call calc_F(time)
   !print *, "ns, 1"
   if (iwavy.eq. 7)then
      if(isbot)then
         ! get pas, air pressure on surface at time step n+1/2
         
         if(iwcontrol==1)then
            pp (:,:,1) = pp(:,:,1)+padd
         endif
         
         call get_pas(time)
         !print *, "ns, 2" 
         do ith =1, nth
            time_hos=time-(nth-ith)*dt/nth
            call hos_wave_3d(eta_hos,vps_hos,dt/nth,pa_hos)
         enddo
         !print *, "ns, 3"
         call bottom_hos_les(time)
      endif
   else
      call bottom2(time)
   end if
   !print *, "ns, 4"
   ! transfer data to the upper cpus
   call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
   call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
   call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
   call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
   call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
   call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
   call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
   call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)

   ! plyunote: commented by me
   !call top_fs_les_hos
   
   call nl_coef

   if(istop)then
      ut=u(:,:,xsz(3))
      vt=v(:,:,xsz(3))
      wt=0
   endif

   call press_g(time)

   !----------------
   ! correction step
   !----------------   
   call correction_us
   call bc_les
   
   !----------------
   ! Check validity
   !----------------
   call volume_lnr2(time)

   call dealiasxy(u(:,:,1:xsz(3)))
   call dealiasxy(v(:,:,1:xsz(3)))
   call dealiasxy(w(:,:,1:xsz(3)))

   if(istop)then 
      call dealiasxy(u(:,:,xsz(3)+1))
      call dealiasxy(v(:,:,xsz(3)+1))
      call dealiasxy(w(:,:,xsz(3)+1))
   end if 

   call update_ghost(u, level)
   call update_ghost(v, level)
   call update_ghost(w, level)

   deallocate(etatmp, ettmp, etat, ett, padd)

 end subroutine nlfs_rk2_nn_les

 
 !checked
 !---------------------------------------
 !input:time
 !input from param: ipa, timep, tcp, rdgl
 !input: pp, pas0
 !output: pas
 !---------------------------------------
 subroutine get_pas(time)
   use decomp
   use param, only: ipa, timep, tcp, rdgl, nxs, nys, np1
   ! use navier, only: pp,pas,pas0
   use navier, only: pp
   use hos_param, only : pa_hos, pa0_hos, vps_hos
   implicit none

   real(wp), intent(IN) :: time
   !integer :: i,j
   real(wp), dimension(nxs, nys/np1) :: temp
   real(wp) :: time0,fexp
   
   if(ipa.eq.0)then
      pa_hos =0
   elseif(ipa.eq.1)then
      call padding(pp(:,:,1),pa_hos)
      temp = pa_hos
      pa_hos = 2.0_wp*temp-pa0_hos
      pa0_hos=temp
      if (time < timep)then
         time0=0
      else
         time0=time-timep
      endif
      fexp=exp(-tcp*time0**2.0_wp)
      pa_hos=rdgl*pa_hos*(1.0_wp-fexp)
   else
      print *, "invalid ipa"
   endif
       
 endsubroutine get_pas

 !checked
 !--------------------------------------
 !input: hh
 !output: uzfs, vzfs, w
 !output: eta0, eta
 !--------------------------------------
 subroutine top_fs_les_hos
   use decomp
   use grid, only: eta, eta0, hh
   use navier, only: uzfs,vzfs,w
   implicit none

   eta0 = 0
   eta = eta0 +hh

   if(istop)then
      uzfs =0
      vzfs =0
      w(:,:,xsz(3)-1)=0
   endif
 endsubroutine top_fs_les_hos

 !checked
 !--------------------------------------
 !input: ub, vb, wb
 !output: u, v, w
 !--------------------------------------
 subroutine bc_les
   use decomp
   use navier, only: u,v,w,ub,vb,wb
   implicit none

   if(isbot)then
      u(:,:,1)=ub
      v(:,:,1)=vb
      w(:,:,1)=wb
   endif

   if(istop)then
      u(:,:,xsz(3)+1)=u(:,:,xsz(3)-1)
      v(:,:,xsz(3)+1)=v(:,:,xsz(3)-1)
      u(:,:,xsz(3))=u(:,:,xsz(3)-1)
      v(:,:,xsz(3))=v(:,:,xsz(3)-1)
      w(:,:,xsz(3)-1)=0
      w(:,:,xsz(3))=-w(:,:,xsz(3)-2)
   endif
   
 endsubroutine bc_les

 !checked
 !--------------------------------------
 !input: fs
 !input: nxs, nys, np1
 !output: f
 !--------------------------------------
! subroutine filtering(fs,f)
!   use decomp
!   use param, only: nxs,nys,np1
!   use MPI
!   use fft
!   implicit none
!   include "fftw3.f"
  
!   real(wp), dimension(nxs,nys/np1), intent(IN) :: fs
!   real(wp), dimension(xsz(1),xsz(2)), intent(OUT) :: f 
!   real(wp), allocatable, dimension(:,:) :: temp1
!   real(wp), allocatable, dimension(:,:) :: temp2
!   real(wp), allocatable, dimension(:,:) :: temp3
!   real(wp), allocatable, dimension(:,:) :: temp4
!   real(wp), allocatable, dimension(:,:) :: temp5
!   real(wp), allocatable, dimension(:,:) :: temp6
!   integer :: ierror
!   integer*8 :: plan1, plan2
!   real(wp), dimension(nxs) :: in1  
!   real(wp), dimension(nxs+2) :: out1
!   real(wp), dimension(nys) :: in2  
!   real(wp), dimension(nys+2) :: out2
  
!   integer :: i,j
  
!   allocate(temp1(nxs,nys))
!   allocate(temp2(nxs,nys))
!   allocate(temp3(nx_global+2,nys))
!   allocate(temp4(nys,nx_global+2))
!   allocate(temp5(ny_global,nx_global))
!   allocate(temp6(ysz(1),ysz(2)))
  
!   print*, 'temp1'
!   temp1(:,myid1*nys/np1+1:(myid1+1)*nys/np1)=fs(:,:)
!   print*, 'temp2'
  
!   call MPI_Allreduce(temp1, temp2, nxs*nys, real_type, &
!        MPI_SUM, MPI_COMM_2D_COL, ierror)
  
!   call dfftw_plan_dft_r2c_1d(plan1,nxs,in1,out1,FFTW_ESTIMATE)
!   call dfftw_plan_dft_r2c_1d(plan2,nys,in2,out2,FFTW_ESTIMATE)
  
!   do j =1,nys
!      call dfftw_execute_dft_r2c(plan1, temp2(1:nxs,j), out1)
!      do i=1,nx_global+2
!         temp3(i,j)=out1(i)
!      enddo
!   enddo
  
!   !transpose temp3->temp4
!   temp4 = transpose(temp3)

!   do j =1,nx_global
!      call dfftw_execute_dft_r2c(plan2, temp4(1:nys,j), out2)
!      do i =1,ny_global
!         temp5(i,j)=out2(i)
!      enddo
!   enddo
!   temp5 = temp5/nxs/nys

!   ! temp5->local size array
!   temp6 = temp5(yst(1):yend(1),yst(2):yend(2))

!   call fft_c2r_xy(temp6,f)
  
!   call dfftw_destroy_plan(plan1)
!   call dfftw_destroy_plan(plan2)

! end subroutine filtering
subroutine filtering(nxs, nys, np1, fs,f)
   use decomp
   use fft_hos
   use fft
   use MPI
   ! use param, only: nxs,nys,np1

   implicit none

   integer, intent(IN) :: nxs, nys, np1
   real(wp), dimension(nxs,nys/np1), intent(IN) :: fs
   real(wp), dimension(xsz(1),xsz(2)), intent(OUT) :: f 

   real(wp), dimension(nxs, nys/np1) :: ft1
   real(wp), dimension(xsz(1), nys/np1) :: ft2
   real(wp), dimension(ysz(1), ysz(2)) :: ft3
   real(wp), dimension(nys, xsz(1)/np1) :: temp

   integer :: j

   ft1=fs
   call fft_for_x_hos(ft1)

   do j=1, nys/np1
      ft2(1:xsz(1),j)=ft1(1:xsz(1),j)
   end do 

   call transpose_2d(ft2, temp, xsz(1), nys/np1)
   call fft_for_x_hos(temp, nys, xsz(1)/np1)

   ! print*, ysz(2), xsz(1)/np1

   do j=1,ysz(2) 
      ft3(1:ysz(1),j)=temp(1:ysz(1),j)
   end do 

   call fft_c2r_xy(ft3,f)

end subroutine filtering
!==========================================
 !checked
 !--------------------------------------
 !input: f
 !input: nxs, nys, np1
 !output: fs
 !--------------------------------------
! subroutine padding(f,fs)
!   use decomp
!   use param, only: nxs,nys,np1
!   use MPI
!   use fft
!   implicit none
!   include "fftw3.f"

!   real(wp), dimension(xsz(1),xsz(2)), intent(IN) :: f
!   real(wp), dimension(nxs,nys/np1), intent(OUT) :: fs

!   real(wp), dimension(xsz(1),xsz(2)) :: temp0
!   real(wp), allocatable, dimension(:,:) :: temp1
!   real(wp), allocatable, dimension(:,:) :: temp2
!   real(wp), allocatable, dimension(:,:) :: temp3
!   real(wp), allocatable, dimension(:,:) :: temp4
!   real(wp), allocatable, dimension(:,:) :: temp5
!   real(wp), allocatable, dimension(:,:) :: temp6
!   real(wp), allocatable, dimension(:,:) :: temp7
!   real(wp), allocatable, dimension(:,:) :: temp8

!   integer :: ierror
!   integer*8 :: plan1, plan2

!   real(wp), dimension(nys+2) :: in1
!   real(wp), dimension(nys) ::  out1
!   real(wp), dimension(nxs+2) :: in2
!   real(wp), dimension(nxs) :: out2
  
!   integer :: i,j
  
!   allocate(temp1(ysz(1),ysz(2)))
!   allocate(temp2(ysz(1),ysz(2)))
!   allocate(temp3(ny_global,nx_global))
!   allocate(temp4(nys+2,nx_global))
!   allocate(temp5(nys,nx_global))
!   allocate(temp6(nx_global,nys))  
!   allocate(temp7(nxs+2,nys))
!   allocate(temp8(nxs,nys))

!   call dfftw_plan_dft_r2c_1d(plan1,nys+2,in1,out1,FFTW_ESTIMATE)
!   call dfftw_plan_dft_r2c_1d(plan2,nxs+2,in2,out2,FFTW_ESTIMATE)

!   temp0 = f
!   call fft_r2c_xy(temp0,temp1)

!   temp1 = temp1/nx_global/ny_global

!   temp2(yst(1):yend(1),yst(2):yend(2))=temp1(:,:)
  
!   call MPI_Allreduce(temp2, temp3, ny_global*nx_global, real_type, &
!        MPI_SUM, MPI_COMM_2D_ROW, ierror)

!   do j =1, nx_global
!      do i =1, ny_global
!         temp4(i,j)=temp3(i,j)
!      enddo
!      do i =ny_global+1,nys+2
!         temp4(i,j)=0.
!      enddo
!      call dfftw_plan_dft_c2r_1d(plan1,nys+2,temp4(1:nys+2,j),temp5(1:nys,j),FFTW_ESTIMATE)
!   enddo

!   temp6=transpose(temp5)
  
!   do j =1, nys
!      do i =1, nx_global
!         temp7(i,j)=temp6(i,j)
!      enddo
!      do i =nx_global+1,nxs+2
!         temp7(i,j)=0.
!      enddo
!      call dfftw_plan_dft_c2r_1d(plan2,nxs+2,temp7(1:nxs+2,j),temp8(1:nxs,j),FFTW_ESTIMATE)
!   enddo

!   fs(:,1:nys/np1)=temp8(:,(myid1*nys/np1+1):(myid1+1)*nys/np1)

!   call dfftw_destroy_plan(plan1)
!   call dfftw_destroy_plan(plan2)
    
! end subroutine padding

subroutine padding(f,fs)
  use decomp
  use param, only: nxs,nys,np1
  use MPI
  use fft
  use fft_hos

  implicit none

  real(wp), intent(IN) :: f(xsz(1),xsz(2))
  real(wp), intent(OUT) :: fs(nxs,nys/np1)

  real(wp), dimension(:,:), allocatable :: ft, ft1, ft2, ft3, ft4
  ! real(wp), dimension(xsz(1),xsz(2)) :: ft
  ! real(wp), dimension(nxs,xsz(2)) :: ft1
  ! real(wp), dimension(ysz(1),nxs/np1) :: ft2
  ! real(wp), dimension(nys,nxs/np1) :: ft3
  ! real(wp), dimension(nxs,nys/np1) :: ft4

  integer :: j

  allocate(ft(xsz(1),xsz(2))) 
  allocate(ft1(nxs,xsz(2)))
  allocate(ft2(ysz(1),nxs/np1))
  allocate(ft3(nys,nxs/np1))
  allocate(ft4(nxs,nys/np1))

  ft=f; ft1=0; ft2=0; ft3=0; ft4=0

  call fft_r2c_x(ft)
  ft=ft/xsz(1)

  do j=1, xsz(2)
     ft1(1:xsz(1),j)=ft(1:xsz(1),j)
     ft1(xsz(1)+1:nxs,j)=0
  end do 


  call transpose_2d(ft1, ft2, nxs, xsz(2))

  call fft_for_x_hos(ft2, ysz(1), nxs/np1)

  do j=1, nxs/np1
     ft3(1:ysz(1),j)=ft2(1:ysz(1),j)
     ft3(ysz(1)+1:nys, j)=0
  end do 
  
  call fft_bac_x_hos(ft3, nys, nxs/np1)
  call transpose_2d(ft3, ft4, nys, nxs/np1)
  call fft_bac_x_hos(ft4, fs)

  deallocate(ft) 
  deallocate(ft1)
  deallocate(ft2)
  deallocate(ft3)
  deallocate(ft4)

end subroutine padding
!===================================

subroutine input_hos_par
  use decomp
  use param
  use hos_param
  use constants

  implicit none

   myid_hos = myid1
   nxhos = nxs
   nyhos = nys
   ncpu_hos = np1
   npw = 3
   pex_hos = pex
   pey_hos = pey

   dx_hos = twopi / pex_hos / nxhos
   dy_hos = twopi / pey_hos / nyhos
   dt_hos = dt/ nth

   fr2_hos = fr2
   
end subroutine input_hos_par

subroutine form_drag_k(fpk, betak, aa, fpt)

  use MPI
  use decomp, only : myid, isbot, xst, xend, xsz, nx_global, ny_global, nz_global
  use param, only : pex, usbot
  use navier
  use utils
  use fft_hos
  use spectral_hos
  use grid, only : hh
  use constants, only : twopi

  implicit none

  real(wp), intent(out) :: fpt
  real(wp), intent(out), dimension(nx_global) :: fpk, betak, aa

  real(wp), allocatable, dimension(:) :: fp0, beta0, betak0
  real(wp), allocatable, dimension(:) :: a, b, thetae, thetap, dtheta, tmp_alpha
  real(wp), allocatable, dimension(:,:) :: eta_surf, p1, hhall, fpk0, aa0
  real(wp), allocatable, dimension(:,:) :: ex, ey
  real(wp), allocatable, dimension(:,:,:) :: ppall

  integer :: i,j,k,root0, ierror

  allocate(fp0(nx_global))
  allocate(beta0(nx_global))
  allocate(betak0(nx_global))
  allocate(fpk0(nx_global,ny_global))
  allocate(aa0(nx_global,ny_global))
  allocate(eta_surf(nx_global,1))
  allocate(hhall(nx_global,ny_global))
  allocate(p1(nx_global,1))
  allocate(a(nx_global))
  allocate(b(nx_global))
  allocate(thetae(nx_global))
  allocate(thetap(nx_global))
  allocate(dtheta(nx_global))
  allocate(tmp_alpha(nx_global))
  allocate(ppall(nx_global,ny_global,nz_global))
  allocate(ex(nx_global,ny_global))
  allocate(ey(nx_global,ny_global))

  fp0 = 0
  fpk = 0
  fpk0 = 0
  beta0=0
  betak0=0
  betak=0
  aa0=0
  aa=0

  root0 = 0

  call gather_2d_xy(hh(1:xsz(1),1:xsz(2)),hhall,root0)
  call gather_3d_xyz(pp(1:xsz(1),1:xsz(2),1:xsz(3)),ppall,root0)

  if (myid == root0) then
     do j = 1, ny_global
        eta_surf(:,1) = -hhall(:,j)
        p1(:,1) = ppall(:,j,1)

        call get_amp_phase(eta_surf,a,thetae,nx_global)       
        call get_amp_phase(p1,b,thetap,nx_global)

        do i=2,nx_global/2+1
           if(abs(a(i)) > 1.e-18) then
              tmp_alpha(i)=b(i)/a(i)/(i-1.0)
           else
              tmp_alpha(i)=0.
           endif
           dtheta(i)=thetap(i)-thetae(i)
           if(dtheta(i) > twopi) then
              dtheta(i)=dtheta(i)-twopi
           elseif(dtheta(i) < 0.) then
              dtheta(i)=dtheta(i)+twopi
           endif
        enddo
        
        do i=2,nx_global/2+1
           beta0(i)=-tmp_alpha(i)/usbot**2*sin(dtheta(i))
           fp0(i)=beta0(i)*0.5*(a(i)*(i-1)*pex)**2*usbot**2
        enddo


        !-------------------------------------
        !     REMOVE ZERO MODE
        !     PERFORM SPANWISE AVERAGING
        !-------------------------------------
        
        do i = 1, nx_global/2
           fpk0(i,j) = fpk0(i,j) + fp0(i+1) 
           aa0(i,j) = aa0(i,j) + a(i+1) 
        end do
     end do

     do i = 1, nx_global/2/3*2-1
        aa(i) = sum(aa0(i,:)) / ny_global
        fpk(i) = sum(fpk0(i,:)) / ny_global
        betak(i) = fpk(i) / (0.5*(aa(i)*i*pex)**2*usbot**2)
     end do

     call pdfx_hos(hhall, ex, pex, nx_global, ny_global)

     fpt = 0
     do j = 1, ny_global
        do i = 1, nx_global 
           fpt = fpt + ppall(i,j,1) * ex(i,j) 
        end do
     end do
     fpt = fpt / nx_global / ny_global

  end if

  deallocate(fp0, beta0, betak0, fpk0, aa0, ex, ey)
  deallocate(eta_surf, p1, a, b, thetae, thetap, dtheta, tmp_alpha)

end subroutine form_drag_k

subroutine turb_analysis!(au_wv,aw_wv, phi_u,phi_w)

  use MPI
  use decomp, only : myid, isbot, xst, xend, xsz, nx_global, ny_global, nz_global
  use param, only : pex, usbot
  use navier
  use utils
  use fft_hos
  use spectral_hos
  use grid, only : hh
  use constants, only : twopi

  implicit none

!  real(wp), intent(out), dimension(:) :: au_wv,aw_wv, phi_u,phi_w

  real(wp), allocatable, dimension(:) :: a, amp_u, amp_w, thetae, thetau,thetaw,thetap, dtheta, tmp_alpha
  real(wp), allocatable, dimension(:,:) :: eta_surf, p1, hhall, u1, w1
  real(wp), allocatable, dimension(:,:,:) :: ppall, uall, wall

  real(wp) dtheu(nz_global),dthew(nz_global), au(nz_global), aw(nz_global)

  integer :: i,j,k,root0, ierror

  allocate(eta_surf(nx_global,1))
  allocate(hhall(nx_global,ny_global))
  allocate(u1(nx_global,1))
  allocate(w1(nx_global,1))
  allocate(a(nx_global))
  allocate(amp_u(nx_global))
  allocate(amp_w(nx_global))
  allocate(thetae(nx_global))
  allocate(thetap(nx_global))
  allocate(thetau(nx_global))
  allocate(thetaw(nx_global))
  allocate(tmp_alpha(nx_global))
  allocate(ppall(nx_global,ny_global,nz_global))
  allocate(uall(nx_global,ny_global,nz_global))
  allocate(wall(nx_global,ny_global,nz_global))

  root0 = 0

  call gather_2d_xy(hh(1:xsz(1),1:xsz(2)),hhall,root0)
  call gather_3d_xyz(u(1:xsz(1),1:xsz(2),1:xsz(3)),uall,root0)
  call gather_3d_xyz(w(1:xsz(1),1:xsz(2),1:xsz(3)),wall,root0)
  call gather_3d_xyz(pp(1:xsz(1),1:xsz(2),1:xsz(3)),ppall,root0)

  if (myid == root0) then
     dtheu = 0
     dthew = 0
     au = 0
     aw = 0
     do k = 1, nz_global
        
        do j = 1, ny_global
           eta_surf(:,1) = -hhall(:,j)
           u1(:,1) = uall(:,j,k)
           w1(:,1) = wall(:,j,k)

           call get_amp_phase(eta_surf,a,thetae,nx_global)          
           call get_amp_phase(u1,amp_u,thetau,nx_global)
           call get_amp_phase(w1,amp_w,thetaw,nx_global)
           
           dtheu(k) = dtheu(k) + (thetau(5) - thetae(5)) / ny_global
           dthew(k) = dthew(k) + (thetaw(5) - thetae(5)) / ny_global
           au(k) = au(k) + amp_u(5) / ny_global
           aw(k) = aw(k) + amp_w(5) / ny_global
        end do        
        write(296,'(25e12.4)') k*1.0,dtheu(k),dthew(k),au(k),aw(k)
     end do

  end if

!  deallocate(hhall,uall,wall,ppall)
!  deallocate(eta_surf, u1, w1, a, thetae, thetap, dtheta, tmp_alpha)

end subroutine turb_analysis


subroutine get_beta_miles(wvn, beta)

  use decomp
  use param
  use hos_param
  use constants

  implicit none
  
  real(wp), intent(out) :: beta
  real(wp), intent(in) :: wvn
  
  real(wp) zc, zetac, kappa

  kappa = 0.41
  zc = z0 * exp((kappa / usbot ) * sqrt(1 / fr2 / wvn))
  zetac = zc * wvn
  if (zetac > 0.281) then
     beta = 2 * log(0.281  / zetac)
  else
     beta = (pi * zetac / kappa**2) * (log(0.281  / zetac))**4 + &
          2 * log(0.281  / zetac)
  end if

end subroutine get_beta_miles

