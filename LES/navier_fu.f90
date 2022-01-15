!all checked
!input & output checked
!========================================================
   !> @brief Calculate FU, FV, FW.
   !
   !> @param[in]  time    time
   !> @param[in]  dt      time step
   !> @param[in]  Re      Reynolds number
   !> @param[in]  bforce  linear forcing coefficient
   !> @param[in]  tf1,tf2 shear force in x and y direction
   !========================================================

subroutine calc_F(time)
  
  use param, only : pex, pey 
  use constants, only : TWOPI
  use spectral

  !> added by plyu
  use turbine_model
  use discontinuity_smooth, only: analyze_limit3d_from_3d, tbn_lim_x, &
    tbn_lim_y, lim_dir

  implicit none
  
  real(wp), intent(IN) :: time
      
  real(wp) :: deltax, deltay
  integer :: i, j
  
  !added by plyu, 12/14/2016
  !calculate wind turbine force
  !print *, 'plyudebug, iturbine,nxwt,nywt=', iturbine, nxwt, nywt 
  if (iturbine .eq. 1) then
    call discloca_v3
    call veldisc_v2(u, level, time)
    call wind_turbine_force(wtforce, level, fturbinex)      
  else if(iturbine .eq. 3) then
    call rotor_model_acl(time, u, v, w, wtforce, wtforce_y, wtforce_z, &
      fturbinex, fturbiney, level, flag_de)
  else if(iturbine .eq. 5) then
    call rotor_model_acs(time, u, v, w, wtforce, wtforce_y, wtforce_z, &
      fturbinex, fturbiney, level, flag_de)
  else if(iturbine .eq. 7) then
    call rotor_model_admr(time, u, v, w, wtforce, wtforce_y, wtforce_z, &
      fturbinex, fturbiney, level, flag_de)
  endif

  !> added by plyu: identify locations of discontinuity (assume that turbine
  !  force causes this discontinuity)
  if (ids>0 .and. (iturbine.eq.3 .or. iturbine.eq.5 .or. iturbine.eq.7)) then
    lim_dir = 1 !> indicate it is analyzing x direction
    call analyze_limit3d_from_3d (flag_de, tbn_lim_x)
    call transpose_xy(flag_de, bufy)
    lim_dir = 2 !> indicate it is analyzing y direction
    call analyze_limit3d_from_3d (bufy, tbn_lim_y)
    lim_dir = 0
  endif
  
  !defined in navier.f/navier_init:
  !input: u, v, w (navier_init)
  !output: ux, uy, uz_u (navier_init)
  call gradu(u, ux, uy, uz_u)!; print *, 'after gradu'
  call gradu(v, vx, vy, vz_u)!; print *, 'after gradv'
  call gradw(w, wx, wy, wz_w)!; print *, 'after gradw'
    
  deltax = 2.0_wp * TWOPI/pex/nx_global
  deltay = 2.0_wp * TWOPI/pey/ny_global
  
  uf = u; vf = v; wf = w

  ! if(myid==0)print*,'u1', w(1,1,:)

  !common.f/solver_common module
  !mfilt: filter type indicator, initialized in navier_les/les_init
  !mfilt=1 for now
  call les_filter(uf(:,:,1:xsz(3)), mfilt, deltax, deltay)
  call les_filter(vf(:,:,1:xsz(3)), mfilt, deltax, deltay)
  call les_filter(wf(:,:,1:xsz(3)), mfilt, deltax, deltay)
  
  call update_ghost(uf, level)
  call update_ghost(vf, level)
  call update_ghost(wf, level)

  if (istop) then
     call les_filter(uf(:,:,xsz(3)+1), mfilt, deltax, deltay)
     call les_filter(vf(:,:,xsz(3)+1), mfilt, deltax, deltay)
     call les_filter(wf(:,:,xsz(3)+1), mfilt, deltax, deltay)
  end if

   ! if(myid==0)print*,'uf', uf(1,1,:)
   ! if(myid==2)print*,'uf', uf(1,1,:)
  ! if(myid==2)print*,'uf2', wf(1,1,:)

  ! do j=1, xsz(2)
  !    do i=1, xsz(1)
  !       read(104,*)u(i,j,32)
  !    end do 
  ! end do 
   ! if(myid==0) print*,'u1', u(1:xsz(1),4,32)
   ! call dealiasxy(u(:,:,32))
   ! if(myid==0) print*,'u2', u(1:xsz(1),4,32)

  !navier_fu.f
  !apply wall model
  !output: tauwx, tauwy
  if(isbot)then
     !if (iwavy .eq. 1 .and. dabs(hka) .le. 1e-6) then
     !  !> plyunote: for planar bottom, try this
     !  call wall_model_planar_average
     !else
       call wall_model_v3
     !endif
  endif

  !navier_les.f
  !input: u, v, w
  !output: s11~s33, s11w~s33w
  call get_strain(u, v, w)!; print *, 'after get_strain'

  !input: u, v, w, uf, vf, wf, deltax, deltay
  !output: nut
  call get_nut(u, v, w, uf, vf, wf, deltax, deltay)

  !output: nutw
  call get_nutw(u, v, w, uf, vf, wf, deltax, deltay)

   ! print*, nutw(1,1,1:xsz(3))
   ! if(myid==0)print*,'nutw', nutw(1,1,:)
   ! if(myid==1)print*,'nutw', nutw(1,1,:)
  nut = nut + 1.0_wp/re
  nutw = nutw +1.0_wp/re

  ! print*, 1.0_wp/re
  ! print*, 'nut', nut(1,1,1:nz)
  ! print*, 'nutw', nutw(1,1,1:nz)
  ! print*, re

  !input: s11~s33, s11w~s33w, nut, nutw
  !output: t11~t33, t11w~t33w
  call get_SGS_stress
  if (myid==0) print *, 'SGS stress computed.'

  !navier_fu.f
  !input: tauwx, tauwy

  !input&output: hu, hv, hw
  !input: u, v, w
  !input: ux, uy, uz_u
  !input: vx, vy, vz_u
  !intput: wx, wy, wz_w
  !input: nut
  !input: t11, t12, t13w, t22, t23w
  !input: wtforce (navier_init)
  !output: fu, fv, fw
  !output(wall model): t13w(:,:,1), t23w(:,:,1)
  !output: t13w(:,:,xsz(3)),t23w(:,:,xsz(3))
  !output: t11(:,:,xsz(3)+1), t12(:,:,xsz(3)+1)
  !output: t22(:,:,xsz(3)+1)
  
  call fun_u_nn_cb(fu,hu,time,tauwx)!; print *, 'after fun_u'
  call fun_v_nn_cb(fv,hv,time,tauwy)!; print *, 'after fun_v'
  call fun_w_nn_cb(fw,hw,time)!; print *, 'after fun_w'
  ! if(myid==0)print*, fw(1,1,:)
  ! if(myid==1)print*, fw(1,1,:)
  
  if(myid==0) print*, 'prediction!'
  call velocity_prediction(fu,fv,fw)

  !> added by plyu, parameters refer to module turbine_model
  !print *, 2, ti
  if (iturbine .gt. 0 .and. time>0 .and. (InletRelaxation .eq. 1 .or. &
    InletRelaxation >= 3)) call Add_Velocity_Relaxation_Block
  
end subroutine calc_F

subroutine coeff_stepfun_linear(x_in, x_out)
  implicit none
  real(wp) :: x_in, x_out
  
  x_out = x_in
end subroutine coeff_stepfun_linear

subroutine coeff_stepfun_smoothheaviside(x_in, x_out)
  implicit none
  real(wp) :: x_in, x_out
  
  !> plyunote: search this on google to see the plot:
  !! y=0.5*(1.0+(x-0.5)/0.5+1.0/3.1415926*sin((x-0.5)*3.1415926/0.5))
  x_out = 0.5*(1.0+(x_in-0.5)/0.5+1.0/3.1415926*sin((x_in-0.5)*3.1415926/0.5))
end subroutine coeff_stepfun_smoothheaviside

subroutine coeff_stepfun(x_in, x_out, i_stepfun)
  implicit none
  real(wp) :: x_in, x_out
  integer :: i_stepfun

  if (i_stepfun .eq. 1) then
    call coeff_stepfun_linear(x_in, x_out)
  elseif (i_stepfun .eq. 2) then
    call coeff_stepfun_smoothheaviside(x_in, x_out)
  else
    print *, "coeff_stepfun: i_stepfun not implemented."
  endif
end subroutine coeff_stepfun

subroutine Add_Velocity_Relaxation_Block
  use param
  use decomp
  use spectral
  use turbine_model
  !use mpi
  implicit none

  !real(wp), allocatable, dimension(:,:) :: 
  real(wp) :: var1, var2, var3
  real(wp), allocatable, dimension(:,:,:,:) :: veltar, velsec
  real(wp), allocatable, dimension(:) :: vel_send

  integer :: i, j, k, imin, imax
  integer :: ti_group, ti_infile, v_i, data_pos
  character(len=64) :: filen, groupname
  real(wp) :: z
  real(kind=4) :: bin_temp
  real(wp), allocatable, dimension(:) :: c_smooth
  
  integer :: ierr, recl_size, file_size

  if (restop .gt. 1.e-6) then
    print *, "plyudebug: This case is not yet included in Add_Inlet_Forcing"
  endif

  !> IR_bw: block width. it has been moved to input file control.dat
  !IR_bw = 9  

  !> plyunote: read from binaries
  !if (ti == ti_first + 1) then 
    allocate(veltar(IR_bw, xsz(2), 1-level:xsz(3)+level, 3))
    allocate(velsec(IR_bw, ny_global, nz_global, 3))
    allocate(vel_send(IR_bw*ny_global*nz_global*3))
    veltar(:,:,:,:)=0.0_wp; velsec(:,:,:,:)=0.0_wp; vel_send(:) = 0.0_wp

    allocate(c_smooth(IR_bw))
  !endif
  
  !> bw is pre-assigned here. It might also be read from user input
  !! here 1 is to downgrade to slice mode;
  !! and 7 consists of a 3-grid central band and two 2-grid wing band,
  !! and 9 consists of a 3-grid central band and two 3-grid wing band
  if (IR_bw == 1) then
    c_smooth = (/1.0/)
  elseif (IR_bw == 7) then
    c_smooth = (/0.25,0.75,1.0,1.0,1.0,0.75,0.25/)
  elseif (IR_bw == 9) then
    !c_smooth = (/0.166667, 0.5, 0.833333, 1.0, 1.0, 1.0, 0.833333, 0.5, 0.166667/)
    c_smooth = (/0.166667, 0.5, 0.833333, 1.0, 1.0, 1.0, 1.0, 0.75, 0.25/)
  endif

  ti_group = (ti+IR_ti_shift-1)/1000+1
  ti_infile = ti - (ti_group-1)*1000

  write(filen, '(a,i0.4,a,i0.6,a4)') '../inflowdata/bin_', IR_dataindex, &
  '/inlet_', ti_group, '.bin'
  
  if (myid .eq. 0) then
    open(16001,FILE=filen,STATUS='OLD',ACTION='READ',form='unformatted',&
      access='direct',recl=4,iostat=ierr)
    if (ierr /= 0) then
      print *, "File not found:", filen 
    endif
    
    !inquire(16001, recl=recl_size, size=file_size)
    !print *, 'recl=', recl_size, 'file_size=',file_size, &
    !  ', expected nt=1000, bw=',bw,', ny=',&
    !  ny_global, ', nz=', nz_global, ', nvec=3', &
    !  ', total=recl*nt*bw*ny*nz*nvec=', &
    !  recl_size*1000*bw*ny_global*nz_global*3

    do i = 1, IR_bw
      do j = 1, ny_global
        do k = 1, nz_global
          do v_i = 1, 3
            data_pos = v_i+3*(k-1+nz_global*(j-1+ny_global*(i-1+IR_bw*(ti_infile-1))))
            !print *, 'plyudebug, ti_infile=',ti_infile,", j=",j,&
            !  ',k=',k,',v_i=',v_i,',data_pos=',data_pos
            read(16001, rec=data_pos, iostat=ierr) bin_temp
            !if (ierr .ne. 0) then
            !  print *, filen, i, j, k, v_i, data_pos
            !endif
            velsec(i,j,k,v_i) = bin_temp
          enddo
        enddo
      enddo
    enddo

    close(16001)
  endif

  vel_send = reshape(velsec, (/IR_bw*ny_global*nz_global*3/))
  
  !> plyunote: velsec is target velocity in global coordinate (ny,nz)
  call mpi_bcast(vel_send, IR_bw*ny_global*nz_global*3, mpi_double_precision, 0, &
    mpi_comm_2d_cart, ierr)
  velsec = reshape(vel_send, (/IR_bw, ny_global, nz_global, 3/))
  
  if (isbot .and. istop) then
    veltar(:,:,2:xsz(3),:) = velsec(:, xst(2):xst(2)+xsz(2)-1, &
      xst(3)+1:xst(3)+xsz(3)-1, 1:3)
  elseif (isbot ) then
    veltar(:,:,2:xsz(3)+level,:) = velsec(:, xst(2):xst(2)+xsz(2)-1, &
      xst(3)+1:xst(3)+xsz(3)-1+level, 1:3)
  elseif(istop) then
    veltar(:,:,1-level:xsz(3),:) = velsec(:, xst(2):xst(2)+xsz(2)-1, &
      xst(3)-level:xst(3)+xsz(3)-1, 1:3)
  else
    veltar(:,:,:,:) = velsec(:, xst(2):xst(2)+xsz(2)-1, &
      xst(3)-level:xst(3)+xsz(3)-1+level, 1:3)
  endif
  
  !> Get du_relax(xsz(1),xsz(2),1-level:xsz(3)+level,3)
  imin = 1; imax = IR_bw
  if (xst(1)<=imax .and. (xst(1)+xsz(1)-1)>=imin) then
    do k = 1-level, xsz(3)+level
      z = zz(k) * hbar
      if (isbot .and. (k==1 .or. k==0)) then

      elseif(istop .and. k==xsz(3)+level) then

      ! for wave cases, z might less than 0, hard to use log law.
      elseif(z>0.0_wp .and. zz(k)<1.0_wp) then
        !print *, 'myid, k_local, k_global, dz, z=',myid,k,xst(3)+k-1,dz(k),z
        do i = 1, IR_bw
          do j = 1, xsz(2)
            u(i,j,k) = veltar(i,j,k,1)*c_smooth(i)+u(i,j,k)*(1.0-c_smooth(i))
            v(i,j,k) = veltar(i,j,k,2)*c_smooth(i)+v(i,j,k)*(1.0-c_smooth(i))
            w(i,j,k) = veltar(i,j,k,3)*c_smooth(i)+w(i,j,k)*(1.0-c_smooth(i))
          enddo
        enddo
      endif
    enddo
      
  endif
  
  deallocate(veltar, velsec, vel_send)
  deallocate(c_smooth)
     
  !call dealiasxy(u(:,:,1:xsz(3)))
  !call dealiasxy(v(:,:,1:xsz(3)))
  !call dealiasxy(w(:,:,1:xsz(3)))  

  !if(istop)then 
  !   call dealiasxy(u(:,:,xsz(3)+1))
  !   call dealiasxy(v(:,:,xsz(3)+1))
  !   call dealiasxy(w(:,:,xsz(3)+1))  
  !end if 

  ! if(myid==2)print*,'xsz(3)',size(u,3) 
  call MPI_BARRIER(mpi_comm_world, ierr)
  call update_ghost(u, level)
  call update_ghost(v, level)
  call update_ghost(w, level)
  call MPI_BARRIER(mpi_comm_world, ierr)

end subroutine Add_Velocity_Relaxation_Block

!subroutine Add_Velocity_Relaxation_Slice
!  use param
!  use decomp
!  use spectral
!  use turbine_model
!  !use mpi
!  implicit none
!
!  !real(wp), allocatable, dimension(:,:) :: 
!  real(wp) :: var1, var2, var3
!  real(wp), allocatable, dimension(:,:,:) :: veltar, velsec
!  real(wp), allocatable, dimension(:) :: vel_send
!  integer :: it ! ti_shift has been moved to turbine_model parameters
!
!  integer :: i, j, k, imin, imax
!  integer :: ti_group, ti_infile, v_i, data_pos
!  character(len=64) :: filen, groupname
!  real(wp) :: z
!  real(kind=4) :: bin_temp
!  
!  integer :: ierr
!
!  if (restop .gt. 1.e-6) then
!    print *, "plyudebug: This case is not yet included in Add_Inlet_Forcing"
!  endif
!
!  !> plyunote: read from binaries
!  allocate(veltar(xsz(2), 1-level:xsz(3)+level, 3))
!  allocate(velsec(ny_global, nz_global, 3))
!  allocate(vel_send(ny_global*nz_global*3))
!  veltar(:,:,:)=0.0_wp; velsec(:,:,:)=0.0_wp; vel_send(:) = 0.0_wp
!
!  ti_group = (ti+IR_ti_shift-1)/1000+1
!  ti_infile = ti - (ti_group-1)*1000
!
!  write(filen, '(a,i0.4,a,i0.6,a4)') '../inflowdata/bin_', IR_dataindex, &
!  '/inlet_', ti_group, '.bin'
!
!  if (myid .eq. 0) then
!    open(16001,FILE=filen,STATUS='OLD',ACTION='READ',form='unformatted',&
!      access='direct',recl=4)
!
!    do j = 1, ny_global
!      do k = 1, nz_global
!        do v_i = 1, 3
!          data_pos = v_i+3*(k-1+nz_global*(j-1+ny_global*(ti_infile-1)))
!          !print *, 'plyudebug, ti_infile=',ti_infile,", j=",j,&
!          !  ',k=',k,',v_i=',v_i,',data_pos=',data_pos
!          read(16001,rec=data_pos) bin_temp
!          velsec(j,k,v_i) = bin_temp
!        enddo
!      enddo
!    enddo
!
!    close(16001)
!  endif
!
!  vel_send = reshape(velsec, (/nz_global*ny_global*3/))
!
!  !> plyunote: velsec is target velocity in global coordinate (ny,nz)
!  call mpi_bcast(vel_send, nz_global*ny_global*3, mpi_double_precision, 0, &
!    mpi_comm_2d_cart, ierr)
!  velsec = reshape(vel_send, (/ny_global, nz_global, 3/))
!
!  if (isbot ) then
!    veltar(:,2:xsz(3)+level,:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!      xst(3)+1:xst(3)+xsz(3)-1+level, 1:3)
!  elseif(istop) then
!    veltar(:,1-level:xsz(3),:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!      xst(3)-level:xst(3)+xsz(3)-1, 1:3)
!  else
!    veltar(:,:,:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!      xst(3)-level:xst(3)+xsz(3)-1+level, 1:3)
!  endif
!
!  !> Get du_relax(xsz(1),xsz(2),1-level:xsz(3)+level,3)
!  imin = 1; imax = 1
!  if (xst(1)<=imax .and. (xst(1)+xsz(1)-1)>=imin) then
!    do k = 1-level, xsz(3)+level
!      z = zz(k) * hbar
!      if (isbot .and. (k==1 .or. k==0)) then
!
!      elseif(istop .and. k==xsz(3)+level) then
!
!      ! for wave cases, z might less than 0, hard to use log law.
!      elseif(z>0.0_wp) then
!        !print *, 'myid, k_local, k_global, dz, z=',myid,k,xst(3)+k-1,dz(k),z
!        i = 1
!        do j = 1, xsz(2)
!          u(i,j,k) = veltar(j,k,1)
!          v(i,j,k) = veltar(j,k,2)
!          w(i,j,k) = veltar(j,k,3)
!        enddo
!
!      endif
!    enddo
!      
!  endif
!    
!  deallocate(veltar, velsec, vel_send)
!     
!  !call dealiasxy(u(:,:,1:xsz(3)))
!  !call dealiasxy(v(:,:,1:xsz(3)))
!  !call dealiasxy(w(:,:,1:xsz(3)))  
!
!  !if(istop)then 
!  !   call dealiasxy(u(:,:,xsz(3)+1))
!  !   call dealiasxy(v(:,:,xsz(3)+1))
!  !   call dealiasxy(w(:,:,xsz(3)+1))  
!  !end if 
!
!  ! if(myid==2)print*,'xsz(3)',size(u,3) 
!  call MPI_BARRIER(mpi_comm_world, ierr)
!  call update_ghost(u, level)
!  call update_ghost(v, level)
!  call update_ghost(w, level)
!  call MPI_BARRIER(mpi_comm_world, ierr)
!
!end subroutine Add_Velocity_Relaxation_Slice
!
!
!!> added by plyu, to restore velocity at upstream to certain value
!!! Controled by parameter InletRelaxation:
!!! 1: A log-law mean profile + ZERO turbulence
!!! 4: Import a turbulent infow database generated from Toni's code (Mann's method) 
!!! 5 or 3: Import a turbulent inflow database generated by precursor simulations. 
!!!         The difference between 5 and 3 is the source folder location.
!!! 6: Similar with 5, but we are reading a block region (9 grids streamwise),
!!!    instead of a one-grid thickness slice in case 3.
!subroutine Add_Velocity_Relaxation_old
!  use param
!  use decomp
!  use spectral
!  use turbine_model
!  !use mpi
!  implicit none
!  !real(wp), dimension(:,:,:) :: u(xsz(1), xsz(2), 1-level:xsz(3)+level)
!  !real(wp), dimension(:,:,:) :: v(xsz(1), xsz(2), 1-level:xsz(3)+level)
!  !real(wp), dimension(:,:,:) :: w(xsz(1), xsz(2), 1-level:xsz(3)+level)
!  !real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), 1-level:xsz(3)+level)
!  !real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), 1-level:xsz(3)+level)
!  !real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), 1-level:xsz(3)+level)
!  !integer :: level
!
!  !> u_tgt: target velocity
!  real(wp) :: r, xi, f_max, utgt_fmax, u_fmax, urelax_fmax, ir_vol, ir_turb_vol
!  real(wp) :: ir_rms, dvol, dturb, ir_c_rms, ir_c1, ir_usum
!  real(wp), allocatable, dimension(:) :: rvec, su
!  real(wp) :: u_relax(3), temp_out(3), temp_in(3)
!  real(wp) :: IR_len, IR_ws1, IR_ws2, temp_coeff
!  !real(wp), allocatable, dimension(:,:) :: f_relax, f_temp
!  
!  !real(wp), allocatable, dimension(:,:) :: 
!  real(wp) :: var1, var2, var3, var4, var5, var6
!  integer :: int1, int2, int3
!  real(wp), allocatable, dimension(:,:,:) :: veltar, velsec
!  real(wp), allocatable, dimension(:,:,:,:) :: veltar_b, velsec_b
!  integer :: it ! ti_shift has been moved to turbine_model parameters
!
!  real(wp), allocatable, dimension(:,:,:,:) :: du_relax, du_temp
!  real(wp) :: zlbot, z, zs, u1, u2
!  integer :: i, j, k, j1, imin, imax, ismooth, i_fmax, j_fmax, k_fmax
!  integer :: ir_num, DynamicIRCoeff
!  character(len=64) :: filen, groupname
!
!  real(wp), allocatable, dimension(:) :: vel_send
!  integer :: ierr
!
!  f_max = 0.0_wp; i_fmax=1; j_fmax=1; k_fmax=1
!  u_fmax=0.0_wp; utgt_fmax=0.0_wp; urelax_fmax=0.0_wp
!
!  ir_num = 0; ir_vol = 0.0_wp; ir_turb_vol = 0.0_wp; ir_rms = 0.0_wp
!  ir_usum = 0.0_wp
! 
!  !nsmooth = 2  
!  DynamicIRCoeff = 0
!   
!  allocate(su(1-level:xsz(3)+level))
!  allocate(rvec(xsz(1)))
!  allocate(veltar(xsz(2), 1-level:xsz(3)+level, 3))
!  allocate(du_relax(xsz(1),xsz(2),1-level:xsz(3)+level,3))
!  allocate(du_temp(xsz(1),xsz(2),1-level:xsz(3)+level,3))
!
!  su(:) = 0.0_wp
!  du_relax(:,:,:,:) = 0.0_wp
!  veltar(:,:,:) = 0.0_wp 
!
!  imin = floor(IR_start/(xl/nx))+1
!  imax = floor(IR_end/(xl/nx))+1
!  if (InletRelaxation .eq. 7) then
!    imin = 1
!    imax = 1
!  endif
!
!  if (InletRelaxation .eq. 6 .and. (imax-imin)>9) then
!    print *, "InletRelaxation error, width(", imax-imin,") should no be larger than",&
!      " inflow database limit (9 grids)"
!  endif
!
!  !print *, 'IR_1'    
!  if (restop .le. 1.e-6) then
!    zlbot = hbar * resbot/(resbot+restop)
!    !zlsbot = 1.0_wp / resbot
!    !zlstop = 0.0_wp
!    usbot = 1.0_wp / (2.5_wp * log(hbar/z0))
!    !ustop = 0.0_wp
!    !print *, 'IR_2'
!
!    !> Get su(xsz(3))
!    !utop = usbot * (2.5_wp*log(hbar/z0))
!    do k = 1-level, xsz(3)+level
!      z = zz(k) * hbar
!      ! for wave cases, z might less than 0, hard to use log law.
!      if (z .le. zlbot .and. z>0.0) then
!        zs = z/zlsbot
!        u1 = usbot * zs
!        u2 = usbot * (2.5_wp * log(z/z0))
!        if (u1 .lt. u2) su(k) = u1
!        if (u1 .ge. u2) su(k) = u2
!        if (u2 .lt. 0.0_wp) su(k) = u1
!      endif        
!    enddo
!
!    !> Get rvec(xsz(1)), it is the weight of original velocity in synthetic vel.
!    rvec(:) = 0.0_wp
!    IR_len = IR_end - IR_start
!    IR_ws1 = IR_start + 0.2 * IR_len
!    IR_ws2 = IR_start + 0.8 * IR_len
!
!    if (xst(1)<=imax .and. (xst(1)+xsz(1)-1)>=imin) then
!      do i = 1, xsz(1)
!        xi = (xst(1)+i-1)*(xl/nx)
!        
!        !> plyunote: original simple version, the weight linear decreases
!        !if (xi>=IR_start .and. xi <= IR_end) then
!        !  rvec(i) = (IR_end - xi)/(IR_end - IR_start)
!        !endif
!
!        !> plyunote: 3-segment version, decrease-stable-increase
!        if (xi >= IR_start .and. xi <= IR_ws1) then
!          !> plyunote: in first segment, rvec decreases from 1.0 to 0.0
!          temp_coeff = (xi - IR_start) / (IR_ws1 - IR_start) !> scaling
!          call coeff_stepfun(temp_coeff, rvec(i), 2)
!          rvec(i) = 1.0 - rvec(i)
!        elseif (xi >= IR_ws1 .and. xi <= IR_ws2) then
!          !> plyunote: in second segment, rvec is constantly 0.0
!          rvec(i) = 0.0
!        elseif (xi >= IR_ws2 .and. xi <= IR_end) then
!          !> plyunote: in third segment, rvec increases from 0.0 to 1.0
!          temp_coeff = (xi - IR_ws2) / (IR_end - IR_ws2)
!          call coeff_stepfun(temp_coeff, rvec(i), 2)
!        
!        !> plyunote: following lines are used for debug.
!        !else
!        !    if(myid.eq.0) then
!        !      print *, "Velocity_Relaxation error: rvec not decided."
!        !      print *, "Debug info: i=",i,", imin=",imin,", imax=",imax, &
!        !        ", xi=",xi,", IR_start=", IR_start,", IR_end=",IR_end, &
!        !        ", IR_ws1=", IR_ws1, ", IR_ws2=", IR_ws2
!        !    endif
!        endif
!      enddo
!      if((ti-ti_first)<100) then
!        rvec(:) = 1.0*(1.0-(ti-ti_first)/100.0)+rvec(:) *(ti-ti_first)/100.0 
!      endif
!    endif
!    
!    if(InletRelaxation .eq. 1) then
!      do k = 1-level, xsz(3)+level
!        veltar(:, k, 1) = su(k)
!        veltar(:, k, 2:3) = 0.0_wp
!      enddo
!    else if (InletRelaxation .eq.3 .or. InletRelaxation .eq. 5) then
!      allocate(velsec(ny_global, nz_global, 3))
!      allocate(vel_send(ny_global*nz_global*3))
!      velsec(:,:,:) = 0.0_wp; vel_send(:) = 0.0_wp
!      
!      !IR_ti_shift=35000
!      if (InletRelaxation .eq. 3) then
!        write(filen, '(a,i0.10,a4)') '../admr02/output_inlet/inlet_',ti+IR_ti_shift,'.dat'
!      elseif(InletRelaxation .eq. 5) then
!        write(groupname, *) (ti+IR_ti_shift-1)/1000
!        write(filen, '(a,i0.4,a,a,a,i0.10,a4)') '../inflowdata/case_', IR_dataindex, &
!          '/group_', trim(adjustl(groupname)), '/inlet_',ti+IR_ti_shift,'.dat'
!      else
!        print *, "InletRelaxation type not implemented"
!      endif
!
!      open(16001,FILE=filen,STATUS='OLD',ACTION='READ')
!      read(16001,*)
!      read(16001,*)
!      do k = 1, nz_global
!        do j = 1, ny_global
!          read(16001,*) var1, var2, var3, var4, var5, var6
!          velsec(j,k,1:3)=(/var4, var5, var6/)
!        enddo
!      enddo
!      close(16001)
!      if (isbot ) then
!        veltar(:,2:xsz(3)+level,:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!        xst(3)+1:xst(3)+xsz(3)-1+level, 1:3)
!      elseif(istop) then
!        veltar(:,1-level:xsz(3),:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!        xst(3)-level:xst(3)+xsz(3)-1, 1:3)
!      else
!        veltar(:,:,:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!        xst(3)-level:xst(3)+xsz(3)-1+level, 1:3)
!      endif
!    else if (InletRelaxation .eq. 6) then
!      allocate(veltar_b(9, xsz(2), 1-level:xsz(3)+level, 3))
!      allocate(velsec_b(9, ny_global, nz_global, 3))
!      allocate(vel_send(9*ny_global*nz_global*3))
!      veltar_b(:,:,:,:)=0.0_wp; velsec_b(:,:,:,:)=0.0_wp; vel_send(:) = 0.0_wp
!      
!      !IR_ti_shift=212000
!      !if(InletRelaxation .eq. 6) then
!        write(groupname, *) (ti+IR_ti_shift-1)/1000+1
!        write(filen, '(a,i0.4,a,a,a,i0.10,a4)') '../inflowdata/case_', IR_dataindex, &
!          '/group_', trim(adjustl(groupname)), '/inlet_',ti+IR_ti_shift,'.dat'
!      !else
!      !  print *, "InletRelaxation type not implemented"
!      !endif
!
!      if (myid .eq. 0) then
!        open(16001,FILE=filen,STATUS='OLD',ACTION='READ')
!        read(16001,*)
!        read(16001,*)
!        do k = 1, nz_global
!          do j = 1, ny_global
!            do i = 1, 9
!              read(16001,*) var1, var2, var3, var4, var5, var6
!              velsec_b(i,j,k,1:3)=(/var4, var5, var6/) * IR_scale
!              !> plyunote: IR_scale should be around 1.0, if it is different
!              !from 1.0, based on Taylor frozen theory, the different advection
!              !speed will require you to do a interpolation of the inflow
!              !database, which is InletRelaxation mode 7.
!            enddo
!          enddo
!        enddo
!        close(16001)
!      endif
!      
!      vel_send = reshape(velsec_b, (/nz_global*ny_global*9*3/))
!      
!      !> plyunote: velsec is target velocity in global coordinate (ny,nz)
!      call mpi_bcast(vel_send, nz_global*ny_global*9*3, mpi_double_precision, 0, &
!        mpi_comm_2d_cart, ierr)
!      velsec_b = reshape(vel_send, (/9, ny_global, nz_global, 3/))
!
!      if (isbot ) then
!        veltar_b(:,:,2:xsz(3)+level,:) = velsec_b(1:9,xst(2):xst(2)+xsz(2)-1, &
!          xst(3)+1:xst(3)+xsz(3)-1+level, 1:3)
!      elseif(istop) then
!        veltar_b(:,:,1-level:xsz(3),:) = velsec_b(1:9,xst(2):xst(2)+xsz(2)-1, &
!          xst(3)-level:xst(3)+xsz(3)-1, 1:3)
!      else
!        veltar_b(:,:,:,:) = velsec_b(1:9,xst(2):xst(2)+xsz(2)-1, &
!          xst(3)-level:xst(3)+xsz(3)-1+level, 1:3)
!      endif
!    else if (InletRelaxation .eq. 4) then
!      allocate(velsec(ny_global, nz_global, 3))
!      allocate(vel_send(ny_global*nz_global*3))
!      velsec(:,:,:) = 0.0_wp; vel_send(:) = 0.0_wp
!
!      if (myid.eq.0) then
!        IR_ti_shift = ((ti-1) / 100) *100 + 1
!        !write(filen, '(a,i0.4,a,i0.6,a)') '../inflowdata/case_', IR_dataindex, &
!        !  '/inflow_', ti_shift, '_dt=0.0845204.txt'
!        write(filen, '(a,i0.4,a,i0.6,a)') '../inflowdata/case_', IR_dataindex, &
!          '/inflow_', IR_ti_shift, '.txt'
!        !if(myid.eq.0) then
!        print*, 'read inlet file', filen
!        !endif
!        open(16002,file=filen,status='old',action='read')
!        do it = 1, (ti - IR_ti_shift)
!          if (ti>IR_ti_shift) then
!            do k = 1, nz_global+1
!              do j = 1, ny_global+1
!                read(16002, *)
!                !read(16002,*) int1, int2, int3, var4, var5, var6
!                !if(myid.eq.0) then
!                !  print*, 'skip line:', it, k, j, int1, int2, int3, var4, var5, var6
!                !endif
!              enddo
!            enddo
!          endif
!        enddo
!
!        do j = 1, ny_global+1
!          read(16002,*)
!        enddo
!        do k = 2, nz_global
!          do j = 1, ny_global
!            read(16002,*) int1, int2, int3, var4, var5, var6
!            !if (myid.eq.0) then
!            !  print *, 'ti=',ti,', read inlet debug:',&
!            !    int1, int2, int3, var4, var5, var6
!            !endif
!            velsec(j,k,1:3) = (/var6, var4, var5/) * IR_scale 
!          enddo
!          read(16002,*)
!        enddo
!        close(16002)
!      endif
!
!      vel_send = reshape(velsec, (/nz_global*ny_global*3/))
!      
!      !> plyunote: velsec is target velocity in global coordinate (ny,nz)
!      call mpi_bcast(vel_send, nz_global*ny_global*3, mpi_double_precision, 0, &
!        mpi_comm_2d_cart, ierr)
!      velsec = reshape(vel_send, (/ny_global, nz_global, 3/))
!        
!      !> plyunote: veltar(xsz(1), xsz(2), xsz(3)+2*level) is target velocity
!      if (isbot ) then
!        veltar(:,2:xsz(3)+level,:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!        xst(3)+1:xst(3)+xsz(3)-1+level, 1:3)
!      elseif(istop) then
!        veltar(:,1-level:xsz(3),:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!        xst(3)-level:xst(3)+xsz(3)-1, 1:3)
!      else
!        veltar(:,:,:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!        xst(3)-level:xst(3)+xsz(3)-1+level, 1:3)
!      endif
!    else if (InletRelaxation .eq. 7) then
!      !> plyunote: read from precursor simulation, but only read one slice and apply
!      !!   to a relaxation zone of single-grid width (actually it's the first plane.)
!      allocate(velsec(ny_global, nz_global, 3))
!      allocate(vel_send(ny_global*nz_global*3))
!      veltar(:,:,:)=0.0_wp; velsec(:,:,:)=0.0_wp; vel_send(:) = 0.0_wp
!      
!      write(groupname, *) (ti+IR_ti_shift-1)/1000+1
!      write(filen, '(a,i0.4,a,a,a,i0.10,a4)') '../inflowdata/case_', IR_dataindex, &
!        '/group_', trim(adjustl(groupname)), '/inlet_',ti+IR_ti_shift,'.dat'
!
!      if (myid .eq. 0) then
!        open(16001,FILE=filen,STATUS='OLD',ACTION='READ')
!        read(16001,*)
!        read(16001,*)
!        !> plyunote: the following iteration is wrong. x iteration should be included.
!        do k = 1, nz_global
!          do j = 1, ny_global
!            read(16001,*) var1, var2, var3, var4, var5, var6
!            velsec(j,k,1:3)=(/var4, var5, var6/) * IR_scale
!            !> plyunote: IR_scale should be around 1.0, if it is different
!            !from 1.0, based on Taylor frozen theory, the different advection
!            !speed will require you to do a interpolation of the inflow
!            !database, which is InletRelaxation mode 7.
!          enddo
!        enddo
!        close(16001)
!      endif
!      
!      vel_send = reshape(velsec, (/nz_global*ny_global*3/))
!      
!      !> plyunote: velsec is target velocity in global coordinate (ny,nz)
!      call mpi_bcast(vel_send, nz_global*ny_global*3, mpi_double_precision, 0, &
!        mpi_comm_2d_cart, ierr)
!      velsec = reshape(vel_send, (/ny_global, nz_global, 3/))
!
!      if (isbot ) then
!        veltar(:,2:xsz(3)+level,:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!          xst(3)+1:xst(3)+xsz(3)-1+level, 1:3)
!      elseif(istop) then
!        veltar(:,1-level:xsz(3),:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!          xst(3)-level:xst(3)+xsz(3)-1, 1:3)
!      else
!        veltar(:,:,:) = velsec(xst(2):xst(2)+xsz(2)-1, &
!          xst(3)-level:xst(3)+xsz(3)-1+level, 1:3)
!      endif
!    endif
!
!    !> Get du_relax(xsz(1),xsz(2),1-level:xsz(3)+level,3)
!    
!    if (xst(1)<=imax .and. (xst(1)+xsz(1)-1)>=imin) then
!      do k = 1-level, xsz(3)+level
!        z = zz(k) * hbar
!        if (isbot .and. (k==1 .or. k==0)) then
!
!        elseif(istop .and. k==xsz(3)+level) then
!
!        ! for wave cases, z might less than 0, hard to use log law.
!        elseif(z>0.0_wp) then
!          !print *, 'myid, k_local, k_global, dz, z=',myid,k,xst(3)+k-1,dz(k),z
!          if (InletRelaxation .eq. 7) then
!            i = 1
!            do j = 1, xsz(2)
!              u_relax(1) = (1.0_wp - rvec(i))*veltar(j,k,1)+rvec(i)*u(i,j,k)
!              u_relax(2) = (1.0_wp - rvec(i))*veltar(j,k,2)+rvec(i)*v(i,j,k)
!              u_relax(3) = (1.0_wp - rvec(i))*veltar(j,k,3)+rvec(i)*w(i,j,k)
!              u(i,j,k) = u_relax(1)
!              v(i,j,k) = u_relax(2)
!              w(i,j,k) = u_relax(3)
!            enddo
!          else
!            do j = 1, xsz(2)
!              do i = 1, xsz(1)
!                xi = (xst(1)+i-1)*(xl/nx)
!                if (xi>=IR_start .and. xi <= IR_end) then
!                  if (InletRelaxation .eq. 6) then
!                    u_relax(1) = (1.0_wp - rvec(i))*veltar_b(i-imin+1,j,k,1)+rvec(i)*u(i,j,k)
!                    u_relax(2) = (1.0_wp - rvec(i))*veltar_b(i-imin+1,j,k,2)+rvec(i)*v(i,j,k)
!                    u_relax(3) = (1.0_wp - rvec(i))*veltar_b(i-imin+1,j,k,3)+rvec(i)*w(i,j,k)
!                  else
!                    u_relax(1) = (1.0_wp - rvec(i))*veltar(j,k,1)+rvec(i)*u(i,j,k)
!                    u_relax(2) = (1.0_wp - rvec(i))*veltar(j,k,2)+rvec(i)*v(i,j,k)
!                    u_relax(3) = (1.0_wp - rvec(i))*veltar(j,k,3)+rvec(i)*w(i,j,k)
!                  endif
!
!                  u(i,j,k) = u_relax(1)
!                  v(i,j,k) = u_relax(2)
!                  w(i,j,k) = u_relax(3)
!                               
!                endif
!              enddo
!            enddo
!          endif
!        endif
!      enddo
!        
!    endif
!  else
!    print *, "plyudebug: This case is not yet included in Add_Inlet_Forcing"
!  endif
!    
!  deallocate(su, rvec)
!  deallocate(du_relax, du_temp)
!  deallocate(veltar)
!  if(InletRelaxation .eq. 3 .or. InletRelaxation .eq. 4 .or. InletRelaxation .eq.7)  then
!    deallocate(velsec, vel_send)
!  endif
!  if (InletRelaxation .eq. 6) then
!    deallocate(velsec_b, veltar_b, vel_send)
!  endif
!      
!  !call dealiasxy(u(:,:,1:xsz(3)))
!  !call dealiasxy(v(:,:,1:xsz(3)))
!  !call dealiasxy(w(:,:,1:xsz(3)))  
!
!  !if(istop)then 
!  !   call dealiasxy(u(:,:,xsz(3)+1))
!  !   call dealiasxy(v(:,:,xsz(3)+1))
!  !   call dealiasxy(w(:,:,xsz(3)+1))  
!  !end if 
!
!  ! if(myid==2)print*,'xsz(3)',size(u,3) 
!  call MPI_BARRIER(mpi_comm_world, ierr)
!  call update_ghost(u, level)
!  call update_ghost(v, level)
!  call update_ghost(w, level)
!  call MPI_BARRIER(mpi_comm_world, ierr)
!
!    !do j=1, xsz(2)
!    !  if (abs(f_relax(j))>f_max) then
!    !    f_max = abs(f_relax(j))
!    !    i_fmax = xst(1)+i-1; j_fmax=xst(2)+j-1; k_fmax=xst(3)+j-1
!    !    u_fmax = u(i,j,k); utgt_fmax=u_tgt(j); urelax_fmax=u_relax(j)
!    !  endif
!    !enddo
!    !
!
!    !print *, 'plyudebug, cpu ', myid, 'i, j, k=',i_fmax,j_fmax,k_fmax,&
!
!end subroutine Add_Velocity_Relaxation_old


!========================================================
   !> @brief Calculates the explicit nonlinear terms in
   !>  u momentum equation.
   !
   !> Calculates the explicit nonlinear terms FU and HU.
   !> \f[
   !>  FU^n = \frac{1}{2}(3HU^n - HU^{n-1})+\frac{1}{2 Re}\nabla^2 u^n,
   !> \f]
   !> where
   !> \f[
   !>   HU^n=-\zeta_t \frac{\partial u}{\partial \zeta}
   !>        -\frac{\partial (uu)}{\partial x}
   !>        -\frac{\partial (uv)}{\partial y}
   !>        -\frac{\partial (uw)}{\partial z}
   !>        + f.
   !> \f]
   !> \f$f\f$ is the \f$x\f$ component of the body foce.
   !! Here the body force is the linear forcing term and/or adverse 
   !! pressure gradient.
   !
   !> This subroutine applies to a mesh with constant grid spacing at the
   !> bottom and half-grid at the surface. Shear free boundary conditions
   !> are applied at the bottom.
   !
   !> See printed code on 1/12/2015.
   !
   !> @param[in]  time    time
   !> @param[in]  dt      time step
   !> @param[in]  Re      Reynolds number
   !> @param[in]  bforce  optional, linear forcing coefficient
   !> @param[in]  tf      optional, shear force at the surface
   !========================================================
!checked
!input: u, v, w
!input: ux, uy, uz_u

!input: tauwx, nut
!input: t11, t12, t13w
!input: wtforce (navier_init)
!input&output: hu
!output: fu, fv, fw
!output: t13w(:,:,1),t13w(:,:,xsz(3))
!output: t23w(:,:,1),t23w(:,:,xsz(3))
!output: t11(:,:,xsz(3)+1),t12(:,:,xsz(3)+1)

!internal global variables: tmp_x3, tmp_x4, tmp_x6
subroutine fun_u_nn_cb(fu,hu,time,tauwx)
  use grid, only : dz,  zetat, zetax
  use spectral
  use utils
  use param, only : pex, iturbine, timeturb, bforce, dt
  
  implicit none
  
  real(wp), intent(IN) :: time
  real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(INOUT) :: hu
  real(wp), dimension(:,:), intent(IN) :: tauwx
  real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: fu
  
  integer :: i,j,k  
  real(wp), allocatable, dimension(:,:) :: t13wx

  allocate(t13wx(xsz(1),xsz(2)))
  ! interpolate w at u,v points
  call update_ghost(w,level)
  do k =1, xsz(3)
     tmp_x3(:,:,k) = (w(:,:,k)+w(:,:,k-1))/2.0_wp
  enddo
  
  if(istop) tmp_x3(:,:,xsz(3))= w(:,:,xsz(3)-1)
  
  ! tmp_x4 is the convection term in u direction
  ! tmp_x4 = -u*ux-v*uy-tmp_x3*uz_u
  tmp_x4=0

  tmp_x7=u*ux
  call dealiasxy(tmp_x7(:,:,1:xsz(3)))
  tmp_x4=tmp_x4-tmp_x7

  tmp_x7=v*uy
  call dealiasxy(tmp_x7(:,:,1:xsz(3)))
  tmp_x4=tmp_x4-tmp_x7

  tmp_x7=tmp_x3*uz_u
  call dealiasxy(tmp_x7(:,:,1:xsz(3)))
  tmp_x4=tmp_x4-tmp_x7

  ! print*, u(1,1,1:xsz(3))*ux(1,1,1:xsz(3))

  tmp_x3=0  
  !apply wall model  
  if(isbot)then
     call pdfx(w(:,:,1), t13wx, pex)
     t13w(:,:,1)=tauwx(:,:)-nut(:,:,1)*(t13wx(:,:)+zetax(:,:,1)*(w(:,:,2)-w(:,:,1))/dz(2))
     call dealiasxy(t13w(:,:,1))
  endif
  
  ! print*, t13w(1,1,1)
  !apply top boundary condition
  if(istop)then
     t11(:,:,xsz(3)+1)=t11(:,:,xsz(3)-1)
     t12(:,:,xsz(3)+1)=t12(:,:,xsz(3)-1)
     t13w(:,:,xsz(3))=-t13w(:,:,xsz(3)-2)         
  endif
  
  call div_tau(t11, t12,  t13w, tmp_x6(:,:,1:xsz(3)))
      
  ! print*, 'tmp_x6', tmp_x6(1,1,1:xsz(3))

  tmp_x4 = tmp_x4 - tmp_x6
  
  ! compute du/dzeta
  call calc_uzeta(u, u_zeta(:,:,1:), level)
  call update_ghost(u_zeta,level)

  if(istop)then
     u_zeta(:,:,xsz(3)-1) = (u(:,:,xsz(3)+1) - u(:,:,xsz(3)-2))/(2*dz(xsz(3)-2))
     u_zeta(:,:,xsz(3)) = (u(:,:,xsz(3)+1) - u(:,:,xsz(3)-1))/(dz(xsz(3)-2))
  endif
  
  ! tmp_x4(:,:,1:xsz(3)) = tmp_x4(:,:,1:xsz(3))-zetat(:,:,1:xsz(3))*u_zeta(:,:,1:xsz(3))
  tmp_x7(:,:,1:xsz(3)) = zetat(:,:,1:xsz(3))*u_zeta(:,:,1:xsz(3))
  call dealiasxy(tmp_x7(:,:,1:xsz(3)))
  tmp_x4(:,:,1:xsz(3)) = tmp_x4(:,:,1:xsz(3))-tmp_x7(:,:,1:xsz(3))
  
  ! print*,'tmp_x4', tmp_x4(1,1,1:xsz(3))


  if (abs(time-dt) .le. 1e-6) hu= tmp_x4(:,:,1:xsz(3))
  
  if(iturbine.le.0.or.time.le.timeturb)then
     fu = (3.0_wp*tmp_x4(:,:,1:xsz(3))-hu)/2.0_wp+bforce
  else
     fu = (3.0_wp*tmp_x4(:,:,1:xsz(3))-hu)/2.0_wp+bforce + fturbinex &
          +wtforce*(1.0_wp-exp(-1.0_wp/(timeturb**2+1.0e-9_wp)*(time-timeturb)**2))
  endif

!  do i = 1, xsz(1)
!    do j = 1, xsz(2)
!      do k = 1, xsz(3)
!        if (abs(wtforce(i,j,k))>1e-7) then
!          print *, i, j, k, xst(1)+i-1, xst(2)+j-1, xst(3)+k-1, &
!            wtforce(i,j,k), fu(i,j,k)
!        endif
!      enddo
!    enddo
!  enddo
  
  hu= tmp_x4(:,:,1:xsz(3))

  ! print*, 'fu', fu(1,1,1:xsz(3))

  tmp_x3=0
  tmp_x4=0
  tmp_x6=0
  tmp_x7=0

  deallocate(t13wx)
end subroutine fun_u_nn_cb

   !========================================================
   !> @brief Calculates the explicit nonlinear terms in
   !>  v momentum equation.
   !
   !> Calculates the explicit nonlinear terms FV and HV.
   !> \f[
   !>  FV^n = \frac{1}{2}(3HV^n - HV^{n-1})+\frac{1}{2 Re}\nabla^2 v^n,
   !> \f]
   !> where
   !> \f[
   !>   HV^n=-\zeta_t \frac{\partial v}{\partial \zeta}
   !>        -\frac{\partial (vu)}{\partial x}
   !>        -\frac{\partial (vv)}{\partial y}
   !>        -\frac{\partial (vw)}{\partial z}
   !>        + f.
   !> \f]
   !> \f$f\f$ is the \f$y\f$ component of the body foce.
   !! Here the body force is the linear forcing term and/or adverse 
   !! pressure gradient.
   !
   !> This subroutine applies to a mesh with constant grid spacing at the
   !> bottom and half-grid at the surface. Shear free boundary conditions
   !> are applied at the bottom.
   !
   !> See printed code on 1/12/2015.
   !
   !> @param[in]  time    time
   !> @param[in]  dt      time step
   !> @param[in]  Re      Reynolds number
   !> @param[in]  bforce  optional, linear forcing coefficient
   !> @param[in]  tf      optional, shear force at the surface
   !========================================================
!checked
!used global tmp variables: tmp_x3, tmp_x4, tmp_x6
subroutine fun_v_nn_cb(fv,hv,time,tauwy)
  use grid, only : dz, zetat, zetay
  use spectral
  use utils
  use param, only: pex, dt
  implicit none
  
  real(wp), intent(IN) :: time
  real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: fv
  real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(INOUT) :: hv
  real(wp), dimension(:,:), intent(IN) :: tauwy
  integer :: k
  
  real(wp), allocatable, dimension(:,:) :: t23wy
  
  allocate(t23wy(xsz(1),xsz(2)))
  ! interpolate w at u,v points
  call update_ghost(w,level)
  do k =1, xsz(3)
     tmp_x3(:,:,k) = (w(:,:,k)+w(:,:,k-1))/2.0_wp
  enddo
  
  if(istop) tmp_x3(:,:,xsz(3))= w(:,:,xsz(3)-1)  
  
  !tmp_x4 is the convection term in u direction
  ! tmp_x4 = -u*vx-v*vy-tmp_x3*vz_u
  tmp_x4=0

  tmp_x7=u*vx
  call dealiasxy(tmp_x7(:,:,1:xsz(3)))
  tmp_x4=tmp_x4-tmp_x7

  tmp_x7=v*vy
  call dealiasxy(tmp_x7(:,:,1:xsz(3)))
  tmp_x4=tmp_x4-tmp_x7

  tmp_x7=tmp_x3*vz_u
  call dealiasxy(tmp_x7(:,:,1:xsz(3)))
  tmp_x4=tmp_x4-tmp_x7
  !apply wall model
    
  if(isbot)then
     call pdfy_x(w(:,:,1), t23wy, pex)
     t23w(:,:,1)=tauwy-nut(:,:,1)*(t23wy+zetay(:,:,1)*(w(:,:,2)-w(:,:,1))/dz(2))
     call dealiasxy(t23w(:,:,1))
  endif
  
  !apply top boundary condition  
  if(istop)then
     t12(:,:,xsz(3)+1)=t12(:,:,xsz(3)-1)
     t22(:,:,xsz(3)+1)=t22(:,:,xsz(3)-1)
     t23w(:,:,xsz(3))=-t23w(:,:,xsz(3)-2)         
  endif
  
  call div_tau(t12, t22, t23w, tmp_x6(:,:,1:xsz(3)))

  ! why do not add wall model stress in fun_v
  if(isbot)tmp_x6(:,:,2)=0
  tmp_x4 = tmp_x4 - tmp_x6
   
  ! compute dv/dzeta
  call calc_uzeta(v, u_zeta(:,:,1:), level)
  call update_ghost(u_zeta,level)

  if(istop)then
     u_zeta(:,:,xsz(3)-1) = (v(:,:,xsz(3)+1) - v(:,:,xsz(3)-2))/(2*dz(xsz(3)-2))
  endif

  ! tmp_x4(:,:,1:xsz(3)) = tmp_x4(:,:,1:xsz(3))-zetat(:,:,1:xsz(3))*u_zeta(:,:,1:xsz(3))
  tmp_x7(:,:,1:xsz(3)) = zetat(:,:,1:xsz(3))*u_zeta(:,:,1:xsz(3))
  call dealiasxy(tmp_x7(:,:,1:xsz(3)))
  tmp_x4(:,:,1:xsz(3)) = tmp_x4(:,:,1:xsz(3))-tmp_x7(:,:,1:xsz(3))

  if (abs(time-dt).le. 1e-6) hv= tmp_x4(:,:,1:xsz(3))

  !> modified by plyu
  !fv = (3.0_wp*tmp_x4(:,:,1:xsz(3))-hv)/2.0_wp
  if(iturbine.le.0.or.time.le.timeturb)then
     fv = (3.0_wp*tmp_x4(:,:,1:xsz(3))-hv)/2.0_wp
  else
     fv = (3.0_wp*tmp_x4(:,:,1:xsz(3))-hv)/2.0_wp + fturbiney &
          +wtforce_y*(1.0_wp-exp(-1.0_wp/(timeturb**2+1.0e-9_wp)*(time-timeturb)**2))
  endif
 
  hv = tmp_x4(:,:,1:xsz(3))

  tmp_x3=0
  tmp_x4=0
  tmp_x6=0
  tmp_x7=0

  deallocate(t23wy)
  ! print*, 'fv', fv(1,1,1:xsz(3))
end subroutine fun_v_nn_cb

   !========================================================
   !> @brief Calculates the explicit nonlinear terms in
   !>  w momentum equation.
   !
   !> Calculates the explicit nonlinear terms FW and HW.
   !> \f[
   !>  FW^n = \frac{1}{2}(3HW^n - HW^{n-1})+\frac{1}{2 Re}\nabla^2 w^n,
   !> \f]
   !> where
   !> \f[
   !>   HW^n=-\zeta_t \frac{\partial w}{\partial \zeta}
   !>        -\frac{\partial (wu)}{\partial x}
   !>        -\frac{\partial (wv)}{\partial y}
   !>        -\frac{\partial (ww)}{\partial z}
   !>        + f_z.
   !> \f]
   !> \f$f_z\f$ is the \f$z\f$ component of the body foce. Here the body force is
   !> the linear forcing term.
   !
   !> This subroutine applies to a mesh with constant grid spacing at the
   !> bottom and half-grid at the surface. Shear free boundary conditions
   !> are applied at the bottom.
   !
   !> See printed code on 1/12/2015.
   !
   !> @param[in]  time    time
   !> @param[in]  dt      time step
   !> @param[in]  Re      Reynolds number
   !> @param[in]  bforce  optional, linear forcing coefficient
!========================================================

!checked
!used global tmp variables: tmp_x6, tmp_x5, tmp_x4, tmp_x2, tmp_x3, 
subroutine fun_w_nn_cb(fw,hw,time)
      use grid, only :   zetat_w
      use spectral
      use utils
      use param, only: dt

      implicit none

      real(wp), intent(IN) :: time
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(INOUT) :: hw
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: fw

      integer :: k
      
      !apply top boundary condition

      if(istop)then
         t13w(:,:,xsz(3))=-t13w(:,:,xsz(3)-2)
         t23w(:,:,xsz(3))=-t23w(:,:,xsz(3)-2)         
      endif
      
      call div_tau_w(t13, t23, t33, t13w, t23w, tmp_x6(:,:,1:xsz(3)))
      
      tmp_x4 = tmp_x4 - tmp_x6

      call update_ghost(u,level)
      call update_ghost(v,level)
      
      ! interpolate u,v at w points
      do k =1, xsz(3)
         tmp_x2(:,:,k) = (u(:,:,k)+u(:,:,k+1))/2.0_wp
         tmp_x3(:,:,k) = (v(:,:,k)+v(:,:,k+1))/2.0_wp
      enddo

      !tmp_x4 is the convection term in u direction
      ! tmp_x4 = tmp_x4-tmp_x2*wx-tmp_x3*wy-w*wz_w
      tmp_x7=tmp_x2*wx
      call dealiasxy(tmp_x7(:,:,1:xsz(3)))
      tmp_x4=tmp_x4-tmp_x7

      tmp_x7=tmp_x3*wy
      call dealiasxy(tmp_x7(:,:,1:xsz(3)))
      tmp_x4=tmp_x4-tmp_x7

      tmp_x7=w*wz_w
      call dealiasxy(tmp_x7(:,:,1:xsz(3)))
      tmp_x4=tmp_x4-tmp_x7

      ! compute dw/dzeta
      call calc_wzeta(w, u_zeta(:,:,1:), level)
      call update_ghost(u_zeta,level)

      ! tmp_x4(:,:,1:xsz(3)) = tmp_x4(:,:,1:xsz(3))-zetat_w(:,:,1:xsz(3))*u_zeta(:,:,1:xsz(3))

      tmp_x7(:,:,1:xsz(3)) = zetat_w(:,:,1:xsz(3))*u_zeta(:,:,1:xsz(3))
      call dealiasxy(tmp_x7(:,:,1:xsz(3)))
      tmp_x4(:,:,1:xsz(3)) = tmp_x4(:,:,1:xsz(3))-tmp_x7(:,:,1:xsz(3))

      if (abs(time-dt).le.1e-6) hw= tmp_x4(:,:,1:xsz(3))

      !> modified by plyu
      !fw = (3.0_wp*tmp_x4(:,:,1:xsz(3))-hw)/2.0_wp
      if(iturbine.le.0.or.time.le.timeturb)then
        fw = (3.0_wp*tmp_x4(:,:,1:xsz(3))-hw)/2.0_wp
      else
        fw = (3.0_wp*tmp_x4(:,:,1:xsz(3))-hw)/2.0_wp &
          +wtforce_z*(1.0_wp-exp(-1.0_wp/(timeturb**2+1.0e-9_wp)*(time-timeturb)**2))
      endif

      hw = tmp_x4(:,:,1:xsz(3))

      ! print*, 'fw',fw(1,1,1:xsz(3))

      tmp_x2=0
      tmp_x3=0
      tmp_x4=0
      tmp_x5=0
      tmp_x6=0
      tmp_x7=0
      
   end subroutine fun_w_nn_cb

   !checked
   !input: fu, fv, fw, dt
   !intput & output: u, v, w
   
   subroutine velocity_prediction(fu,fv,fw)
     use spectral
     use param, only: dt
     use mpi
     implicit none
     real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(IN) :: fu,fv,fw
     integer :: ierr
      integer :: i,j,k
     
      ! if(myid==2)print*,'xsz(3)',size(u,3) 
      if(isbot .and. istop)then 
         u(:,:,2:xsz(3)-1) = u(:,:,2:xsz(3)-1) + dt*fu(:,:,2:xsz(3)-1)
         v(:,:,2:xsz(3)-1) = v(:,:,2:xsz(3)-1) + dt*fv(:,:,2:xsz(3)-1)
         w(:,:,2:xsz(3)-2) = w(:,:,2:xsz(3)-2) + dt*fw(:,:,2:xsz(3)-2)
      else if(isbot)then
         u(:,:,2:xsz(3)) = u(:,:,2:xsz(3)) + dt*fu(:,:,2:xsz(3))
         v(:,:,2:xsz(3)) = v(:,:,2:xsz(3)) + dt*fv(:,:,2:xsz(3))
         w(:,:,2:xsz(3)) = w(:,:,2:xsz(3)) + dt*fw(:,:,2:xsz(3))
      else if(istop)then
         u(:,:,1:xsz(3)-1) = u(:,:,1:xsz(3)-1) + dt*fu(:,:,1:xsz(3)-1)
         v(:,:,1:xsz(3)-1) = v(:,:,1:xsz(3)-1) + dt*fv(:,:,1:xsz(3)-1)
         w(:,:,1:xsz(3)-2) = w(:,:,1:xsz(3)-2) + dt*fw(:,:,1:xsz(3)-2)
      else
         u(:,:,1:xsz(3)) = u(:,:,1:xsz(3)) + dt*fu(:,:,1:xsz(3))
         v(:,:,1:xsz(3)) = v(:,:,1:xsz(3)) + dt*fv(:,:,1:xsz(3))
         w(:,:,1:xsz(3)) = w(:,:,1:xsz(3)) + dt*fw(:,:,1:xsz(3))        
      endif

      ! print*, u(1,1,1:xsz(3)+1)

      call dealiasxy(u(:,:,1:xsz(3)))
      call dealiasxy(v(:,:,1:xsz(3)))
      call dealiasxy(w(:,:,1:xsz(3)))  

      if(istop)then 
         call dealiasxy(u(:,:,xsz(3)+1))
         call dealiasxy(v(:,:,xsz(3)+1))
         call dealiasxy(w(:,:,xsz(3)+1))  
      end if 

     ! do i = 1, xsz(1)
     !   do j = 1, xsz(2)
     !     do k = 1, xsz(3)
     !       if (abs(wtforce(i,j,k))>1e-7) then
     !         print *,'a,', i, j, k, xst(1)+i-1, xst(2)+j-1, xst(3)+k-1, &
     !           wtforce(i,j,k), fu(i,j,k), u(i,j,k)
     !       endif
     !     enddo
     !   enddo
     ! enddo


      ! if(myid==2)print*,'xsz(3)',size(u,3) 
      call MPI_BARRIER(mpi_comm_world, ierr)
      call update_ghost(u, level)
      call update_ghost(v, level)
      call update_ghost(w, level)
      call MPI_BARRIER(mpi_comm_world, ierr)
      ! print*, 'end'
      
      !do i = 1, xsz(1)
      !  do j = 1, xsz(2)
      !    do k = 1, xsz(3)
      !      if (abs(wtforce(i,j,k))>1e-7) then
      !        print *,'b,', i, j, k, xst(1)+i-1, xst(2)+j-1, xst(3)+k-1, &
      !          wtforce(i,j,k), fu(i,j,k), u(i,j,k)
      !      endif
      !    enddo
      !  enddo
      !enddo

      ! if(myid==0) print*, u(1,1,:)
      ! if(myid==1) print*, u(1,1,:)
 endsubroutine velocity_prediction

 !checked
 !input: hx, hy, eta, hh, zz, 
 !input: uf, vf, wf
 !input: z0, hbar
 !output: tauwx, tauwy
 subroutine wall_model_v3
   use decomp
   use grid, only : hx, hy, eta, hh, zz
   !use navier, only : uf, vf, wf, tauwx,tauwy
   use param, only: z0, hbar
   use spectral
   implicit none

   real(wp) :: utal, utals, vtal, vtals
   real(wp) :: t1x,t1z,t2y,t2z
   real(wp) :: z, vel, tt, kappa
   integer :: i,j

   DATA KAPPA/0.4_wp/

   ! call dealiasxy(uf(:,:,1:xsz(3)))
   ! if(myid==0) print*, uf(:,1,1)

      do j =1, xsz(2)
         do i =1, xsz(1)
            t1x=1.0_wp/sqrt(1+hx(i,j)**2)
            t1z=-hx(i,j)/sqrt(1+hx(i,j)**2)
            t2y=1.0_wp/sqrt(1+hy(i,j)**2)
            t2z=-hy(i,j)/sqrt(1+hy(i,j)**2)
            utal=uf(i,j,2)*t1x+wf(i,j,2)*t1z
            utals=uf(i,j,1)*t1x+wf(i,j,1)*t1z
            vtal=vf(i,j,2)*t2y+wf(i,j,2)*t2z
            vtals=vf(i,j,1)*t2y+wf(i,j,1)*t2z
            z=zz(2)*(hbar+eta(i,j))-hh(i,j)
            vel=sqrt((utal-utals)**2+(vtal-vtals)**2)
! cc            vel=sqrt(utal**2+vtal**2)
            tt=-(kappa/log((z+hh(i,j))/z0))**2*vel
            tauwx(i,j)=tt*(utal-utals)
            tauwy(i,j)=tt*(vtal-vtals)
! cc            tauwx(i,j)=tt*utal
! cc            tauwy(i,j)=tt*vtal

        
         ! if(i==1 .and. j==1.and.myid==0)then
            ! print*, "wall model",log((z+hh(i,j))/z0)
            ! print*, "wall model",log((z+hh(i,j)))
            ! print*, "wall model",z+hh(i,j)
            ! print*, "wall model",tauwx(1,1)
            ! print*, "wall model",tauwy(1,1)
         ! end if 

         enddo
      enddo

      ! if(myid==0)print*,'1', tauwx(:,4)
      do j=1, xsz(2)
         do i=1, xsz(1)
            ! read(103,*)tauwx(i,j)
         end do 
      end do 
      ! print*, 'sum',sum(tauwx(1:xsz(1),4))/xsz(1)
      ! if(myid==0)print*,'1', tauwx(:,4)
   call dealiasxy(tauwx)
   ! stop
   call dealiasxy(tauwy)
      ! if(myid==0)print*,'1', tauwx(:,4)

 endsubroutine wall_model_v3

 !> added by plyu, for plane channel inflow generation
 subroutine wall_model_planar_average
   use decomp
   use grid, only : hx, hy, eta, hh, zz
   !use navier, only : uf, vf, wf, tauwx,tauwy
   use param, only: z0, hbar
   use spectral
   use mpi
   implicit none

   real(wp) :: utal, utals, vtal, vtals
   real(wp) :: t1x,t1z,t2y,t2z
   real(wp) :: z, vel, tt, kappa
   integer :: i,j

   real(wp) :: velsum_l, velsum_g, areasum_l, areasum_g
   real(wp) :: dsend(2), drecv(2), velm_p
   integer :: ierror

   DATA KAPPA/0.4_wp/

   ! call dealiasxy(uf(:,:,1:xsz(3)))
   ! if(myid==0) print*, uf(:,1,1)

      velsum_l = 0.0; areasum_l = 0.0; areasum_l = 0.0; areasum_g = 0.0
      do j =1, xsz(2)
         do i =1, xsz(1)
            t1x=1.0_wp/sqrt(1+hx(i,j)**2)
            t1z=-hx(i,j)/sqrt(1+hx(i,j)**2)
            t2y=1.0_wp/sqrt(1+hy(i,j)**2)
            t2z=-hy(i,j)/sqrt(1+hy(i,j)**2)
            utal=uf(i,j,2)*t1x+wf(i,j,2)*t1z
            utals=uf(i,j,1)*t1x+wf(i,j,1)*t1z
            vtal=vf(i,j,2)*t2y+wf(i,j,2)*t2z
            vtals=vf(i,j,1)*t2y+wf(i,j,1)*t2z
            z=zz(2)*(hbar+eta(i,j))-hh(i,j)
            vel=sqrt((utal-utals)**2+(vtal-vtals)**2)
! cc            vel=sqrt(utal**2+vtal**2)
            velsum_l = velsum_l + vel*sqrt(hx(i,j)**2+hy(i,j)**2+1)
            areasum_l = areasum_l + sqrt(hx(i,j)**2+hy(i,j)**2+1)
        enddo
      enddo

      dsend(1) = velsum_l; dsend(2) = areasum_l; drecv = 0.0
      call MPI_Allreduce(dsend, drecv, 2, mpi_double_precision, &
            MPI_SUM, MPI_COMM_2D_COL, ierror)
      velsum_g = drecv(1); areasum_g = drecv(2)
      velm_p = velsum_g / areasum_g
      if(myid.eq.0) print *, "velm_bottom=", velm_p

      do j =1, xsz(2)
         do i =1, xsz(1)
            tt=-(kappa/log((z+hh(i,j))/z0))**2*velm_p
            tauwx(i,j)=tt*(utal-utals)
            tauwy(i,j)=tt*(vtal-vtals)
! cc            tauwx(i,j)=tt*utal
! cc            tauwy(i,j)=tt*vtal

        
         ! if(i==1 .and. j==1.and.myid==0)then
            ! print*, "wall model",log((z+hh(i,j))/z0)
            ! print*, "wall model",log((z+hh(i,j)))
            ! print*, "wall model",z+hh(i,j)
            ! print*, "wall model",tauwx(1,1)
            ! print*, "wall model",tauwy(1,1)
         ! end if 

         enddo
      enddo

      ! if(myid==0)print*,'1', tauwx(:,4)
      do j=1, xsz(2)
         do i=1, xsz(1)
            ! read(103,*)tauwx(i,j)
         end do 
      end do 
      ! print*, 'sum',sum(tauwx(1:xsz(1),4))/xsz(1)
      ! if(myid==0)print*,'1', tauwx(:,4)
   call dealiasxy(tauwx)
   ! stop
   call dealiasxy(tauwy)
      ! if(myid==0)print*,'1', tauwx(:,4)

 end subroutine wall_model_planar_average
