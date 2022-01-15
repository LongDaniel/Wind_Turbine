module post_proc
  
  use hos_param, only : wp, nxhos, nyhos, ncpu_hos, pex_hos, pey_hos
  use fft_hos

  public :: post_les_hos, get_1dspectrum
  public :: post_les_hos_wave_coherent_analysis, post_les_hos_k_omg_spectrum

  !> added by plyu
  public :: post_les_hos_turbine, post_tke_budget, post_vortex_dynamics

contains

  
  subroutine post_les_hos
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
    
    use hdf_io
    !end
    
    implicit none
    
    integer :: ioutc, ioutd
    real(wp) :: time
    
    integer :: ierror, ierr
    
    double precision :: t1, t2
    
    !HOS variables
    integer icon
    
    integer(INT64) :: id
    integer :: i, j, k, it, root0
    
    real(wp), allocatable, dimension(:) :: fpk, betak, ampfd, epskx
    
    real(wp) beta_miles
    
    character (len=100) :: fileid, fgrid, fparam_les, fdata_les, faux, fparam_hos, fdata_hos

    real(wp) c_phase, u_lambda, u_half_lambda, wage1, wage2, wage3
    real(wp) clen, cspeed, ctime, x, y, fpt

    real(wp), allocatable, dimension(:,:) :: etaall, hhall, tmp
    real(wp) ftn, rek, dwk, disspm, l_kom_ave

    real(wp), allocatable, dimension(:,:,:) :: uupall, wwpall, velall, ppall, etasig, pp_wave
    real(wp), allocatable, dimension(:,:,:,:) :: ppsig
    real(wp), allocatable, dimension(:,:) :: dtheta, uxy, vxy, wxy
    real(wp), allocatable, dimension(:,:) :: ek11, ek22, ek33
    real(wp), allocatable, dimension(:) :: dissp, ekm, l_kom

    real(wp), allocatable, dimension(:) :: q1, q2, q3, q4, um, vm, wm, uup, vvp, wwp, zzall, uw_shear
   
    integer ini_it, nk
    
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

    allocate(fpk(nx_global))
    allocate(betak(nx_global))
    allocate(ampfd(nx_global))
    allocate(epskx(nx_global))
    allocate(dissp(nz_global))
    allocate(l_kom(nz_global))

    allocate(etaall(nx_global,ny_global))
    allocate(hhall(nx_global,ny_global))
    allocate(tmp(nx_global,ny_global))

    allocate(uupall(nx_global,ny_global,nz_global))
    allocate(wwpall(nx_global,ny_global,nz_global))
    allocate(velall(nx_global,ny_global,nz_global))
    allocate(ppall(nx_global,ny_global,nz_global))

    allocate(uxy(nx_global,ny_global))
    allocate(vxy(nx_global,ny_global))
    allocate(wxy(nx_global,ny_global))
    
    nk = nx_global / 3
    allocate(ek11(nk,nz_global))
    allocate(ek22(nk,nz_global))
    allocate(ek33(nk,nz_global))
    allocate(ekm(nx_global))

    allocate(q1(nz_global))
    allocate(q2(nz_global))
    allocate(q3(nz_global))
    allocate(q4(nz_global))
    allocate(um(nz_global), uup(nz_global))
    allocate(vm(nz_global), vvp(nz_global))
    allocate(wm(nz_global), wwp(nz_global))
    allocate(uw_shear(nz_global))
    allocate(zzall(nz_global))

    read(14,*) clen, cspeed, ctime

    if (myid == 0) then
       open(66,file="wind_input.dat")
       open(68,file="turb_spectrum.dat")
       open(70,file="mean_flux.dat")
       open(72,file="q_analysis.dat")
    end if

    ini_it = 600
    do it = ini_it + 1, ini_it + ntime
       write(fileid,*) it
       fileid = trim(adjustl(fileid))
       
       fparam_les = trim("restart_param.dat"//fileid)
       fdata_les = trim("restart.h5"//fileid)
       faux =  trim("restart_aux.h5"//fileid)
       fgrid = trim("grid.h5"//fileid)
       
       fparam_hos = trim("restart_param_hos.dat"//fileid)
       fdata_hos = trim("restart_hos.h5"//fileid)
       
       call reread(time, ioutd, ioutc, fgrid, fparam_les, fdata_les, faux)
       call grid_gen
       
       !need to add hos initialization subroutines below
       if(myid==0) print *, 'HOS initialization started!'
       
       if (it == ini_it + 1) then
          call hos_init
       end if
       call read_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos, fparam_hos, fdata_hos)
       
       call wavenum(wvn_hos)
       call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
       call zeta(eta_hos,zp_hos)
       tmp=1.0_wp
       call boundvp(vps_hos,r_hos,zp_hos)
       call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
       call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
       
       call bottom_hos_les(time)
       
       if(myid==0) print *, 'HOS initialized!'
       
       ! transfer data to the upper cpus
       call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       
       eta=hh+0
       
       call init_p_poisson
       call nl_coef
       
       !test starts here!
       if (isbot) then
          tmp = 0
          tmp(xst(1):xend(1), xst(2):xend(2)) = eta(:,:)
          call MPI_Reduce(tmp, etaall, nx_global*ny_global, &
             mpi_double_precision, MPI_SUM, 0, MPI_COMM_2D_COL, ierror)
          epskx = 0
          call fft_for_x_hos(etaall,tmp,nx_global,ny_global)
          do i = 3, nx_global, 2
             do j = 1, ny_global
                rek = sqrt(((i-1.0)/2*pex)**2+((j-1.0)/2*pey)**2)
                ftn = tmp(i,j)**2 + tmp(i+1,j)**2
                epskx((i-1)/2) = epskx((i-1)/2) + 2*ftn*rek**2/ny_global
             end do
          end do
          do i = 1, nx_global
             epskx(i) = sqrt(2*epskx(i))
          end do
       end if
       
       call form_drag_k(fpk, betak, ampfd, fpt)
       
       root0 = 0
       
       call gather_2d_xy(hh(1:xsz(1),1:xsz(2)),etaall,root0)
       
       !     call turb_analysis
       if (myid == 0) then
          write(66,*) " VARIABLES = wvn,c,wage1, wage2, wage3,betak, beta_plg, beta_psm,beta_m" &
               // ",gamma, gam_don99, gam_don06, ka"
          write(66,"(A,E16.4,A,I5)") ' ZONE T="',time,'" I=', nx_global/2/3*2-1
          write(66,"(A,E16.4,A)") ' AUXDATA time="',time,'"'
          write(66,"(A,E16.4,A)") ' AUXDATA clen="',clen,'"'
          write(66,"(A,E16.4,A)") ' AUXDATA ctime="',ctime,'"'
          write(66,"(A,E16.4,A)") ' AUXDATA form_drag="',fpt,'"'
          
          !        print *,"usbot=",usbot,"z0=",z0, "fr2=", fr2
          do j = 1, nx_global/2/3*2-1
             call get_beta_miles(j*pex,beta_miles)
             c_phase = sqrt(1/fr2/j/pex)
             wage1 = c_phase / usbot
             u_lambda = 2.5*usbot*log(twopi/20/pex/z0)
             u_half_lambda = 2.5*usbot*log(pi/20/pex/z0)
             
             wage2 = c_phase / u_lambda
             wage3 = c_phase / u_half_lambda
             write(66,'(25e12.4)') j * pex, c_phase, wage1, wage2, wage3, betak(j), 16.0, 48.0,beta_miles &
                  , betak(j) / wage1**2, 0.28*(1/wage3-1)*abs(1/wage3-1), 0.17*(1/wage3-1)*abs(1/wage3-1),&
                  epskx(j)
          end do
!          write(67,'(25e12.4)') it*1.0, fpt
          print *,"finish step ",it
       end if
       
       !test ends here!
       
       call gather_3d_xyz(u(1:xsz(1),1:xsz(2),1:xsz(3)),velall,root0)
       do k = 1, nz_global
          call get_1dspectrum(velall(:,:,k),ek11(:,k),dwk,nk)
       end do
       uupall = velall

       call gather_3d_xyz(v(1:xsz(1),1:xsz(2),1:xsz(3)),velall,root0)
       do k = 1, nz_global
          call get_1dspectrum(velall(:,:,k),ek22(:,k),dwk,nk)
       end do

       call gather_3d_xyz(w(1:xsz(1),1:xsz(2),1:xsz(3)),velall,root0)
       do k = 1, nz_global
          call get_1dspectrum(velall(:,:,k),ek33(:,k),dwk,nk)
       end do

       call get_mean_flux(um,vm,wm,uup,vvp,wwp,zzall,uupall,wwpall)
       call get_quadrant(q1,q2,q3,q4,uupall,wwpall)

       do k = 1, nz_global
          uw_shear(k) = sum(uupall(:,:,k)*wwpall(:,:,k)) / nx_global / ny_global
       end do

       if (myid == 0) then

          write(70,*) " VARIABLES = z,um,vm,wm,uf,vf,wf,q1,q2,q3,q4,uw"
          write(70,"(A,E16.4,A,I5)") ' ZONE T="',time,'" I=', nz_global
          write(70,"(A,E16.4,A)") ' AUXDATA time="',time,'"'
          write(70,"(A,E16.4,A)") ' AUXDATA ustar="',usbot,'"'
          write(70,"(A,E16.4,A)") ' AUXDATA hbar="',hbar,'"'
          write(70,"(A,E16.4,A)") ' AUXDATA z0="',z0,'"'

          do k = 1, nz_global
             write(70,'(25e12.4)')  zzall(k)*hbar, um(k), vm(k), wm(k), &
               uup(k), vvp(k), wwp(k), q1(k), q2(k), q3(k), q4(k), uw_shear(k)
          end do

          write(72,*) " VARIABLES = uf_nw,wf_nw,uf_out,wf_out"
          write(72,"(A,E16.4,A,I5)") ' ZONE T="',time,'" I=', nx_global * ny_global
          write(72,"(A,E16.4,A)") ' AUXDATA time="',time,'"'
          write(72,"(A,E16.4,A)") ' AUXDATA ustar="',usbot,'"'

          k = 80
          do j = 1, ny_global
             do i = 1, nx_global
                write(72,'(25e12.4)') uupall(i,j,5), wwpall(i,j,5),uupall(i,j,k), wwpall(i,j,k)
             end do
          end do
             

          do k = 1, nz_global / 2
             dissp(k) = 0
             do j = 1, nk
                dissp(k) = dissp(k) + (2*usbot/resbot)*(j * dwk)**2 &
                     * (ek11(j,k)+ek22(j,k)+ek33(j,k)) * dwk
             end do

             l_kom(k) = ((usbot/resbot)**3/dissp(k))**0.25
          end do

          ekm = 0
          do j = 1, nk
             do k = 1, nz_global / 2
                ekm(j) = ekm(j) + (ek11(j,k) + ek22(j,k) + ek33(j,k)) / (nz_global / 2)
             end do
          end do
          disspm = 0
          do j = 1, nk
             disspm = disspm + (2*usbot/resbot)*(j * dwk)**2 &
                  * ekm(j) * dwk
          end do
          l_kom_ave = ((usbot/resbot)**3/disspm)**0.25

          write(68,*) " VARIABLES = k,ek11_nw,ek22_nw,ek33_nw,ek11_out,ek22_out,ek33_out, ekm"
          write(68,"(A,E16.4,A,I5)") ' ZONE T="',time,'" I=', nk
          write(68,"(A,E16.4,A)") ' AUXDATA time="',time,'"'
          write(68,"(A,E16.4,A)") ' AUXDATA disp_nw="',dissp(5),'"'
          write(68,"(A,E16.4,A)") ' AUXDATA l_kom_nw="',l_kom(5),'"'
          write(68,"(A,E16.4,A)") ' AUXDATA disp_out="',dissp(80),'"'
          write(68,"(A,E16.4,A)") ' AUXDATA l_kom_out="',l_kom(80),'"'
          write(68,"(A,E16.4,A)") ' AUXDATA disp="',disspm,'"'
          write(68,"(A,E16.4,A)") ' AUXDATA l_kom="',l_kom_ave,'"'
          write(68,"(A,E16.4,A)") ' AUXDATA nu="',usbot/resbot,'"'

          k = 80
          do j = 1, nk
             write(68,'(25e12.4)') j * dwk, ek11(j,5), ek22(j,5), ek33(j,5), ek11(j,k), &
                  ek22(j,k), ek33(j,k), ekm(j)
          end do
       end if
       
       ! kfilt = 0
    end do

    if (myid == 0) then
       close(66)
       close(68)
    end if

    
    deallocate(fpk,betak,ampfd, epskx)
    
    deallocate(etaall,hhall,tmp)
    
    deallocate(uupall,wwpall,velall)

    deallocate(uxy, vxy, wxy)

    deallocate(ek11, ek22, ek33, dissp, ekm, l_kom)

    deallocate(q1,q2,q3,q4)

    deallocate(um,vm,wm,uup,vvp,wwp)

    call fft_finalize
    call decomp_finalize
    !HOS
    call fft_finalize_hos
    
    call mpi_finalize(ierror)
  end subroutine post_les_hos
  
  subroutine comp2cart (uin, uout, dir)
    use param
    use decomp
    use grid, only : eta, hh, zz, zw, hx, hy, her, level
    use mpi
    use lib_array
    implicit none
    real(wp), dimension(:,:,1-level:), intent(in) :: uin
    real(wp), dimension(:,:,1-level:), intent(out) :: uout
    integer :: dir ! dir should be 1, 2, or 3

    real(wp), dimension(:), allocatable :: ucolin, ucolout
    integer :: i, j, k
    real(wp), dimension(:), allocatable :: cartzcol, cartzwcol, zzi, zwi
    real(wp), dimension(:), allocatable :: zin, zout
    integer :: iprint

    iprint = 0
    allocate(ucolin(1-level:xsz(3)+level))
    allocate(ucolout(1-level:xsz(3)+level))
    allocate(zin(1-level:xsz(3)+level))
    allocate(zout(1-level:xsz(3)+level))
    ucolin = 0.0; ucolout = 0.0; zin = 0.0; zout = 0.0
    if (dir.ne.3) then
      allocate(cartzcol(1-level:xsz(3)+level))
      allocate(zzi(1-level:xsz(3)+level))
      cartzcol = 0.0; zzi = zz * hbar; zout = zzi
    else
      allocate(cartzwcol(1-level:xsz(3)+level))
      allocate(zwi(1-level:xsz(3)+level))
      cartzwcol = 0.0; zwi = zw * hbar; zout = zwi
    endif
    
    do i = 1, xsz(1)
      do j = 1, xsz(2)
        ucolin = uin(i,j,1-level:xsz(3)+level)
        
        !! interpolation. Fill with -1000.0 if out of range
        do k = 1-level, xsz(3)+level
          if (dir .ne. 3) then
            cartzcol(k) = zz(k) * (hbar+hh(i,j)) - hh(i,j)
            zin = cartzcol
          else
            cartzwcol(k) = zw(k) * (hbar+hh(i,j)) - hh(i,j)
            zin = cartzwcol
          endif
        enddo
        ucolout = interp1d(zin, ucolin, zout, .false., -1000.0_wp)

        do k = 1, xsz(3)
          if (zout(k)<(-hh(i,j))) then
            !> this should happen only in bottom subdomains.
            !> use the wave surface velocity directly. u,v,w overlapped here.
            ucolout(k) = ucolin(1)
          else if (zout(k)<zin(2)) then
            !> it need log law interpolation here. let's do it later.
            exit  !> replace it later.
          else if (zout(k)>=zin(2)) then
            exit
          endif
        enddo

        do k = 1, xsz(3)
          if (ucolout(k)<-10.0 .and. dir.eq.1) then
            iprint = 1
          endif
        enddo
        
        if (iprint > 0) then
          print *, 'myid,',myid,'i,j,',i,j,',xst(1:3)=',xst(1:3)
          print *, 'zin,', zin
          print *, 'ucolin,', ucolin
          print *, 'zout,',zout 
          print *, 'ucolout,',ucolout
        endif
        iprint = 0

        !> plyunote: we hope zout(1:xsz(3)) in range of zin(1-level:xsz(3)+level)
        !> for small wave amplitude case, this might be true
        !> for large wave amplitude case, or if dz is very small, here need
        !  interpolation in global zzall.

       ! if (i.eq.1 .and. j.eq.1) then
       !   print *, 'myid,', myid, ' zin,', zin
       !   print *, 'myid,', myid, ' uin,', ucolin
       !   print *, 'myid,', myid, ' zout,', zout
       !   print *, 'myid,', myid, ' uout,', ucolout
       ! endif
        
        ! consider move this MPI comm out of i,j loop
        !call update_ghost(ucolout, level)       

        !! if wave surface > 0, first grid above surface use log interpolation

        !! if 0<z<wave surface, then use wave surface vel for those points below
        !! surface

        uout(i,j,1:xsz(3)) = ucolout(1:xsz(3))
      enddo
    enddo

    call update_ghost(uout, level)

    deallocate(ucolin, ucolout, zin, zout)
    if (dir .ne. 3) then
      deallocate(cartzcol, zzi)
    else 
      deallocate(cartzwcol, zwi)
    endif
  end subroutine comp2cart
  
  subroutine comp2cart_global (uin, zbasic, hho, dir)
    use param
    use decomp
    !use grid, only : eta, hh, zz, zw, hx, hy, her, level
    use mpi
    use lib_array
    implicit none
    real(wp), dimension(nx_global,ny_global,nz_global), intent(inout) :: uin
    real(wp), dimension(nz_global), intent(in) :: zbasic
    real(wp), dimension(nx_global, ny_global), intent(in) :: hho
    !real(wp), dimension(:,:,:), intent(out) :: uout
    integer :: dir ! dir should be 1, 2, or 3, 4, 5

    real(wp), dimension(nz_global) :: ucolin, ucolout, zin, zout
    integer :: i, j, k
    !real(wp), dimension(:), allocatable :: cartzcol, cartzwcol, zzi, zwi
    integer :: iprint

    iprint = 0
    !allocate(ucolin(1:nz_global))
    !allocate(ucolout(1:nz_global))
    !allocate(zin(1-level:xsz(3)+level))
    !allocate(zout(1-level:xsz(3)+level))
    ucolin = 0.0; ucolout = 0.0; !zin = 0.0; zout = 0.0
    !if (dir.ne.3) then
      !allocate(cartzcol(1-level:xsz(3)+level))
      !allocate(zzi(1-level:xsz(3)+level))
      !cartzcol = 0.0; zzi = zz * hbar; zout = zzi
    !else
    !  allocate(cartzwcol(1-level:xsz(3)+level))
    !  allocate(zwi(1-level:xsz(3)+level))
    !  cartzwcol = 0.0; zwi = zw * hbar; zout = zwi
    !endif
    
    do i = 1, nx_global
      do j = 1, ny_global
        ucolin = uin(i,j,:)
        
        !! interpolation. Fill with -1000.0 if out of range
        do k = 1, nz_global
          zin = zbasic(k) * (hbar+hho(i,j)) - hho(i,j)
          zout = zbasic(k) * hbar
        enddo
        ucolout = interp1d(zin, ucolin, zout, .false., -1000.0_wp)

        do k = 1, nz_global 
          if (zout(k)<(-hho(i,j))) then
            !> this should happen only in bottom subdomains.
            !> use the wave surface velocity directly. u,v,w overlapped here.
            ucolout(k) = ucolin(1)
          else if (zout(k)<zin(2)) then
            !> it need log law interpolation here. let's do it later.
            exit  !> replace it later.
          else if (zout(k)>=zin(2)) then
            exit
          endif
        enddo

        do k = 1, nz_global
          if (ucolout(k)<-10.0 .and. dir .eq. 1) then
            iprint = 1
          endif
        enddo
        
        if (iprint > 0 .and. myid.eq.0) then
          print *, 'myid,',myid,',tag=', dir,',i,j,',i,j
          print *, 'zin,', zin
          print *, 'ucolin,', ucolin
          print *, 'zout,',zout 
          print *, 'ucolout,',ucolout
        endif
        iprint = 0

        !> plyunote: we hope zout(1:xsz(3)) in range of zin(1-level:xsz(3)+level)
        !> for small wave amplitude case, this might be true
        !> for large wave amplitude case, or if dz is very small, here need
        !  interpolation in global zzall.

       ! if (i.eq.1 .and. j.eq.1) then
       !   print *, 'myid,', myid, ' zin,', zin
       !   print *, 'myid,', myid, ' uin,', ucolin
       !   print *, 'myid,', myid, ' zout,', zout
       !   print *, 'myid,', myid, ' uout,', ucolout
       ! endif
        
        ! consider move this MPI comm out of i,j loop
        !call update_ghost(ucolout, level)       

        !! if wave surface > 0, first grid above surface use log interpolation

        !! if 0<z<wave surface, then use wave surface vel for those points below
        !! surface

        uin(i,j,1:nz_global) = ucolout(1:nz_global)
      enddo
    enddo

    !call update_ghost(uout, level)

    !deallocate(ucolin, ucolout)
    !if (dir .ne. 3) then
    !  deallocate(cartzcol, zzi)
    !else 
    !  deallocate(cartzwcol, zwi)
    !endif
  end subroutine comp2cart_global

!  diff_central(x, y):
!    x0 = x[:-2]
!    x1 = x[1:-1]
!    x2 = x[2:]
!    y0 = y[:-2]
!    y1 = y[1:-1]
!    y2 = y[2:]
!    f = (x2 - x1)/(x2 - x0)
!    f1 = (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0)
!    f2 = x.copy()
!    f2[1:-1] = f1
!    f2[0] = f1[0]
!    f2[-1] = f1[-1]
!    return f2   

  subroutine one_side_diff_3elem(x_, y_, dydx_, i)
    implicit none
    real(wp), intent(in) :: x_(3), y_(3)
    integer, intent(in) :: i
    real(wp), intent(out) :: dydx_
    real(wp) :: dx1, dx2

    dx1 = x_(2)-x_(1); dx2 = x_(3)-x_(2)
    if (i.eq.1) then
      dydx_ = (y_(2)-y_(1))*(dx1+dx2)/(dx1*dx2)+(y_(1)-y_(3))*dx1/dx2/(dx1+dx2)
    elseif (i.eq.3) then
      dydx_ = (y_(3)-y_(2))*(dx1+dx2)/dx1/dx2-(y_(3)-y_(1))*dx2/dx1/(dx1+dx2)
    else
      print *, 'Not implemented for this i.(one_side_diff_3elem)'
      stop
    endif

  end subroutine one_side_diff_3elem

  subroutine center_diff(x_, y_, dydx_)
    implicit none
    real(wp), intent(in) :: x_(:), y_(:)
    real(wp), intent(out) :: dydx_(:)

    integer :: n, i
    real(wp), dimension(:), allocatable :: f_

    n = size(x_)
    allocate(f_(2:n-1))

    f_(:) = (x_(3:)-x_(2:n-1))/(x_(3:)-x_(1:n-2))
    do i = 2, n-1
      dydx_(i) = (1.-f_(i))*(y_(i+1)-y_(i))/(x_(i+1)-x_(i)) &
        + f_(i)*(y_(i)-y_(i-1))/(x_(i)-x_(i-1))
    enddo

    !!> lazy as me
    !dydx_(1) = (y_(2)-y_(1))/(x_(2)-x_(1))
    !dydx_(n) = (y_(n)-y_(n-1))/(x_(n)-x_(n-1))
    call one_side_diff_3elem(x_(1:3), y_(1:3), dydx_(1), 1)
    call one_side_diff_3elem(x_(n-2:n), y_(n-2:n), dydx_(n), 3)

    deallocate(f_)
  end subroutine center_diff

  subroutine calc_budget(z_, dpdx_, um_, wm_, upwpm_, uppwppm_, wtfxm_, &
    taub_, mb_, mkeb_, dsgsm_, dsgsum_, nutm_, nacfxm_, fname_)
    use decomp, only : myid
    implicit none
    real(wp), dimension(:), intent(in):: z_, um_, wm_, upwpm_, uppwppm_, wtfxm_
    real(wp), dimension(:), intent(in):: dsgsm_, dsgsum_, nutm_, nacfxm_
    real(wp), intent(in) :: dpdx_
    character(len=*), intent(in) :: fname_ 

    !> mb: momentum budget
    !> mkeb: mean kinetic energy budget
    !real(wp), intent(out) :: taub_(nz_global, 3), mb_(nz_global, 5), &
    !  mkeb_(nz_global, 6)
    real(wp), intent(out) :: taub_(:,:), mb_(:,:), mkeb_(:,:)

    integer :: k, n
    
    n = size(z_)

    taub_(:,2) = -upwpm_
    taub_(:,3) = -uppwppm_
    taub_(:,1) = -upwpm_ - uppwppm_

    !> momentum term 1: sum of following terms
    !> momentum term 2: pressure gradient
    !> momentum term 3: reynolds stress
    !> momentum term 4: dispersive shear stress
    !> momentum term 5: wind turbine force
    !> momentum term 6: SubGridScale (SGS) stress
    mb_(:,2) = dpdx_
    call center_diff(z_, -upwpm_, mb_(:,3))
    call center_diff(z_, -uppwppm_, mb_(:,4))
    mb_(:,5) = wtfxm_(:)
    mb_(:,6) = dsgsm_(:)
    call center_diff(z_, um_, mb_(:,7))
    mb_(:,7) = -wm_(:) * mb_(:,7)
    mb_(:,8) = nacfxm_
    mb_(:,1) = mb_(:,2)+mb_(:,3)+mb_(:,4)+mb_(:,5)+mb_(:,6)+mb_(:,8)

    !> mean kinetic energy part 1: work done by pressure
    mkeb_(:,3) = dpdx_ * um_
    call center_diff(z_, taub_(:,1)*um_, mkeb_(:,4))
    call center_diff(z_, um_, mkeb_(:,5))
    mkeb_(:,5) = mkeb_(:,5) * taub_(:,1)
    mkeb_(:,6) = -wtfxm_ * um_
    mkeb_(:,7) = dsgsum_(:)
    mkeb_(:,8) = mb_(:,7) * um_
    mkeb_(:,9) = -nacfxm_ * um_
    mkeb_(:,1) = mkeb_(:,3) + mkeb_(:,4)
    mkeb_(:,2) = mkeb_(:,5) + mkeb_(:,6) + mkeb_(:,7) + mkeb_(:,9)
    
    if (myid.eq.0) then
      !open(1004, file="mean_field_1d.dat")
      print *, 'Writing output to ', trim(fname_)
      open(1004, file=trim(fname_))
      write(1004,*) 'variables=z, um, wm, nutm', &
        ', taub_sum, taub_upwp, taub_uppwpp', &
        ', mb_sum, mb_dpdx, mb_upwp, mb_uppwpp, mb_wtfx, mb_nacfx, mb_sgs, mb_conv', &
        ', mkeb_sum1, mkeb_sum2, mkeb_dpdx, mkeb_phie, mkeb_epsie, mkeb_wt', &
        ', mkeb_nacfx, mkeb_sgsu, mkeb_conv'
      write(1004,*) 'ZONE I=1 J=1 K=', n, ' DATAPACKING=POINT'
      do k = 1, n
        write(1004, *) z_(k), um_(k), wm_(k), nutm_(k), taub_(k,1), taub_(k,2),&
          taub_(k,3), mb_(k,1), mb_(k,2), mb_(k,3), mb_(k,4), mb_(k,5),mb_(k,6),&
          mb_(k,7), mb_(k,8), mkeb_(k,1), mkeb_(k,2), mkeb_(k,3), mkeb_(k,4), &
          mkeb_(k,5), mkeb_(k,6), mkeb_(k,7), mkeb_(k,8), mkeb_(k,9)
      enddo
      close(1004)
    endif
  end subroutine calc_budget


  subroutine reread_rotor_force (time, it, level)
    use param, only : iturbine
    use navier, only : wtforce, wtforce_y, wtforce_z, u, flag_de, fturbinex,&
      fturbiney   
    use turbine_model
    implicit none
    real(wp) :: time
    integer :: it, level

    integer :: ibi, n_elmt_1, i, j, nb, n1e, n2e, n3e

   !> we use output file 'line_XXXXXX_XXX_nf.dat' from Export_LineLocation
    if (iturbine .eq. 1) then
      call discloca_v3
      call veldisc_v2(u, level, time)
      call wind_turbine_force(wtforce, level, fturbinex)
    elseif (iturbine .eq. 3) then
      call Read_LineLocation(wtm, NumberOfTurbines, it)
      !> need these variables:
      !! Pre_process: ibm(ibi)%i_min 
      !! rotor_Rot: ibm(ibi)%dA(l), ibm(ibi)%cent_x(l)
      !! collect_grid, update_zgrid: cartx, carty, cartz

      ! collect_grid has been executed in read_turbine_control
      call update_zgrid

      do ibi = 1, NumberOfTurbines
        n_elmt_1 = wtm(ibi)%n_elmt / num_blade
        do nb = 1, num_blade
          do j = 1, n_elmt_1
            i = (nb-1) * n_elmt_1 + j
            n1e = wtm(ibi)%nv1(i); n2e = wtm(ibi)%nv2(i)
            
            !> the lines until dA are not necessary. Because the value doesn't
            !! change with timesteps
            !dr_vec(1) = wtm(ibi)%x_bp(n2e) - wtm(ibi)%x_bp(n1e)
            !dr_vec(2) = wtm(ibi)%y_bp(n2e) - wtm(ibi)%y_bp(n1e)
            !dr_vec(3) = wtm(ibi)%z_bp(n2e) - wtm(ibi)%z_bp(n1e)
            !wtm(ibi)%dA(i) = norm2(dr_vec)
            
            wtm(ibi)%cent_x(i) = (wtm(ibi)%x_bp(n1e)+wtm(ibi)%x_bp(n2e))/2.
            wtm(ibi)%cent_y(i) = (wtm(ibi)%y_bp(n1e)+wtm(ibi)%y_bp(n2e))/2.
            wtm(ibi)%cent_z(i) = (wtm(ibi)%z_bp(n1e)+wtm(ibi)%z_bp(n2e))/2.
          enddo
        enddo
      enddo

      call Pre_process(wtm, NumberOfTurbines, 1)
      
      call Calc_F_eul(wtforce, wtforce_y, wtforce_z, fturbinex, fturbiney, level, &
        wtm, fsi_wt, NumberOfTurbines, 1.0_wp, 1, flag_de) 
      if(myid.eq.0) print *, 'bforce=', bforce, ', fturb_x=', fturbinex 
    
    elseif (iturbine .eq. 5 .or. iturbine .eq. 7) then
      call Read_SurfaceLocation(wtm, NumberOfTurbines, it)
      call update_zgrid

      do ibi = 1, NumberOfTurbines
        do i = 1, wtm(ibi)%n_elmt
          n1e = wtm(ibi)%nv1(i); n2e = wtm(ibi)%nv2(i); n3e = wtm(ibi)%nv3(i)
          
          !> the lines until dA are not necessary. Because the value doesn't
          !! change with timesteps
          !ds12(1) = wtm(ibi)%x_bp(n2e) - wtm(ibi)%x_bp(n1e)
          !ds12(2) = wtm(ibi)%y_bp(n2e) - wtm(ibi)%y_bp(n1e)
          !ds12(3) = wtm(ibi)%z_bp(n2e) - wtm(ibi)%z_bp(n1e)

          !ds13(1) = wtm(ibi)%x_bp(n3e) - wtm(ibi)%x_bp(n1e)
          !ds13(2) = wtm(ibi)%y_bp(n3e) - wtm(ibi)%y_bp(n1e)
          !ds13(3) = wtm(ibi)%z_bp(n3e) - wtm(ibi)%z_bp(n1e)

          !call crossx(ds12, ds13, tmparray)
          !dr = norm2(tmparray)
          !wtm(ibi)%dA(i) = dr/2.0_wp
          
          wtm(ibi)%cent_x(i) = (wtm(ibi)%x_bp(n1e)+wtm(ibi)%x_bp(n2e)+wtm(ibi)%x_bp(n3e))/3.0_wp
          wtm(ibi)%cent_y(i) = (wtm(ibi)%y_bp(n1e)+wtm(ibi)%y_bp(n2e)+wtm(ibi)%y_bp(n3e))/3.0_wp
          wtm(ibi)%cent_z(i) = (wtm(ibi)%z_bp(n1e)+wtm(ibi)%z_bp(n2e)+wtm(ibi)%z_bp(n3e))/3.0_wp
        enddo
      enddo
      
      call Pre_process(wtm, NumberOfTurbines, 1)
      
      call Calc_F_eul(wtforce, wtforce_y, wtforce_z, fturbinex, fturbiney, level, &
        wtm, fsi_wt, NumberOfTurbines, 1.0_wp, 1, flag_de) 
      if(myid.eq.0) print *, 'bforce=', bforce, ', fturb_x=', fturbinex 
    else
      print *, 'Forces for this turbine mode have not been implemented yet.'
    endif
  end subroutine reread_rotor_force

  subroutine reread_nacelle_force(it, fx, fy, fz, level)
    use param, only : iturbine
    use decomp
    use turbine_model
    use navier, only : flag_de, fturbinex, fturbiney
    implicit none
    real(wp), dimension(:,:,:) :: fx, fy, fz
    integer :: it, level

    integer :: ibi, i, n1e, n2e, n3e

    if (iturbine .ne. 0 .and. nacelle_model .eq. 1) then
      call Read_NacelleLocation(ibm_nac, NumberOfNacelle, it)
      ! call update_zgrid

      if (fsitype .eq. 1) then
        do ibi = 1, NumberOfNacelle
          do i = 1, ibm_nac(ibi)%n_elmt
            n1e = ibm_nac(ibi)%nv1(i); n2e = ibm_nac(ibi)%nv2(i); n3e = ibm_nac(ibi)%nv3(i)
            ibm_nac(ibi)%cent_x(i) = (ibm_nac(ibi)%x_bp(n1e)+ibm_nac(ibi)%x_bp(n2e) &
              +ibm_nac(ibi)%x_bp(n3e))/3.0_wp
            ibm_nac(ibi)%cent_y(i) = (ibm_nac(ibi)%y_bp(n1e)+ibm_nac(ibi)%y_bp(n2e) &
              +ibm_nac(ibi)%y_bp(n3e))/3.0_wp
            ibm_nac(ibi)%cent_z(i) = (ibm_nac(ibi)%z_bp(n1e)+ibm_nac(ibi)%z_bp(n2e) &
              +ibm_nac(ibi)%z_bp(n3e))/3.0_wp
          enddo

          !if(myid.eq.0) print *, "plyudebug: read_nacelle_param, Pre_process 1"
          call Pre_process(ibm_nac, NumberOfNacelle, 1)

          !if(myid.eq.0) print *, "plyudebug: read_nacelle_param, Coordinates_IP"
          call Coordinates_IP(ibm_nac, NumberOfNacelle)

          !> Pre_process_IP is integrated into Pre_process
          !if(myid.eq.0) print *, "plyudebug: read_nacelle_param, Pre_process 2"
          call Pre_process(ibm_nac, NumberOfNacelle, 2)
        enddo
      endif

      ibi = 1; i = 1
      n1e = ibm_nac(ibi)%nv1(i); n2e = ibm_nac(ibi)%nv2(i);
      n3e = ibm_nac(ibi)%nv3(i)
      !print*, "bp=",ibm_nac(1)%x_bp(n1e), ibm_nac(1)%x_bp(n2e), ibm_nac(1)%x_bp(n3e), &
      !  ", centx=",ibm_nac(1)%cent_x(1), ", i_min=", ibm_nac(1)%i_min(1), &
      !  ", i_max=", ibm_nac(1)%i_max(1), ", F_lagr_x=", ibm_nac(1)%F_lagr_x(1)

      fx = 0.0; fy = 0.0; fz = 0.0
      call Calc_F_eul(fx(:,:,1:xsz(3)), fy(:,:,1:xsz(3)), fz(:,:,1:xsz(3)), &
        fturbinex, fturbiney, level, ibm_nac, fsi_nac, NumberOfNacelle, 1.0_wp, 2, flag_de)
      call update_ghost(fx, level)
    endif

  end subroutine reread_nacelle_force

!> plyunote:
! 6. calculate budget of mean momentum and mean kinetic energy
! 5. calculate surface effective friction velocity
! 4. zw interpolate to xypoint

! 3. deallocate
! 2. consider interpolation scheme for z < eta
! 1. add z to h5 output

  subroutine post_les_hos_turbine
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

    !> added by plyu
    use turbine_model
    use solver_common, only : tmp_x6
    use lib_array, only: is_nan
    use discontinuity_smooth, only: analyze_limit3d_from_3d, tbn_lim_x, &
        tbn_lim_y, lim_dir 
    
    use hdf_io
    !end
    
    implicit none
    
    integer :: ioutc, ioutd
    real(wp) :: time
    
    integer :: ierror, ierr
    
    double precision :: t1, t2
    
    !HOS variables
    integer icon
    
    integer(INT64) :: id
    integer :: i, j, k, it, root0
    
    real(wp), allocatable, dimension(:) :: fpk, betak, ampfd, epskx
    
    real(wp) beta_miles
    
    character (len=100) :: fileid, fgrid, fparam_les, fdata_les, faux, fparam_hos, fdata_hos

    real(wp) c_phase, u_lambda, u_half_lambda, wage1, wage2, wage3
    real(wp) clen, cspeed, ctime, x, y, fpt

    real(wp), allocatable, dimension(:,:) :: etaall, hhall, tmp
    real(wp) ftn, rek, dwk, disspm, l_kom_ave

    real(wp), allocatable, dimension(:,:,:) :: uupall, wwpall, velall, ppall, etasig, pp_wave
    real(wp), allocatable, dimension(:,:,:,:) :: ppsig
    real(wp), allocatable, dimension(:,:) :: dtheta, uxy, vxy, wxy
    real(wp), allocatable, dimension(:,:) :: ek11, ek22, ek33
    real(wp), allocatable, dimension(:,:) :: ek11_tm, ek22_tm, ek33_tm
    real(wp), allocatable, dimension(:,:) :: eki11, eki22, eki33
    real(wp), allocatable, dimension(:,:) :: eki11_tm, eki22_tm, eki33_tm
    real(wp), allocatable, dimension(:) :: dissp, ekm, l_kom, ekim

    real(wp), allocatable, dimension(:) :: q1, q2, q3, q4, um, vm, wm, uup, vvp, wwp, zzall, uw_shear
   
    !> variables added by plyu
    integer :: idp ! counter for debug
    integer :: ini_it, end_it, nk, n_it, rem_it
    integer :: rankrow, rankcol, sizerow, sizecol
    real(wp) :: deltax, deltay
    real(wp), allocatable, dimension(:) :: um_l, vm_l, wm_l 
    real(wp), allocatable, dimension(:) :: upwpm_l, uppwppm_l
    real(wp), allocatable, dimension(:) :: wtfxm_l
    real(wp), allocatable, dimension(:) :: upwpm, uppwppm
    real(wp), allocatable, dimension(:,:) :: taubm, mbm, mkebm

    real(wp), allocatable, dimension(:,:,:) :: u_tm, v_tm, w_tm, pp_tm
    real(wp), allocatable, dimension(:,:,:) :: uu_tm, vv_tm, ww_tm, uw_tm, tke_tm
    real(wp), allocatable, dimension(:,:,:) :: templocal3d, temp_de_1 
    real(wp) :: pp_ref

    real(wp), allocatable, dimension(:) :: wtfxm, wtfym, wtfzm 
    real(wp), allocatable, dimension(:,:,:) :: wtfx_tm, wtfy_tm, wtfz_tm
    !> u^\prime
    !real(wp), allocatable, dimension(:,:,:) :: up1_tm, vp1_tm, wp1_tm
    !> \bar{u}^{\prime\prime}
    !real(wp), allocatable, dimension(:,:,:) :: up2_tm, vp2_tm, wp2_tm

    real(wp) :: volflux(3)

    type(HDFObj) :: writer
    integer :: ih5s

    !> ui, vi, wi are velocity components described in cartesian grid, which are
    !! interpolated from boundary-fited grid
    integer :: i_interp
    real(wp), allocatable, dimension(:,:,:) :: ui, vi, wi, ppi
    real(wp), allocatable, dimension(:,:,:) :: ui_tm, vi_tm, wi_tm, ppi_tm
    real(wp), allocatable, dimension(:,:,:) :: uui_tm, vvi_tm, wwi_tm, uwi_tm, tkei_tm
    real(wp), allocatable, dimension(:) :: uim, vim, wim, upwpim, uppwppim
    real(wp), allocatable, dimension(:,:) :: taubim, mbim, mkebim
    real(wp), allocatable, dimension(:,:,:) :: wtfxi, wtfxi_tm 
    real(wp), allocatable, dimension(:) :: wtfxim
    !> u^\prime
    !real(wp), allocatable, dimension(:,:,:) :: up1i_tm, vp1i_tm, wp1i_tm
    !> \bar{u}^{\prime\prime}
    !real(wp), allocatable, dimension(:,:,:) :: up2i_tm, vp2i_tm, wp2i_tm

    !> local variables for turbine model
    integer :: ibi, n_elmt_1, nb, n1e, n2e, n3e
    real(wp), dimension(3) :: dr_vec, ds12, ds13, tmparray
    real(wp) :: dr
    
    integer :: ii, iii, jjj, itemp, j_1d1, ktemp
    integer :: post_n_1d1, post_n_1d2, post_n_2d1, post_n_2d2, post_n_2d3
    integer, allocatable, dimension(:) :: post_i_1d1, post_i_1d2
    integer, allocatable, dimension(:) :: post_i_2d1, post_i_2d2, post_i_2d3
    character (len=64) :: fout
    real(wp) :: tb_xc, tb_yc, tb_zc, tb_d, xtemp, ytemp, ztemp
    real(wp), allocatable, dimension(:) :: tempz, zwall
    real(wp), allocatable, dimension(:) :: coord_1d1, coord_1d2 
    real(wp), allocatable, dimension(:) :: coord_2d1, coord_2d2, coord_2d3 
    !> coord_1d1_g is the grid coord nearest to coord_1d1
    real(wp), allocatable, dimension(:) :: coord_1d1_g, coord_1d2_g 
    real(wp), allocatable, dimension(:) :: coord_2d1_g, coord_2d2_g, coord_2d3_g
    real(wp), allocatable, dimension(:,:) :: tempw_1d1
    real(wp), allocatable, dimension(:,:,:) :: temp3d
    real(wp) :: utopm
    
    !> comment this line when turbine_model is used
    !real(wp), allocatable, dimension(:,:) :: eo, hho, hxo, hyo
    !real(wp), allocatable, dimension(:,:) :: hxo 

    real(wp), allocatable, dimension(:,:) :: tauwxo, tauwyo, pso
    real(wp), allocatable, dimension(:,:) :: tau_vis, tau_pre
    real(wp), allocatable, dimension(:,:) :: tau_viso, tau_preo 
    real(wp) :: tau_vis_m, tau_pre_m, tau_all_m, utau_eff, utau_eff_tm
    
    integer :: nm_eta
    real(wp) :: eta_nm(5)
    
    real(wp), allocatable, dimension(:,:) :: eta_sp, eta_theta 

    !> variable for SGS
    real(wp), allocatable, dimension(:,:) :: t13wx
    real(wp), allocatable, dimension(:,:,:) :: dsgs_tm, dsgsu_tm, nut_tm
    real(wp), allocatable, dimension(:,:,:) :: dsgsi_tm, dsgsui_tm, nuti_tm
    real(wp), allocatable, dimension(:) :: dsgsm, dsgsum, nutm
    real(wp), allocatable, dimension(:) :: dsgsim, dsgsuim, nutim
    real(wp), allocatable, dimension(:) :: dsgsm_l, dsgsum_l, nutm_l

    !> variable for nacelle model force
    real(wp), allocatable, dimension(:,:,:) :: nacfx, nacfy, nacfz
    real(wp), allocatable, dimension(:,:,:) :: nacfxi, nacfx_tm, nacfxi_tm
    real(wp), allocatable, dimension(:) :: nacfxm, nacfxim, nacfxm_l

    !> u_1d1: velocity along a vertical line, col1 is index of coordinates,
    !         col2 is component index which indicates it is u, v or w,
    !         col3 is index of section in list of post-process.
    !> um_1d1: mean velocity
    !> up_1d1: u prime, namely the fluctuating part of velocity
    !> uu_1d1: u'v', and so on. col2 and col3 are its indexes.
    real(wp), allocatable, dimension(:,:,:) :: u_1d1, ui_1d1
    real(wp), allocatable, dimension(:,:,:) :: um_1d1, uim_1d1
    real(wp), allocatable, dimension(:,:,:) :: up_1d1
    real(wp), allocatable, dimension(:,:,:,:) :: uu_1d1

    !> u_1d2: velocity along a streamwise horizontal line
    !real(wp), allocatable, dimension(:,:,:) :: u_1d2
    !real(wp), allocatable, dimension(:,:,:) :: um_1d2
    !real(wp), allocatable, dimension(:,:,:) :: up_1d2
    !real(wp), allocatable, dimension(:,:,:,:) :: uu_1d2
    
    !> u_2d1: velocity in a vertical plane which is perpendicular to x-axis, 
    !         col1,2 is index of coordinates,
    !         col3 is component index which indicates it is u, v or w,
    !         col4 is index of section in list of post-process.
    !> um_2d1: mean velocity
    !> up_2d1: u prime, namely the fluctuating part of velocity
    !> uu_2d1: u'v', and so on. col2 and col3 are its indexes.
    real(wp), allocatable, dimension(:,:,:,:) :: u_2d1, ui_2d1
    real(wp), allocatable, dimension(:,:,:,:) :: um_2d1, uim_2d1
    real(wp), allocatable, dimension(:,:,:,:) :: up_2d1
    real(wp), allocatable, dimension(:,:,:,:,:) :: uu_2d1
    
    !> u_2d2: velocity in a horizontal plane perpendicular to z-axis
    real(wp), allocatable, dimension(:,:,:,:) :: u_2d2, ui_2d2
    real(wp), allocatable, dimension(:,:,:,:) :: um_2d2, uim_2d2
    real(wp), allocatable, dimension(:,:,:,:) :: up_2d2
    real(wp), allocatable, dimension(:,:,:,:,:) :: uu_2d2
    
    !> u_2d3: velocity in a vertical plane perpendicular to y-axis
    real(wp), allocatable, dimension(:,:,:,:) :: u_2d3, ui_2d3
    real(wp), allocatable, dimension(:,:,:,:) :: um_2d3, uim_2d3
    real(wp), allocatable, dimension(:,:,:,:) :: up_2d3
    real(wp), allocatable, dimension(:,:,:,:,:) :: uu_2d3

    !> vortex dynamics: it is implemented in another  individual subroutine
    !real(wp), allocatable, dimension(:,:,:,:,:) :: d_u_x
    !real(wp), allocatable, dimension(:,:,:,:) :: vor, lap_vor

    logical :: inan
   
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
    
    call grid_gen
     
    icon = 16
    
    if (myid .eq. 0) then
      print *, "********* This is les-hos-turbine post-processing script *******"
    endif

    !allocate(fpk(nx_global))
    !allocate(betak(nx_global))
    !allocate(ampfd(nx_global))
    !allocate(epskx(nx_global))
    !allocate(dissp(nz_global))
    !allocate(l_kom(nz_global))

    !allocate(etaall(nx_global,ny_global))
    !allocate(hhall(nx_global,ny_global))
    allocate(tmp(nx_global,ny_global))

    !allocate(uupall(nx_global,ny_global,nz_global))
    !allocate(wwpall(nx_global,ny_global,nz_global))
    !allocate(velall(nx_global,ny_global,nz_global))
    !allocate(ppall(nx_global,ny_global,nz_global))

    !allocate(uxy(nx_global,ny_global))
    !allocate(vxy(nx_global,ny_global))
    !allocate(wxy(nx_global,ny_global))
    !
    nk = nx_global / 3
    allocate(ek11(nk,nz_global), ek11_tm(nk,nz_global))
    allocate(ek22(nk,nz_global), ek22_tm(nk,nz_global))
    allocate(ek33(nk,nz_global), ek33_tm(nk,nz_global))
    allocate(ekm(nx_global))
    ek11_tm = 0.0_wp; ek22_tm = 0.0_wp; ek33_tm = 0.0_wp
    if (iwavy .eq. 7) then
      allocate(eki11(nk,nz_global), eki11_tm(nk,nz_global))
      allocate(eki22(nk,nz_global), eki22_tm(nk,nz_global))
      allocate(eki33(nk,nz_global), eki33_tm(nk,nz_global))
      allocate(ekim(nx_global))
      eki11_tm = 0.0_wp; eki22_tm = 0.0_wp; eki33_tm = 0.0_wp
    endif

    !allocate(q1(nz_global))
    !allocate(q2(nz_global))
    !allocate(q3(nz_global))
    !allocate(q4(nz_global))

    allocate(um(nz_global), vm(nz_global), wm(nz_global))
    allocate(um_l(nz_global), vm_l(nz_global), wm_l(nz_global))
    allocate(upwpm_l(nz_global), uppwppm_l(nz_global))
    allocate(upwpm(nz_global), uppwppm(nz_global))
    allocate(wtfxm(nz_global))
    allocate(wtfxm_l(nz_global))
    allocate(taubm(nz_global,3))
    allocate(mbm(nz_global,8))
    allocate(mkebm(nz_global,9))

    um = 0.0; vm = 0.0; wm = 0.0
    um_l = 0.0; vm_l = 0.0; wm_l = 0.0
    upwpm_l = 0.0; uppwppm_l = 0.0
    upwpm = 0.0; uppwppm = 0.0; 
    wtfxm = 0.0; wtfxm_l = 0.0
    taubm = 0.0; mbm = 0.0; mkebm = 0.0
    
    if (iwavy .eq. 7) then
      allocate(uim(nz_global), vim(nz_global), wim(nz_global))
      allocate(upwpim(nz_global), uppwppim(nz_global))
      allocate(wtfxim(nz_global))
      allocate(taubim(nz_global,3))
      allocate(mbim(nz_global,8))
      allocate(mkebim(nz_global,9))
      uim = 0.0; vim = 0.0; wim = 0.0
      upwpim = 0.0; uppwppim = 0.0
      wtfxim = 0.0
      taubim = 0.0; mbim = 0.0; mkebim = 0.0
    endif
    !allocate(um(nz_global), uup(nz_global))
    !allocate(vm(nz_global), vvp(nz_global))
    !allocate(wm(nz_global), wwp(nz_global))
    !allocate(uw_shear(nz_global))
    !allocate(zzall(nz_global))

    !read(14,*) clen, cspeed, ctime

    level = 1
    call mpi_comm_rank(mpi_comm_2d_row, rankrow, ierr) ! rank in z
    call mpi_comm_rank(mpi_comm_2d_col, rankcol, ierr) ! rank in y
    call mpi_comm_size(mpi_comm_2d_row, sizerow, ierr)
    call mpi_comm_size(mpi_comm_2d_col, sizecol, ierr)
    if(myid.eq.0) print *, 'sizerow=',sizerow,', sizecol=',sizecol
    
    ! Gather zz coordinates
    allocate(tempz(nz_global), zzall(nz_global), zwall(nz_global))
    tempz = 0.0; zzall = 0.0; zwall = 0.0 
    
    tempz(xst(3):xend(3)) = zz(1:xsz(3))
    call MPI_Allreduce(tempz, zzall, nz_global, mpi_double_precision, &
         MPI_SUM, MPI_COMM_2D_ROW, ierror)
    tempz = 0.0
    tempz(xst(3):xend(3)) = zw(1:xsz(3))
    call MPI_Allreduce(tempz, zwall, nz_global, mpi_double_precision, &
         MPI_SUM, MPI_COMM_2D_ROW, ierror)
    
    allocate(tauwxo(nx,ny), tauwyo(nx,ny), pso(nx,ny))
    allocate(tau_vis(xsz(1),xsz(2)), tau_pre(xsz(1),xsz(2)))
    !allocate(hxo(nx,ny))
    allocate(tau_viso(nx,ny), tau_preo(nx,ny))
    utau_eff_tm = 0.0_wp
    tauwxo(:,:) = 0.0_wp; tauwyo(:,:) = 0.0_wp; pso(:,:) = 0.0_wp
    tau_viso(:,:) = 0.0_wp; tau_preo(:,:) = 0.0_wp

    allocate(eta_sp(nx,ny), eta_theta(nx,ny))
    eta_sp(:,:) = 0.0; eta_theta(:,:) = 0.0
    
    allocate(t13wx(xsz(1),xsz(2)))
    allocate(dsgs_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(dsgsu_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(nut_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(dsgsm(nz_global), dsgsum(nz_global), nutm(nz_global))
    allocate(dsgsm_l(nz_global), dsgsum_l(nz_global), nutm_l(nz_global))
    dsgs_tm = 0.0; dsgsu_tm = 0.0; nut_tm = 0.0
    dsgsm(:) = 0.0; dsgsum(:) = 0.0; nutm(:) = 0.0

    if (iwavy .eq. 7) then
      allocate(dsgsi_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(dsgsui_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(nuti_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(dsgsim(nz_global), dsgsuim(nz_global), nutim(nz_global))
      dsgsi_tm = 0.0; dsgsui_tm = 0.0; nuti_tm = 0.0
      dsgsim(:) = 0.0; dsgsuim(:) = 0.0; nutim(:) = 0.0
    endif
    
    allocate(nacfx(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(nacfy(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(nacfz(xsz(1), xsz(2), 1-level:xsz(3)+level))
    nacfx = 0.0; nacfy = 0.0; nacfz = 0.0

    allocate(nacfx_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(nacfxm(nz_global))
    nacfx_tm = 0.0; nacfxm = 0.0 
    
    allocate(nacfxm_l(nz_global))
    nacfxm_l = 0.0

    if (iwavy .eq. 7) then
      allocate(nacfxi(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(nacfxi_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(nacfxim(nz_global))
      nacfxi = 0.0; nacfxi_tm = 0.0; nacfxim = 0.0 
    endif

    !allocate(d_u_x(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    !allocate(vor(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    !allocate(lap_vor(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    !d_u_x = 0.0; vor = 0.0;! lap_vor = 0.0

    !allocate(d_u_xm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    !allocate(vorm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    !allocate(lap_vorm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    !d_u_xm = 0.0; vorm = 0.0; lap_vorm = 0.0
    
    !allocate(d_u_xi(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    !allocate(vori(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    !allocate(lap_vori(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))

    !allocate(d_u_xim(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    !allocate(vorim(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    !allocate(lap_vorim(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))

    volflux = 0.0
    
    if (myid == 0) then
       open(66,file="wind_input.dat")
       open(68,file="turb_spectrum.dat")
       open(70,file="mean_flux.dat")
       open(72,file="q_analysis.dat")
    end if

   
    if (myid.eq.0) then
      do i = 1, nz_global
        print *, "i, zzall, zwall", i, zzall(i), zwall(i)
      enddo
    endif
    
    root0 = 0
    ini_it = 2600


    open(1001, file="post_control.inp")
    read(1001, *) ini_it, end_it
    !read(1001, *) post_n_1d1, post_n_1d2, post_n_2d1, post_n_2d2
    read(1001, *) tb_xc, tb_yc, tb_zc
    read(1001, *) tb_d
    close(1001)
    
    ! interpolate velocity from boundary fitted grid to cartesian grid
    if (iwavy .eq. 7) then
      allocate(ui(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(vi(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(wi(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(ppi(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(wtfxi(xsz(1), xsz(2), 1-level:xsz(3)+level))
      
      ui = 0.0; vi = 0.0; wi = 0.0; ppi = 0.0
      wtfxi = 0.0
      !call comp2cart(u,ui,1); call comp2cart(v,vi,2); call comp2cart(w,wi,3)
      !call comp2cart(pp,ppi,1)
    endif

    allocate(u_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(v_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(w_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(pp_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(uw_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(tke_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(templocal3d(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(temp_de_1(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(wtfx_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))

    u_tm(:,:,:) = 0.0; v_tm(:,:,:) = 0.0; w_tm(:,:,:) = 0.0;
    pp_tm(:,:,:) = 0.0; uw_tm(:,:,:) = 0.0; tke_tm(:,:,:) = 0.0
    templocal3d(:,:,:) = 0.0; temp_de_1(:,:,:) = 0.0
    wtfx_tm(:,:,:) = 0.0

    if (iwavy .eq. 7) then
      allocate(ui_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(vi_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(wi_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(ppi_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(uwi_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(tkei_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      allocate(wtfxi_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
      ui_tm(:,:,:) = 0.0; vi_tm(:,:,:) = 0.0; wi_tm(:,:,:) = 0.0
      ppi_tm(:,:,:) = 0.0; uwi_tm(:,:,:) = 0.0; tkei_tm(:,:,:) = 0.0
      wtfxi_tm = 0.0
    endif

    !> comment this line when turbine_model is used
    !allocate(eo(nx_global, ny_global), hho(nx_global, ny_global))

    if(mod(end_it-ini_it, noutd).ne.0) then
      print *, 'end_it-ini_it and noutd are not dividable'
      stop      
    endif

    open(1002, file="post_coord_input.inp")
    
    read(1002, *) post_n_1d1
    allocate(post_i_1d1(post_n_1d1))
    allocate(coord_1d1(post_n_1d1),coord_1d1_g(post_n_1d1))
    j_1d1 = nint(tb_yc/(yl/ny))

    do ii = 1, post_n_1d1
      read(1002, *) coord_1d1(ii)
      !tb_post_x(ii) = tb_xc + tb_d * tb_post_x(ii)
      post_i_1d1(ii) = nint(coord_1d1(ii) / (xl/nx))
      if (post_i_1d1(ii)<1) then
        post_i_1d1(ii) = 1
      else if (post_i_1d1(ii)>nx_global) then
        post_i_1d1(ii) = nx_global
      endif
      coord_1d1_g(ii) = post_i_1d1(ii) * (xl/nx)
      if(myid .eq. 0) then
        print *, 'coord_1d1: ii, ix, x_g, x=', ii, post_i_1d1(ii), &
          coord_1d1_g(ii), coord_1d1(ii)
      endif
    enddo
    
    read(1002, *) post_n_2d1
    allocate(post_i_2d1(post_n_2d1))
    allocate(coord_2d1(post_n_2d1), coord_2d1_g(post_n_2d1))
    !j_2d1 = nint(tb_yc/(yl/ny))

    do ii = 1, post_n_2d1
      read(1002, *) coord_2d1(ii)
      !tb_post_x(ii) = tb_xc + tb_d * tb_post_x(ii)
      post_i_2d1(ii) = nint(coord_2d1(ii) / (xl/nx))
      if (post_i_2d1(ii)<1) then
        post_i_2d1(ii) = 1
      else if (post_i_2d1(ii)>nx_global) then
        post_i_2d1(ii) = nx_global
      endif
      coord_2d1_g(ii) = post_i_2d1(ii) * (xl/nx)
      if(myid .eq. 0) then
        print *, 'coord_2d1: ii, ix, x_g, x=', ii, post_i_2d1(ii), &
          coord_2d1_g(ii), coord_2d1(ii)
      endif
    enddo
    
    read(1002, *) post_n_2d2
    allocate(post_i_2d2(post_n_2d2))
    allocate(coord_2d2(post_n_2d2), coord_2d2_g(post_n_2d2))

    do ii = 1, post_n_2d2
      read(1002, *) coord_2d2(ii)
      !tb_post_x(ii) = tb_xc + tb_d * tb_post_x(ii)
      post_i_2d2(ii) = 0 
      
      do k = 1, nz_global-1
        if (zzall(k)*hbar <= coord_2d2(ii) .and. zzall(k+1)*hbar > coord_2d2(ii)) then
          post_i_2d2(ii) = k
        endif
      enddo
      
      if (post_i_2d2(ii)<1) then
        post_i_2d2(ii) = 1
      else if (post_i_2d2(ii)>nz_global) then
        post_i_2d2(ii) = nz_global
      endif
          
      coord_2d2_g(ii) = zzall(post_i_2d2(ii))*hbar

      !j_1d1 = nint(tb_yc/(yl/ny))
      if(myid .eq. 0) then
        print *, 'coord_2d2: ii, ix, x_g, x=', ii, post_i_2d2(ii), &
          coord_2d2_g(ii), coord_2d2(ii)
      endif
    enddo
    
    read(1002, *) post_n_2d3
    allocate(post_i_2d3(post_n_2d3))
    allocate(coord_2d3(post_n_2d3), coord_2d3_g(post_n_2d3))
    !j_2d3 = nint(tb_yc/(yl/ny))

    do ii = 1, post_n_2d3
      read(1002, *) coord_2d3(ii)
      !tb_post_x(ii) = tb_xc + tb_d * tb_post_x(ii)
      post_i_2d3(ii) = nint(coord_2d3(ii) / (yl/ny))
      if (post_i_2d3(ii)<1) then
        post_i_2d3(ii) = 1
      else if (post_i_2d3(ii)>ny_global) then
        post_i_2d3(ii) = ny_global
      endif
      coord_2d3_g(ii) = post_i_2d3(ii) * (yl/ny)
      if(myid .eq. 0) then
        print *, 'coord_2d3: ii, ix, x_g, x=', ii, post_i_2d3(ii), &
          coord_2d3_g(ii), coord_2d3(ii)
      endif
    enddo
    
    close(1002)

    if(myid.eq.0) print *, 'Global nx, ny, nz=', nx_global, ny_global, nz_global
    
    allocate(temp3d(nx_global, ny_global, nz_global))
    ! recall the definition
    !> u_1d1: velocity along a vertical line, col1 is index of coordinates,
    !         col2 is component index which indicates it is u, v or w,
    !         col3 is index of section in list of post-process.
    !> um_1d1: mean velocity
    !> up_1d1: u prime, namely the fluctuating part of velocity
    !> uu_1d1: u'v', and so on. col2 and col3 are its indexes.
    !real(wp), allocatable, dimension(:,:,:) :: u_1d1
    !real(wp), allocatable, dimension(:,:,:) :: um_1d1
    !real(wp), allocatable, dimension(:,:,:) :: up_1d1
    !real(wp), allocatable, dimension(:,:,:,:) :: uu_1d1

    allocate(u_1d1(nz_global, 3, post_n_1d1), um_1d1(nz_global, 3, post_n_1d1))
    !allocate(up_1d1(nz_global, 3, post_n_1d1), uu_1d1(nz_global, 3, 3, post_n_1d1))
    
    allocate(u_2d1(ny_global, nz_global, 3, post_n_2d1))
    allocate(um_2d1(ny_global, nz_global, 3, post_n_2d1))
    !allocate(up_2d1(ny_global, nz_global, 3, post_n_2d1))
    !allocate(uu_2d1(ny_global, nz_global, 3, 3, post_n_2d1))
    
    allocate(u_2d2(nx_global, ny_global, 3, post_n_2d2))
    allocate(um_2d2(nx_global, ny_global, 3, post_n_2d2))
    !allocate(up_2d2(nx_global, ny_global, 3, post_n_2d2))
    !allocate(uu_2d2(nx_global, ny_global, 3, 3, post_n_2d2))
    
    allocate(u_2d3(nx_global, nz_global, 3, post_n_2d3))
    allocate(um_2d3(nx_global, nz_global, 3, post_n_2d3))
    
    u_1d1 = 0.0_wp; um_1d1 = 0.0_wp
    u_2d1 = 0.0_wp; um_2d1 = 0.0_wp
    u_2d2 = 0.0_wp; um_2d2 = 0.0_wp
    u_2d3 = 0.0_wp; um_2d3 = 0.0_wp
    
    if (iwavy .eq. 7) then
      allocate(ui_1d1(nz_global, 3, post_n_1d1), uim_1d1(nz_global, 3, post_n_1d1))
      !allocate(up_1d1(nz_global, 3, post_n_1d1), uu_1d1(nz_global, 3, 3, post_n_1d1))
      
      allocate(ui_2d1(ny_global, nz_global, 3, post_n_2d1))
      allocate(uim_2d1(ny_global, nz_global, 3, post_n_2d1))
      !allocate(up_2d1(ny_global, nz_global, 3, post_n_2d1))
      !allocate(uu_2d1(ny_global, nz_global, 3, 3, post_n_2d1))
      
      allocate(ui_2d2(nx_global, ny_global, 3, post_n_2d2))
      allocate(uim_2d2(nx_global, ny_global, 3, post_n_2d2))
      !allocate(up_2d2(nx_global, ny_global, 3, post_n_2d2))
      !allocate(uu_2d2(nx_global, ny_global, 3, 3, post_n_2d2))
      
      allocate(ui_2d3(nx_global, nz_global, 3, post_n_2d3))
      allocate(uim_2d3(nx_global, nz_global, 3, post_n_2d3))
      !allocate(up_2d3(nx_global, nz_global, 3, post_n_2d3))
      !allocate(uu_2d3(nx_global, nz_global, 3, 3, post_n_2d3))
      
      ui_1d1 = 0.0_wp; uim_1d1 = 0.0_wp
      ui_2d1 = 0.0_wp; uim_2d1 = 0.0_wp
      ui_2d2 = 0.0_wp; uim_2d2 = 0.0_wp
      ui_2d3 = 0.0_wp; uim_2d3 = 0.0_wp
    endif
    
    
   ! n_it = 0; time = 0.0;
   ! rem_it = ini_it - nint(real(ini_it)/real(noutd))*noutd
   ! 
   ! !> run for 1st time to get mean field
   ! do ii = nint(real(ini_it)/real(noutd)), &
   !   nint(real(end_it)/real(noutd))
   !   n_it = n_it + 1
   !   it = ii * noutd + rem_it
   !   write(fileid,'(a4,I0.10)') 'DAT_',it
   !   if (myid.eq.0) print *, 'Post analysis output: ', fileid

   !   !fdat = trim(fileid//".dat")
   !   !fh5 = trim(fileid//".h5")
   !   call reread_output(it, fileid)
   !   if(myid.eq.0) print *, 'time=', time
   !   call update_ghost(u, level); call update_ghost(v, level); 
   !   call update_ghost(w, level); call update_ghost(pp, level);
   !   if (it == ini_it ) then
   !     call hos_init
   !     !if(myid.eq.0) print *, 'hos_init called'
   !   endif

   !   if (iwavy .eq. 7) then
   !     write(fileid, '(i0.8)') it 
   !     fileid = trim(adjustl(fileid))
   !     fparam_hos = trim("restart_param_hos.dat"//fileid)
   !     fdata_hos = trim("restart_hos.h5"//fileid)
   !     call read_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos,&
   !       fparam_hos, fdata_hos)
   !     !if(myid.eq.0) print *, 'hos restarted'
   !  
   !     !> variables in following lines allocated in hos_init 
   !     call wavenum(wvn_hos)
   !     call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
   !     call zeta(eta_hos,zp_hos)
   !     tmp=1.0_wp
   !     call boundvp(vps_hos,r_hos,zp_hos)
   !     call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
   !     call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
   !     !if(myid.eq.0) print *, 'hos internal computation finished'
   !     
   !     if (isbot) call bottom_hos_les(time)
   !     !if(myid.eq.0) print *, 'bottom_hos_les called. Got ubs, eta, ub, hh,&
   !     !  & ht,hx'
   !   endif

   !   if(myid==0) print *, 'HOS initialized!'

   !   ! transfer data to the upper cpus
   !   call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
   !   call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
   !   call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
   !   call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
   !   call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
   !   call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
   !   call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
   !   call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
   !   if(myid.eq.0) print *, 'HOS data transfered to upper cpus'

   !   call init_p_poisson
   !   
   !   !> nl_coef: input eta, hx, hy
   !   !>          output ex, ey, ehx, ehy, reh, her, exr, zetax,...
   !   call nl_coef
   !   
   !   !> gather eta and hh
   !   call gather_2d_xy(eta, eo)
   !   call gather_2d_xy(hh, hho)

   !   !> filtered u
   !   deltax = 2.0_wp * TWOPI/pex/nx_global
   !   deltay = 2.0_wp * TWOPI/pey/ny_global
   !   
   !   !uf = u; vf = v; wf = w
   !   !!if(myid.eq.0) print *, 'mfilt=', mfilt
   !   !call les_filter(uf(:,:,1:xsz(3)), mfilt, deltax, deltay)
   !   !call les_filter(vf(:,:,1:xsz(3)), mfilt, deltax, deltay)
   !   !call les_filter(wf(:,:,1:xsz(3)), mfilt, deltax, deltay)
   !   !call update_ghost(uf, level)
   !   !call update_ghost(vf, level)
   !   !call update_ghost(wf, level)

   !   if (iwavy.eq.7) then
   !     call comp2cart(u,ui,1); call comp2cart(v,vi,2); call comp2cart(w,wi,3)
   !     call comp2cart(pp,ppi,1)
   !   endif
   !   
   !   u_tm = u_tm + u; v_tm = v_tm + v; w_tm = w_tm + w
   !   if (iwavy .eq. 7) then
   !     ui_tm = ui_tm + ui; vi_tm = vi_tm + vi; wi_tm = wi_tm + wi
   !   endif
   ! enddo

    !> read turbine model parameters
    write(fileid,'(a4,I0.10)') 'DAT_', ini_it
    if (myid .eq. 0) then
      open(92, file=trim(fileid)//".dat")
      read(92, *) time, nx_global, ny_global, nz_global, xl, yl, zl, hbar
      close(92)
    endif
    call MPI_BCAST(time,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
    call MPI_BCAST(nx_global,1,mpi_integer,0,mpi_comm_2d_cart,ierr)
    call MPI_BCAST(ny_global,1,mpi_integer,0,mpi_comm_2d_cart,ierr)
    call MPI_BCAST(nz_global,1,mpi_integer,0,mpi_comm_2d_cart,ierr)
    call MPI_BCAST(xl,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
    call MPI_BCAST(yl,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
    call MPI_BCAST(zl,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
    call MPI_BCAST(hbar,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
    
    call read_turbine_model_param(iturbine, time)
    !> plyunote: zz(-1)
    zz(1-level) = 2.0*zz(1) - zz(1+level)
    zw(1-level) = 2.0*zw(1) - zw(1+level)
    zz(xsz(3)+level) = 2.0 * zz(xsz(3)) - zz(xsz(3)-1)
    zw(xsz(3)+level) = 2.0 * zw(xsz(3)) - zw(xsz(3)-1)

    !> run for 2nd time to calc 2nd moments
    n_it = 0; time = 0.0;
    rem_it = ini_it - nint(real(ini_it)/real(noutd))*noutd
    
    do ii = nint(real(ini_it)/real(noutd)), &
      nint(real(end_it)/real(noutd))
      n_it = n_it + 1
      it = ii * noutd + rem_it
      write(fileid,'(a4,I0.10)') 'DAT_',it
      if (myid.eq.0) print *, 'Post analysis output: ', fileid

      !fdat = trim(fileid//".dat")
      !fh5 = trim(fileid//".h5")
      call reread_output(it, fileid, time)
      if(myid.eq.0) print *, 'time=', time

      !if (rankcol.eq.0) then
      !  print *, 'rank_z=',rankrow, ', u1=',u(1,1,:)
      !endif      
      
      call update_ghost(u, level); call update_ghost(v, level); 
      call update_ghost(w, level); call update_ghost(pp, level);
      
      !if (rankcol.eq.0) then
      !  print *, 'rank_z=',rankrow, ', u2=',u(1,1,:)
      !endif      
      
      !> plyunote: don't need to recalculate zz
      !call grid_gen

      !print *, "mark, 1" 
      if (it == ini_it ) then
        call hos_init
        !if(myid.eq.0) print *, 'hos_init called'
      endif

      !print *, "mark, 2"
      if (iwavy .eq. 7 .and. isbot) then
        write(fileid, '(i0.8)') it 
        fileid = trim(adjustl(fileid))
        fparam_hos = trim("restart_param_hos.dat"//fileid)
        fdata_hos = trim("restart_hos.h5"//fileid)
        call read_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos,&
          fparam_hos, fdata_hos)
        !if(myid.eq.0) print *, 'hos restarted'
        !print *, "eta_hos=", eta_hos(1:5,1)
     
        !do i = 1, size(eta_hos, 1)
        !  do j = 1, size(eta_hos, 2)
        !    if (isnan(eta_hos(i,j))) print *, "nan in eta_hos_in_post_main:", i, j
        !  enddo
        !enddo
        
        !> variables in following lines allocated in hos_init 
        call wavenum(wvn_hos)
        call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
        call zeta(eta_hos,zp_hos)
        tmp=1.0_wp
        call boundvp(vps_hos,r_hos,zp_hos)
        call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
        call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
        !if(myid.eq.0) print *, 'hos internal computation finished'
        
        call bottom_hos_les(time)
        !if(myid.eq.0) print *, 'bottom_hos_les called. Got ubs, eta, ub, hh,&
        !  & ht,hx'
      endif

      !if(myid==0) print *, 'HOS initialized!'

      ! transfer data to the upper cpus
      !print *, "mark, 3"
      call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.1"
      call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.2"
      call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.3"
      call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.4"
      call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.5"
      call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.6"
      call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.7"
      call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.8"
      !if(myid.eq.0) print *, 'HOS data transfered to upper cpus'

      !> we didn't calculate pressure here. Might be needed to consider spatial
      ! distribution and average.
      !call init_p_poisson
      
      !> nl_coef: input eta, hx, hy
      !>          output ex, ey, ehx, ehy, reh, her, exr, zetax,...
      call nl_coef
      
      !> gather eta and hh
      !print *, "mark, 4"
      call gather_2d_xy(eta, eo)
      call gather_2d_xy(hh, hho)

      !> check problematic data
      !print *, "mark, 5"
      do i = 1, xsz(1)
        do j = 1, xsz(2)
          inan = is_nan(eta(i,j))
          if (inan) print *, "NaN in eta"
          inan = is_nan(hh(i,j))
          if (inan) print *, "NaN in hh"
          do k = 1, xsz(3)
            inan = is_nan(u(i,j,k))
            if (inan) print *, "NaN in u"
            inan = is_nan(v(i,j,k))
            if (inan) print *, "NaN in v"
            inan = is_nan(w(i,j,k))
            if (inan) print *, "NaN in w"
          enddo
        enddo
      enddo
      do k = 1, xsz(3)
        inan = is_nan(zz(k))
        if (inan) print *, "NaN in zz"
      enddo


      !> calculate volume flux
      !print *, "mark, 6"
      call calc_vol_flux(u(:,:,1:xsz(3)), volflux)

      !> filtered u
      deltax = 2.0_wp * TWOPI/pex/nx_global
      deltay = 2.0_wp * TWOPI/pey/ny_global
      
      uf = u; vf = v; wf = w
      !if(myid.eq.0) print *, 'mfilt=', mfilt
      call les_filter(uf(:,:,1:xsz(3)), mfilt, deltax, deltay)
      call les_filter(vf(:,:,1:xsz(3)), mfilt, deltax, deltay)
      call les_filter(wf(:,:,1:xsz(3)), mfilt, deltax, deltay)
      call update_ghost(uf, level)
      call update_ghost(vf, level)
      call update_ghost(wf, level)

      !print *, "mark, 7"
      if (istop) then
        !> ignore this les_filter step of xsz(3)+1 in navier_fu.f90
      endif

      if (isbot) then
        call wall_model_v3

        call spec_theta_kx_ky(eo,eta_sp,eta_theta)
        
        tauwxo = 0.0; tauwyo = 0.0; pso = 0.0
        call gather_2d_xy(tauwx, tauwxo)
        call gather_2d_xy(tauwy, tauwyo)
        call gather_2d_xy(pp(:,:,1), pso)
        
        tau_vis = -tauwx * xl/nx * yl/ny 
        tau_pre = pp(:,:,1) * hx * xl/nx * yl/ny
        call gather_2d_xy(hx, hxo)
        call gather_2d_xy(tau_vis, tau_viso)
        call gather_2d_xy(tau_pre, tau_preo)
        tau_vis_m = sum(tau_viso)/xl/yl
        tau_pre_m = sum(tau_pre)/xl/yl
        tau_all_m = tau_vis_m + tau_pre_m
        utau_eff = sqrt(tau_all_m)
        utau_eff_tm = utau_eff_tm + utau_eff
        if (myid .eq. 0) then
          print *, 'Surface stress, viscous part=', tau_vis_m,  &
            ', pressure part=', tau_pre_m, ', sum=', tau_all_m, &
            ', effective utau=', utau_eff
        
          fout=''
          write(fout, '(a,I0.10,a4)') 'SURFACE_STRESS_',it,'.DAT'
          open(1003, file=fout)
          write(1003,*) "VARIABLES = X, Y, Z, TAUWX, TAUWY, PPS, HX, SP, THETA"
          write(1003,*) 'ZONE T="', it*dt, '" I=', nx_global,' J=',&
            ny_global, ' K=1', ' SOLUTIONTIME=', it*dt
          
          do j = 1, ny_global
            do i = 1, nx_global
              xtemp = i * (xl/nx); ytemp = j * (yl/ny)
              ztemp = zzall(1) * (hbar+hho(i, j)) - hho(i,j)    
              write(1003, '(9E18.8E3)') xtemp, ytemp, ztemp,&
                tauwxo(i,j), tauwyo(i,j), pso(i,j), hxo(i,j), &
                eta_sp(i,j), eta_theta(i,j)
            enddo
          enddo

          close(1003)
        endif
      endif

      !> Get SGS stress and its contribution to x momentum

      !navier_les.f
      !input: u, v, w
      !output: s11~s33, s11w~s33w
      s11 = 0.0; s12 = 0.0
      call get_strain(u, v, w)
      !if (isbot) then
      !  print *, "s11=",s11(1:5,1,1), ", s12=",s12(1:5,1,1)
      !endif

      !input: u, v, w, uf, vf, wf, deltax, deltay
      !output: nut
      call get_nut(u, v, w, uf, vf, wf, deltax, deltay)
      nut = nut + 1.0_wp/re
      
      !input: s11~s33, s11w~s33w, nut, nutw
      !output: t11~t33, t11w~t33w
      t11 = 0.0; t12 = 0.0
      call get_SGS_stress

      !apply wall model  
      if(isbot)then
         call pdfx(w(:,:,1), t13wx, pex)
         t13w(:,:,1)=tauwx(:,:)-nut(:,:,1)*(t13wx(:,:)+zetax(:,:,1)*(w(:,:,2)-w(:,:,1))/dz(2))
         call dealiasxy(t13w(:,:,1))

         !print *, "t13w debug: tauwx,", tauwx(1:5,1), "; nut,", nut(1:5,1,1),&
         !  "; t13wx,", t13wx(1:5,1), "; zetax,", zetax(1:5,1,1),"; w,",w(1:5,1,1)
      endif

      ! print*, t13w(1,1,1)
      !apply top boundary condition
      if(istop)then
         t11(:,:,xsz(3)+1)=t11(:,:,xsz(3)-1)
         t12(:,:,xsz(3)+1)=t12(:,:,xsz(3)-1)
         t13w(:,:,xsz(3))=-t13w(:,:,xsz(3)-2)
      endif

      tmp_x6(1:xsz(1),1:xsz(2),1-level:xsz(3)+level) = 0.0
      !if (isbot) then
      !  print *, "t11=",t11(1:5,1,1),", t12=",t12(1:5,1,1),", t13w=",t13w(1:5,1,1)
      !endif
      call div_tau(t11, t12,  t13w, tmp_x6(:,:,1:xsz(3)))
      !if (isbot) then
      !  print *, "tmp_x6:", tmp_x6(1:5,1,1)
      !endif
      call update_ghost(tmp_x6, level)

      nut_tm = nut_tm + nut
      dsgs_tm = dsgs_tm + tmp_x6
      dsgsu_tm = dsgsu_tm + tmp_x6*u

      !> calculate moments of eta
      nm_eta = 5 ! 5 at most
      if (nm_eta>5 .and. myid.eq.0) print *,'nm_eta value incorrect.'
      eta_nm(:) = 0.0
      do j = 1, ny_global
        do i = 1, nx_global
          do k = 1, nm_eta
            eta_nm(k) = eta_nm(k) + eo(i,j)**k
          enddo
        enddo
      enddo
      eta_nm(1:nm_eta) = eta_nm(1:nm_eta)/nx_global/ny_global
      if(myid.eq.0) then
        print *, 'eta_nm =', eta_nm(1:nm_eta)
      endif

      !> calculate turbine force
      !> as we didn't save the rotation angle, so we cannot use following
      !> steps to calc force
     ! if (iturbine .eq. 1) then
     !   call discloca_v3
     !   call veldisc_v2(u, level)
     !   call wind_turbine_force(wtforce, level, fturbinex)      
     ! else if(iturbine .eq. 3) then
     !   call rotor_model_acl(time, u, v, w, wtforce, wtforce_y, wtforce_z, &
     !     fturbinex, fturbiney, level)
     ! else if(iturbine .eq. 5) then
     !   call rotor_model_acs(time, u, v, w, wtforce, wtforce_y, wtforce_z, &
     !     fturbinex, fturbiney, level)
     ! else if(iturbine .eq. 7) then
     !   call rotor_model_admr(time, u, v, w, wtforce, wtforce_y, wtforce_z, &
     !     fturbinex, fturbiney, level)
     ! endif

     if (iturbine .ne. 0) then
       call reread_rotor_force (time, it, level)
       if (nacelle_model .eq. 1) then
         call reread_nacelle_force (it, nacfx, nacfy, nacfz, level)
         
         nacfx_tm = nacfx_tm + nacfx
         if (iwavy .eq. 7) then
           call comp2cart(nacfx, nacfxi, 1)
           nacfxi_tm = nacfxi_tm + nacfxi
         endif
       endif
     endif

      !!> nacfx has been dealiased, wtforce might be not dealiased.
      !if (nacelle_model .ne. 0) then
      !  call dealiasxy(wtforce(:,:,1:xsz(3)))
      !endif

      !> analyze discontinuity infor.
      if (iturbine .ge. 3 .and. (ids .eq. 1 .or. idsd .eq. 1)) then
        lim_dir = 1 !> indicate it is analyzing x direction
        call analyze_limit3d_from_3d (flag_de, tbn_lim_x)
        call transpose_xy(flag_de, bufy)
        lim_dir = 2 !> indicate it is analyzing y direction
        call analyze_limit3d_from_3d (bufy, tbn_lim_y)
        lim_dir = 0
      endif
      
      templocal3d(:,:,1:xsz(3)) = wtforce(:,:,1:xsz(3))
      call update_ghost(templocal3d, level)
      wtfx_tm = wtfx_tm + templocal3d
      
      if (iwavy.eq.7) then
        ! i_interp=2 not completed.
        i_interp = 1

        if (i_interp .eq. 1) then
          ! plyunote: do interpolation in each processor. a ghost cell might be not
          ! enough.
          if (myid.eq.0 .and. it.eq.ini_it) then
            print *, 'By selecting i_interp as 1, the comp2cart will run in &
              &local domains in each processor. The 1 layer ghost cell can &
              &sometimes provide enough buffering height for interpolation, &
              &while sometimes not. If you found result abnormal, &
              &please consider set the division number in z direction as 1. &
              &Or, change i_interp to 2, which will be slower.'
          endif
          call comp2cart(u,ui,1); call comp2cart(v,vi,2); call comp2cart(w,wi,3)
          call comp2cart(pp,ppi,1)
          call comp2cart(templocal3d, wtfxi, 1)
          
          call comp2cart(nut, templocal3d, 1)
          nuti_tm = nuti_tm + templocal3d
          call comp2cart(tmp_x6, templocal3d, 1)
          dsgsi_tm = dsgsi_tm + templocal3d
          call comp2cart(tmp_x6*u, templocal3d, 1)
          dsgsui_tm = dsgsui_tm + templocal3d
        else if (i_interp .eq. 2) then
          ! plyunote: do interpolation in one core
          if (myid.eq.0 .and. it.eq.ini_it) then
            print *, 'By selecting i_interp as 2 the comp2cart() will run in &
            &global mode. Every processor will gather data of entire domain. &
            &This will make sure our interpolation always in the right range. &
            &But it is really slow.'
          endif

          call gather_3d_xyz(u(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d)
          call comp2cart_global (temp3d, zzall, hho, 1)
          ui(1:xsz(1),1:xsz(2),1:xsz(3)) = &
            temp3d(xst(1):xend(1),xst(2):xend(2),xst(3):xend(3))
          call update_ghost(ui, level)
          
          call gather_3d_xyz(v(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d)
          call comp2cart_global (temp3d, zzall, hho, 2)
          vi(1:xsz(1),1:xsz(2),1:xsz(3)) = &
            temp3d(xst(1):xend(1),xst(2):xend(2),xst(3):xend(3))
          call update_ghost(vi, level)
          
          call gather_3d_xyz(w(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d)
          call comp2cart_global (temp3d, zwall, hho, 3)
          wi(1:xsz(1),1:xsz(2),1:xsz(3)) = &
            temp3d(xst(1):xend(1),xst(2):xend(2),xst(3):xend(3))
          call update_ghost(wi, level)
          
          call gather_3d_xyz(pp(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d)
          call comp2cart_global (temp3d, zzall, hho, 4)
          ppi(1:xsz(1),1:xsz(2),1:xsz(3)) = &
            temp3d(xst(1):xend(1),xst(2):xend(2),xst(3):xend(3))
          call update_ghost(ppi, level)
          
          call gather_3d_xyz(templocal3d(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d)
          call comp2cart_global (temp3d, zzall, hho, 5)
          wtfxi(1:xsz(1),1:xsz(2),1:xsz(3)) = &
            temp3d(xst(1):xend(1),xst(2):xend(2),xst(3):xend(3))
          call update_ghost(wtfxi, level)
        else
          print *, 'Error value for i_interp.'
        endif
      endif
      
      !if (rankcol.eq.0) then
      !  print *, 'rank_z=',rankrow, ', u3=',u(1,1,:)
      !endif      

      !print *, 1, (ui(1:5,1,1:5))
            
      !> interpolate w to u,v,p grid.
      do k = 1, xsz(3)
        if ((xst(3)+k-1)==1) then
          templocal3d(:,:,k) = w(:,:,k)
        else
          templocal3d(:,:,k) = 0.5*w(:,:,k)+0.5*w(:,:,k-1)
        endif
      enddo
      call update_ghost(templocal3d, level)

      !!> Vorticity
      !call gradu(u, d_u_x(:,:,:,1,1), d_u_x(:,:,:,1,2), d_u_x(:,:,:,1,3))
      !call gradu(v, d_u_x(:,:,:,2,1), d_u_x(:,:,:,2,2), d_u_x(:,:,:,2,3))
      !call gradu(templocal3d, d_u_x(:,:,:,3,1), d_u_x(:,:,:,2,2), d_u_x(:,:,:,3,3))   
      !vor(:,:,:,1) = d_u_x(:,:,:,3,2) - d_u_x(:,:,:,2,3)
      !vor(:,:,:,2) = d_u_x(:,:,:,1,3) - d_u_x(:,:,:,3,1)
      !vor(:,:,:,3) = d_u_x(:,:,:,2,1) - d_u_x(:,:,:,1,2)
      
      u_tm = u_tm + u; v_tm = v_tm + v; w_tm = w_tm + templocal3d;
      pp_tm = pp_tm + pp
      temp_de_1 = u * templocal3d
      call dealiasxy(temp_de_1(:,:,1:xsz(3)))
      uw_tm = uw_tm + temp_de_1
      
      temp_de_1 = u * u 
      call dealiasxy(temp_de_1(:,:,1:xsz(3)))
      tke_tm = tke_tm + temp_de_1 
      temp_de_1 = v * v 
      call dealiasxy(temp_de_1(:,:,1:xsz(3)))
      tke_tm = tke_tm + temp_de_1 
      temp_de_1 = templocal3d * templocal3d 
      call dealiasxy(temp_de_1(:,:,1:xsz(3)))
      tke_tm = tke_tm + temp_de_1 
      
      !wtfx_tm(:,:,1:xsz(3)) = wtfx_tm(:,:,1:xsz(3)) + wtforce(:,:,1:xsz(3))
      if (iwavy .eq. 7) then
        do k = 1, xsz(3)
          if ((xst(3)+k-1)==1) then
            templocal3d(:,:,k) = wi(:,:,k)
          else
            templocal3d(:,:,k) = 0.5*wi(:,:,k)+0.5*wi(:,:,k-1)
          endif
        enddo
        call update_ghost(templocal3d, level)
        
        ui_tm = ui_tm + ui; vi_tm = vi_tm + vi; wi_tm = wi_tm + templocal3d
        temp_de_1 = ui * templocal3d; ppi_tm = ppi_tm + ppi
        call dealiasxy(temp_de_1(:,:,1:xsz(3)))
        uwi_tm = uwi_tm + temp_de_1
        
        temp_de_1 = ui * ui 
        call dealiasxy(temp_de_1(:,:,1:xsz(3)))
        tkei_tm = tkei_tm + temp_de_1 
        temp_de_1 = vi * vi 
        call dealiasxy(temp_de_1(:,:,1:xsz(3)))
        tkei_tm = tkei_tm + temp_de_1 
        temp_de_1 = templocal3d * templocal3d 
        call dealiasxy(temp_de_1(:,:,1:xsz(3)))
        tkei_tm = tkei_tm + temp_de_1 

        wtfxi_tm = wtfxi_tm + wtfxi
      endif
      
      !up1_tm = u - u_tm; vp1_tm = v - v_tm; wp1_tm = w - w_tm 
      !> plyunote: And up1_tm initialization. And up1*wp1. Stop here. Because we have
      !!           better solution
      
      !> gather instantaneous velocity: u
      call gather_3d_xyz(u(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d,root0)
      do i = 1, post_n_1d1
        u_1d1(1:nz_global, 1, i) = temp3d(post_i_1d1(i), j_1d1, 1:nz_global)
      enddo
      do i = 1, post_n_2d1
        u_2d1(1:ny_global, 1:nz_global, 1, i) = temp3d(post_i_2d1(i), 1:ny_global, 1:nz_global)
      enddo
      do i = 1, post_n_2d2
        u_2d2(1:nx_global, 1:ny_global, 1, i) = temp3d(1:nx_global, 1:ny_global, post_i_2d2(i))
      enddo
      do i = 1, post_n_2d3
        u_2d3(1:nx_global, 1:nz_global, 1, i) = temp3d(1:nx_global, post_i_2d3(i), 1:nz_global)
      enddo

      do k = 1, nz_global
        call get_1dspectrum(temp3d(:,:,k), ek11(:,k), dwk, nk)
      end do

      !> output average top velocity
      utopm = sum(temp3d(:,:,nz))/nx/ny
      if(myid.eq.0) print *, 'UTOPM = ', utopm

      ! a small test of hh, eta, zz, zz*hbar, cartz, u, ui, v, vi, w, wi
      !do k = 1, xsz(3)
      !  print *, myid, k, xst(3)+k-1, hh(1,1), eta(1,1), &
      !    zz(k), zz(k)*hbar, zz(k) * (hbar+hh(1,1)) - hh(1,1),&
      !    u(1,1,k), ui(1,1,k), v(1,1,k), vi(1,1,k), w(1,1,k),wi(1,1,k)
      !enddo
      
      if (iwavy .eq. 7) then
        call gather_3d_xyz(ui(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d,root0)
        do i = 1, post_n_1d1
          ui_1d1(1:nz_global, 1, i) = temp3d(post_i_1d1(i), j_1d1, 1:nz_global)
        enddo
        do i = 1, post_n_2d1
          ui_2d1(1:ny_global, 1:nz_global, 1, i) = temp3d(post_i_2d1(i), 1:ny_global, 1:nz_global)
        enddo
        do i = 1, post_n_2d2
          ui_2d2(1:nx_global, 1:ny_global, 1, i) = temp3d(1:nx_global, 1:ny_global, post_i_2d2(i))
        enddo
        do i = 1, post_n_2d3
          ui_2d3(1:nx_global, 1:nz_global, 1, i) = temp3d(1:nx_global, post_i_2d3(i), 1:nz_global)
        enddo
        
        do k = 1, nz_global
          call get_1dspectrum(temp3d(:,:,k), eki11(:,k), dwk, nk)
        end do
      endif
      
      !> gather v
      call gather_3d_xyz(v(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d,root0)
      do i = 1, post_n_1d1
        u_1d1(1:nz_global, 2, i) = temp3d(post_i_1d1(i), j_1d1, 1:nz_global)
      enddo
      do i = 1, post_n_2d1
        u_2d1(1:ny_global, 1:nz_global, 2, i) = temp3d(post_i_2d1(i), 1:ny_global, 1:nz_global)
      enddo
      do i = 1, post_n_2d2
        u_2d2(1:nx_global, 1:ny_global, 2, i) = temp3d(1:nx_global, 1:ny_global, post_i_2d2(i))
      enddo
      do i = 1, post_n_2d3
        u_2d3(1:nx_global, 1:nz_global, 2, i) = temp3d(1:nx_global, post_i_2d3(i), 1:nz_global)
      enddo

      do k = 1, nz_global
        call get_1dspectrum(temp3d(:,:,k), ek22(:,k), dwk, nk)
      end do

      if (iwavy .eq. 7) then
        call gather_3d_xyz(vi(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d,root0)
        do i = 1, post_n_1d1
          ui_1d1(1:nz_global, 2, i) = temp3d(post_i_1d1(i), j_1d1, 1:nz_global)
        enddo
        do i = 1, post_n_2d1
          ui_2d1(1:ny_global, 1:nz_global, 2, i) = temp3d(post_i_2d1(i), 1:ny_global, 1:nz_global)
        enddo
        do i = 1, post_n_2d2
          ui_2d2(1:nx_global, 1:ny_global, 2, i) = temp3d(1:nx_global, 1:ny_global, post_i_2d2(i))
        enddo
        do i = 1, post_n_2d3
          ui_2d3(1:nx_global, 1:nz_global, 2, i) = temp3d(1:nx_global, post_i_2d3(i), 1:nz_global)
        enddo
        
        do k = 1, nz_global
          call get_1dspectrum(temp3d(:,:,k), eki22(:,k), dwk, nk)
        end do
      endif

      !> gather w, considering staggered grid      
      call gather_3d_xyz(w(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d,root0)
      do i = 1, post_n_1d1
        u_1d1(1, 3, i) = temp3d(post_i_1d1(i), j_1d1, 1)
        u_1d1(2:(nz_global-1), 3, i) = 0.5_wp * temp3d(post_i_1d1(i), j_1d1, 2:(nz_global-1)) &
          + 0.5_wp * temp3d(post_i_1d1(i), j_1d1, 1:(nz_global-2))
        u_1d1(nz_global, 3, i) = temp3d(post_i_1d1(i), j_1d1, nz_global-1)
        !> k=nz, the third formula can be merged into the second.
      enddo
      do i = 1, post_n_2d1
        u_2d1(1:ny_global, 1, 3, i) = temp3d(post_i_2d1(i), 1:ny_global, 1)
        u_2d1(1:ny_global, 2:(nz_global-1), 3, i) = &
          temp3d(post_i_2d1(i), 1:ny_global, 2:(nz_global-1)) * 0.5_wp &
          + 0.5_wp * temp3d(post_i_2d1(i), 1:ny_global, 1:(nz_global-2))
        u_2d1(1:ny_global, nz_global, 3, i) = temp3d(post_i_2d1(i), 1:ny_global, nz_global-1)
      enddo
      do i = 1, post_n_2d2
        if (post_i_2d2(i)==1) then
          u_2d2(1:nx_global, 1:ny_global, 3, i) = temp3d(1:nx_global, 1:ny_global, 1)
        elseif (post_i_2d2(i)>1) then
          u_2d2(1:nx_global, 1:ny_global, 3, i) = temp3d(1:nx_global, 1:ny_global, post_i_2d2(i))&
            * 0.5_wp + 0.5_wp * temp3d(1:nx_global, 1:ny_global, post_i_2d2(i)-1)
        endif
      enddo
      do i = 1, post_n_2d3
        u_2d3(1:nx_global, 1, 3, i) = temp3d(1:nx_global, post_i_2d3(i), 1)
        u_2d3(1:nx_global, 2:nz_global, 3, i) = temp3d(1:nx_global, post_i_2d3(i), 2:nz_global)&
          * 0.5_wp + 0.5_wp * temp3d(1:nx_global, post_i_2d3(i), 1:(nz_global-1))
      enddo
            
      do k = 1, nz_global
        call get_1dspectrum(temp3d(:,:,k), ek33(:,k), dwk, nk)
      end do

      if (iwavy .eq. 7) then
        call gather_3d_xyz(wi(1:xsz(1),1:xsz(2),1:xsz(3)),temp3d,root0)
        do i = 1, post_n_1d1
          ui_1d1(1, 3, i) = temp3d(post_i_1d1(i), j_1d1, 1)
          ui_1d1(2:nz_global, 3, i) = 0.5_wp * temp3d(post_i_1d1(i), j_1d1, 2:nz_global) &
            + 0.5_wp * temp3d(post_i_1d1(i), j_1d1, 1:(nz_global-1))
        enddo
        do i = 1, post_n_2d1
          ui_2d1(1:ny_global, 1, 3, i) = temp3d(post_i_2d1(i), 1:ny_global, 1)
          ui_2d1(1:ny_global, 2:nz_global, 3, i) = temp3d(post_i_2d1(i), 1:ny_global, 2:nz_global)&
            * 0.5_wp + 0.5_wp * temp3d(post_i_2d1(i), 1:ny_global, 1:(nz_global-1))
        enddo
        do i = 1, post_n_2d2
          if (post_i_2d2(i)>1) then
            ui_2d2(1:nx_global, 1:ny_global, 3, i) = temp3d(1:nx_global, 1:ny_global, post_i_2d2(i))&
              * 0.5_wp + 0.5_wp * temp3d(1:nx_global, 1:ny_global, post_i_2d2(i)-1)
          endif
        enddo
        do i = 1, post_n_2d3
          ui_2d3(1:nx_global, 1, 3, i) = temp3d(1:nx_global, post_i_2d3(i), 1)
          ui_2d3(1:nx_global, 2:nz_global, 3, i) = temp3d(1:nx_global, post_i_2d3(i), 2:nz_global)&
            * 0.5_wp + 0.5_wp * temp3d(1:nx_global, post_i_2d3(i), 1:(nz_global-1))
        enddo
        
        do k = 1, nz_global
          call get_1dspectrum(temp3d(:,:,k), eki33(:,k), dwk, nk)
        end do
      endif
 
      !> sum for calculating mean velocity
      um_1d1 = um_1d1 + u_1d1
      um_2d1 = um_2d1 + u_2d1
      um_2d2 = um_2d2 + u_2d2
      um_2d3 = um_2d3 + u_2d3

      ek11_tm = ek11_tm + ek11
      ek22_tm = ek22_tm + ek22
      ek33_tm = ek33_tm + ek33

      if (iwavy .eq. 7) then
        uim_1d1 = uim_1d1 + ui_1d1
        uim_2d1 = uim_2d1 + ui_2d1
        uim_2d2 = uim_2d2 + ui_2d2
        uim_2d3 = uim_2d3 + ui_2d3
      
        eki11_tm = eki11_tm + eki11
        eki22_tm = eki22_tm + eki22
        eki33_tm = eki33_tm + eki33
      endif

      !> save instantaneous velocity to file
      if (myid .eq. 0) then
        do i =1, post_n_1d1 
          fout=''
          write(fout,'(a,I0.4)') 'mkdir -p POST_U_1D1_', i
          call system(fout)
          
          fout=''
          write(fout, '(a,I0.4,a,I0.10,a1,I0.4,a4)') './POST_U_1D1_',i,'/POST_U_1D1_',it,'_',i,'.DAT'
          open(1003, file=fout)
          write(1003,*) "VARIABLES = X, Y, Z, U, V, W"
          write(1003,*) 'ZONE T="', it*dt, '" I=1 J=1 K=', nz_global, &
            ' SOLUTIONTIME=', it*dt
          
          xtemp = coord_1d1_g(i)
          !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
          ytemp = j_1d1 * (yl/ny)
          iii = post_i_1d1(i)
          jjj = j_1d1
          do k = 1, nz_global      
            ztemp = zzall(k) * (hbar+hho(iii, jjj)) - hho(iii,jjj)
            write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
              u_1d1(k,1,i), u_1d1(k,2,i), u_1d1(k,3,i)
          enddo
          close(1003)
        enddo
       
        if (iwavy .eq. 7) then
          do i =1, post_n_1d1 
            fout=''
            write(fout,'(a,I0.4)') 'mkdir -p POST_UI_1D1_', i
            call system(fout)
            
            fout=''
            write(fout, '(a,I0.4,a,I0.10,a1,I0.4,a4)') './POST_UI_1D1_',i,'/POST_UI_1D1_',it,'_',i,'.DAT'
            open(1003, file=fout)
            write(1003,*) "VARIABLES = X, Y, Z, U, V, W"
            write(1003,*) 'ZONE T="', it*dt, '" I=1 J=1 K=', nz_global, &
              ' SOLUTIONTIME=', it*dt
            
            xtemp = coord_1d1_g(i)
            !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
            ytemp = j_1d1 * (yl/ny)
            do k = 1, nz_global      
              ztemp = zzall(k)*hbar
              write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
                ui_1d1(k,1,i), ui_1d1(k,2,i), ui_1d1(k,3,i)
            enddo
            close(1003)
          enddo
        endif
        
        do i =1, post_n_2d1 
          fout=''
          write(fout,'(a,I0.4)') 'mkdir -p POST_U_2D1_', i
          call system(fout)
          
          fout=''
          write(fout, '(a,I0.4,a,I0.10,a1,I0.4,a4)') './POST_U_2D1_',i,'/POST_U_2D1_',it,'_',i,'.DAT'
          open(1003, file=fout)
          write(1003,*) "VARIABLES = X, Y, Z, U, V, W"
          write(1003,*) 'ZONE T="', it*dt, '" I=1 J=', ny_global,' K=', nz_global, &
            ' SOLUTIONTIME=', it*dt
          
          xtemp = coord_2d1_g(i)
          !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
          iii = post_i_2d1(i)
          do k = 1, nz_global      
            do j = 1, ny_global
              ytemp = j * (yl/ny)
              ztemp = zzall(k) * (hbar+hho(iii, j)) - hho(iii,j)    
              write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
                u_2d1(j,k,1,i), u_2d1(j,k,2,i), u_2d1(j,k,3,i)
            enddo
          enddo
          close(1003)
        enddo
        
        if (iwavy .eq. 7) then
          do i =1, post_n_2d1 
            fout=''
            write(fout,'(a,I0.4)') 'mkdir -p POST_UI_2D1_', i
            call system(fout)
            
            fout=''
            write(fout, '(a,I0.4,a,I0.10,a1,I0.4,a4)') './POST_UI_2D1_',i,'/POST_UI_2D1_',it,'_',i,'.DAT'
            open(1003, file=fout)
            write(1003,*) "VARIABLES = X, Y, Z, U, V, W"
            write(1003,*) 'ZONE T="', it*dt, '" I=1 J=', ny_global,' K=', nz_global, &
              ' SOLUTIONTIME=', it*dt
            
            xtemp = coord_2d1_g(i)
            !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
            do k = 1, nz_global      
              ztemp = zzall(k)*hbar
              do j = 1, ny_global
                ytemp = j * (yl/ny)
                write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
                  ui_2d1(j,k,1,i), ui_2d1(j,k,2,i), ui_2d1(j,k,3,i)
              enddo
            enddo
            close(1003)
          enddo
        endif
         
        do i =1, post_n_2d2 
          fout=''
          write(fout,'(a,I0.4)') 'mkdir -p POST_U_2D2_', i
          call system(fout)
          
          fout=''
          write(fout, '(a,I0.4,a,I0.10,a1,I0.4,a4)') './POST_U_2D2_',i,'/POST_U_2D2_',it,'_',i,'.DAT'
          open(1003, file=fout)
          write(1003,*) "VARIABLES = X, Y, Z, U, V, W"
          write(1003,*) 'ZONE T="', it*dt, '" I=', nx_global, ' J=', ny_global,' K=1', &
            ' SOLUTIONTIME=', it*dt
          
          !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
          do j = 1, ny_global      
            ytemp = j * (yl/ny)
            do k = 1, nx_global
              xtemp = k * (xl/nx)
              ztemp = zzall(post_i_2d2(i)) * (hbar+hho(k, j)) - hho(k,j)    
              write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
                u_2d2(k,j,1,i), u_2d2(k,j,2,i), u_2d2(k,j,3,i)
            enddo
          enddo
          close(1003)
        enddo
        
        if (iwavy .eq. 7) then
          do i =1, post_n_2d2 
            fout=''
            write(fout,'(a,I0.4)') 'mkdir -p POST_UI_2D2_', i
            call system(fout)
            
            fout=''
            write(fout, '(a,I0.4,a,I0.10,a1,I0.4,a4)') './POST_UI_2D2_',i,'/POST_UI_2D2_',it,'_',i,'.DAT'
            open(1003, file=fout)
            write(1003,*) "VARIABLES = X, Y, Z, U, V, W"
            write(1003,*) 'ZONE T="', it*dt, '" I=', nx_global, ' J=', ny_global,' K=1', &
              ' SOLUTIONTIME=', it*dt
            
            ztemp = coord_2d2_g(i)
            !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
            do j = 1, ny_global      
              ytemp = j * (yl/ny)
              do k = 1, nx_global
                xtemp = k * (xl/nx)
                write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
                  ui_2d2(k,j,1,i), ui_2d2(k,j,2,i), ui_2d2(k,j,3,i)
              enddo
            enddo
            close(1003)
          enddo
        endif

        do i =1, post_n_2d3 
          fout=''
          write(fout,'(a,I0.4)') 'mkdir -p POST_U_2D3_', i
          call system(fout)
          
          fout=''
          write(fout, '(a,I0.4,a,I0.10,a1,I0.4,a4)') './POST_U_2D3_',i,'/POST_U_2D3_',it,'_',i,'.DAT'
          open(1003, file=fout)
          write(1003,*) "VARIABLES = X, Y, Z, U, V, W"
          write(1003,*) 'ZONE T="', it*dt, '" I=', nx_global, ' J=1 K=', nz_global, &
            ' SOLUTIONTIME=', it*dt
          
          ytemp = coord_2d3_g(i)
          jjj = post_i_2d3(i)
          !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
          do k = 1, nz_global      
            do j = 1, nx_global
              xtemp = j * (xl/nx)
              ztemp = zzall(k) * (hbar+hho(j, jjj)) - hho(j, jjj)    
              write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
                u_2d3(j,k,1,i), u_2d3(j,k,2,i), u_2d3(j,k,3,i)
            enddo
          enddo
          close(1003)
        enddo
        
        if (iwavy .eq. 7) then
          do i =1, post_n_2d3 
            fout=''
            write(fout,'(a,I0.4)') 'mkdir -p POST_UI_2D3_', i
            call system(fout)
            
            fout=''
            write(fout, '(a,I0.4,a,I0.10,a1,I0.4,a4)') './POST_UI_2D3_',i,'/POST_UI_2D3_',it,'_',i,'.DAT'
            open(1003, file=fout)
            write(1003,*) "VARIABLES = X, Y, Z, U, V, W"
            write(1003,*) 'ZONE T="', it*dt, '" I=', nx_global, ' J=1 K=', nz_global, &
              ' SOLUTIONTIME=', it*dt
            
            ytemp = coord_2d3_g(i)
            !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
            do k = 1, nz_global      
              ztemp = zzall(k)*hbar
              do j = 1, nx_global
                xtemp = j * (xl/nx)
                write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
                  u_2d3(j,k,1,i), u_2d3(j,k,2,i), u_2d3(j,k,3,i)
              enddo
            enddo
            close(1003)
          enddo
        endif

      endif

      ! 1d spectrum
      ekm = 0
      do j = 1, nk
        do k = 1, nz_global / 2
          ekm(j) = ekm(j) + (ek11(j,k) + ek22(j,k) + ek33(j,k)) / (nz_global / 2)
        end do
      end do
      if (iwavy .eq. 7) then
        ekim = 0
        do j = 1, nk
          do k = 1, nz_global / 2
            ekim(j) = ekim(j) + (eki11(j,k) + eki22(j,k) + eki33(j,k)) / (nz_global / 2)
          end do
        end do
      endif

      if (myid.eq.0) then
        fout = ''
        write(fout, '(a,I0.10,a)') 'spectrum_', it, '.dat'
        open(1003, file=fout) 
        if (iwavy .eq. 7) then
          write(1003, *) "VARIABLES = X, Z, EK11, EK22, EK33, EKI11, EKI22, EKI33"
        else
          write(1003, *) "VARIABLES = X, Z, EK11, EK22, EK33"
        endif
        write(1003,*) 'ZONE T="', it*dt, '" I=', nk, ' J=1 K=', nz_global, &
          ' SOLUTIONTIME=', it*dt
        do k = 1, nz_global
          ztemp = zzall(k) * hbar
          do j = 1, nk
            xtemp = j  
            if (iwavy .eq. 7) then
              write(1003, *) xtemp, ztemp, ek11(j,k), ek22(j,k), ek33(j,k), &
                eki11(j,k), eki22(j,k), eki33(j,k)
            else
              write(1003, *) xtemp, ztemp, ek11(j,k), ek22(j,k), ek33(j,k)
            end if 
          enddo
        enddo
        close(1003)
      end if
    enddo

    utau_eff_tm = utau_eff_tm / n_it
    if(myid.eq.0) print *, 'Time averaged utau_eff = ', utau_eff_tm
    
    !> plyunote: let's do spatial and time average
    um_1d1 = um_1d1 / n_it
    um_2d1 = um_2d1 / n_it
    um_2d2 = um_2d2 / n_it
    um_2d3 = um_2d3 / n_it

    ek11_tm = ek11_tm/n_it; ek22_tm = ek22_tm/n_it; ek33_tm = ek33_tm/n_it
    if (iwavy .eq. 7) then
      eki11_tm = eki11_tm/n_it; eki22_tm = eki22_tm/n_it; eki33_tm = eki33_tm/n_it
    endif
    
    call IO_Init
    
    !> u_tm is time average of each grid (nx*ny*nz), it needn't mpi comm
    u_tm = u_tm / n_it; v_tm = v_tm / n_it; w_tm = w_tm / n_it
    pp_tm = pp_tm / n_it; uw_tm = uw_tm / n_it
    tke_tm = tke_tm / n_it
    
    temp_de_1 = u_tm ** 2
    call dealiasxy(temp_de_1(:,:,1:xsz(3)))
    tke_tm = tke_tm - temp_de_1
    temp_de_1 = v_tm ** 2
    call dealiasxy(temp_de_1(:,:,1:xsz(3)))
    tke_tm = tke_tm - temp_de_1
    temp_de_1 = w_tm ** 2
    call dealiasxy(temp_de_1(:,:,1:xsz(3)))
    tke_tm = tke_tm - temp_de_1
    !tke_tm = tke_tm - u_tm**2 - v_tm**2 - w_tm**2
    
    wtfx_tm = wtfx_tm / n_it
    nut_tm = nut_tm / n_it; dsgs_tm = dsgs_tm / n_it; dsgsu_tm = dsgsu_tm / n_it
    nacfx_tm = nacfx_tm / n_it
    if (iwavy .eq. 7) then
      ui_tm = ui_tm / n_it; vi_tm = vi_tm / n_it; wi_tm = wi_tm / n_it
      ppi_tm = ppi_tm / n_it; uwi_tm = uwi_tm / n_it
      tkei_tm = tkei_tm / n_it
      
      temp_de_1 = ui_tm ** 2
      call dealiasxy(temp_de_1(:,:,1:xsz(3)))
      tkei_tm = tkei_tm - temp_de_1
      temp_de_1 = vi_tm ** 2
      call dealiasxy(temp_de_1(:,:,1:xsz(3)))
      tkei_tm = tkei_tm - temp_de_1
      temp_de_1 = wi_tm ** 2
      call dealiasxy(temp_de_1(:,:,1:xsz(3)))
      tkei_tm = tkei_tm - temp_de_1
      !tkei_tm = tkei_tm - ui_tm**2 - vi_tm**2 - wi_tm**2
      
      wtfxi_tm = wtfxi_tm / n_it
      nuti_tm = nuti_tm / n_it; dsgsi_tm = dsgsi_tm / n_it
      dsgsui_tm = dsgsui_tm / n_it
      nacfxi_tm = nacfxi_tm / n_it
    endif

    !> dealiasing for tke and turbine_force
    call dealiasxy(tke_tm(:,:,1:xsz(3)))
    if (iwavy .eq. 7) then
      call dealiasxy(tkei_tm(:,:,1:xsz(3)))
    endif

    !> um is horizontal spatial average of u_tm, it need mpi communication
    
    !if (rankcol .eq. 0) then
    !  print *, 'rank_z=', rankrow, ', level=',level, ',size(u)=', &
    !    size(u,3),', size(u_tm)=', size(u_tm,3),', size(um_l)=',  &
    !    size(um_l), ', u_tm=',u_tm(1,1,:)
    !endif
    !call MPI_BARRIER(mpi_comm_world, ierr)

    !um = 0.0; vm = 0.0; wm = 0.0
    um_l = 0.0; !vm_l = 0.0;
    wm_l = 0.0
    upwpm_l = 0.0; wtfxm_l = 0.0
    nutm_l = 0.0; dsgsm_l = 0.0; dsgsum_l = 0.0; nacfxm_l = 0.0
    do k = 1, xsz(3)
      um_l(xst(3)+k-1) = sum(u_tm(:,:,k)) / xsz(1) / xsz(2)
      !vm_l(xst(3)+k-1) = sum(v_tm(:,:,k)) / xsz(1) / xsz(2)
      wm_l(xst(3)+k-1) = sum(w_tm(:,:,k)) / xsz(1) / xsz(2)
      upwpm_l(xst(3)+k-1) = sum(uw_tm(:,:,k)-u_tm(:,:,k)*w_tm(:,:,k)) / xsz(1) / xsz(2)
      wtfxm_l(xst(3)+k-1) = sum(wtfx_tm(:,:,k)) / xsz(1) / xsz(2)
      nutm_l(xst(3)+k-1) = sum(nut_tm(:,:,k)) / xsz(1) / xsz(2)
      dsgsm_l(xst(3)+k-1) = sum(dsgs_tm(:,:,k)) / xsz(1) / xsz(2)
      dsgsum_l(xst(3)+k-1) = sum(dsgsu_tm(:,:,k)) / xsz(1) / xsz(2)
      nacfxm_l(xst(3)+k-1) = sum(nacfx_tm(:,:,k)) / xsz(1) / xsz(2)
    enddo
    !if (rankcol .eq. 0) print *, 'rank_in_z=',rankrow,', xst=', xst(1:3),&
    !  ', xsz=', xsz(1:3), ', um=', um_l(1+xst(3)-1:xsz(3)+xst(3)-1)
    call MPI_Allreduce(um_l, um, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_ROW, ierror)
    call MPI_Allreduce(wm_l, wm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_ROW, ierror)
    call MPI_Allreduce(upwpm_l, upwpm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_ROW, ierror)
    call MPI_Allreduce(wtfxm_l, wtfxm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_ROW, ierror)
    call MPI_Allreduce(nutm_l, nutm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_ROW, ierror)
    call MPI_Allreduce(dsgsm_l, dsgsm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_ROW, ierror)
    call MPI_Allreduce(dsgsum_l, dsgsum, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_ROW, ierror)
    call MPI_Allreduce(nacfxm_l, nacfxm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_ROW, ierror)
    
    !call MPI_BARRIER(mpi_comm_world, ierr)
    !if (rankcol .eq. 0) print *, 'after MPI comm in z, rank_in_z=', &
    !  rankrow, ', um=', um(1:nz_global)
    
    !call MPI_BARRIER(mpi_comm_world, ierr)
    !if (rankrow .eq. 0) print *, 'after MPI comm in z, rank_in_y=', &
    !  rankcol, ', um at bottom domain is ', um(1+xst(3)-1:xsz(3)+xst(3)-1) 
    um_l = um; !vm_l = vm;
    wm_l = wm; upwpm_l = upwpm
    wtfxm_l = wtfxm
    nutm_l = nutm; dsgsm_l = dsgsm; dsgsum_l = dsgsum; nacfxm_l = nacfxm
    call MPI_Allreduce(um_l, um, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_COL, ierror)
    call MPI_Allreduce(wm_l, wm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_COL, ierror)
    call MPI_Allreduce(upwpm_l, upwpm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_COL, ierror)
    call MPI_Allreduce(wtfxm_l, wtfxm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_COL, ierror)
    call MPI_Allreduce(nutm_l, nutm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_COL, ierror)
    call MPI_Allreduce(dsgsm_l, dsgsm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_COL, ierror)
    call MPI_Allreduce(dsgsum_l, dsgsum, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_COL, ierror)
    call MPI_Allreduce(nacfxm_l, nacfxm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_COL, ierror)
    um = um / sizecol; wm = wm / sizecol; upwpm = upwpm / sizecol
    wtfxm = wtfxm / sizecol
    nutm = nutm / sizecol; dsgsm = dsgsm / sizecol; dsgsum = dsgsum / sizecol
    nacfxm = nacfxm / sizecol
    !if (rankrow .eq. 0) print *, 'after MPI comm in y, rank_in_y=', &
    !  rankcol, ', um at bottom domain is ', um(1+xst(3)-1:xsz(3)+xst(3)-1) 

    uppwppm_l(:) = 0.0; uppwppm(:) = 0.0
    do k = 1, xsz(3)
      uppwppm_l(xst(3)+k-1) = sum((u_tm(:,:,k) - &
        RESHAPE((/(um(xst(3)+k-1),i=1,xsz(1)*xsz(2))/),(/xsz(1),xsz(2)/)))&
        *(w_tm(:,:,k)-&
        RESHAPE((/(wm(xst(3)+k-1),i=1,xsz(1)*xsz(2))/),(/xsz(1),xsz(2)/)))&
        ) / xsz(1) / xsz(2)
    enddo
    call MPI_Allreduce(uppwppm_l, uppwppm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_ROW, ierror)
    uppwppm_l = uppwppm
    call MPI_Allreduce(uppwppm_l, uppwppm, nz_global, mpi_double_precision, &
      MPI_SUM, MPI_COMM_2D_COL, ierror)
    uppwppm = uppwppm / sizecol

    do i = 1, xsz(1)
      do j = 1, xsz(2)
        do k = 1, xsz(3)
          !templocal3d(i,j,k) = zz(k)*(hbar+hh(i,j))-hh(i,j)
          templocal3d(i,j,k) = zz(k) * hbar
        enddo
      enddo
    enddo

    !print *, "dsgsm:", dsgsm(1:5)
    
    call calc_budget(zzall*hbar, bforce+fturbinex, um, wm, upwpm, uppwppm, wtfxm, &
      taubm, mbm, mkebm, dsgsm, dsgsum, nutm, nacfxm, "mean_field_1d.dat")
    
    if(myid.eq.0) print *, "Writing output to mean_field_3d.h5"
    call writerOpen("mean_field_3d.h5", mpi_comm_2d_cart, writer)
    call write3d_xyz(writer, "z", templocal3d(:,:,1:xsz(3)))
    call write3d_xyz(writer, "u", u_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "v", v_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "w", w_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "pp", pp_tm(:,:,1:xsz(3)))
    
    temp_de_1 = uw_tm - u_tm * w_tm
    call dealiasxy(temp_de_1(:,:,1:xsz(3)))
    call write3d_xyz(writer, "uw", temp_de_1(:,:,1:xsz(3)))

    call write3d_xyz(writer, "tke", tke_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "wtfx", wtfx_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "nut", nut_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "dsgs", dsgs_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "dsgsu", dsgsu_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "nacfx", nacfx_tm(:,:,1:xsz(3)))
    call writerClose(writer)


    if (iwavy.eq.7) then
      !> repeat the code for UI 
      !if (rankcol .eq. 0) then
      !  print *, 'rank_z=', rankrow, ', level=',level, ',size(u)=', &
      !    size(u,3),', size(u_tm)=', size(u_tm,3),', size(um_l)=',  &
      !    size(um_l), ', u_tm=',u_tm(1,1,:)
      !endif
      !call MPI_BARRIER(mpi_comm_world, ierr)

      !um = 0.0; vm = 0.0; wm = 0.0
      um_l = 0.0; !vm_l = 0.0; 
      wm_l = 0.0
      upwpm_l = 0.0; wtfxm_l = 0.0
      nutm_l = 0.0; dsgsm_l = 0.0; dsgsum_l = 0.0; nacfxm_l = 0.0
      do k = 1, xsz(3)
        um_l(xst(3)+k-1) = sum(ui_tm(:,:,k)) / xsz(1) / xsz(2)
        !vm_l(xst(3)+k-1) = sum(vi_tm(:,:,k)) / xsz(1) / xsz(2)
        wm_l(xst(3)+k-1) = sum(wi_tm(:,:,k)) / xsz(1) / xsz(2)
        upwpm_l(xst(3)+k-1) = sum(uwi_tm(:,:,k)-ui_tm(:,:,k)*wi_tm(:,:,k)) / xsz(1) / xsz(2)
        wtfxm_l(xst(3)+k-1) = sum(wtfxi_tm(:,:,k)) / xsz(1) / xsz(2)
        nutm_l(xst(3)+k-1) = sum(nuti_tm(:,:,k)) / xsz(1) / xsz(2)
        dsgsm_l(xst(3)+k-1) = sum(dsgsi_tm(:,:,k)) / xsz(1) / xsz(2)
        dsgsum_l(xst(3)+k-1) = sum(dsgsui_tm(:,:,k)) / xsz(1) / xsz(2)
        nacfxm_l(xst(3)+k-1) = sum(nacfxi_tm(:,:,k)) / xsz(1) / xsz(2)
      enddo
      !if (rankcol .eq. 0) print *, 'rank_in_z=',rankrow,', xst=', xst(1:3),&
      !  ', xsz=', xsz(1:3), ', um=', um_l(1+xst(3)-1:xsz(3)+xst(3)-1)
      call MPI_Allreduce(um_l, uim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_ROW, ierror)
      call MPI_Allreduce(wm_l, wim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_ROW, ierror)
      call MPI_Allreduce(upwpm_l, upwpim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_ROW, ierror)
      call MPI_Allreduce(wtfxm_l, wtfxim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_ROW, ierror)
      call MPI_Allreduce(nutm_l, nutim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_ROW, ierror)
      call MPI_Allreduce(dsgsm_l, dsgsim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_ROW, ierror)
      call MPI_Allreduce(dsgsum_l, dsgsuim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_ROW, ierror)
      call MPI_Allreduce(nacfxm_l, nacfxim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_ROW, ierror)
      
      !call MPI_BARRIER(mpi_comm_world, ierr)
      !if (rankcol .eq. 0) print *, 'after MPI comm in z, rank_in_z=', &
      !  rankrow, ', um=', um(1:nz_global)
      
     ! call MPI_BARRIER(mpi_comm_world, ierr)
      !if (rankrow .eq. 0) print *, 'after MPI comm in z, rank_in_y=', &
      !  rankcol, ', um at bottom domain is ', um(1+xst(3)-1:xsz(3)+xst(3)-1) 
      um_l = uim; !vm_l = vim;
      wm_l = wim; upwpm_l = upwpim
      wtfxm_l = wtfxim
      nutm_l = nutim; dsgsm_l = dsgsim; dsgsum_l = dsgsuim; nacfxm_l = nacfxim
      call MPI_Allreduce(um_l, uim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_COL, ierror)
      call MPI_Allreduce(wm_l, wim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_COL, ierror)
      call MPI_Allreduce(upwpm_l, upwpim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_COL, ierror)
      call MPI_Allreduce(wtfxm_l, wtfxim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_COL, ierror)
      call MPI_Allreduce(nutm_l, nutim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_COL, ierror)
      call MPI_Allreduce(dsgsm_l, dsgsim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_COL, ierror)
      call MPI_Allreduce(dsgsum_l, dsgsuim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_COL, ierror)
      call MPI_Allreduce(nacfxm_l, nacfxim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_COL, ierror)
      uim = uim / sizecol; wim = wim / sizecol; upwpim = upwpim / sizecol
      wtfxim = wtfxim / sizecol
      nutim = nutim / sizecol; dsgsim = dsgsim / sizecol; dsgsuim = dsgsuim / sizecol
      nacfxim = nacfxim / sizecol
      !if (rankrow .eq. 0) print *, 'after MPI comm in y, rank_in_y=', &
      !  rankcol, ', um at bottom domain is ', um(1+xst(3)-1:xsz(3)+xst(3)-1) 

      uppwppm_l(:) = 0.0; uppwppim(:) = 0.0
      do k = 1, xsz(3)
        uppwppm_l(xst(3)+k-1) = sum((ui_tm(:,:,k) - &
          RESHAPE((/(uim(xst(3)+k-1),i=1,xsz(1)*xsz(2))/),(/xsz(1),xsz(2)/)))&
          *(wi_tm(:,:,k)-&
          RESHAPE((/(wim(xst(3)+k-1),i=1,xsz(1)*xsz(2))/),(/xsz(1),xsz(2)/)))&
          ) / xsz(1) / xsz(2)
      enddo
      call MPI_Allreduce(uppwppm_l, uppwppim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_ROW, ierror)
      uppwppm_l = uppwppim
      call MPI_Allreduce(uppwppm_l, uppwppim, nz_global, mpi_double_precision, &
        MPI_SUM, MPI_COMM_2D_COL, ierror)
      uppwppim = uppwppim / sizecol
    
      call calc_budget(zzall*hbar, bforce+fturbinex, uim, wim, upwpim, &
        uppwppim, wtfxim, taubim, mbim, mkebim, dsgsim, dsgsuim, &
        nutim, nacfxim, "mean_field_1di.dat")

      if (myid.eq.0) print *, "Writing output to mean_field_3di.h5"
      call writerOpen("mean_field_3di.h5", mpi_comm_2d_cart, writer)
      call write3d_xyz(writer, "z", templocal3d(:,:,1:xsz(3)))
      call write3d_xyz(writer, "u", ui_tm(:,:,1:xsz(3)))
      call write3d_xyz(writer, "v", vi_tm(:,:,1:xsz(3)))
      call write3d_xyz(writer, "w", wi_tm(:,:,1:xsz(3)))
      call write3d_xyz(writer, "pp", wi_tm(:,:,1:xsz(3)))
      
      temp_de_1 = uwi_tm - ui_tm * wi_tm
      call dealiasxy(temp_de_1(:,:,1:xsz(3)))
      call write3d_xyz(writer, "uw", temp_de_1(:,:,1:xsz(3)))
      call write3d_xyz(writer, "tke", tkei_tm(:,:,1:xsz(3)))
      call write3d_xyz(writer, "wtfx", wtfxi_tm(:,:,1:xsz(3)))
      call write3d_xyz(writer, "nut", nuti_tm(:,:,1:xsz(3)))
      call write3d_xyz(writer, "dsgs", dsgsi_tm(:,:,1:xsz(3)))
      call write3d_xyz(writer, "dsgsu", dsgsui_tm(:,:,1:xsz(3)))
      call write3d_xyz(writer, "nacfx", nacfxi_tm(:,:,1:xsz(3)))
      
      call writerClose(writer)

    endif

    call IO_Finalize

    !> save mean velocity to file
    if (myid .eq. 0) then
      !call system('mkdir -p POST_U_1D1')
      do i =1, post_n_1d1 
        write(fout, '(a,I0.4,a4)') 'POST_UM_1D1_',i,'.DAT'
        open(1003, file=fout)
        write(1003,*) "VARIABLES = X, Y, Z, UM, VM, WM"
        write(1003,*) 'ZONE T="', it*dt, '" I=1 J=1 K=', nz_global, &
          ' SOLUTIONTIME=', it*dt
        
        xtemp = coord_1d1_g(i)
        !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
        ytemp = j_1d1 * (yl/ny)
        do k = 1, nz_global      
          ztemp = zzall(k) * hbar
          write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
            um_1d1(k,1,i), um_1d1(k,2,i), um_1d1(k,3,i)
        enddo
        close(1003)
      enddo
      
      !call system('mkdir -p POST_U_2D1')
      do i =1, post_n_2d1 
        write(fout, '(a,I0.4,a4)') 'POST_UM_2D1_',i,'.DAT'
        open(1003, file=fout)
        write(1003,*) "VARIABLES = X, Y, Z, UM, VM, WM"
        write(1003,*) 'ZONE T="', it*dt, '" I=1 J=', ny_global,' K=', nz_global, &
          ' SOLUTIONTIME=', it*dt
        
        xtemp = coord_2d1_g(i)
        !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
        do k = 1, nz_global      
          ztemp = zzall(k) * hbar
          do j = 1, ny_global
            ytemp = j * (yl/ny)
            write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
              um_2d1(j,k,1,i), um_2d1(j,k,2,i), um_2d1(j,k,3,i)
          enddo
        enddo
        close(1003)
      enddo
      
      !call system('mkdir -p POST_U_2D2')
      do i =1, post_n_2d2 
        write(fout, '(a,I0.4,a4)') 'POST_UM_2D2_',i,'.DAT'
        open(1003, file=fout)
        write(1003,*) "VARIABLES = X, Y, Z, UM, VM, WM"
        write(1003,*) 'ZONE T="', it*dt, '" I=', nx_global, ' J=', ny_global,' K=1', &
          ' SOLUTIONTIME=', it*dt
        
        ztemp = coord_2d2_g(i)
        !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
        do j = 1, ny_global      
          ytemp = j * (yl/ny)
          do k = 1, nx_global
            xtemp = k * (xl/nx)
            write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
              um_2d2(k,j,1,i), um_2d2(k,j,2,i), um_2d2(k,j,3,i)
          enddo
        enddo
        close(1003)
      enddo

      !call system('mkdir -p POST_U_2D3')
      do i =1, post_n_2d3 
        write(fout, '(a,I0.4,a4)') 'POST_UM_2D3_',i,'.DAT'
        open(1003, file=fout)
        write(1003,*) "VARIABLES = X, Y, Z, UM, VM, WM"
        write(1003,*) 'ZONE T="', it*dt, '" I=', nx_global, ' J=1 K=', nz_global, &
          ' SOLUTIONTIME=', it*dt
        
        ytemp = coord_2d3_g(i)
        !print *, 'plyudebug, post_i,', ii, i, tb_post_ix(i), xl/nx, xtemp
        do k = 1, nz_global      
          ztemp = zzall(k) * hbar
          do j = 1, nx_global
            xtemp = j * (xl/nx)
            write(1003, '(6E18.8E3)') xtemp, ytemp, ztemp,&
              um_2d3(j,k,1,i), um_2d3(j,k,2,i), um_2d3(j,k,3,i)
          enddo
        enddo
        close(1003)
      enddo

      fout = ''
      write(fout, '(a)') 'spectrum_mean.dat'
      open(1003, file=fout) 
      if (iwavy .eq. 7) then
        write(1003, *) "VARIABLES = X, Z, EK11, EK22, EK33, EKI11, EKI22, EKI33"
      else
        write(1003, *) "VARIABLES = X, Z, EK11, EK22, EK33"
      endif
      write(1003,*) 'ZONE T="', it*dt, '" I=', nk, ' J=1 K=', nz_global, &
        ' SOLUTIONTIME=', it*dt
      do k = 1, nz_global
        ztemp = zzall(k) * hbar
        do j = 1, nk
          xtemp = j 
          if (iwavy .eq. 7) then
            write(1003, *) xtemp, ztemp, ek11_tm(j,k), ek22_tm(j,k), ek33_tm(j,k), &
              eki11_tm(j,k), eki22_tm(j,k), eki33_tm(j,k)
          else
            write(1003, *) xtemp, ztemp, ek11_tm(j,k), ek22_tm(j,k), ek33_tm(j,k)
          end if 
        enddo
      enddo
      close(1003)

    endif
       
    deallocate(tmp)
    deallocate(um, vm, wm, um_l, vm_l, wm_l)
    deallocate(upwpm_l, uppwppm_l, upwpm, uppwppm)
    deallocate(wtfxm, wtfxm_l)
    if (iwavy .eq. 7) then
      deallocate(uim, vim, wim, upwpim, uppwppim)
      deallocate(wtfxim)
    endif

    deallocate(tempz, zzall, zwall)
    deallocate(eo, hho) 
    deallocate(tauwxo, tauwyo, pso)
    
    deallocate(t13wx)
    deallocate(dsgs_tm, dsgsu_tm, nut_tm)
    deallocate(dsgsm, dsgsum, nutm)
    deallocate(dsgsm_l, dsgsum_l, nutm_l)

    if (iwavy .eq. 7) then
      deallocate(dsgsi_tm, dsgsui_tm, nuti_tm)
      deallocate(dsgsim, dsgsuim, nutim)
    endif
    
    deallocate(temp3d)
    deallocate(u_1d1, um_1d1)
    deallocate(u_2d1)
    deallocate(um_2d1)
    deallocate(u_2d2)
    deallocate(um_2d2)
    deallocate(u_2d3)
    deallocate(um_2d3)
    
    if (iwavy .eq. 7) then
      deallocate(ui, vi, wi, ppi)
      deallocate(ui_1d1, uim_1d1)
      deallocate(ui_2d1)
      deallocate(uim_2d1)
      deallocate(ui_2d2)
      deallocate(uim_2d2)
      deallocate(ui_2d3)
      deallocate(uim_2d3)
    endif    

    deallocate(u_tm, v_tm, w_tm, pp_tm, uw_tm, tke_tm)
    deallocate(templocal3d, temp_de_1)
    deallocate(wtfx_tm, nacfx_tm)
    deallocate(taubm, mbm, mkebm)
    deallocate(nacfxm, nacfxm_l)
    if (iwavy .eq. 7) then
      deallocate(ui_tm, vi_tm, wi_tm, ppi_tm, uwi_tm, tkei_tm)
      deallocate(wtfxi_tm, nacfxi_tm)
      deallocate(taubim, mbim, mkebim)
      deallocate(nacfxim)
    endif

    deallocate(ek11, ek22, ek33, ek11_tm, ek22_tm, ek33_tm)
    deallocate(ekm)
    if (iwavy .eq. 7) then
      deallocate(eki11, eki22, eki33, eki11_tm, eki22_tm, eki33_tm)
      deallocate(ekim)
    endif

    deallocate(post_i_1d1)
    deallocate(coord_1d1, coord_1d1_g)
    deallocate(post_i_2d1)
    deallocate(coord_2d1, coord_2d1_g)
    deallocate(post_i_2d2)
    deallocate(coord_2d2, coord_2d2_g)
    deallocate(post_i_2d3)
    deallocate(coord_2d3, coord_2d3_g)    
    
    !ini_it = 600
    !do it = ini_it + 1, ini_it + ntime
    !   write(fileid,*) it
    !   fileid = trim(adjustl(fileid))
    !   
    !   fparam_les = trim("restart_param.dat"//fileid)
    !   fdata_les = trim("restart.h5"//fileid)
    !   faux =  trim("restart_aux.h5"//fileid)
    !   fgrid = trim("grid.h5"//fileid)
    !   
    !   fparam_hos = trim("restart_param_hos.dat"//fileid)
    !   fdata_hos = trim("restart_hos.h5"//fileid)
    !   
    !   call reread(time, ioutd, ioutc, fgrid, fparam_les, fdata_les, faux)
    !   call grid_gen
    !   
    !   !need to add hos initialization subroutines below
    !   if(myid==0) print *, 'HOS initialization started!'
    !   
    !   if (it == ini_it + 1) then
    !      call hos_init
    !   end if
    !   call read_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos, fparam_hos, fdata_hos)
    !   
    !   call wavenum(wvn_hos)
    !   call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
    !   call zeta(eta_hos,zp_hos)
    !   tmp=1.0_wp
    !   call boundvp(vps_hos,r_hos,zp_hos)
    !   call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
    !   call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
    !   
    !   call bottom_hos_les(time)
    !   
    !   if(myid==0) print *, 'HOS initialized!'
    !   
    !   ! transfer data to the upper cpus
    !   call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
    !   call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
    !   call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
    !   call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
    !   call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
    !   call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
    !   call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
    !   call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
    !   
    !   eta=hh+0
    !   
    !   call init_p_poisson
    !   call nl_coef
    !   
    !   !test starts here!
    !   if (isbot) then
    !      tmp = 0
    !      tmp(xst(1):xend(1), xst(2):xend(2)) = eta(:,:)
    !      call MPI_Reduce(tmp, etaall, nx_global*ny_global, &
    !         mpi_double_precision, MPI_SUM, 0, MPI_COMM_2D_COL, ierror)
    !      epskx = 0
    !      call fft_for_x_hos(etaall,tmp,nx_global,ny_global)
    !      do i = 3, nx_global, 2
    !         do j = 1, ny_global
    !            rek = sqrt(((i-1.0)/2*pex)**2+((j-1.0)/2*pey)**2)
    !            ftn = tmp(i,j)**2 + tmp(i+1,j)**2
    !            epskx((i-1)/2) = epskx((i-1)/2) + 2*ftn*rek**2/ny_global
    !         end do
    !      end do
    !      do i = 1, nx_global
    !         epskx(i) = sqrt(2*epskx(i))
    !      end do
    !   end if
    !   
    !   call form_drag_k(fpk, betak, ampfd, fpt)
    !   
    !   root0 = 0
    !   
    !   call gather_2d_xy(hh(1:xsz(1),1:xsz(2)),etaall,root0)
    !   
    !   !     call turb_analysis
    !   if (myid == 0) then
    !      write(66,*) " VARIABLES = wvn,c,wage1, wage2, wage3,betak, beta_plg, beta_psm,beta_m" &
    !           // ",gamma, gam_don99, gam_don06, ka"
    !      write(66,"(A,E16.4,A,I5)") ' ZONE T="',time,'" I=', nx_global/2/3*2-1
    !      write(66,"(A,E16.4,A)") ' AUXDATA time="',time,'"'
    !      write(66,"(A,E16.4,A)") ' AUXDATA clen="',clen,'"'
    !      write(66,"(A,E16.4,A)") ' AUXDATA ctime="',ctime,'"'
    !      write(66,"(A,E16.4,A)") ' AUXDATA form_drag="',fpt,'"'
    !      
    !      !        print *,"usbot=",usbot,"z0=",z0, "fr2=", fr2
    !      do j = 1, nx_global/2/3*2-1
    !         call get_beta_miles(j*pex,beta_miles)
    !         c_phase = sqrt(1/fr2/j/pex)
    !         wage1 = c_phase / usbot
    !         u_lambda = 2.5*usbot*log(twopi/20/pex/z0)
    !         u_half_lambda = 2.5*usbot*log(pi/20/pex/z0)
    !         
    !         wage2 = c_phase / u_lambda
    !         wage3 = c_phase / u_half_lambda
    !         write(66,'(25e12.4)') j * pex, c_phase, wage1, wage2, wage3, betak(j), 16.0, 48.0,beta_miles &
    !              , betak(j) / wage1**2, 0.28*(1/wage3-1)*abs(1/wage3-1), 0.17*(1/wage3-1)*abs(1/wage3-1),&
    !              epskx(j)
    !      end do
!   !       write(67,'(25e12.4)') it*1.0, fpt
    !      print *,"finish step ",it
    !   end if
    !   
    !   !test ends here!
    !   
    !   call gather_3d_xyz(u(1:xsz(1),1:xsz(2),1:xsz(3)),velall,root0)
    !   do k = 1, nz_global
    !      call get_1dspectrum(velall(:,:,k),ek11(:,k),dwk,nk)
    !   end do
    !   uupall = velall

    !   call gather_3d_xyz(v(1:xsz(1),1:xsz(2),1:xsz(3)),velall,root0)
    !   do k = 1, nz_global
    !      call get_1dspectrum(velall(:,:,k),ek22(:,k),dwk,nk)
    !   end do

    !   call gather_3d_xyz(w(1:xsz(1),1:xsz(2),1:xsz(3)),velall,root0)
    !   do k = 1, nz_global
    !      call get_1dspectrum(velall(:,:,k),ek33(:,k),dwk,nk)
    !   end do

    !   call get_mean_flux(um,vm,wm,uup,vvp,wwp,zzall,uupall,wwpall)
    !   call get_quadrant(q1,q2,q3,q4,uupall,wwpall)

    !   do k = 1, nz_global
    !      uw_shear(k) = sum(uupall(:,:,k)*wwpall(:,:,k)) / nx_global / ny_global
    !   end do

    !   if (myid == 0) then

    !      write(70,*) " VARIABLES = z,um,vm,wm,uf,vf,wf,q1,q2,q3,q4,uw"
    !      write(70,"(A,E16.4,A,I5)") ' ZONE T="',time,'" I=', nz_global
    !      write(70,"(A,E16.4,A)") ' AUXDATA time="',time,'"'
    !      write(70,"(A,E16.4,A)") ' AUXDATA ustar="',usbot,'"'
    !      write(70,"(A,E16.4,A)") ' AUXDATA hbar="',hbar,'"'
    !      write(70,"(A,E16.4,A)") ' AUXDATA z0="',z0,'"'

    !      do k = 1, nz_global
    !         write(70,'(25e12.4)')  zzall(k)*hbar, um(k), vm(k), wm(k), uup(k), vvp(k), wwp(k), q1(k), q2(k), q3(k), q4(k), uw_shear(k)
    !      end do

    !      write(72,*) " VARIABLES = uf_nw,wf_nw,uf_out,wf_out"
    !      write(72,"(A,E16.4,A,I5)") ' ZONE T="',time,'" I=', nx_global * ny_global
    !      write(72,"(A,E16.4,A)") ' AUXDATA time="',time,'"'
    !      write(72,"(A,E16.4,A)") ' AUXDATA ustar="',usbot,'"'

    !      k = 80
    !      do j = 1, ny_global
    !         do i = 1, nx_global
    !            write(72,'(25e12.4)') uupall(i,j,5), wwpall(i,j,5),uupall(i,j,k), wwpall(i,j,k)
    !         end do
    !      end do
    !         

    !      do k = 1, nz_global / 2
    !         dissp(k) = 0
    !         do j = 1, nk
    !            dissp(k) = dissp(k) + (2*usbot/resbot)*(j * dwk)**2 &
    !                 * (ek11(j,k)+ek22(j,k)+ek33(j,k)) * dwk
    !         end do

    !         l_kom(k) = ((usbot/resbot)**3/dissp(k))**0.25
    !      end do

    !      ekm = 0
    !      do j = 1, nk
    !         do k = 1, nz_global / 2
    !            ekm(j) = ekm(j) + (ek11(j,k) + ek22(j,k) + ek33(j,k)) / (nz_global / 2)
    !         end do
    !      end do
    !      disspm = 0
    !      do j = 1, nk
    !         disspm = disspm + (2*usbot/resbot)*(j * dwk)**2 &
    !              * ekm(j) * dwk
    !      end do
    !      l_kom_ave = ((usbot/resbot)**3/disspm)**0.25

    !      write(68,*) " VARIABLES = k,ek11_nw,ek22_nw,ek33_nw,ek11_out,ek22_out,ek33_out, ekm"
    !      write(68,"(A,E16.4,A,I5)") ' ZONE T="',time,'" I=', nk
    !      write(68,"(A,E16.4,A)") ' AUXDATA time="',time,'"'
    !      write(68,"(A,E16.4,A)") ' AUXDATA disp_nw="',dissp(5),'"'
    !      write(68,"(A,E16.4,A)") ' AUXDATA l_kom_nw="',l_kom(5),'"'
    !      write(68,"(A,E16.4,A)") ' AUXDATA disp_out="',dissp(80),'"'
    !      write(68,"(A,E16.4,A)") ' AUXDATA l_kom_out="',l_kom(80),'"'
    !      write(68,"(A,E16.4,A)") ' AUXDATA disp="',disspm,'"'
    !      write(68,"(A,E16.4,A)") ' AUXDATA l_kom="',l_kom_ave,'"'
    !      write(68,"(A,E16.4,A)") ' AUXDATA nu="',usbot/resbot,'"'

    !      k = 80
    !      do j = 1, nk
    !         write(68,'(25e12.4)') j * dwk, ek11(j,5), ek22(j,5), ek33(j,5), ek11(j,k), &
    !              ek22(j,k), ek33(j,k), ekm(j)
    !      end do
    !   end if
    !   
    !   ! kfilt = 0
    !end do

    !if (myid == 0) then
    !   close(66)
    !   close(68)
    !end if

    !
    !deallocate(fpk,betak,ampfd, epskx)
    !
    !deallocate(etaall,hhall,tmp)
    !
    !deallocate(uupall,wwpall,velall)

    !deallocate(uxy, vxy, wxy)

    !deallocate(ek11, ek22, ek33, dissp, ekm, l_kom)

    !deallocate(q1,q2,q3,q4)

    !deallocate(um,vm,wm,uup,vvp,wwp)

    call fft_finalize
    call decomp_finalize
    !HOS
    call fft_finalize_hos
   
    if (myid.eq.0) print *, "Post process end successfully." 
    call mpi_finalize(ierror)
  end subroutine post_les_hos_turbine

  subroutine reread_mean(filename)
    use hdf_io
    use iso_fortran_env, only : INT64
    use param
    use navier
    use decomp, only : myid, xsz, MPI_COMM_2D_CART
    use mpi
    implicit none

    character (len=*), intent(in) :: filename      

    type(HDFObj) :: reader

    call IO_Init

    call readerOpen(trim(filename), MPI_COMM_2D_CART, reader)
    call read3D_xyz(reader, "u", u(:,:,1:xsz(3)))
    call read3D_xyz(reader, "v", v(:,:,1:xsz(3)))
    call read3D_xyz(reader, "w", w(:,:,1:xsz(3)))
    call read3D_xyz(reader, "pp", pp(:,:,1:xsz(3)))
    call readerClose(reader)

    call IO_Finalize
  end subroutine reread_mean

  !> added by plyu
  subroutine post_tke_budget
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

    !> added by plyu
    use turbine_model
    use solver_common, only : tmp_x6
    use lib_array, only: is_nan
    use discontinuity_smooth, only: analyze_limit3d_from_3d, tbn_lim_x, &
        tbn_lim_y, lim_dir 
    
    use hdf_io
    !end
    
    implicit none
    
    integer :: ioutc, ioutd
    real(wp) :: time
    
    integer :: ierror, ierr
    
    double precision :: t1, t2
    
    !HOS variables
    integer icon
    
    integer(INT64) :: id
    integer :: i, j, k, it, root0
    
    real(wp), allocatable, dimension(:) :: fpk, betak, ampfd, epskx
    
    real(wp) beta_miles
    
    character (len=100) :: fileid, fgrid, fparam_les, fdata_les, faux, fparam_hos, fdata_hos

    real(wp) c_phase, u_lambda, u_half_lambda, wage1, wage2, wage3
    real(wp) clen, cspeed, ctime, x, y, fpt

    real(wp), allocatable, dimension(:,:) :: etaall, hhall, tmp
    real(wp) ftn, rek, dwk, disspm, l_kom_ave

    real(wp), allocatable, dimension(:,:,:) :: uupall, wwpall, velall, ppall, etasig, pp_wave
    real(wp), allocatable, dimension(:,:,:,:) :: ppsig
    real(wp), allocatable, dimension(:,:) :: dtheta, uxy, vxy, wxy
    real(wp), allocatable, dimension(:,:) :: ek11, ek22, ek33
    real(wp), allocatable, dimension(:,:) :: ek11_tm, ek22_tm, ek33_tm
    real(wp), allocatable, dimension(:,:) :: eki11, eki22, eki33
    real(wp), allocatable, dimension(:,:) :: eki11_tm, eki22_tm, eki33_tm
    real(wp), allocatable, dimension(:) :: dissp, ekm, l_kom, ekim

    real(wp), allocatable, dimension(:) :: q1, q2, q3, q4, um, vm, wm, uup, vvp, wwp, zzall, uw_shear
   
    !> variables added by plyu
    integer :: ii
    real(wp), allocatable, dimension(:,:,:) :: h5_3dbuffer
    
    integer :: idp ! counter for debug
    integer :: ini_it, end_it, nk, n_it, rem_it
    integer :: rankrow, rankcol, sizerow, sizecol
    real(wp) :: deltax, deltay
    
    real(wp), allocatable, dimension(:,:,:) :: templocal3d 

    real(wp) :: volflux(3)

    type(HDFObj) :: writer
    integer :: ih5s

    !> local variables for turbine model
    integer :: ibi, n_elmt_1, nb, n1e, n2e, n3e
    real(wp), dimension(3) :: dr_vec, ds12, ds13, tmparray
    real(wp) :: dr
    
    real(wp) :: tb_xc, tb_yc, tb_zc, tb_d, xtemp, ytemp, ztemp
    real(wp), allocatable, dimension(:) :: tempz, zwall
    real(wp), allocatable, dimension(:,:,:) :: temp3d
    real(wp) :: utopm
    
    !> comment this line when turbine_model is used
    !real(wp), allocatable, dimension(:,:) :: eo, hho, hxo, hyo
    !real(wp), allocatable, dimension(:,:) :: hxo 

    !> variable for SGS
    real(wp), allocatable, dimension(:,:) :: t13wx

    !> variable for nacelle model force
    real(wp), allocatable, dimension(:,:,:) :: nacfx, nacfy, nacfz

    !>vortex dynamics: it is implemented in another  individual subroutine
    !real(wp), allocatable, dimension(:,:,:,:,:) :: d_u_x, d_u_x_tm, vorj_duidxj_tm
    !real(wp), allocatable, dimension(:,:,:,:) :: vor, vor_tm
    !real(wp), allocatable, dimension(:,:,:) :: q, q_tm

    !> t.k.e budget
    !! comprehend variable name sijp2_tm as temporal mean of s_{ij}^\prime^2
    real(wp), allocatable, dimension(:,:,:,:) :: u_tm, temp_4d, uip_p_tm, uip_ujp_ujp_tm, ujp_sijp_tm
    real(wp), allocatable, dimension(:,:,:) :: pp_tm, prod_tm, dissp_tm, tke_tm
    real(wp), allocatable, dimension(:,:,:) :: tran_turb_tm, tran_pres_tm, tran_visc_tm, adv_tm
    real(wp), allocatable, dimension(:,:,:,:,:) :: d_u_x, sijp2_tm
    real(wp), allocatable, dimension(:,:,:,:,:) :: uip_ujp_tm
    real(wp), allocatable, dimension(:,:,:) :: de3d

    logical :: inan

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
    
    call grid_gen
     
    icon = 16

    if (myid .eq. 0) then
      print *, "********* This is TKE post-processing script *******"
    endif

    !> plyunote: zz(-1)
    zz(1-level) = 2.0*zz(1) - zz(1+level)
    zw(1-level) = 2.0*zw(1) - zw(1+level)
    zz(xsz(3)+level) = 2.0 * zz(xsz(3)) - zz(xsz(3)-1)
    zw(xsz(3)+level) = 2.0 * zw(xsz(3)) - zw(xsz(3)-1)
    allocate(tmp(nx_global,ny_global))

    level = 1
    call mpi_comm_rank(mpi_comm_2d_row, rankrow, ierr) ! rank in z
    call mpi_comm_rank(mpi_comm_2d_col, rankcol, ierr) ! rank in y
    call mpi_comm_size(mpi_comm_2d_row, sizerow, ierr)
    call mpi_comm_size(mpi_comm_2d_col, sizecol, ierr)
    if(myid.eq.0) print *, 'sizerow=',sizerow,', sizecol=',sizecol
    
    ! Gather zz coordinates
    allocate(tempz(nz_global), zzall(nz_global), zwall(nz_global))
    tempz = 0.0; zzall = 0.0; zwall = 0.0 
    
    tempz(xst(3):xend(3)) = zz(1:xsz(3))
    call MPI_Allreduce(tempz, zzall, nz_global, mpi_double_precision, &
         MPI_SUM, MPI_COMM_2D_ROW, ierror)
    tempz = 0.0
    tempz(xst(3):xend(3)) = zw(1:xsz(3))
    call MPI_Allreduce(tempz, zwall, nz_global, mpi_double_precision, &
         MPI_SUM, MPI_COMM_2D_ROW, ierror)
    
    allocate(t13wx(xsz(1),xsz(2)))
    
    allocate(nacfx(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(nacfy(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(nacfz(xsz(1), xsz(2), 1-level:xsz(3)+level))
    nacfx = 0.0; nacfy = 0.0; nacfz = 0.0

    !allocate(h5_3dbuffer(xsz(1), xsz(2), xsz(3)))
    !h5_3dbuffer = 0.0_wp
    
    allocate(u_tm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    allocate(temp_4d(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    allocate(uip_p_tm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    allocate(pp_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(prod_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(dissp_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(tke_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(tran_turb_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(tran_pres_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(tran_visc_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(adv_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(d_u_x(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    allocate(sijp2_tm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    allocate(uip_ujp_tm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    allocate(uip_ujp_ujp_tm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    allocate(ujp_sijp_tm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    allocate(de3d(xsz(1), xsz(2), 1-level:xsz(3)+level))

    u_tm(:,:,:,:) = 0.0_wp; temp_4d(:,:,:,:) = 0.0_wp
    uip_p_tm(:,:,:,:) = 0.0_wp
    pp_tm(:,:,:) = 0.0_wp; prod_tm(:,:,:) = 0.0_wp
    dissp_tm(:,:,:) = 0.0_wp; tke_tm(:,:,:) = 0.0_wp
    tran_turb_tm(:,:,:) = 0.0_wp; tran_visc_tm(:,:,:) = 0.0_wp
    tran_pres_tm(:,:,:) = 0.0_wp; adv_tm(:,:,:) = 0.0_wp
    d_u_x(:,:,:,:,:) = 0.0_wp;  sijp2_tm(:,:,:,:,:) = 0.0_wp
    uip_ujp_tm(:,:,:,:,:) = 0.0_wp; uip_ujp_ujp_tm(:,:,:,:) = 0.0_wp
    ujp_sijp_tm(:,:,:,:) = 0.0_wp; de3d(:,:,:) = 0.0_wp
  
    !call collect_grid
    call read_turbine_model_param(iturbine, time)
    
    volflux = 0.0
   
    if (myid.eq.0) then
      do i = 1, nz_global
        print *, "i, zzall, zwall", i, zzall(i), zwall(i)
      enddo
    endif
    
    root0 = 0
    ini_it = 2600

    open(1001, file="post_control.inp")
    read(1001, *) ini_it, end_it
    !read(1001, *) post_n_1d1, post_n_1d2, post_n_2d1, post_n_2d2
    read(1001, *) tb_xc, tb_yc, tb_zc
    read(1001, *) tb_d
    close(1001)
    
    allocate(templocal3d(xsz(1), xsz(2), 1-level:xsz(3)+level))

    templocal3d(:,:,:) = 0.0

    !> comment this line when turbine_model is used
    !allocate(eo(nx_global, ny_global), hho(nx_global, ny_global))

    if(mod(end_it-ini_it, noutd).ne.0) then
      print *, 'end_it-ini_it and noutd are not dividable'
      stop      
    endif
    
    if(myid.eq.0) print *, 'Global nx, ny, nz=', nx_global, ny_global, nz_global
    
    !allocate(temp3d(nx_global, ny_global, nz_global))

    !> plyunote: zz(-1)
    zz(1-level) = 2.0*zz(1) - zz(1+level)
    zw(1-level) = 2.0*zw(1) - zw(1+level)
    zz(xsz(3)+level) = 2.0 * zz(xsz(3)) - zz(xsz(3)-1)
    zw(xsz(3)+level) = 2.0 * zw(xsz(3)) - zw(xsz(3)-1)
    
    n_it = 0; time = 0.0;
    rem_it = ini_it - nint(real(ini_it)/real(noutd))*noutd
    
    if (myid.eq.0) call system('mkdir -p MEAN')

    call reread_mean("mean_field_3d.h5")
    !call dealiasxy(u(:,:,1:xsz(3)))
    !call dealiasxy(v(:,:,1:xsz(3)))
    !call dealiasxy(w(:,:,1:xsz(3)))
    !call dealiasxy(pp(:,:,1:xsz(3)))
    call update_ghost(u, level); call update_ghost(v, level); 
    call update_ghost(w, level); call update_ghost(pp, level);
    u_tm(:,:,:,1) = u(:,:,:); u_tm(:,:,:,2) = v(:,:,:)
    u_tm(:,:,:,3) = w(:,:,:); pp_tm(:,:,:) = pp(:,:,:)

    do ii = nint(real(ini_it)/real(noutd)), &
      nint(real(end_it)/real(noutd))
      n_it = n_it + 1
      it = ii * noutd + rem_it
      write(fileid,'(a4,I0.10)') 'DAT_',it
      if (myid.eq.0) print *, 'Post analysis output: ', fileid

      call reread_output(it, fileid, time)
      if(myid.eq.0) print *, 'time=', time

      !if (rankcol.eq.0) then
      !  print *, 'rank_z=',rankrow, ', u1=',u(1,1,:)
      !endif      
      
      call update_ghost(u, level); call update_ghost(v, level); 
      call update_ghost(w, level); call update_ghost(pp, level);
      
      !if (rankcol.eq.0) then
      !  print *, 'rank_z=',rankrow, ', u2=',u(1,1,:)
      !endif      
      
      !print *, "mark, 1" 
      if (it == ini_it ) then
        call hos_init
        !if(myid.eq.0) print *, 'hos_init called'
      endif

      !print *, "mark, 2"
      if (iwavy .eq. 7 .and. isbot) then
        write(fileid, '(i0.8)') it 
        fileid = trim(adjustl(fileid))
        fparam_hos = trim("restart_param_hos.dat"//fileid)
        fdata_hos = trim("restart_hos.h5"//fileid)
        call read_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos,&
          fparam_hos, fdata_hos)
        !if(myid.eq.0) print *, 'hos restarted'
        !print *, "eta_hos=", eta_hos(1:5,1)
     
        !do i = 1, size(eta_hos, 1)
        !  do j = 1, size(eta_hos, 2)
        !    if (isnan(eta_hos(i,j))) print *, "nan in eta_hos_in_post_main:", i, j
        !  enddo
        !enddo
        
        !> variables in following lines allocated in hos_init 
        call wavenum(wvn_hos)
        call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
        call zeta(eta_hos,zp_hos)
        tmp=1.0_wp
        call boundvp(vps_hos,r_hos,zp_hos)
        call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
        call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
        !if(myid.eq.0) print *, 'hos internal computation finished'
        
        call bottom_hos_les(time)
        !if(myid.eq.0) print *, 'bottom_hos_les called. Got ubs, eta, ub, hh,&
        !  & ht,hx'
      endif

      !if(myid==0) print *, 'HOS initialized!'

      ! transfer data to the upper cpus
      !print *, "mark, 3"
      call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.1"
      call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.2"
      call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.3"
      call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.4"
      call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.5"
      call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.6"
      call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.7"
      call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.8"
      !if(myid.eq.0) print *, 'HOS data transfered to upper cpus'

      !> we didn't calculate pressure here. Might be needed to consider spatial
      ! distribution and average.
      !call init_p_poisson
      
      !> nl_coef: input eta, hx, hy
      !>          output ex, ey, ehx, ehy, reh, her, exr, zetax,...
      call nl_coef
      
      !> gather eta and hh
      !print *, "mark, 4"
      call gather_2d_xy(eta, eo)
      call gather_2d_xy(hh, hho)

      !> check problematic data
      !print *, "mark, 5"
      do i = 1, xsz(1)
        do j = 1, xsz(2)
          inan = is_nan(eta(i,j))
          if (inan) print *, "NaN in eta"
          inan = is_nan(hh(i,j))
          if (inan) print *, "NaN in hh"
          do k = 1, xsz(3)
            inan = is_nan(u(i,j,k))
            if (inan) print *, "NaN in u"
            inan = is_nan(v(i,j,k))
            if (inan) print *, "NaN in v"
            inan = is_nan(w(i,j,k))
            if (inan) print *, "NaN in w"
          enddo
        enddo
      enddo
      do k = 1, xsz(3)
        inan = is_nan(zz(k))
        if (inan) print *, "NaN in zz"
      enddo

      !> filtered u
      deltax = 2.0_wp * TWOPI/pex/nx_global
      deltay = 2.0_wp * TWOPI/pey/ny_global
      
      uf = u; vf = v; wf = w
      !if(myid.eq.0) print *, 'mfilt=', mfilt
      call les_filter(uf(:,:,1:xsz(3)), mfilt, deltax, deltay)
      call les_filter(vf(:,:,1:xsz(3)), mfilt, deltax, deltay)
      call les_filter(wf(:,:,1:xsz(3)), mfilt, deltax, deltay)
      call update_ghost(uf, level)
      call update_ghost(vf, level)
      call update_ghost(wf, level)

      !print *, "mark, 7"
      if (istop) then
        !> ignore this les_filter step of xsz(3)+1 in navier_fu.f90
      endif

      if (isbot) then
        call wall_model_v3

        call gather_2d_xy(hx, hxo) 
      endif

      !> Get SGS stress and its contribution to x momentum

      !navier_les.f
      !input: u, v, w
      !output: s11~s33, s11w~s33w
      s11 = 0.0; s12 = 0.0
      call get_strain(u, v, w)
      !if (isbot) then
      !  print *, "s11=",s11(1:5,1,1), ", s12=",s12(1:5,1,1)
      !endif

      !input: u, v, w, uf, vf, wf, deltax, deltay
      !output: nut
      call get_nut(u, v, w, uf, vf, wf, deltax, deltay)
      nut = nut + 1.0_wp/re
      
      !> get turbine force and analyze discontinuity infor.
      if (iturbine .ne. 0) then
        call reread_rotor_force (time, it, level)
        if (nacelle_model .eq. 1) then
          call reread_nacelle_force (it, nacfx, nacfy, nacfz, level)
        endif
      endif
      if (iturbine .ge. 3 .and. (ids .eq. 1 .or. idsd .eq. 1)) then
        lim_dir = 1 !> indicate it is analyzing x direction
        call analyze_limit3d_from_3d (flag_de, tbn_lim_x)
        call transpose_xy(flag_de, bufy)
        lim_dir = 2 !> indicate it is analyzing y direction
        call analyze_limit3d_from_3d (bufy, tbn_lim_y)
        lim_dir = 0
      endif

      !> interpolate w to u,v,p grid.
      do k = 1, xsz(3)
        if ((xst(3)+k-1)==1) then
          templocal3d(:,:,k) = w(:,:,k)
        else
          templocal3d(:,:,k) = 0.5*w(:,:,k)+0.5*w(:,:,k-1)
        endif
      enddo
      call update_ghost(templocal3d, level)

      temp_4d(:,:,:,1) = u - u_tm(:,:,:,1)
      temp_4d(:,:,:,2) = v - u_tm(:,:,:,2)
      temp_4d(:,:,:,3) = templocal3d - u_tm(:,:,:,3)
      
      do i = 1, 3
        de3d = 0.5_wp * temp_4d(:,:,:,i)**2
        call dealiasxy(de3d(:,:,1:xsz(3)))
        call update_ghost(de3d, level)
        tke_tm = tke_tm + de3d
      enddo

      do i = 1, 3
        do j = 1, 3
          de3d = temp_4d(:,:,:,i) * temp_4d(:,:,:,j)
          call dealiasxy(de3d(:,:,1:xsz(3)))
          call update_ghost(de3d, level)
          uip_ujp_tm(:,:,:,i,j) = uip_ujp_tm(:,:,:,i,j) + de3d
        enddo
        call gradu(temp_4d(:,:,:,i), d_u_x(:,:,:,i,1), d_u_x(:,:,:,i,2), d_u_x(:,:,:,i,3))
      enddo

      do i = 1, 3
        do j = 1, 3
          de3d = 0.25_wp * (d_u_x(:,:,:,i,j) + d_u_x(:,:,:,j,i))**2
          call dealiasxy(de3d(:,:,1:xsz(3)))
          de3d = -2.0_wp * nut * de3d
          call dealiasxy(de3d(:,:,1:xsz(3)))
          call update_ghost(de3d, level)
          sijp2_tm(:,:,:,i,j) = sijp2_tm(:,:,:,i,j) + de3d 
        enddo
      enddo

      do i = 1, 3
        de3d = temp_4d(:,:,:,i) * (pp - pp_tm)
        call dealiasxy(de3d(:,:,1:xsz(3)))
        call update_ghost(de3d, level)
        uip_p_tm(:,:,:,i) = uip_p_tm(:,:,:,i) + de3d
        do j = 1, 3
          de3d = temp_4d(:,:,:,i) * temp_4d(:,:,:,j)
          call dealiasxy(de3d(:,:,1:xsz(3)))
          de3d = de3d * temp_4d(:,:,:,j)
          call dealiasxy(de3d(:,:,1:xsz(3)))
          call update_ghost(de3d, level)
          uip_ujp_ujp_tm(:,:,:,i) = uip_ujp_ujp_tm(:,:,:,i) + de3d
          
          de3d = temp_4d(:,:,:,j) * 0.5_wp * (d_u_x(:,:,:,i,j) + d_u_x(:,:,:,j,i))
          call dealiasxy(de3d(:,:,1:xsz(3)))
          de3d = de3d * nut
          call dealiasxy(de3d(:,:,1:xsz(3)))
          call update_ghost(de3d, level)
          ujp_sijp_tm(:,:,:,i) = ujp_sijp_tm(:,:,:,i) + de3d
        enddo
      enddo
    enddo

    tke_tm = tke_tm / n_it
    uip_ujp_tm = uip_ujp_tm / n_it
    sijp2_tm = sijp2_tm / n_it
    uip_p_tm = uip_p_tm / n_it
    uip_ujp_ujp_tm = uip_ujp_ujp_tm / n_it
    ujp_sijp_tm = ujp_sijp_tm / n_it

    !> plyunote: doing dealiasing is too late here. Later modify it.
    !call dealiasxy(tke_tm(:,:,1:xsz(3)))
    !do i = 1, 3
    !  call dealiasxy(uip_p_tm(:,:,1:xsz(3),i))
    !  call dealiasxy(uip_ujp_ujp_tm(:,:,1:xsz(3),i))
    !  call dealiasxy(ujp_sijp_tm(:,:,1:xsz(3),i))
    !  do j = 1, 3
    !    call dealiasxy(uip_ujp_tm(:,:,1:xsz(3),i, j))
    !    call dealiasxy(sijp2_tm(:,:,1:xsz(3),i,j))
    !  enddo
    !enddo

    call gradu(tke_tm(:,:,:), d_u_x(:,:,:,1,1), d_u_x(:,:,:,1,2), d_u_x(:,:,:,1,3))
    do i = 1, 3
      de3d = u_tm(:,:,:,i) * d_u_x(:,:,:,1,i)
      call dealiasxy(de3d(:,:,1:xsz(3)))
      call update_ghost(de3d, level)
      adv_tm = adv_tm - de3d
    enddo

    prod_tm(:,:,:) = 0.0_wp

    do i =1, 3
      call gradu(u_tm(:,:,:,i), d_u_x(:,:,:,i,1), d_u_x(:,:,:,i,2), d_u_x(:,:,:,i,3))
    enddo

    do i = 1, 3
      do j = 1, 3
        de3d = uip_ujp_tm(:,:,:,i,j) * d_u_x(:,:,:,i,j)
        call dealiasxy(de3d(:,:,1:xsz(3)))
        prod_tm(:,:,:) = prod_tm(:,:,:) - de3d
        dissp_tm(:,:,:) = dissp_tm(:,:,:) + sijp2_tm(:,:,:,i,j)
      enddo
    enddo

    !call dealiasxy(prod_tm(:,:,1:xsz(3)))

    call calc_divu_center(tran_turb_tm(:,:,1:xsz(3)), uip_ujp_ujp_tm(:,:,:,1), &
      uip_ujp_ujp_tm(:,:,:,2), uip_ujp_ujp_tm(:,:,:,3)) 
    call calc_divu_center(tran_pres_tm(:,:,1:xsz(3)), uip_p_tm(:,:,:,1), &
      uip_p_tm(:,:,:,2), uip_p_tm(:,:,:,3)) 
    call calc_divu_center(tran_visc_tm(:,:,1:xsz(3)), ujp_sijp_tm(:,:,:,1), &
      ujp_sijp_tm(:,:,:,2), ujp_sijp_tm(:,:,:,3)) 

    tran_turb_tm = tran_turb_tm * (-0.5_wp)
    tran_pres_tm = tran_pres_tm * (-1.0_wp)
    tran_visc_tm = tran_visc_tm * (2.0_wp)

    !call update_ghost(tran_turb_tm, level)
    !call update_ghost(tran_pres_tm, level)
    !call update_ghost(tran_visc_tm, level)

    do i = 1, xsz(1)
      do j = 1, xsz(2)
        do k = 1, xsz(3)
          !templocal3d(i,j,k) = zz(k)*(hbar+hh(i,j))-hh(i,j)
          templocal3d(i,j,k) = zz(k) * hbar
        enddo
      enddo
    enddo

    if (myid.eq.0) print *, "Writing output to mean_tke_3d.h5"
    call IO_Init
    call writerOpen("mean_tke_3d.h5", mpi_comm_2d_cart, writer)
    call write3d_xyz(writer, "z", templocal3d(:,:,1:xsz(3)))
    call write3d_xyz(writer, "TKE", tke_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "Advection", adv_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "Production", prod_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "Dissipation", dissp_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "TransportByTurbulence", tran_turb_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "TransportByPressure", tran_pres_tm(:,:,1:xsz(3)))
    call write3d_xyz(writer, "TransportByViscousStress", tran_visc_tm(:,:,1:xsz(3)))
    call writerClose(writer)
    call IO_Finalize  

    deallocate(tmp)
    deallocate(tempz, zzall, zwall)
    deallocate(t13wx)

    deallocate(templocal3d)
    deallocate(u_tm, temp_4d, uip_p_tm, uip_ujp_ujp_tm, ujp_sijp_tm)
    deallocate(pp_tm, prod_tm, dissp_tm, tke_tm, adv_tm)
    deallocate(tran_turb_tm, tran_visc_tm, tran_pres_tm)
    deallocate(d_u_x, sijp2_tm)
    deallocate(uip_ujp_tm)
    deallocate(de3d)

    call fft_finalize
    call decomp_finalize
    !HOS
    call fft_finalize_hos
   
    if (myid.eq.0) print *, "Post process end successfully." 
    call mpi_finalize(ierror)

  end subroutine post_tke_budget

  !> added by plyu
  subroutine post_vortex_dynamics
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

    !> added by plyu
    use turbine_model
    use solver_common, only : tmp_x6
    use lib_array, only: is_nan 
    use discontinuity_smooth, only: analyze_limit3d_from_3d, tbn_lim_x, &
        tbn_lim_y, lim_dir 
    
    use hdf_io
    !end
    
    implicit none
    
    integer :: ioutc, ioutd
    real(wp) :: time
    
    integer :: ierror, ierr
    
    double precision :: t1, t2
    
    !HOS variables
    integer icon
    
    integer(INT64) :: id
    integer :: i, j, k, it, root0
    
    real(wp), allocatable, dimension(:) :: fpk, betak, ampfd, epskx
    
    real(wp) beta_miles
    
    character (len=100) :: fileid, fgrid, fparam_les, fdata_les, faux, fparam_hos, fdata_hos

    real(wp) c_phase, u_lambda, u_half_lambda, wage1, wage2, wage3
    real(wp) clen, cspeed, ctime, x, y, fpt

    real(wp), allocatable, dimension(:,:) :: etaall, hhall, tmp
    real(wp) ftn, rek, dwk, disspm, l_kom_ave

    real(wp), allocatable, dimension(:,:,:) :: uupall, wwpall, velall, ppall, etasig, pp_wave
    real(wp), allocatable, dimension(:,:,:,:) :: ppsig
    real(wp), allocatable, dimension(:,:) :: dtheta, uxy, vxy, wxy
    real(wp), allocatable, dimension(:,:) :: ek11, ek22, ek33
    real(wp), allocatable, dimension(:,:) :: ek11_tm, ek22_tm, ek33_tm
    real(wp), allocatable, dimension(:,:) :: eki11, eki22, eki33
    real(wp), allocatable, dimension(:,:) :: eki11_tm, eki22_tm, eki33_tm
    real(wp), allocatable, dimension(:) :: dissp, ekm, l_kom, ekim

    real(wp), allocatable, dimension(:) :: q1, q2, q3, q4, um, vm, wm, uup, vvp, wwp, zzall, uw_shear
   
    !> variables added by plyu
    integer :: ii
    real(wp), allocatable, dimension(:,:,:) :: h5_3dbuffer
    
    integer :: idp ! counter for debug
    integer :: ini_it, end_it, nk, n_it, rem_it
    integer :: rankrow, rankcol, sizerow, sizecol
    real(wp) :: deltax, deltay
    
    real(wp), allocatable, dimension(:,:,:) :: templocal3d 

    real(wp) :: volflux(3)

    type(HDFObj) :: writer
    integer :: ih5s

    !> local variables for turbine model
    integer :: ibi, n_elmt_1, nb, n1e, n2e, n3e
    real(wp), dimension(3) :: dr_vec, ds12, ds13, tmparray
    real(wp) :: dr
    
    real(wp) :: tb_xc, tb_yc, tb_zc, tb_d, xtemp, ytemp, ztemp
    real(wp), allocatable, dimension(:) :: tempz, zwall
    real(wp), allocatable, dimension(:,:,:) :: temp3d
    real(wp) :: utopm
    
    !> comment this line when turbine_model is used
    !real(wp), allocatable, dimension(:,:) :: eo, hho, hxo, hyo
    !real(wp), allocatable, dimension(:,:) :: hxo 

    !> variable for SGS
    real(wp), allocatable, dimension(:,:) :: t13wx

    !> variable for nacelle model force
    real(wp), allocatable, dimension(:,:,:) :: nacfx, nacfy, nacfz

    !>vortex dynamics: it is implemented in another  individual subroutine
    real(wp), allocatable, dimension(:,:,:,:,:) :: d_u_x, d_u_x_tm, vorj_duidxj_tm
    real(wp), allocatable, dimension(:,:,:,:) :: vor, vor_tm
    real(wp), allocatable, dimension(:,:,:) :: q, q_tm


    logical :: inan

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
    
    call grid_gen
     
    icon = 16

    if (myid .eq. 0) then
      print *, "********* This is vortex-dynamics post-processing script *******"
    endif

    !> plyunote: zz(-1)
    zz(1-level) = 2.0*zz(1) - zz(1+level)
    zw(1-level) = 2.0*zw(1) - zw(1+level)
    zz(xsz(3)+level) = 2.0 * zz(xsz(3)) - zz(xsz(3)-1)
    zw(xsz(3)+level) = 2.0 * zw(xsz(3)) - zw(xsz(3)-1)
    allocate(tmp(nx_global,ny_global))

    level = 1
    call mpi_comm_rank(mpi_comm_2d_row, rankrow, ierr) ! rank in z
    call mpi_comm_rank(mpi_comm_2d_col, rankcol, ierr) ! rank in y
    call mpi_comm_size(mpi_comm_2d_row, sizerow, ierr)
    call mpi_comm_size(mpi_comm_2d_col, sizecol, ierr)
    if(myid.eq.0) print *, 'sizerow=',sizerow,', sizecol=',sizecol
    
    ! Gather zz coordinates
    allocate(tempz(nz_global), zzall(nz_global), zwall(nz_global))
    tempz = 0.0; zzall = 0.0; zwall = 0.0 
    
    tempz(xst(3):xend(3)) = zz(1:xsz(3))
    call MPI_Allreduce(tempz, zzall, nz_global, mpi_double_precision, &
         MPI_SUM, MPI_COMM_2D_ROW, ierror)
    tempz = 0.0
    tempz(xst(3):xend(3)) = zw(1:xsz(3))
    call MPI_Allreduce(tempz, zwall, nz_global, mpi_double_precision, &
         MPI_SUM, MPI_COMM_2D_ROW, ierror)
    
    allocate(t13wx(xsz(1),xsz(2)))
    
    allocate(nacfx(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(nacfy(xsz(1), xsz(2), 1-level:xsz(3)+level))
    allocate(nacfz(xsz(1), xsz(2), 1-level:xsz(3)+level))
    nacfx = 0.0; nacfy = 0.0; nacfz = 0.0

    !allocate(h5_3dbuffer(xsz(1), xsz(2), xsz(3)))
    !h5_3dbuffer = 0.0_wp
    
    allocate(d_u_x(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    allocate(vor(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    allocate(q(xsz(1), xsz(2), 1-level:xsz(3)+level))
    d_u_x(:,:,:,:,:) = 0.0_wp; vor(:,:,:,:) = 0.0_wp; q(:,:,:) = 0.0_wp

    allocate(d_u_x_tm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    allocate(vor_tm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3))
    allocate(vorj_duidxj_tm(xsz(1), xsz(2), 1-level:xsz(3)+level, 3, 3))
    allocate(q_tm(xsz(1), xsz(2), 1-level:xsz(3)+level))
    d_u_x_tm(:,:,:,:,:) = 0.0_wp; q_tm(:,:,:) = 0.0_wp
    vor_tm(:,:,:,:) = 0.0_wp; vorj_duidxj_tm(:,:,:,:,:) = 0.0_wp
    
    !call collect_grid
    call read_turbine_model_param(iturbine, time)
    
    volflux = 0.0
   
    if (myid.eq.0) then
      do i = 1, nz_global
        print *, "i, zzall, zwall", i, zzall(i), zwall(i)
      enddo
    endif
    
    root0 = 0
    ini_it = 2600

    open(1001, file="post_control.inp")
    read(1001, *) ini_it, end_it
    !read(1001, *) post_n_1d1, post_n_1d2, post_n_2d1, post_n_2d2
    read(1001, *) tb_xc, tb_yc, tb_zc
    read(1001, *) tb_d
    close(1001)
    
    allocate(templocal3d(xsz(1), xsz(2), 1-level:xsz(3)+level))

    templocal3d(:,:,:) = 0.0

    !> comment this line when turbine_model is used
    !allocate(eo(nx_global, ny_global), hho(nx_global, ny_global))

    if(mod(end_it-ini_it, noutd).ne.0) then
      print *, 'end_it-ini_it and noutd are not dividable'
      stop      
    endif
    
    if(myid.eq.0) print *, 'Global nx, ny, nz=', nx_global, ny_global, nz_global
    
    !allocate(temp3d(nx_global, ny_global, nz_global))

    n_it = 0; time = 0.0;
    rem_it = ini_it - nint(real(ini_it)/real(noutd))*noutd
    
    if (myid.eq.0) call system('mkdir -p VORTEX')
    if (myid.eq.0) call system('mkdir -p MEAN')

    do ii = nint(real(ini_it)/real(noutd)), &
      nint(real(end_it)/real(noutd))
      n_it = n_it + 1
      it = ii * noutd + rem_it
      write(fileid,'(a4,I0.10)') 'DAT_',it
      if (myid.eq.0) print *, 'Post analysis output: ', fileid

      call reread_output(it, fileid, time)
      if(myid.eq.0) print *, 'time=', time

      !if (rankcol.eq.0) then
      !  print *, 'rank_z=',rankrow, ', u1=',u(1,1,:)
      !endif      
      
      call update_ghost(u, level); call update_ghost(v, level); 
      call update_ghost(w, level); call update_ghost(pp, level);
      
      !if (rankcol.eq.0) then
      !  print *, 'rank_z=',rankrow, ', u2=',u(1,1,:)
      !endif      
      
      !print *, "mark, 1" 
      if (it == ini_it ) then
        call hos_init
        !if(myid.eq.0) print *, 'hos_init called'
      endif

      !print *, "mark, 2"
      if (iwavy .eq. 7 .and. isbot) then
        write(fileid, '(i0.8)') it 
        fileid = trim(adjustl(fileid))
        fparam_hos = trim("restart_param_hos.dat"//fileid)
        fdata_hos = trim("restart_hos.h5"//fileid)
        call read_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos,&
          fparam_hos, fdata_hos)
        !if(myid.eq.0) print *, 'hos restarted'
        !print *, "eta_hos=", eta_hos(1:5,1)
     
        !do i = 1, size(eta_hos, 1)
        !  do j = 1, size(eta_hos, 2)
        !    if (isnan(eta_hos(i,j))) print *, "nan in eta_hos_in_post_main:", i, j
        !  enddo
        !enddo
        
        !> variables in following lines allocated in hos_init 
        call wavenum(wvn_hos)
        call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
        call zeta(eta_hos,zp_hos)
        tmp=1.0_wp
        call boundvp(vps_hos,r_hos,zp_hos)
        call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
        call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
        !if(myid.eq.0) print *, 'hos internal computation finished'
        
        call bottom_hos_les(time)
        !if(myid.eq.0) print *, 'bottom_hos_les called. Got ubs, eta, ub, hh,&
        !  & ht,hx'
      endif

      !if(myid==0) print *, 'HOS initialized!'

      ! transfer data to the upper cpus
      !print *, "mark, 3"
      call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.1"
      call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.2"
      call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.3"
      call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.4"
      call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.5"
      call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.6"
      call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.7"
      call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row,ierr)
      !print *, "mark, 3.8"
      !if(myid.eq.0) print *, 'HOS data transfered to upper cpus'

      !> we didn't calculate pressure here. Might be needed to consider spatial
      ! distribution and average.
      !call init_p_poisson
      
      !> nl_coef: input eta, hx, hy
      !>          output ex, ey, ehx, ehy, reh, her, exr, zetax,...
      call nl_coef
      
      !> gather eta and hh
      !print *, "mark, 4"
      call gather_2d_xy(eta, eo)
      call gather_2d_xy(hh, hho)

      !> check problematic data
      !print *, "mark, 5"
      do i = 1, xsz(1)
        do j = 1, xsz(2)
          inan = is_nan(eta(i,j))
          if (inan) print *, "NaN in eta"
          inan = is_nan(hh(i,j))
          if (inan) print *, "NaN in hh"
          do k = 1, xsz(3)
            inan = is_nan(u(i,j,k))
            if (inan) print *, "NaN in u"
            inan = is_nan(v(i,j,k))
            if (inan) print *, "NaN in v"
            inan = is_nan(w(i,j,k))
            if (inan) print *, "NaN in w"
          enddo
        enddo
      enddo
      do k = 1, xsz(3)
        inan = is_nan(zz(k))
        if (inan) print *, "NaN in zz"
      enddo

      !> get turbine force and analyze discontinuity infor.
     if (iturbine .ne. 0) then
       call reread_rotor_force (time, it, level)
       if (nacelle_model .eq. 1) then
         call reread_nacelle_force (it, nacfx, nacfy, nacfz, level)
       endif
     endif
      if (iturbine .ge. 3 .and. (ids .eq. 1 .or. idsd .eq. 1)) then
        lim_dir = 1 !> indicate it is analyzing x direction
        call analyze_limit3d_from_3d (flag_de, tbn_lim_x)
        call transpose_xy(flag_de, bufy)
        lim_dir = 2 !> indicate it is analyzing y direction
        call analyze_limit3d_from_3d (bufy, tbn_lim_y)
        lim_dir = 0
      endif

     !> calculate volume flux
      !print *, "mark, 6"
      call calc_vol_flux(u(:,:,1:xsz(3)), volflux)

      !> filtered u
      deltax = 2.0_wp * TWOPI/pex/nx_global
      deltay = 2.0_wp * TWOPI/pey/ny_global
      
      !> interpolate w to u,v,p grid.
      do k = 1, xsz(3)
        if ((xst(3)+k-1)==1) then
          templocal3d(:,:,k) = w(:,:,k)
        else
          templocal3d(:,:,k) = 0.5*w(:,:,k)+0.5*w(:,:,k-1)
        endif
      enddo
      call update_ghost(templocal3d, level)

      !> Vorticity
      call gradu(u, d_u_x(:,:,:,1,1), d_u_x(:,:,:,1,2), d_u_x(:,:,:,1,3))
      call gradu(v, d_u_x(:,:,:,2,1), d_u_x(:,:,:,2,2), d_u_x(:,:,:,2,3))
      call gradu(templocal3d, d_u_x(:,:,:,3,1), d_u_x(:,:,:,3,2), d_u_x(:,:,:,3,3))   
      vor(:,:,:,1) = d_u_x(:,:,:,3,2) - d_u_x(:,:,:,2,3)
      vor(:,:,:,2) = d_u_x(:,:,:,1,3) - d_u_x(:,:,:,3,1)
      vor(:,:,:,3) = d_u_x(:,:,:,2,1) - d_u_x(:,:,:,1,2)

      q(:,:,:) = 0.25*(vor(:,:,:,1)**2+vor(:,:,:,2)**2+vor(:,:,:,3)**2) &
        - 0.5*(d_u_x(:,:,:,1,1)**2 + d_u_x(:,:,:,2,2)**2 + d_u_x(:,:,:,3,3)**2) &
        - 0.25*(d_u_x(:,:,:,1,2)+d_u_x(:,:,:,2,1))**2 &
        - 0.25*(d_u_x(:,:,:,1,3)+d_u_x(:,:,:,3,1))**2 &
        - 0.25*(d_u_x(:,:,:,3,2)+d_u_x(:,:,:,2,3))**2 

      d_u_x_tm = d_u_x_tm + d_u_x
      vor_tm = vor_tm + vor
      do i = 1, 3
        do j = 1, 3
          vorj_duidxj_tm(:,:,:,i,j) = vorj_duidxj_tm(:,:,:,i,j) + vor(:,:,:,j) &
            * d_u_x(:,:,:,i,j)
        enddo
      enddo
      q_tm = q_tm + q

      do i = 1, xsz(1)
        do j = 1, xsz(2)
          do k = 1, xsz(3)
            templocal3d(i,j,k) = zz(k)*(hbar+hh(i,j))-hh(i,j)
            !templocal3d(i,j,k) = zz(k) * hbar
          enddo
        enddo
      enddo
    
      call IO_Init
      write(fileid,'(a11,I0.10,a3)') 'VORTEX/VOR_',it,'.h5'
      if (myid.eq.0) print *, "Writing output to ", fileid
      call writerOpen(fileid, mpi_comm_2d_cart, writer)
      !h5_3dbuffer(:,:,:) = templocal3d(:,:,1:xsz(3))
      call write3D_xyz(writer, "z", templocal3d(:,:,1:xsz(3)))
      !call write3d_xyz(writer, "u", u_tm(:,:,1:xsz(3)))
      !call write3d_xyz(writer, "v", v_tm(:,:,1:xsz(3)))
      !call write3d_xyz(writer, "w", w_tm(:,:,1:xsz(3)))
      call write3d_xyz(writer, "VOR_X", vor(:,:,1:xsz(3),1))
      call write3d_xyz(writer, "VOR_Y", vor(:,:,1:xsz(3),2))
      call write3d_xyz(writer, "VOR_Z", vor(:,:,1:xsz(3),3))
      call write3d_xyz(writer, "Q", q(:,:,1:xsz(3)))
      call writerClose(writer)
      call IO_Finalize

      
      !if (myid.eq.0) then
      !  !> plyunote: use ascii dat format to output 3D data is irreasonable.
      !  write(fileid,'(a,I0.10,a)') 'VORTEX/VOR_',it,'.dat'
      !  if (myid.eq.0) print *, "Writing output to ", fileid
      !  open(1003, file=fileid, action='write') 
      !  write(1003, *) "VARIABLES = X, Y, Z, VOR_X, VOR_Y, VOR_Z, Q"
      !  write(1003,*) 'ZONE T="', it*dt, '" I=', nx_global, ' J=', ny_global, &
      !    ' K=', nz_global, ' SOLUTIONTIME=', it*dt
      !  do k = 1, nz_global
      !    ztemp = zz(k) * (hbar+hh(i,j))-hh(i,j)
      !    do j = 1, ny_global
      !      do i = 1, nx_global
      !        write(1003, *) cartx(i), carty(j), ztemp, vor(i,j,k,1), &
      !          vor(i,j,k,2), vor(i,j,k,3), q(i,j,k)
      !      enddo
      !    enddo
      !  enddo
      !  close(1003)
      !endif

    enddo

    d_u_x_tm = d_u_x_tm / n_it
    vor_tm = vor_tm / n_it
    vorj_duidxj_tm = vorj_duidxj_tm / n_it
    q_tm = q_tm / n_it

    do i = 1, xsz(1)
      do j = 1, xsz(2)
        do k = 1, xsz(3)
          !templocal3d(i,j,k) = zz(k)*(hbar+hh(i,j))-hh(i,j)
          templocal3d(i,j,k) = zz(k) * hbar
        enddo
      enddo
    enddo

    call write_h5_3d_one(templocal3d(:,:,1:xsz(3)), "z", "z_tm.h5")
    call write_h5_3d_one(vor_tm(:,:,1:xsz(3),1), "VOR_X", "vor_1_tm.h5")
    call write_h5_3d_one(vor_tm(:,:,1:xsz(3),2), "VOR_Y", "vor_2_tm.h5")
    call write_h5_3d_one(vor_tm(:,:,1:xsz(3),3), "VOR_Z", "vor_3_tm.h5")
    call write_h5_3d_one(vorj_duidxj_tm(:,:,1:xsz(3),1,1), "VORX_DUDX", "vorj_duidxj_1_1_tm.h5")
    call write_h5_3d_one(vorj_duidxj_tm(:,:,1:xsz(3),2,1), "VORX_DVDX", "vorj_duidxj_2_1_tm.h5")
    call write_h5_3d_one(vorj_duidxj_tm(:,:,1:xsz(3),3,1), "VORX_DWDX", "vorj_duidxj_3_1_tm.h5")
    call write_h5_3d_one(vorj_duidxj_tm(:,:,1:xsz(3),1,2), "VORY_DUDY", "vorj_duidxj_1_2_tm.h5")
    call write_h5_3d_one(vorj_duidxj_tm(:,:,1:xsz(3),2,2), "VORY_DVDY", "vorj_duidxj_2_2_tm.h5")
    call write_h5_3d_one(vorj_duidxj_tm(:,:,1:xsz(3),3,2), "VORY_DWDY", "vorj_duidxj_3_2_tm.h5")
    call write_h5_3d_one(vorj_duidxj_tm(:,:,1:xsz(3),1,3), "VORZ_DUDZ", "vorj_duidxj_1_3_tm.h5")
    call write_h5_3d_one(vorj_duidxj_tm(:,:,1:xsz(3),2,3), "VORZ_DVDZ", "vorj_duidxj_2_3_tm.h5")
    call write_h5_3d_one(vorj_duidxj_tm(:,:,1:xsz(3),3,3), "VORZ_DWDZ", "vorj_duidxj_3_3_tm.h5")
    call write_h5_3d_one(d_u_x_tm(:,:,1:xsz(3),1,1), "DUDX", "duidxj_1_1_tm.h5")
    call write_h5_3d_one(d_u_x_tm(:,:,1:xsz(3),1,2), "DUDY", "duidxj_1_2_tm.h5")
    call write_h5_3d_one(d_u_x_tm(:,:,1:xsz(3),1,3), "DUDZ", "duidxj_1_3_tm.h5")
    call write_h5_3d_one(d_u_x_tm(:,:,1:xsz(3),2,1), "DVDX", "duidxj_2_1_tm.h5")
    call write_h5_3d_one(d_u_x_tm(:,:,1:xsz(3),2,2), "DVDY", "duidxj_2_2_tm.h5")
    call write_h5_3d_one(d_u_x_tm(:,:,1:xsz(3),2,3), "DVDZ", "duidxj_2_3_tm.h5")
    call write_h5_3d_one(d_u_x_tm(:,:,1:xsz(3),3,1), "DWDX", "duidxj_3_1_tm.h5")
    call write_h5_3d_one(d_u_x_tm(:,:,1:xsz(3),3,2), "DWDY", "duidxj_3_2_tm.h5")
    call write_h5_3d_one(d_u_x_tm(:,:,1:xsz(3),3,3), "DWDZ", "duidxj_3_3_tm.h5")
    call write_h5_3d_one(q_tm(:,:,1:xsz(3)), "QM2", "q2_tm.h5")
    
    q_tm(:,:,:) = 0.25*(vor_tm(:,:,:,1)**2+vor_tm(:,:,:,2)**2+vor_tm(:,:,:,3)**2) &
      - 0.5*(d_u_x_tm(:,:,:,1,1)**2 + d_u_x_tm(:,:,:,2,2)**2 + d_u_x_tm(:,:,:,3,3)**2) &
      - 0.25*(d_u_x_tm(:,:,:,1,2)+d_u_x_tm(:,:,:,2,1))**2 &
      - 0.25*(d_u_x_tm(:,:,:,1,3)+d_u_x_tm(:,:,:,3,1))**2 &
      - 0.25*(d_u_x_tm(:,:,:,3,2)+d_u_x_tm(:,:,:,2,3))**2 
    call write_h5_3d_one(q_tm(:,:,1:xsz(3)), "Q", "q_tm.h5")
    
    deallocate(templocal3d)

    deallocate(vor, vor_tm)
    deallocate(q, q_tm)
    deallocate(vorj_duidxj_tm, d_u_x, d_u_x_tm)

    deallocate(tmp)
    deallocate(tempz, zzall, zwall)
    deallocate(t13wx)

    call fft_finalize
    call decomp_finalize
    !HOS
    call fft_finalize_hos
   
    if (myid.eq.0) print *, "Post process end successfully." 
    call mpi_finalize(ierror)
  end subroutine post_vortex_dynamics

  subroutine write_h5_3d_one(var, varname, filename)
    use hdf_io
    use decomp
    use mpi

    implicit none

    real(wp), intent(in) :: var(1:xsz(1), 1:xsz(2), 1:xsz(3))
    character (len=*), intent(in)  :: varname, filename
    character(len=128) :: varname2, filename2

    type(HDFObj) :: writer

    write(filename2, *) filename
    write(varname2, *) varname
    call IO_Init
    if (myid.eq.0) print *, "Writing output to ", trim(adjustl(filename2))
    call writerOpen(trim(adjustl(filename2)), mpi_comm_2d_cart, writer)
    call write3D_xyz(writer, trim(adjustl(varname2)), var)
    call writerClose(writer)
    call IO_Finalize  

  end subroutine write_h5_3d_one

  subroutine post_les_hos_wave_coherent_analysis
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

    use hdf_io
    !end

    implicit none

    integer :: ioutc, ioutd
    real(wp) :: time

    integer :: ierror, ierr

    double precision :: t1, t2

    !HOS variables
    integer icon

    integer(INT64) :: id
    integer :: i, j, k, it, root0

    character (len=100) :: fileid, fgrid, fparam_les, fdata_les, faux, fparam_hos, fdata_hos

    real(wp) c_phase, u_lambda, u_half_lambda, wage1, wage2, wage3
    real(wp) clen, cspeed, ctime, x, y, fpt

    real(wp), allocatable, dimension(:,:) :: etaall, hhall, tmp
    real(wp) ftn, rek

    real(wp), allocatable, dimension(:,:,:) :: ppall, etasig, pp_wave
    real(wp), allocatable, dimension(:,:,:,:) :: ppsig
    real(wp), allocatable, dimension(:,:) :: dtheta

    integer ini_it, nk

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

    allocate(etaall(nx_global,ny_global))
    allocate(hhall(nx_global,ny_global))
    allocate(tmp(nx_global,ny_global))

    allocate(ppall(nx_global,ny_global,nz_global))

    allocate(etasig(nx_global,ny_global,ntime))
    allocate(pp_wave(ntime,nx_global,ny_global))
    allocate(ppsig(ntime,nx_global,ny_global,6))
    allocate(dtheta(nx_global,ny_global))

    read(14,*) clen, cspeed, ctime

    if (myid == 0) then
!       open(72,file="q_analysis.dat")
    end if

    ini_it = 600
    do it = ini_it + 1, ini_it + ntime
       write(fileid,*) it
       fileid = trim(adjustl(fileid))

       fparam_les = trim("restart_param.dat"//fileid)
       fdata_les = trim("restart.h5"//fileid)
       faux =  trim("restart_aux.h5"//fileid)
       fgrid = trim("grid.h5"//fileid)

       fparam_hos = trim("restart_param_hos.dat"//fileid)
       fdata_hos = trim("restart_hos.h5"//fileid)

       call reread(time, ioutd, ioutc, fgrid, fparam_les, fdata_les, faux)
       call grid_gen

       !need to add hos initialization subroutines below
       if(myid==0) print *, 'HOS initialization started!'

       if (it == ini_it + 1) then
          call hos_init
       end if
       call read_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos, fparam_hos, fdata_hos)

       call wavenum(wvn_hos)
       call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
       call zeta(eta_hos,zp_hos)
       tmp=1.0_wp
       call boundvp(vps_hos,r_hos,zp_hos)
       call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
       call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

       call bottom_hos_les(time)

       if(myid==0) print *, 'HOS initialized!'

       ! transfer data to the upper cpus
       call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)

       eta=hh+0

       call init_p_poisson
       call nl_coef

        !test starts here!
       root0 = 0

       call gather_2d_xy(hh(1:xsz(1),1:xsz(2)),etaall,root0)
       etasig(:,:,it-ini_it) = etaall * (-1.0)
       call gather_3d_xyz(pp(1:xsz(1),1:xsz(2),1:xsz(3)),ppall,root0)
       do k = 1, 6
          do j = 1, ny_global
            do i = 1, nx_global
                ppsig(it-ini_it,i,j,k) = ppall(i,j,k)
             end do
          end do
       end do

!       call turb_analysis
       if(myid==0) print *,"finish step ",it
       !test ends here!

    end do

    if (myid == 0) then

       do it = 1, ntime
          write(936,*) "VARIABLES = x, y, pp, eta"
          write(936,"(A,F20.8,A,I5,A,I5,A,F20.8)") 'ZONE T="',it*1.0,'" I=',nx_global, &
               ' J=',ny_global,' SOLUTIONTIME=', it*1.0
          do j=1, ny_global
             do i=1, nx_global
                x = (i-1)*xl/nx_global
                y = (j-1)*yl/ny_global
!                write(936,'(25e12.4)') x, y, ppsig(it,i,j,2), etasig(i,j,it)
            end do
          end do
          write(939,'(25e12.4)') it*1.0, ppsig(it,1,1,2), etasig(1,1,it)
      end do

      do k = 2, 6
         do j = 1, ny_global
            do i = 1, nx_global
               call get_wave_coherent_linft(ppsig(:,i,j,k),pp_wave(:,i,j),etasig(i,j,:),ntime,dtheta(i,j))
            end do
         end do
         do it = 1, ntime
            write(937,*) "VARIABLES = x, y, pp_wave"
            write(937,"(A,F20.8,A,I5,A,I5,A,F20.8)") 'ZONE T="',it*1.0,'" I=',nx_global, &
                 ' J=',ny_global,' SOLUTIONTIME=', it*1.0
            do j=1, ny_global
               do i=1, nx_global
                  x = (i-1)*xl/nx_global
                  y = (j-1)*yl/ny_global
                  !                write(937,'(25e12.4)') x, y, pp_wave(it,i,j)
               end do
            end do
         end do

         write(938,*) "VARIABLES = x, y, dtheta"
         write(938,"(A,F20.8,A,I5,A,I5,A,F20.8)") 'ZONE T="',1.0,'" I=',nx_global, &
              ' J=',ny_global,' SOLUTIONTIME=', k*1.0
         do j=1, ny_global
            do i=1, nx_global
               x = (i-1)*xl/nx_global
               y = (j-1)*yl/ny_global
               write(938,'(25e12.4)') x, y, dtheta(i,j)
            end do
         end do
      end do
    end if
    if (myid == 0) then
       print *,"finish step ",it
    end if

    deallocate(etaall,hhall,tmp)

    deallocate(ppall,etasig,pp_wave,ppsig,dtheta)

    ! call spectral_finalize
    call fft_finalize
    call decomp_finalize
    !HOS
    call fft_finalize_hos

    call mpi_finalize(ierror)
  end subroutine post_les_hos_wave_coherent_analysis
 
  subroutine post_les_hos_k_omg_spectrum
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

    use hdf_io
    !end

    implicit none

    integer :: ioutc, ioutd
    real(wp) :: time

    integer :: ierror, ierr

    double precision :: t1, t2

    !HOS variables
    integer icon

    integer(INT64) :: id
    integer :: i, j, k, it, root0

    character (len=100) :: fileid, fgrid, fparam_les, fdata_les, faux, fparam_hos, fdata_hos

    real(wp) c_phase, u_lambda, u_half_lambda, wage1, wage2, wage3
    real(wp) clen, cspeed, ctime, x, y, fpt

    real(wp), allocatable, dimension(:,:) :: etaall, hhall, tmp
    real(wp) ftn, rek, tmppet

    real(wp), allocatable, dimension(:) :: um, vm, uup, vvp, wwp, zzall
    real(wp), allocatable, dimension(:,:,:) :: uall, vall, usig_nw, usig_out

    real(wp), allocatable, dimension(:,:) :: e11_k1k2, e11_k1omg
    real(wp), allocatable, dimension(:) :: e11_k1, k1co,k2co,k3co
    real(wp), allocatable, dimension(:,:,:) :: ekkomg, ef3d
    integer nomg

    real(wp), allocatable, dimension(:,:,:) :: ppall, etasig, pp_wave
    real(wp), allocatable, dimension(:,:,:,:) :: ppsig
    real(wp), allocatable, dimension(:,:) :: dtheta

    integer ini_it, nk

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

    allocate(etaall(nx_global,ny_global))
    allocate(hhall(nx_global,ny_global))
    allocate(tmp(nx_global,ny_global))

    allocate(uall(nx_global,ny_global,nz_global))
    allocate(vall(nx_global,ny_global,nz_global))
    allocate(usig_nw(nx_global,ny_global,ntime))
    allocate(usig_out(nx_global,ny_global,ntime))
    allocate(um(nz_global), uup(nz_global))
    allocate(vm(nz_global), vvp(nz_global))
!    allocate(wm(nz_global), wwp(nz_global))

    nomg = ntime
    allocate(zzall(nz_global))
    allocate(e11_k1k2(nx_global,ny_global))
    allocate(e11_k1(nx_global))
    allocate(ef3d(nx_global,ny_global,nomg))
    allocate(k1co(nx_global),k2co(ny_global),k3co(nomg))
    allocate(ekkomg(nx_global, ny_global,nomg))
    allocate(e11_k1omg(nx_global,nomg))

!    allocate(ppall(nx_global,ny_global,nz_global))

!    allocate(etasig(nx_global,ny_global,ntime))
!    allocate(pp_wave(ntime,nx_global,ny_global))
!    allocate(ppsig(ntime,nx_global,ny_global,6))
!    allocate(dtheta(nx_global,ny_global))

    read(14,*) clen, cspeed, ctime

    if (myid == 0) then
!       open(72,file="q_analysis.dat")
    end if

    ini_it = 600
    do it = ini_it + 1, ini_it + ntime
       write(fileid,*) it
       fileid = trim(adjustl(fileid))

       fparam_les = trim("restart_param.dat"//fileid)
       fdata_les = trim("restart.h5"//fileid)
       faux =  trim("restart_aux.h5"//fileid)
       fgrid = trim("grid.h5"//fileid)

       fparam_hos = trim("restart_param_hos.dat"//fileid)
       fdata_hos = trim("restart_hos.h5"//fileid)

       call reread(time, ioutd, ioutc, fgrid, fparam_les, fdata_les, faux)
       call grid_gen

       !need to add hos initialization subroutines below
       if(myid==0) print *, 'HOS initialization started!'

       if (it == ini_it + 1) then
          call hos_init
       end if
       call read_hos(time, ioutd, ioutc, eta_hos, vps_hos, pa0_hos, fparam_hos, fdata_hos)

       call wavenum(wvn_hos)
       call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
       call zeta(eta_hos,zp_hos)
       tmp=1.0_wp
       call boundvp(vps_hos,r_hos,zp_hos)
       call wsurf(w_hos,r_hos,zp_hos,wvn_hos)
       call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)

       call bottom_hos_les(time)

       if(myid==0) print *, 'HOS initialized!'

       ! transfer data to the upper cpus
       call MPI_BCAST(eta,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hh,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(ht,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hxy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hxx,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)
       call MPI_BCAST(hyy,xsz(1)*xsz(2),mpi_double_precision,0,mpi_comm_2d_row, ierr)

       eta=hh+0

       call init_p_poisson
       call nl_coef

        !test starts here!
       root0 = 0
       call gather_3d_xyz(u(1:xsz(1),1:xsz(2),1:xsz(3)),uall,root0)

       call get_mean_flux_horiz(um,vm,uup,vvp,zzall)
      
       usig_nw(:,:,it-ini_it) = uall(:,:,2)
       usig_out(:,:,it-ini_it) = uall(:,:,40)

!       call turb_analysis
       if(myid==0) print *,"finish step ",it
       !test ends here!

    end do

    ! start calculating k-omg spectrum
    if (myid == 0) then

       open(878,file="space_time_raw_nw_out_u.dat",status='REPLACE')
       write(878,'(25e12.4)') usig_nw,usig_out
!       open(878,file="space_time_raw_nw_out_u.dat")
!       read(878,*) usig_nw,usig_out
       close(878)

       open(879,file="space_time_raw_mean_flux.dat",status='REPLACE')
       write(879,'(25e12.4)') um,vm,uup,vvp
!       open(879,file="space_time_raw_mean_flux.dat")
!       read(879,*) um,vm,uup,vvp
       close(879)

       tmppet = twopi / nomg / dt

       do k = 2, 40, 38
          if (myid == 0) print *,k
          if (k == 5) then
             call spec_3d(usig_nw,ef3d,k1co,k2co,k3co,pex, pey, tmppet, nx_global, ny_global, nomg)
          else
             call spec_3d(usig_out,ef3d,k1co,k2co,k3co,pex, pey, tmppet, nx_global, ny_global, nomg)
          end if
          do i = 1, nomg
             call spec_2d2kx(ef3d(:,:,i), e11_k1omg(:,i))
          end do

          write(84,*) "VARIABLES = kx, omg, e11_k1omg_num"
          write(84,"(A,F20.8,A,I5,A,I5)") 'ZONE T="',k*1.0,'" I=',nx_global/3, &
               ' J=',nomg
          write(84,"(A,E16.4,A)") ' AUXDATA um="',um(k),'"'
          write(84,"(A,E16.4,A)") ' AUXDATA z="',zzall(k),'"'

          do j=1, nomg
             do i=1, nx_global/3
                write(84,'(25e12.4)') k1co(i), k3co(j), e11_k1omg(i,j)
             end do
          end do
          
          if (k == 5) then
             call spec_kx_ky(usig_nw(:,:,1),e11_k1k2)
          else
             call spec_kx_ky(usig_out(:,:,1),e11_k1k2)
          end if
          !       call spec_2d2kx(e11_k1k2, e11_k1)       

          call spec_kx_ky_omg_lam(e11_k1k2,nomg,k3co,ekkomg,uup(k),vvp(k),um(k))

          do i = 1, nomg
             call spec_2d2kx(ekkomg(:,:,i), e11_k1omg(:,i))
          end do
          
          write(82,*) "VARIABLES = kx, ky, e11_k1k2_ana"
          write(82,"(A,F20.8,A,I5,A,I5)") 'ZONE T="',k*1.0,'" I=',nx_global/3, &
               ' J=',ny_global/3
          write(82,"(A,E16.4,A)") ' AUXDATA um="',um(k),'"'
          write(82,"(A,E16.4,A)") ' AUXDATA z="',zzall(k),'"'
          write(82,"(A,E16.4,A)") ' AUXDATA ustar="',usbot,'"'
          do j=1, ny_global/3
             do i=1, nx_global/3
                write(82,'(25e12.4)') k1co(i), k2co(j), e11_k1k2(i,j)
             end do
          end do

          write(83,*) "VARIABLES = kx, omg, e11_k1omg"
          write(83,"(A,F20.8,A,I5,A,I5)") 'ZONE T="',k*1.0,'" I=',nx_global/3, &
               ' J=',nomg
          write(83,"(A,E16.4,A)") ' AUXDATA um="',um(k),'"'
          write(83,"(A,E16.4,A)") ' AUXDATA z="',zzall(k),'"'
          write(83,"(A,E16.4,A)") ' AUXDATA ustar="',usbot,'"'
          do j=1, nomg
             do i=1, nx_global/3
                write(83,'(25e12.4)') k1co(i), k3co(j), e11_k1omg(i,j)
             end do
          end do
          
          write(85,*) "VARIABLES = kx, omg"
          write(85,"(A,F20.8,A,I5,A,I5)") 'ZONE T="',k*1.0,'" I=',nx_global/3
          write(85,"(A,E16.4,A)") ' AUXDATA um="',um(k),'"'
          write(85,"(A,E16.4,A)") ' AUXDATA z="',zzall(k),'"'
          write(85,"(A,E16.4,A)") ' AUXDATA ustar="',usbot,'"'
          do j=1, nx_global/3
             write(85,'(25e12.4)') k1co(i), k1co(j)*um(k)
          end do
       end do

    end if
    ! end calculating k-omg spectrum
    
    deallocate(etaall,hhall,tmp)

    deallocate(uall,vall,usig_nw, usig_out)
!    deallocate(ppall,etasig,pp_wave,ppsig,dtheta)

    ! call spectral_finalize
    call fft_finalize
    call decomp_finalize
    !HOS
    call fft_finalize_hos

    call mpi_finalize(ierror)
  end subroutine post_les_hos_k_omg_spectrum
  
  subroutine get_wave_coherent(f,f_wave,eta,n,dtheta)
    
    implicit none

    integer, intent(in) :: n
    real(wp), intent(in), dimension(:) :: f, eta
    real(wp), intent(out), dimension(:) :: f_wave
    real(wp), intent(out) :: dtheta

    integer i, j
    real(wp), allocatable, dimension(:) :: heta, deta, dheta, test
    real(wp) sum1, sum2, amp, tmp, alp, f_ave

    allocate(heta(n))
    allocate(deta(n))
    allocate(dheta(n))
    allocate(test(n))

    alp = 0.00001
    call hilbert(eta,heta,n)

    do i = 1, n - 1
       deta(i) = (eta(i+1)-eta(i))
       dheta(i) = (heta(i+1)-heta(i))
    end do
    deta(n) = deta(n-1)
    dheta(n) = dheta(n-1)

    f_ave = sum(f(1:n))/n
    do j = 1, n
       sum1 = 0
       sum2 = 0
       do i = 1, n
          amp = sqrt(deta(i)**2 + dheta(i)**2)
          tmp = sqrt((eta(i)-eta(j))**2 + (heta(i)-heta(j))**2)
          sum1 = sum1 + f(i) * amp * delta(tmp,alp)
          sum2 = sum2 + amp * delta(tmp,alp)
       end do

       f_wave(j) = sum1 / sum2 - f_ave
    end do

    call fft_for_hos(f_wave,test,n)
    test(n/4+1:n) = 0
    call fft_bac_hos(test,f_wave,n)

    sum1 = 0
    sum2 = 0
    do i = 1, n
       sum1 = sum1 + f_wave(i) * heta(i)
       sum2 = sum2 + f_wave(i) * eta(i)
    end do

    dtheta = atan2(sum1,sum2)

    if (dtheta <= 0) then
       dtheta = dtheta + twopi
    end if

    deallocate(heta,deta,dheta,test)

  end subroutine get_wave_coherent

  subroutine get_wave_coherent_linft(f,f_wave,eta,n,dtheta)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in), dimension(:) :: f, eta
    real(wp), intent(out), dimension(:) :: f_wave
    real(wp), intent(out) :: dtheta

    integer i, j
    real(wp), allocatable, dimension(:) :: heta, deta, dheta, tmp
    real(wp), allocatable, dimension(:,:) :: etak,hetak
    real(wp) sum1, sum2, norm, f_ave

    allocate(heta(n))
    allocate(deta(n))
    allocate(dheta(n))
    allocate(tmp(n))
    allocate(etak(n,n/2))
    allocate(hetak(n,n/2))

    etak = 0
    hetak = 0
    call hilbert(eta,heta,n)

    call fft_for_hos(eta,tmp,n)
    do j = 1, n/2-1
       etak(j*2+1,j) = tmp(j*2+1)
       etak(j*2+2,j) = tmp(j*2+2)
    end do

    do j = 1, n/2-1
       tmp = etak(:,j)
       call fft_bac_hos(tmp,etak(:,j),n)
    end do

    call fft_for_hos(heta,tmp,n)
    do j = 1, n/2-1
       hetak(j*2+1,j) = tmp(j*2+1)
       hetak(j*2+2,j) = tmp(j*2+2)
    end do

    do j = 1, n/2-1
       tmp = hetak(:,j)
       call fft_bac_hos(tmp,hetak(:,j),n)
    end do

    f_wave = 0
    do j = 1, n/2-1
       sum1 = sum(f*etak(:,j))
       sum2 = sum(f*hetak(:,j))
       norm = (sum(etak(:,j)*etak(:,j)) + sum(hetak(:,j)*hetak(:,j)))
       f_wave = f_wave + (sum1 * etak(:,j) + sum2 * hetak(:,j)) / norm
    end do

    sum1 = 0
    sum2 = 0
    do i = 1, n
       sum1 = sum1 + f_wave(i) * heta(i)
       sum2 = sum2 + f_wave(i) * eta(i)
    end do

    dtheta = atan2(sum1,sum2)

    if (dtheta <= 0) then
       dtheta = dtheta + twopi
    end if

    deallocate(heta,deta,dheta,tmp)

  end subroutine get_wave_coherent_linft


 !-------------------------------------------------------------------------

  subroutine hilbert(f,hf,n)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in), dimension(:) :: f
    real(wp), intent(out), dimension(:) :: hf

    integer l
    real(wp), allocatable, dimension(:) :: tmp

    allocate(tmp(n))

    call fft_for_hos(f,tmp,n)

    hf = 0
    do l = 3, n - 1, 2
       hf(l) = tmp(l+1)
       hf(l+1) = -tmp(l)
    end do

    tmp = hf(1:n)

    call fft_bac_hos(tmp,hf,n)

    deallocate(tmp)

  end subroutine hilbert

  !-------------------------------------------------------------------------

  subroutine get_inst(f,hf,iamp,iphs,n)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in), dimension(:) :: f, hf
    real(wp), intent(out), dimension(:) :: iamp, iphs

    integer i

    do i = 1, n
       iamp(i) = sqrt(f(i)**2 + hf(i)**2)
       iphs(i) = atan(hf(i)/f(i))
    end do

  end subroutine get_inst

  !-------------------------------------------------------------------------

  real(wp) function delta(x,alpha)

    implicit none

    real(wp), intent(in) :: x, alpha

    real(wp) tmp
    if (alpha <= 0) then
       print *,"Wrong alpha value!"
       stop
    end if

    tmp = -x**2/alpha**2

    if (tmp <= -600) then
       delta = 0
    else
       delta = (1/alpha/sqrt(0.5*twopi)) * exp(tmp)
    end if

  end function delta

  subroutine get_quadrant(q1,q2,q3,q4,uupall,wwpall)

    use decomp
    use param
    use grid, only : zz
    use navier
    use utils

    implicit none

    real(wp), intent(in), dimension(:,:,:) :: uupall, wwpall
    real(wp), intent(out), dimension(:) :: q1, q2, q3, q4

    integer i,j,k
    integer num1, num2, num3, num4

    q1 = 0
    q2 = 0
    q3 = 0
    q4 = 0
    num1 = 0
    num2 = 0
    num3 = 0
    num4 = 0
    do k = 1, nz_global
       do j = 1, ny_global
          do i = 1, nx_global
             if (uupall(i,j,k) > 0 .and. wwpall(i,j,k) > 0) then
                q1(k) = q1(k) + uupall(i,j,k) * wwpall(i,j,k)
                num1 = num1 + 1
             end if
             if (uupall(i,j,k) < 0 .and. wwpall(i,j,k) > 0) then
                q2(k) = q2(k) + uupall(i,j,k) * wwpall(i,j,k)
                num2 = num2 + 1
             end if
             if (uupall(i,j,k) < 0 .and. wwpall(i,j,k) < 0) then
                q3(k) = q3(k) + uupall(i,j,k) * wwpall(i,j,k)
                num3 = num3 + 1
             end if
             if (uupall(i,j,k) > 0 .and. wwpall(i,j,k) < 0) then
                q4(k) = q4(k) + uupall(i,j,k) * wwpall(i,j,k)
                num4 = num4 + 1
             end if
          end do
       end do
    end do

    q1 = q1 / num1
    q2 = q2 / num2
    q3 = q3 / num3
    q4 = q4 / num4
  end subroutine get_quadrant

  subroutine get_mean_flux_horiz(um,vm,uup,vvp,zzall)

    use decomp
    use param
    use grid, only : zz
    use navier
    use utils
    use MPI
    implicit none

    real(wp), intent(out), dimension(:) :: um,vm,uup,vvp,zzall
!    real(wp), intent(in), dimension(:,:,:) :: uall,wall

    real(wp) :: z, time
    real(wp), dimension(nx_global, nz_global) :: umy, vmy, wmy, wmy1, pmy
    real(wp), dimension(nx_global, nz_global) :: uf2my, vf2my, wf2my
    !      real(wp), dimension(ny_global, nz_global) :: umx, vmx, wmx, wmx1, pmx
    real(wp), dimension(1:xsz(1),1:xsz(2),1:xsz(3)) :: tmpf
    real(wp), dimension(nz_global) :: tmpz
    integer :: i, j, k, k1, k2, ierror

    ! Gather zz coordinates
    tmpz = 0
    tmpz(xst(3):xend(3)) = zz(1:xsz(3))
    call MPI_Allreduce(tmpz, zzall, nz_global, real_type, &
         MPI_SUM, MPI_COMM_2D_ROW, ierror)

    ! Compute average over y
    call get_average_y(u(:,:,1:xsz(3)),  umy)
    call get_average_y(v(:,:,1:xsz(3)),  vmy)
    call get_average_y(w(:,:,1:xsz(3)),  wmy)
    call get_average_y(pp(:,:,1:xsz(3)), pmy)

    ! Interpolate w
    wmy1(:,nz_global) = wmy(:,nz_global-1)
    do k=nz_global-1, 2, -1
       wmy1(:,k) = (wmy(:,k)+wmy(:,k-1))/2
    end do
    wmy1(:,1) = wmy(:,1)

    !-----------------------
    ! Compute plane average
    !-----------------------
    do k=1, nz_global
       do i=2, nx_global
          umy(1,k) = umy(1,k) + umy(i,k)
          vmy(1,k) = vmy(1,k) + vmy(i,k)
          wmy1(1,k) = wmy1(1,k) + wmy1(i,k)
          pmy(1,k) = pmy(1,k) + pmy(i,k)
       end do
       umy(1,k) = umy(1,k)/nx_global
       vmy(1,k) = vmy(1,k)/nx_global
       wmy1(1,k) = wmy1(1,k)/nx_global
       pmy(1,k) = pmy(1,k)/nx_global
    end do

    !---------------------------------------------------
    ! Compute the fluctuation relative to plane average
    !---------------------------------------------------
    do k=1, xsz(3)
       do j=1, xsz(2)
          do i=1, xsz(1)
             tmpf(i,j,k) = (u(i,j,k)-umy(1,k+xst(3)-1))
          end do
       end do
    end do
    tmpf = tmpf**2
    call get_average_y(tmpf, uf2my)

    do k=1, xsz(3)
       do j=1, xsz(2)
          do i=1, xsz(1)
             tmpf(i,j,k) = (v(i,j,k)-vmy(1,k+xst(3)-1))
          end do
       end do
    end do
    tmpf = tmpf**2
    call get_average_y(tmpf, vf2my)

    !---------------------
    ! Write plane average
    !---------------------
    um = 0
    vm = 0
    uf = 0
    vf = 0
    if (myid==0) then
       do k=1, nz_global
          do i=2, nx_global
             uf2my(1,k) = uf2my(1,k) + uf2my(i,k)
             vf2my(1,k) = vf2my(1,k) + vf2my(i,k)
          end do
          uf2my(1,k) = uf2my(1,k)/nx_global
          vf2my(1,k) = vf2my(1,k)/nx_global
       end do

       um = umy(1,:)
       vm = vmy(1,:)
       uup = sqrt(uf2my(1,:))
       vvp = sqrt(vf2my(1,:))
    end if

  end subroutine get_mean_flux_horiz


  subroutine get_mean_flux(um,vm,wm,uup,vvp,wwp,zzall,uupall,wwpall)

    use decomp
    use param
    use grid, only : zz
    use navier
    use utils
    use MPI
    implicit none
    
    real(wp), intent(out), dimension(:) :: um,vm,wm,uup,vvp,wwp,zzall
    real(wp), intent(out), dimension(:,:,:) :: uupall,wwpall

    real(wp) :: z, time
    real(wp), dimension(nx_global, nz_global) :: umy, vmy, wmy, wmy1, pmy
    real(wp), dimension(nx_global, nz_global) :: uf2my, vf2my, wf2my
    !      real(wp), dimension(ny_global, nz_global) :: umx, vmx, wmx, wmx1, pmx
    real(wp), dimension(1:xsz(1),1:xsz(2),1:xsz(3)) :: tmpf
    real(wp), dimension(nz_global) :: tmpz
    integer :: i, j, k, k1, k2, ierror
    
    ! Gather zz coordinates
    tmpz = 0
    tmpz(xst(3):xend(3)) = zz(1:xsz(3))
    call MPI_Allreduce(tmpz, zzall, nz_global, real_type, &
         MPI_SUM, MPI_COMM_2D_ROW, ierror)
    
    ! Compute average over y
    call get_average_y(u(:,:,1:xsz(3)),  umy)
    call get_average_y(v(:,:,1:xsz(3)),  vmy)
    call get_average_y(w(:,:,1:xsz(3)),  wmy)
    call get_average_y(pp(:,:,1:xsz(3)), pmy)
    
    ! Interpolate w
    wmy1(:,nz_global) = wmy(:,nz_global-1)
    do k=nz_global-1, 2, -1
       wmy1(:,k) = (wmy(:,k)+wmy(:,k-1))/2
    end do
    wmy1(:,1) = wmy(:,1)
    
    !-----------------------
    ! Compute plane average
    !-----------------------
    do k=1, nz_global
       do i=2, nx_global
          umy(1,k) = umy(1,k) + umy(i,k)
          vmy(1,k) = vmy(1,k) + vmy(i,k)
          wmy1(1,k) = wmy1(1,k) + wmy1(i,k)
          pmy(1,k) = pmy(1,k) + pmy(i,k)
       end do
       umy(1,k) = umy(1,k)/nx_global
       vmy(1,k) = vmy(1,k)/nx_global
       wmy1(1,k) = wmy1(1,k)/nx_global
       pmy(1,k) = pmy(1,k)/nx_global
    end do
    
    !---------------------------------------------------
    ! Compute the fluctuation relative to plane average
    !---------------------------------------------------
    do k=1, xsz(3)
       do j=1, xsz(2)
          do i=1, xsz(1)
             tmpf(i,j,k) = (u(i,j,k)-umy(1,k+xst(3)-1))
          end do
       end do
    end do
    call gather_3d_xyz(tmpf(1:xsz(1),1:xsz(2),1:xsz(3)),uupall,0)
    tmpf = tmpf**2
    call get_average_y(tmpf, uf2my)
    
    do k=1, xsz(3)
       do j=1, xsz(2)
          do i=1, xsz(1)
             tmpf(i,j,k) = (v(i,j,k)-vmy(1,k+xst(3)-1))
          end do
       end do
    end do
    tmpf = tmpf**2
    call get_average_y(tmpf, vf2my)
    
    k1 = 1; k2 = xsz(3)
    if (istop) k2 = xsz(3)-1
    if (isbot) k1 = 2
    do k=k1, k2
       do j=1, xsz(2)
          do i=1, xsz(1)
             tmpf(i,j,k) = ( (w(i,j,k)  -wmy(1,k+xst(3)-1)) &
                  + (w(i,j,k-1)-wmy(1,k+xst(3)-2)) ) / 2
          end do
       end do
    end do
    if (istop) then
       do j=1, xsz(2)
          tmpf(:,j,xsz(3)) = w(:,j,xsz(3)-1)-wmy(1,nz_global-1)
       end do
    end if
    if (isbot) then
       tmpf(:,:,1) = w(:,:,1) - wmy(1,1)
    end if
    call gather_3d_xyz(tmpf(1:xsz(1),1:xsz(2),1:xsz(3)),wwpall,0)    
    tmpf = tmpf**2
    call get_average_y(tmpf, wf2my)
    
    !---------------------
    ! Write plane average
    !---------------------
    um = 0
    vm = 0
    wm = 0
    uf = 0
    vf = 0
    wf = 0
    if (myid==0) then
       do k=1, nz_global
          do i=2, nx_global
             uf2my(1,k) = uf2my(1,k) + uf2my(i,k)
             vf2my(1,k) = vf2my(1,k) + vf2my(i,k)
             wf2my(1,k) = wf2my(1,k) + wf2my(i,k)
          end do
          uf2my(1,k) = uf2my(1,k)/nx_global
          vf2my(1,k) = vf2my(1,k)/nx_global
          wf2my(1,k) = wf2my(1,k)/nx_global
       end do

       um = umy(1,:)
       vm = vmy(1,:)
       wm = wmy(1,:)
       uup = sqrt(uf2my(1,:))
       vvp = sqrt(vf2my(1,:))
       wwp = sqrt(wf2my(1,:))
    end if
    
  end subroutine get_mean_flux
  
  subroutine get_1dspectrum(f,sk,dwk,nk)

    use decomp
    use param
    use fft_hos
    implicit none
    
    integer, intent(in) :: nk
    real(wp), intent(in), dimension(:,:) :: f
    real(wp), intent(out) :: dwk
    real(wp), intent(out), dimension(:) :: sk
    
    real(wp), allocatable, dimension(:,:) :: tmp
    real(wp) fkn, wkx, wky, wav, con
    integer l,m,n,lp1,mp1,root
    
    allocate(tmp(nx_global,ny_global))
    
    sk = 0
    dwk = sqrt((pex*(nx_global/3))**2+(pey*(ny_global/3))**2) / nk
    call fft_for_xy_hos(f,tmp,nx_global,ny_global)
    do lp1 = 1, nx_global / 3
       l = lp1 * 2 - 1
       wkx=(l-1)/2*pex
       do mp1 = 1, ny_global / 3
          m = mp1 * 2 - 1
          wky=(m-1)/2*pey
          wav=sqrt(wkx**2+wky**2)
          do n=1,nk
             if(wav >= ((n-0.5)*dwk).and.wav < ((n+0.5)*dwk))then
                exit
             endif
          enddo

          fkn=tmp(l,m)**2+tmp(l,m+1)**2+tmp(l+1,m)**2+tmp(l+1,m+1)**2
          con=4.
          if(l.eq.1.or.m.eq.1) con=2.
          if(l.eq.1.and.m.eq.1) con=1.
          if (n <= nk) sk(n)=sk(n)+con*fkn
       enddo
    enddo
    
    do n=1,nk
       sk(n)=sk(n)/dwk
    enddo
    deallocate(tmp)
    
  end subroutine get_1dspectrum

  subroutine spec_2d2kx(ek1k2, ek1)

    use decomp
    use param
    use fft_hos
    
    implicit none
    
    real(wp), intent(in), dimension(:,:) :: ek1k2
    real(wp), intent(out), dimension(:) :: ek1

    integer l,m

    ek1 = 0
    do l = 1, nx_global
       ek1(l) = sum(ek1k2(l,:)) * pey
    end do

  end subroutine spec_2d2kx

  subroutine spec_3d(f,ef3d,k1co,k2co,k3co,dk1, dk2, dk3, n1, n2, n3)

    ! Linear Advection Model

    use decomp
    use param
    use fft_hos

    implicit none

    integer, intent(in) :: n1, n2, n3
    real(wp), intent(in) :: dk1, dk2, dk3
    real(wp), intent(in), dimension(:,:,:) :: f
    real(wp), intent(out), dimension(:,:,:) :: ef3d
    real(wp), intent(out), dimension(:) :: k1co,k2co,k3co
    
    real(wp), allocatable, dimension(:,:,:) :: tmp1,tmp2

    integer i,j,k,l,m,n
    real(wp) a1, a2, a3, a4, con

    allocate(tmp1(n1,n2,n3))
    allocate(tmp2(n2,n3,n1))

    !Begin 3D fft
    do j = 1, n2
       do k = 1, n3
          call fft_for_hos(f(:,j,k),tmp1(:,j,k),n1)
       end do
    end do

    do i = 1, n1
       do j = 1, n2
          do k = 1, n3
             tmp2(j,k,i) = tmp1(i,j,k)
          end do
       end do
    end do

    do i = 1, n1
       call fft_for_xy_hos(tmp2(:,:,i),n2,n3)
    end do

    do i = 1, n1
       do j = 1, n2
          do k = 1, n3
             tmp1(i,j,k) = tmp2(j,k,i)
          end do
       end do
    end do
    !End 3D fft

    k1co = 0
    do l = 1, n1/2
       k1co(l) = (l-1) * dk1 
    end do

    k2co = 0
    do m = 1, n2/2 
       k2co(m) = (m-1) * dk2 
    end do

    k3co = 0
    do n = 1, n3
       k3co(n) = (n-n3/2) * dk3 
    end do


    ef3d = 0
    do l = 1, n1, 2
       do m = 1, n2, 2
          do n = 1, n3, 2
             if (l+m+n > 3 .and. (n-1)/2 <= n3/3) then
                !+k1,+k2,-omg
                a1 = 2 * ((tmp1(l,m,n) - tmp1(l,m+1,n+1) - tmp1(l+1,m,n+1) - tmp1(l+1,m+1,n))**2 &
                     + (tmp1(l+1,m,n) -  tmp1(l+1,m+1,n+1) + tmp1(l,m,n+1) + tmp1(l,m+1,n))**2)**0.5
                !+k1,-k2,-omg
                a2 = 2 * ((tmp1(l,m,n) + tmp1(l,m+1,n+1) - tmp1(l+1,m,n+1) + tmp1(l+1,m+1,n))**2 &
                     + (tmp1(l+1,m,n) +  tmp1(l+1,m+1,n+1) + tmp1(l,m,n+1) - tmp1(l,m+1,n))**2)**0.5
                !+k1,+k2,omg
                a3 = 2 * ((tmp1(l,m,n) + tmp1(l,m+1,n+1) + tmp1(l+1,m,n+1) - tmp1(l+1,m+1,n))**2 &
                     + (tmp1(l+1,m,n) + tmp1(l+1,m+1,n+1) - tmp1(l,m,n+1) + tmp1(l,m+1,n))**2)**0.5
                !+k1,-k2,omg
                a4 = 2 * ((tmp1(l,m,n) - tmp1(l,m+1,n+1) + tmp1(l+1,m,n+1) + tmp1(l+1,m+1,n))**2 &
                     + (tmp1(l+1,m,n) -  tmp1(l+1,m+1,n+1) - tmp1(l,m,n+1) - tmp1(l,m+1,n))**2)**0.5
                
                if (l == 1) then
                   if (m == 1) then
                   end if
                   if (n == 1) then
                   end if
                else 
                   if (m == 1) then
                   end if
                end if
                !+omg
                ef3d((l+1)/2,(m+1)/2,n3/2+(n-1)/2) = ef3d((l+1)/2,(m+1)/2,n3/2+(n-1)/2) + 0.5 * (a3**2 + a4**2) / dk1 / dk2 / dk3
                !-omg
                ef3d((l+1)/2,(m+1)/2,n3/2-(n-1)/2) = ef3d((l+1)/2,(m+1)/2,n3/2-(n-1)/2) + 0.5 * (a1**2 + a2**2) / dk1 / dk2 / dk3

                
                if (n == 1) then
                   ! omg = 0 
                   if (l > 1 .and. m > 1) con = 4.0
                   if (l == 1 .or. m == 1) con = 2.0
!                   ef3d((l+1)/2,(m+1)/2,n3/2) = con * (tmp1(l,m,n3/2)**2 + tmp1(l+1,m+1,n3/2)**2 + tmp1(l+1,m,n3/2)**2 + tmp1(l,m+1,n3/2)**2) / dk1 / dk2 / dk3
                   ef3d((l+1)/2,(m+1)/2,n3/2) = 0.5 * (ef3d((l+1)/2,(m+1)/2,n3/2-1) + ef3d((l+1)/2,(m+1)/2,n3/2+1))
                end if
             end if
          end do
       end do
    end do
  end subroutine spec_3d

  subroutine spec_kx_ky_omg_lam(ek1k2,nomg,allomg,ekkomg,urms,vrms,um)
    
    ! Linear Advection Model 

    use decomp
    use param
    use fft_hos

    implicit none

    integer, intent(in) :: nomg
    real(wp), intent(in) :: urms, vrms, um
    real(wp), intent(in), dimension(:) :: allomg
    real(wp), intent(in), dimension(:,:) :: ek1k2
    real(wp), intent(out), dimension(:,:,:) :: ekkomg

    integer k1,k2,m
    real(wp) wk1, wk2, tmp

    do m = 1, nomg
       do k1 = 0, nx_global/2-1
          wk1 = pex * k1
          do k2 = 0, ny_global/2-1
             wk2 = pey * k2
             tmp = -1 * (allomg(m) - wk1 * um)**2 / 2 / (urms**2*wk1**2+vrms**2*wk2**2)
             if (tmp >= -50) then
                ekkomg(k1+1,k2+1,m) = ek1k2(k1+1,k2+1) * exp(tmp) / sqrt(twopi*(urms**2*wk1**2+vrms**2*wk2**2))
             end if
          end do
       end do
    end do

  end subroutine spec_kx_ky_omg_lam

  subroutine spec_kx_ky(f,ef)

    use decomp
    use param
    use fft_hos

    implicit none

    real(wp), intent(in), dimension(:,:) :: f
    real(wp), intent(out), dimension(:,:) :: ef

    integer l,m, kx, ky
    real(wp) con
    real(wp), allocatable, dimension(:,:) :: tmp

    allocate(tmp(nx_global, ny_global))
    call fft_for_xy_hos(f,tmp,nx_global,ny_global)

    ef = 0
    
    do l = 1, nx_global - 1, 2
       kx = (l - 1) / 2
       do m = 1, ny_global - 1, 2
          ky = (m - 1) / 2
          ef(kx+1,ky+1) = tmp(l,m)**2+tmp(l,m+1)**2+tmp(l+1,m)**2+tmp(l+1,m+1)**2
          con = 4.0
          if (kx == 1 .or. ky == 1) con = 2.0
          if (kx == 1 .and. ky == 1) con = 1.0
          ef(kx+1,ky+1) = ef(kx+1,ky+1) * con / pex / pey
       end do
    end do
    deallocate(tmp)

  end subroutine spec_kx_ky
  
  subroutine spec_theta_kx_ky(f,ef,thetaf)

    use decomp
    use param
    use fft_hos

    implicit none

    real(wp), intent(in), dimension(:,:) :: f
    real(wp), intent(out), dimension(:,:) :: ef, thetaf

    integer l,m, kx, ky
    real(wp) con
    real(wp), allocatable, dimension(:,:) :: tmp

    allocate(tmp(nx_global, ny_global))
    call fft_for_xy_hos(f,tmp,nx_global,ny_global)

    ef = 0
    thetaf = 0.0
    
    do l = 1, nx_global - 1, 2
       kx = (l - 1) / 2
       do m = 1, ny_global - 1, 2
          ky = (m - 1) / 2
          ef(kx+1,ky+1) = tmp(l,m)**2+tmp(l,m+1)**2+tmp(l+1,m)**2+tmp(l+1,m+1)**2
          thetaf(kx+1,ky+1) = atan2(tmp(l,m+1)+tmp(l+1,m), tmp(l,m)-tmp(l+1,m+1))
          con = 4.0
          if (kx == 1 .or. ky == 1) con = 2.0
          if (kx == 1 .and. ky == 1) con = 1.0
          ef(kx+1,ky+1) = ef(kx+1,ky+1) * con / pex / pey
       end do
    end do
    deallocate(tmp)

  end subroutine spec_theta_kx_ky
    
end module post_proc
