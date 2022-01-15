module io
   use decomp
   use hdf_io

   implicit none

   public :: saverestart
   public :: reread
   public :: input_iht
   public :: output_all
   public :: output_xy_average
   public :: output_cut
   public :: output_surface
   public :: output_plane
   public :: output_eta

   interface reread
      module procedure reread_no_name, reread_file_name
   end interface reread
      
   interface saverestart
      module procedure save_no_name, save_file_name
   end interface saverestart

   !> added by plyu
   public :: output_all_2
   public :: reread_output

 contains

   !checked
   subroutine save_no_name(time, ioutd, ioutc)
      use param
      use grid, only : zz, zw, dz, dzw, eta, eta0, hh, ht
      use navier
      use wavecontrol, only : at1
      use turbine_model, only: nxwt, nywt, ud, ud1, ud_ref, ud1_ref,&
        NumberOfTurbines, wtm, fsi_wt, NumberOfNacelle, fsi_nac, &
        fsitype, pp_wt, nacelle_model !> added by plyu
      implicit none

      integer,  intent(IN) :: ioutd, ioutc
      real(wp), intent(IN) :: time
      type(HDFObj) :: writer
      integer :: i, j, ibi ! added by plyu

      call IO_Init

      if (myid.eq.0) then
         open(14000, file="restart_param.dat")
         write(14000,*) time, ioutd, ioutc
         close(14000)
      end if

      !> added by plyu: start
      if (myid .eq. 0 .and. iturbine .eq. 1) then
        open(15010, file='restart_ad.dat')
        do i = 1, nxwt
          do j = 1, nywt
            write(15010,*) ud1((i-1)*nywt+j), ud((i-1)*nywt+j), &
              ud1_ref((i-1)*nywt+j), ud_ref((i-1)*nywt+j)
          enddo
        enddo
        close(15010)
      endif
      
      if (myid .eq. 0 .and. iturbine .eq. 3) then
        open(15010, file='restart_uref.dat')
        do ibi = 1, NumberOfTurbines
          write(15010, *) wtm(ibi)%U1_ref, wtm(ibi)%U_ref, &
            ' # U1_ref, U_ref for turbine ', ibi
        enddo
        close(15010)
      endif

      if (myid .eq. 0 .and. iturbine .gt. 1) then
        open(15012, file='Turbine2.inp', action='write', &
          status = 'replace')
        do ibi = 1, NumberOfTurbines
          !> yaw angle to be modified
          write(15012, *) fsi_wt(ibi)%ang_axis, fsi_wt(ibi)%ang0_yaw, &
            fsi_wt(ibi)%angvel_axis, ' # ang_axis, ang0_yaw, angvel_axis', &
            ' for turbine ', ibi 
        enddo

        close(15012)

        if (fsitype.eq.1) then
          open(15013, file='FSI2_rotor.inp', action='write', status='replace')
          do ibi = 1, NumberOfTurbines
            write(15013,*) fsi_wt(ibi)%x_c, fsi_wt(ibi)%y_c, fsi_wt(ibi)%z_c, &
              fsi_wt(ibi)%nx_tb, fsi_wt(ibi)%ny_tb, fsi_wt(ibi)%nz_tb
          enddo
          close(15013)
          if (nacelle_model .ne. 0) then
            open(15014, file='FSI2_nacelle.inp', action='write', status='replace')
            do ibi = 1, NumberOfNacelle
              write(15014,*) fsi_nac(ibi)%x_c, fsi_nac(ibi)%y_c, &
                fsi_nac(ibi)%z_c, fsi_nac(ibi)%nx_tb, fsi_nac(ibi)%ny_tb, &
                fsi_nac(ibi)%nz_tb
            enddo
            close(15014)
          endif  
        endif
      endif
      !> added by plyu: end

      if (isbot .and. myid1==0 .and. iwcontrol==1) then
         open(14000, file="restart_wave.dat")
         write(14000,*) at1
         close(14000)
      end if

      call writerOpen("grid.h5", MPI_COMM_2D_CART, writer)
      call write1D_z(writer, "zz", zz(1:xsz(3)), 0)
      call write1D_z(writer, "zw", zw(1:xsz(3)), 0)
      call write1D_z(writer, "dz", dz(1:xsz(3)), 0)
      call write1D_z(writer, "dzw", dzw(1:xsz(3)), 0)
      call writerClose(writer)

      call writerOpen("restart.h5", MPI_COMM_2D_CART, writer)
      !call write2D_xy(writer, "uzfs", uzfs, topid)
      !call write2D_xy(writer, "vzfs", vzfs, topid)
      !call write2D_xy(writer, "wzfs", wzfs, topid)
      
      !> plyunote: following terms are zeros, so I comment them out.
      !> plyunote: uzfs are used as freesurface velocity, esp in
      !!           wave-above-water case.
      !call write2D_xy(writer, "uzfs", uzfs, 0)
      !call write2D_xy(writer, "vzfs", vzfs, 0)
      !call write2D_xy(writer, "wzfs", wzfs, 0)
      call write2D_xy(writer, "eta", eta, 0)
      
      !> plyunote: eta0 is zero, so I comment it out.
      !call write2D_xy(writer, "eta0", eta0, 0)
      call write2D_xy(writer, "hh", hh, 0)
      
      call write3D_xyz(writer, "u", u(:,:,1:xsz(3)))
      call write3D_xyz(writer, "v", v(:,:,1:xsz(3)))
      call write3D_xyz(writer, "w", w(:,:,1:xsz(3)))
      call write3D_xyz(writer, "pp", pp(:,:,1:xsz(3)))
      call writerClose(writer)

!      print *,"myid=",myid

      call writerOpen("restart_aux.h5", MPI_COMM_2D_CART, writer)
      call write3D_xyz(writer, "hu", hu(:,:,1:xsz(3)))
      call write3D_xyz(writer, "hv", hv(:,:,1:xsz(3)))
      call write3D_xyz(writer, "hw", hw(:,:,1:xsz(3)))
      call write2D_xy(writer, "ht", ht, 0)
      call write2D_xy(writer, "u1", u(:,:,xsz(3)+1), topid)
      !call write2D_xy(writer, "u1", u(:,:,xsz(3)+1), 0)
      call write2D_xy(writer, "v1", v(:,:,xsz(3)+1), topid)
      !call write2D_xy(writer, "v1", v(:,:,xsz(3)+1), 0)
      call writerClose(writer)

!      print *,"myid=",myid

      call IO_Finalize

    end subroutine save_no_name

!checked
   subroutine reread_no_name(time, ioutd, ioutc)
      use param, only: iwcontrol
      use grid, only : zz, zw, dz, dzw, eta, eta0, hh, ht
      use navier
      use wavecontrol, only : at1
      implicit none

      integer,  intent(OUT) :: ioutd, ioutc
      real(wp), intent(OUT) :: time
      real(wp), dimension(xsz(1),xsz(2)) :: u1, v1
      type(HDFObj) :: writer

      ! plyunote: switch for h5 read scheme
      integer :: ih5s

      call IO_Init

      open(14000, file="restart_param.dat", action="READ", status="OLD")
      read(14000,*) time, ioutd, ioutc
      close(14000)

      if (isbot .and. iwcontrol==1) then
         open(14000, file="restart_wave.dat", action="READ", status="OLD")
         read(14000,*) at1
         close(14000)
      end if

      call readerOpen("grid.h5", MPI_COMM_2D_CART, writer)
      call read1D_z(writer, "zz", zz(1:xsz(3)), -1)
      call read1D_z(writer, "zw", zw(1:xsz(3)), -1)
      call read1D_z(writer, "dz", dz(1:xsz(3)), -1)
      call read1D_z(writer, "dzw", dzw(1:xsz(3)), -1)
      call readerClose(writer)

      ih5s = 2 
      if (ih5s .eq. 1) then
        ! plyunote: old way
        call readerOpen("restart.h5", MPI_COMM_2D_CART, writer)
        call read2D_xy(writer, "uzfs", uzfs, -1)
        call read2D_xy(writer, "vzfs", vzfs, -1)
        call read2D_xy(writer, "wzfs", wzfs, -1)
        call read2D_xy(writer, "eta", eta, -1)
        call read2D_xy(writer, "eta0", eta0, -1)
        call read2D_xy(writer, "hh", hh, -1)
        call read2D_xy(writer, "ht", ht, -1)
        call read3D_xyz(writer, "u", u(:,:,1:xsz(3)))
        call read2D_xy(writer, "u1", u1(:,:), topid)
        !call read2D_xy(writer, "u1", u1(:,:), 0)
        call read3D_xyz(writer, "v", v(:,:,1:xsz(3)))
        call read2D_xy(writer, "v1", v1(:,:), topid)
        !call read2D_xy(writer, "v1", v1(:,:), 0)
        call read3D_xyz(writer, "w", w(:,:,1:xsz(3)))
        call read3D_xyz(writer, "pp", pp(:,:,1:xsz(3)))
        call readerClose(writer)

        call readerOpen("restart_aux.h5", MPI_COMM_2D_CART, writer)
        call read3D_xyz(writer, "hu", hu(:,:,1:xsz(3)))
        call read3D_xyz(writer, "hv", hv(:,:,1:xsz(3)))
        call read3D_xyz(writer, "hw", hw(:,:,1:xsz(3)))
        !call read2D_xy(writer, "ht", ht, -1)
        !call read2D_xy(writer, "u1", u1(:,:), topid)
        !call read2D_xy(writer, "v1", v1(:,:), topid)
        call readerClose(writer)
      else if (ih5s .eq. 2) then      
        ! plyunote: new way
        call readerOpen("restart.h5", MPI_COMM_2D_CART, writer)
        !call read2D_xy(writer, "uzfs", uzfs, -1)
        !call read2D_xy(writer, "vzfs", vzfs, -1)
        !call read2D_xy(writer, "wzfs", wzfs, -1)
        call read2D_xy(writer, "eta", eta, -1)
        !call read2D_xy(writer, "eta0", eta0, -1)
        call read2D_xy(writer, "hh", hh, -1)
        !call read2D_xy(writer, "ht", ht, -1)
        call read3D_xyz(writer, "u", u(:,:,1:xsz(3)))
        !call read2D_xy(writer, "u1", u1(:,:), topid)
        !call read2D_xy(writer, "u1", u1(:,:), 0)
        call read3D_xyz(writer, "v", v(:,:,1:xsz(3)))
        !call read2D_xy(writer, "v1", v1(:,:), topid)
        !call read2D_xy(writer, "v1", v1(:,:), 0)
        call read3D_xyz(writer, "w", w(:,:,1:xsz(3)))
        call read3D_xyz(writer, "pp", pp(:,:,1:xsz(3)))
        call readerClose(writer)

        call readerOpen("restart_aux.h5", MPI_COMM_2D_CART, writer)
        call read3D_xyz(writer, "hu", hu(:,:,1:xsz(3)))
        call read3D_xyz(writer, "hv", hv(:,:,1:xsz(3)))
        call read3D_xyz(writer, "hw", hw(:,:,1:xsz(3)))
        call read2D_xy(writer, "ht", ht, -1)
        call read2D_xy(writer, "u1", u1(:,:), topid)
        call read2D_xy(writer, "v1", v1(:,:), topid)
        call readerClose(writer)
      endif 

      call IO_Finalize

      call update_ghost(zz,1)
      call update_ghost(dz,1)
      call update_ghost(zw,1)
      call update_ghost(dzw,1)
      if (istop) zz(xsz(3)+1) = zz(xsz(3))+dz(xsz(3))
      call update_ghost(u,1)
      call update_ghost(v,1)
      call update_ghost(w,1)
      if (istop) then
         u(:,:,xsz(3)+1) = u1(:,:)
         v(:,:,xsz(3)+1) = v1(:,:)
      end if
      call update_ghost(pp,1)

    end subroutine reread_no_name

   !checked
   subroutine save_file_name(time, ioutd, ioutc, fgrid, fparam, fdata, faux)
      use param
      use grid, only : zz, zw, dz, dzw, eta, eta0, hh, ht
      use navier
      use wavecontrol, only : at1
      implicit none

      character(len=*), intent(IN) :: fgrid, fparam, fdata, faux

      integer,  intent(IN) :: ioutd, ioutc
      real(wp), intent(IN) :: time
      type(HDFObj) :: writer

      call IO_Init

      if (myid.eq.0) then
         open(14000, file=fparam)
         write(14000,*) time, ioutd, ioutc
         close(14000)
      end if

      if (isbot .and. myid1==0 .and. iwcontrol==1) then
         open(14000, file="restart_wave.dat")
         write(14000,*) at1
         close(14000)
      end if

      call writerOpen(fgrid, MPI_COMM_2D_CART, writer)
      call write1D_z(writer, "zz", zz(1:xsz(3)), 0)
      call write1D_z(writer, "zw", zw(1:xsz(3)), 0)
      call write1D_z(writer, "dz", dz(1:xsz(3)), 0)
      call write1D_z(writer, "dzw", dzw(1:xsz(3)), 0)
      call writerClose(writer)

      call writerOpen(fdata, MPI_COMM_2D_CART, writer)
      !call write2D_xy(writer, "uzfs", uzfs, topid)
      !call write2D_xy(writer, "vzfs", vzfs, topid)
      !call write2D_xy(writer, "wzfs", wzfs, topid)
      call write2D_xy(writer, "uzfs", uzfs, 0)
      call write2D_xy(writer, "vzfs", vzfs, 0)
      call write2D_xy(writer, "wzfs", wzfs, 0)
      call write2D_xy(writer, "eta", eta, 0)
      call write2D_xy(writer, "eta0", eta0, 0)
      call write2D_xy(writer, "hh", hh, 0)
      call write2D_xy(writer, "ht", ht, 0)
      call write3D_xyz(writer, "u", u(:,:,1:xsz(3)))
      call write2D_xy(writer, "u1", u(:,:,xsz(3)+1), topid)
      !call write2D_xy(writer, "u1", u(:,:,xsz(3)+1), 0)
      call write3D_xyz(writer, "v", v(:,:,1:xsz(3)))
      call write2D_xy(writer, "v1", v(:,:,xsz(3)+1), topid)
      !call write2D_xy(writer, "v1", v(:,:,xsz(3)+1), 0)
      call write3D_xyz(writer, "w", w(:,:,1:xsz(3)))     
      call write3D_xyz(writer, "pp", pp(:,:,1:xsz(3)))
      call writerClose(writer)

!      print *,"myid=",myid

      call writerOpen(faux, MPI_COMM_2D_CART, writer)
      call write3D_xyz(writer, "hu", hu(:,:,1:xsz(3)))
      call write3D_xyz(writer, "hv", hv(:,:,1:xsz(3)))
      call write3D_xyz(writer, "hw", hw(:,:,1:xsz(3)))
      call writerClose(writer)

!      print *,"myid=",myid

      call IO_Finalize

    end subroutine save_file_name
    
!checked
   subroutine reread_file_name(time, ioutd, ioutc, fgrid, fparam, fdata, faux)
      use param, only: iwcontrol
      use grid, only : zz, zw, dz, dzw, eta, eta0, hh, ht
      use navier
      use wavecontrol, only : at1
      implicit none

      character(len=*), intent(IN) :: fgrid, fparam, fdata, faux

      integer,  intent(OUT) :: ioutd, ioutc
      real(wp), intent(OUT) :: time
      real(wp), dimension(xsz(1),xsz(2)) :: u1, v1
      type(HDFObj) :: writer

      call IO_Init

      open(14000, file=fparam, action="READ", status="OLD")
      read(14000,*) time, ioutd, ioutc
      close(14000)

      if (isbot .and. iwcontrol==1) then
         open(14000, file="restart_wave.dat", action="READ", status="OLD")
         read(14000,*) at1
         close(14000)
      end if
      
      call readerOpen(fgrid, MPI_COMM_2D_CART, writer)
      call read1D_z(writer, "zz", zz(1:xsz(3)), -1)
      call read1D_z(writer, "zw", zw(1:xsz(3)), -1)
      call read1D_z(writer, "dz", dz(1:xsz(3)), -1)
      call read1D_z(writer, "dzw", dzw(1:xsz(3)), -1)
      call readerClose(writer)

      call readerOpen(fdata, MPI_COMM_2D_CART, writer)
      call read2D_xy(writer, "uzfs", uzfs, -1)
      call read2D_xy(writer, "vzfs", vzfs, -1)
      call read2D_xy(writer, "wzfs", wzfs, -1)
      call read2D_xy(writer, "eta", eta, -1)
      call read2D_xy(writer, "eta0", eta0, -1)
      call read2D_xy(writer, "hh", hh, -1)
      call read2D_xy(writer, "ht", ht, -1)
      call read3D_xyz(writer, "u", u(:,:,1:xsz(3)))
      call read2D_xy(writer, "u1", u1(:,:), topid)
      !call read2D_xy(writer, "u1", u1(:,:), 0)
      call read3D_xyz(writer, "v", v(:,:,1:xsz(3)))
      call read2D_xy(writer, "v1", v1(:,:), topid)
      !call read2D_xy(writer, "v1", v1(:,:), 0)
      call read3D_xyz(writer, "w", w(:,:,1:xsz(3)))
      call read3D_xyz(writer, "pp", pp(:,:,1:xsz(3)))
      call readerClose(writer)

      call readerOpen(faux, MPI_COMM_2D_CART, writer)
      call read3D_xyz(writer, "hu", hu(:,:,1:xsz(3)))
      call read3D_xyz(writer, "hv", hv(:,:,1:xsz(3)))
      call read3D_xyz(writer, "hw", hw(:,:,1:xsz(3)))
      call readerClose(writer)

      call IO_Finalize

      call update_ghost(zz,1)
      call update_ghost(dz,1)
      call update_ghost(zw,1)
      call update_ghost(dzw,1)
      if (istop) zz(xsz(3)+1) = zz(xsz(3))+dz(xsz(3))
      call update_ghost(u,1)
      call update_ghost(v,1)
      call update_ghost(w,1)
      if (istop) then
         u(:,:,xsz(3)+1) = u1(:,:)
         v(:,:,xsz(3)+1) = v1(:,:)
      end if
      call update_ghost(pp,1)

    end subroutine reread_file_name
    
!checked
   subroutine input_iht
      use param
      use wavecontrol, only: ak, aa
      use constants, only : PI, TWOPI
      use hos_param

      IMPLICIT NONE            
      !-------------------------
      !     READ INPUT FILE
      !-------------------------

      READ(11,*) ISTART
      READ(11,*) ISCALAR
      READ(11,*) ITURBINE
      READ(11,*) IWCONTROL
      READ(11,*) NP1,NP2
      READ(11,*) PEX,PEY
      READ(11,*) ZL,HBAR
      READ(11,*) z0
      READ(11,*) NX,NY,NZ
      READ(11,*) NXS,NYS
      READ(11,*) DT,NTIME
      READ(11,*) ITMAX,ERLIM
      READ(11,*) RESBOT,RESTOP
      READ(11,*) FR2,RWE
      READ(11,*) NOUTD,NOUTC
      READ(11,*) CLBETA,CLGAMA
      READ(11,*) ARM
      !READ(11,*) AKA
      READ(11,*) HKA,NWAVE,CRAT
      READ(11,*) TIMEWAVY,TCOEF     
      !------------------------
      !     WAVY WALL TYPE
      !------------------------

      READ(11,*) IWAVY

      !     IWAVY=1: SOLID WAVY WALL
      !     IWAVY=2: VERTICALLY MOVING WAVY WALL
      !     IWAVY=3: WATER WAVE SURFACE

      !----------------------------------
      !     PARAMETERS FOR FILTERING
      !----------------------------------

      READ(11,*) ERVFILT,NFILT,IFILT

      !     ERVFILT: CRITERIO FOR VELOCITY FILTERING
      !     NFILT: 1/NFILT MODES WILL REMAIN AFTER FILTERING
      !     IFILT: IFILT=1--FILTERING; IFILT=0--NOT FILTERING

      READ(11,*) IPA,NTH
      READ(11,*) TIMEP,TCP
      READ(11,*) TIMEW
      READ(11,*) RDGL
      READ(11,*) TIMETURB

      !wave control parameter
      READ(11,*) aka, nswave

      !-----END HERE

      !HOS parameter
      READ(11,*) nswavex
      READ(11,*) ustar, U10, USS
      READ(11,*) gamma, fetch, g, phi
    

      !--------------------------
      !     SPATIAL PARAMETER
      !--------------------------

      XL=TWOPI/PEX
      YL=TWOPI/PEY
      !NXMOD=NX-1
      !NYMOD=NY-1
      DX=XL/NX
      DY=YL/NY

      !NXMODS=NXS-1
      !NYMODS=NYS-1
      
      !-----END HERE
      !-------------------------------
      !     COEF FOR SURFACE WAVE
      !-------------------------------
      !     AK -- WAVE NUMBER
      !     AA -- WAVE AMPLITUDE
      !     OMEG -- WAVE ANGLE FREQUENCY
      !     NSWAVE --NUMBER OF SURFACE WAVES IN X DIRECTION
      !     AKA -- WAVE SLOPE
      !     CSW -- PHASE SPEED OF SURFACE WAVE

      AK=PEX*NSWAVE
      AA=AKA/AK

      !C--TEST
      ! IF ( MYID == 0 ) THEN
      !    WRITE(*,*) 'AK=',AK
      !    write(*,*) 'AA=',AA
      !    write(*,*) 'OMEGA=',OMEG
      !    write(*,*) 'FR^2=',FR2, '1/WE=',RWE
      ! END IF
      !C--@@@@@@@@@@@@@@@@@@@@@@@@@

      !-----END HERE

      !-------------------------------
      !     TURBULENCE PARAMETERS
      !-------------------------------

      !      GOTO 10
      
      IF(ABS(RESTOP) .LE. 1.E-6)THEN
         ZLSBOT=1./RESBOT
         ZLSTOP=0.
         USBOT=1./(2.5*LOG(hbar/z0))
         USTOP=0.
!========================================
!WHEN TURN ON THE WIND TURBINE,UNCOMMENT
!THE FOLLOWING LINE                  
         ! USBOT = USBOT*2.675
!========================================         
         RE=RESBOT*(2.5*LOG(hbar/z0))
         BFORCE=USBOT**2/HBAR
         GOTO 10
      ENDIF

      IF(ABS(RESBOT) .LE. 1.E-6)THEN
         ZLSBOT=0.
         ZLSTOP=1./RESTOP
         USBOT=0.
         USTOP=1./(2.5*LOG(RESTOP)+5.0)
         RE=RESTOP*(2.5*LOG(RESTOP)+5.0)
         BFORCE=-USTOP**2/HBAR
         GOTO 10
      ENDIF

      ZLSBOT=1./RESBOT
      ZLSTOP=1./RESTOP
      USBOT=1./((2.5*LOG(RESBOT**2/(RESBOT+RESTOP))+5.0) &
      +(2.5*LOG(RESTOP**2/(RESBOT+RESTOP))+5.0)*RESTOP &
      /RESBOT)
      USTOP=1./((2.5*LOG(RESBOT**2/(RESBOT+RESTOP))+5.0)*RESBOT &
      /RESTOP+(2.5*LOG(RESTOP**2/(RESBOT+RESTOP))+5.0))
      RE=RESBOT*(2.5*LOG(RESBOT**2/(RESBOT+RESTOP))+5.0) &
      +RESTOP*(2.5*LOG(RESTOP**2/(RESBOT+RESTOP))+5.0)
      BFORCE=(USBOT**2-USTOP**2)/HBAR

10    CONTINUE

!     re=  98089451.137971550
!     bforce=2.2475641541888847E-004
!====COEF FOR BOTTOM WAVY WALL
!    HK: WAVE NUMBER      
!    HA: WAVE AMPLITUDE
!    HOMEG: WAVE ANGLE FREQUENCY      
!    VPHASE: WAVE PHASE SPEED
!    NWAVE: NUMNER OF WAVE IN X DIRECTION      

      HK = PEX*NWAVE
      HA = HKA/HK
      VPHASE = CRAT * USBOT
      HOMEG = HK*VPHASE
      
    end subroutine input_iht
    
!checked
   subroutine input_LES(infile)
     use param, only : mfilt, icsc, ILASD, TLASD, TLASD0, IVANDRIEST, cs0, aplus, &
          IWAVEBOT, ZCS0
      implicit none

      character(len=*), intent(IN) :: infile

      logical :: isexist

      inquire(FILE=infile, EXIST=isexist)
      if (isexist) then
         open(UNIT=2, FILE=infile, ACTION='READ')
         read(2,*) mfilt
         read(2,*) ICSC
         read(2,*) ILASD, TLASD, TLASD0
         read(2,*) IVANDRIEST, cs0, aplus
         read(2,*) IWAVEBOT, ZCS0
         close(2)
      endif
      
   end subroutine input_LES

!checked   
   subroutine output_all(time, id)
      use iso_fortran_env, only : INT64
      use param
      use grid, only : zz, zw, dz, dzw, eta, eta0, hh, ht
      use navier
      implicit none

      real(wp), intent(IN) :: time
      integer(INT64), intent(OUT) :: id

      type(HDFObj) :: writer
      character(len=128) :: filename

      call IO_Init
      write(filename,'(I14.14)') int(time*1D8,kind=INT64)
      filename = "DAT"//trim(adjustl(filename))

      if (myid == 0) then
         open(92, file=trim(filename)//".dat")
         write(92, *) time, nx_global, ny_global, nz_global, xl, yl, zl, hbar
         close(92)
      end if

      call writerOpen(trim(filename), MPI_COMM_2D_CART, writer)
      call write1D_z(writer, "zz", zz(1:xsz(3)), 0)
      call write1D_z(writer, "zw", zw(1:xsz(3)), 0)
      call write1D_z(writer, "dz", dz(1:xsz(3)), 0)
      call write1D_z(writer, "dzw", dzw(1:xsz(3)), 0)
      call write2D_xy(writer, "eta", eta, 0)
      call write2D_xy(writer, "eta0", eta0, 0)
      call write2D_xy(writer, "hh", hh, 0)
      call write3D_xyz(writer, "u", u(:,:,1:xsz(3)))
      call write3D_xyz(writer, "v", v(:,:,1:xsz(3)))
      call write3D_xyz(writer, "w", w(:,:,1:xsz(3)))
      call write3D_xyz(writer, "pp", pp(:,:,1:xsz(3)))
      call writerClose(writer)

      call IO_Finalize

      id = int(time*1D8)

   end subroutine output_all

   !> plyunote: modified filename extension, added grid data
   subroutine output_all_2(time, id)
      use iso_fortran_env, only : INT64
      use param
      use grid, only : zz, zw, dz, dzw, eta, eta0, hh, ht
      use navier
      use turbine_model, only: nxwt, nywt, ud, ud1, ud_ref, ud1_ref,&
        NumberOfTurbines, wtm, fsi_wt, NumberOfNacelle, fsi_nac, &
        fsitype, pp_wt, nacelle_model !> added by plyu
      implicit none

      real(wp), intent(IN) :: time
      integer(INT64), intent(OUT) :: id

      type(HDFObj) :: writer
      character(len=128) :: filename

      real(wp), allocatable, dimension(:,:,:) :: gridx, gridy, gridz
      integer :: i, j, k, ibi

      allocate(gridx(xsz(1),xsz(2),xsz(3)))
      allocate(gridy(xsz(1),xsz(2),xsz(3)))
      allocate(gridz(xsz(1),xsz(2),xsz(3)))
      do i = 1,xsz(1)
      do j = 1,xsz(2)
      do k = 1,xsz(3)
        gridx(i,j,k) = (xst(1)+i-1)*xl/nx_global
        gridy(i,j,k) = (xst(2)+j-1)*yl/ny_global
        !z = zzall(k)*(etamy(i)+hbar)-hbar
        gridz(i,j,k) = zz(k)*(hbar+hh(i,j)) - hh(i,j) 
      enddo
      enddo
      enddo      

      call IO_Init
      
      if (myid.eq.0) then
         open(14000, file="restart_param_temp.dat")
         write(14000,*) time, ' 0 0'
         close(14000)
      end if

      !> added by plyu: start
      if (myid .eq. 0 .and. iturbine .eq. 1) then
        open(15010, file='restart_ad_temp.dat')
        do i = 1, nxwt
          do j = 1, nywt
            write(15010,*) ud1((i-1)*nywt+j), ud((i-1)*nywt+j), &
              ud1_ref((i-1)*nywt+j), ud_ref((i-1)*nywt+j)
          enddo
        enddo
        close(15010)
      endif
      
      if (myid .eq. 0 .and. iturbine .eq. 3) then
        open(15010, file='restart_uref_temp.dat')
        do ibi = 1, NumberOfTurbines
          write(15010, *) wtm(ibi)%U1_ref, wtm(ibi)%U_ref
        enddo
        close(15010)
      endif

      if (myid .eq. 0 .and. iturbine .gt. 1) then
        open(15012, file='Turbine2_temp.inp', action='write', &
          status = 'replace')
        do ibi = 1, NumberOfTurbines
          !> yaw angle to be modified
          write(15012, *) fsi_wt(ibi)%ang_axis, fsi_wt(ibi)%ang0_yaw, &
            fsi_wt(ibi)%angvel_axis, " $ ang_axis, ang_yaw, angvel_axis" 
        enddo

        close(15012)

        if (fsitype.eq.1) then
          open(15013, file='FSI2_rotor_temp.inp', action='write', status='replace')
          do ibi = 1, NumberOfTurbines
            write(15013,*) fsi_wt(ibi)%x_c, fsi_wt(ibi)%y_c, fsi_wt(ibi)%z_c, &
              fsi_wt(ibi)%nx_tb, fsi_wt(ibi)%ny_tb, fsi_wt(ibi)%nz_tb
          enddo
          close(15013)
          if (nacelle_model .ne. 0) then
            open(15014, file='FSI2_nacelle_temp.inp', action='write', status='replace')
            do ibi = 1, NumberOfNacelle
              write(15014,*) fsi_nac(ibi)%x_c, fsi_nac(ibi)%y_c, &
                fsi_nac(ibi)%z_c, fsi_nac(ibi)%nx_tb, fsi_nac(ibi)%ny_tb, &
                fsi_nac(ibi)%nz_tb
            enddo
            close(15014)
          endif  
        endif
      endif
      !> added by plyu: end
      
      
      write(filename,'(a4,I0.10)') 'DAT_',nint(time/dt)
      !filename = "DAT"//trim(adjustl(filename))

      if (myid == 0) then
         open(92, file=trim(filename)//".dat")
         write(92, *) time, nx_global, ny_global, nz_global, xl, yl, zl, hbar
         close(92)
      end if

      call writerOpen(trim(filename)//".h5", MPI_COMM_2D_CART, writer)
     ! call write3D_xyz(writer, "x", gridx(:,:,1:xsz(3)))
      !call write3D_xyz(writer, "y", gridy(:,:,1:xsz(3)))
      call write3D_xyz(writer, "z", gridz(:,:,1:xsz(3)))
      !call write1D_z(writer, "zz", zz(1:xsz(3)), 0)
      !call write1D_z(writer, "zw", zw(1:xsz(3)), 0)
      !call write1D_z(writer, "dz", dz(1:xsz(3)), 0)
      !call write1D_z(writer, "dzw", dzw(1:xsz(3)), 0)
      call write2D_xy(writer, "eta", eta, 0)
      !call write2D_xy(writer, "eta0", eta0, 0)
      call write2D_xy(writer, "hh", hh, 0)
      call write3D_xyz(writer, "u", u(:,:,1:xsz(3)))
      call write3D_xyz(writer, "v", v(:,:,1:xsz(3)))
      call write3D_xyz(writer, "w", w(:,:,1:xsz(3)))
      call write3D_xyz(writer, "pp", pp(:,:,1:xsz(3)))
      call writerClose(writer)

      call writerOpen("restart_aux_temp.h5", MPI_COMM_2D_CART, writer)
      call write3D_xyz(writer, "hu", hu(:,:,1:xsz(3)))
      call write3D_xyz(writer, "hv", hv(:,:,1:xsz(3)))
      call write3D_xyz(writer, "hw", hw(:,:,1:xsz(3)))
      call write2D_xy(writer, "ht", ht, 0)
      call write2D_xy(writer, "u1", u(:,:,xsz(3)+1), topid)
      call write2D_xy(writer, "v1", v(:,:,xsz(3)+1), topid)
      call writerClose(writer)

      call IO_Finalize

      id = int(time*1D8)

      deallocate(gridx, gridy, gridz)

   end subroutine output_all_2
   
   !> added by plyu
   subroutine reread_output(ti, filename, time)
      use iso_fortran_env, only : INT64
      use param
      use grid, only : zz, zw, dz, dzw, eta, eta0, hh
      use navier
      use decomp, only : myid
      use mpi
      implicit none

      integer, intent(in) :: ti
      character (len=64), intent(in) :: filename      
      real(wp), intent(out) :: time

      integer(INT64) :: id

      type(HDFObj) :: reader

      !real(wp), allocatable, dimension(:,:,:) :: gridx, gridy, gridz
      integer :: i, j, k, ierr
      !print *, "reread, 1"
      call IO_Init
      !print *, "reread, 2"
      if (myid .eq. 0) then 
        open(92, file=trim(filename)//".dat")
        read(92, *) time, nx_global, ny_global, nz_global, xl, yl, zl, hbar
        close(92)
      endif
      !print *, "reread, 3"
      call MPI_BCAST(time,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
      call MPI_BCAST(nx_global,1,mpi_integer,0,mpi_comm_2d_cart,ierr)
      call MPI_BCAST(ny_global,1,mpi_integer,0,mpi_comm_2d_cart,ierr)
      call MPI_BCAST(nz_global,1,mpi_integer,0,mpi_comm_2d_cart,ierr)
      call MPI_BCAST(xl,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
      call MPI_BCAST(yl,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
      call MPI_BCAST(zl,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
      call MPI_BCAST(hbar,1,mpi_double_precision,0,mpi_comm_2d_cart,ierr)
      !print *, "reread, 4"

      if (nint(time/dt).ne.ti) then
        print *, "Timestep doesn't match: ti, time, dt=",ti, time, dt
        stop 101 
      endif
      
      !print *, xsz(1:3)
      !allocate(gridx(xsz(1),xsz(2),xsz(3)))
      !allocate(gridy(xsz(1),xsz(2),xsz(3)))
      !allocate(gridz(xsz(1),xsz(2),xsz(3)))
      !do i = 1,xsz(1)
      !do j = 1,xsz(2)
      !do k = 1,xsz(3)
      !  gridx(i,j,k) = (xst(1)+i-1)*xl/nx_global
      !  gridy(i,j,k) = (xst(2)+j-1)*yl/ny_global
      !  !z = zzall(k)*(etamy(i)+hbar)-hbar
      !  !gridz(i,j,k) = zz(k)*(hbar+hh(i,j)) - hh(i,j) 
      !enddo
      !enddo
      !enddo    
      !print *, "reread, 5"

      call readerOpen(trim(filename)//".h5", MPI_COMM_2D_CART, reader)
      !print *, "readh5, 1"
     ! call write3D_xyz(writer, "x", gridx(:,:,1:xsz(3)))
      !call write3D_xyz(writer, "y", gridy(:,:,1:xsz(3)))
      
      !call read1D_z(reader, "zz", zz(1:xsz(3)), -1)
      !call read1D_z(reader, "zw", zw(1:xsz(3)), -1)
      !call read1D_z(reader, "dz", dz(1:xsz(3)), -1)
      !call read1D_z(reader, "dzw", dzw(1:xsz(3)), -1)
      call read2D_xy(reader, "eta", eta, -1)
      !print *, "readh5, 2"
      !call read2D_xy(reader, "eta0", eta0, -1)
      call read2D_xy(reader, "hh", hh, -1)
      !print *, "readh5, 3"
      call read3D_xyz(reader, "u", u(:,:,1:xsz(3)))
      !print *, "readh5, 4"
      call read3D_xyz(reader, "v", v(:,:,1:xsz(3)))
      !print *, "readh5, 5"
      !print *, size(w(:,:,1:xsz(3)),1),size(w(:,:,1:xsz(3)),2), xsz(3)
      call read3D_xyz(reader, "w", w(:,:,1:xsz(3)))
      !print *, "readh5, 6"
      call read3D_xyz(reader, "pp", pp(:,:,1:xsz(3)))
      !print *, "readh5, 7"
      !call read3D_xyz(reader, "z", gridz(:,:,1:xsz(3)))
      call readerClose(reader)
      !print *, "readh5, 8"

      call IO_Finalize
      !print *, "readh5, 9"

      !id = int(time*1D8)
   end subroutine reread_output

   subroutine output_inlet(time)
      use decomp
      use param
      use grid, only : zz, hh
      use navier
      use utils
      use MPI
      implicit none

      real(wp), intent(in) :: time

      real(wp) :: x,y,z
      real(wp), dimension(nx_global,ny_global,nz_global) :: uall, vall, wall
      real(wp), dimension(nx_global,ny_global) :: hhall
      real(wp), dimension(nz_global) :: tmpz,zzall

      integer :: i,j,k,root0,ierror, ti_
      character(len=128) :: filen


      ! Gather zz coordinates
      tmpz = 0
      tmpz(xst(3):xend(3)) = zz(1:xsz(3))
      call MPI_Allreduce(tmpz, zzall, nz_global, real_type, &
                         MPI_SUM, MPI_COMM_2D_ROW, ierror)

      root0 = 0
      call gather_2d_xy(hh(1:xsz(1),1:xsz(2)),hhall,root0)
      call gather_3d_xyz(u(1:xsz(1),1:xsz(2),1:xsz(3)),uall,root0)
      call gather_3d_xyz(v(1:xsz(1),1:xsz(2),1:xsz(3)),vall,root0)
      call gather_3d_xyz(w(1:xsz(1),1:xsz(2),1:xsz(3)),wall,root0)

      if (myid == root0) then

         ! Export x cut   
         filen = ''
         ti_ = nint(time/dt)
         call system('mkdir -p output_inlet')
         write(filen, '(a21,i0.10,a4)') './output_inlet/inlet_',ti_,'.dat'
         open(18007,FILE=filen,STATUS='REPLACE') ! add by plyu
         write(18007,*) "VARIABLES = X, Y, Z, U, V, W"
         write(18007,"(A,F20.8,A,A,I6,A,I6,A,F20.8)") 'ZONE T="',time,'" I=9',&
              ' J=', ny_global, ' K=',nz_global,' SOLUTIONTIME=', time
         !i = 1
         !x = 0.0_wp
         do k=1,nz_global
           do j=1, ny_global
             do i = 1, 9
               x = (i-1)*xl/nx_global
               y = (j-1)*yl/ny_global
               !z = zzall(k)*hbar-hbar
               z = zzall(k)*(hbar+hhall(i,j)) - hhall(i,j)
               write(18007,*) x, y, z, uall(i,j,k), vall(i,j,k), wall(i,j,k)
             enddo
           enddo
         enddo
         close(18007) ! add by plyu
         print *, 'Inlet information output: ', filen
         
      end if

    end subroutine output_inlet

   subroutine output_cut(time)
      use decomp
      use param
      use grid, only : zz, hh
      use navier
      use utils
      use MPI
      implicit none

      real(wp), intent(in) :: time

      real(wp) :: x,y,z
      real(wp), dimension(nx_global,ny_global,nz_global) :: uall, ppall
      real(wp), dimension(nx_global,ny_global) :: hhall
      real(wp), dimension(nz_global) :: tmpz,zzall

      integer :: i,j,k,root0,ierror


      ! Gather zz coordinates
      tmpz = 0
      tmpz(xst(3):xend(3)) = zz(1:xsz(3))
      call MPI_Allreduce(tmpz, zzall, nz_global, real_type, &
                         MPI_SUM, MPI_COMM_2D_ROW, ierror)

      root0 = 0
      call gather_2d_xy(hh(1:xsz(1),1:xsz(2)),hhall,root0)
      call gather_3d_xyz(u(1:xsz(1),1:xsz(2),1:xsz(3)),uall,root0)
      call gather_3d_xyz(pp(1:xsz(1),1:xsz(2),1:xsz(3)),ppall,root0)

      if (myid == root0) then

         ! Export x cut     
         open(920,FILE='fort.920.dat',STATUS='REPLACE') ! add by plyu
         write(920,*) "VARIABLES = X, Y, Z, U, pp"
         write(920,"(A,F20.8,A,I5,A,I5,A,F20.8)") 'ZONE T="',time,'" I=',ny_global, &
              ' J=',nz_global,' SOLUTIONTIME=', time
         i = 1
         x = 0.0_wp
         do k=1,nz_global
            do j=1, ny_global
               y = (j-1)*yl/ny_global
               !z = zzall(k)*hbar-hbar
               z = zzall(k)*(hbar+hhall(i,j)) - hhall(i,j)
               write(920,*) x, y, z, uall(i,j,k), ppall(i,j,k)
            end do
         end do
         close(920) ! add by plyu
         
         ! Export y cut ()
         open(910,FILE='fort.910.dat',STATUS='REPLACE') ! add by plyu
         write(910,*) "VARIABLES = X, Y, Z, U, pp"
         write(910,"(A,F20.8,A,I5,A,I5,A,F20.8)") 'ZONE T="',time,'" I=',nx_global, &
              ' J=',nz_global,' SOLUTIONTIME=', time
         j = 1
         y = 0.0_wp
         do k=1,nz_global
            do i=1, nx_global
               x = (i-1)*xl/nx_global
               !z = zzall(k)*(etamy(i)+hbar)-hbar
               z = zzall(k)*(hbar+hhall(i,j)) - hhall(i,j)
               write(910,*) x, y, z, uall(i,j,k), ppall(i,j,k)
               !print *,'plyudebug: i,k,x,z,u,p',&
               !  i,k,x,z,uall(i,j,k),ppall(i,j,k)
            end do
         end do         
         close(910) ! add by plyu
      end if

    end subroutine output_cut

   subroutine output_xy_average(time)
      use decomp
      use param
      use grid, only : zz
      use navier
      use utils
      use MPI
      implicit none

      real(wp) :: z, time
      real(wp), dimension(nx_global, nz_global) :: umy, vmy, wmy, wmy1, pmy
      real(wp), dimension(nx_global, nz_global) :: uf2my, vf2my, wf2my
!      real(wp), dimension(ny_global, nz_global) :: umx, vmx, wmx, wmx1, pmx
      real(wp), dimension(1:xsz(1),1:xsz(2),1:xsz(3)) :: tmpf
      real(wp), dimension(nz_global) :: tmpz, zzall
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
               tmpf(i,j,k) = (u(i,j,k)-umy(1,k+xst(3)-1))**2
            end do
         end do
      end do
      call get_average_y(tmpf, uf2my)

      do k=1, xsz(3)
         do j=1, xsz(2)
            do i=1, xsz(1)
               tmpf(i,j,k) = (v(i,j,k)-vmy(1,k+xst(3)-1))**2
            end do
         end do
      end do
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
      tmpf = tmpf**2
      call get_average_y(tmpf, wf2my)

      !---------------------
      ! Write plane average
      !---------------------
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
         open(930,FILE='fort.930.dat',STATUS='REPLACE') ! add by plyu
         write(930,*) "VARIABLES = Z, UM, VM, WM, PM, UF, VF, WF"
         write(930,"(A,F20.8,A,I5,A,F20.8)") 'ZONE T="',time,'" I=',nz_global, &
                                           ' SOLUTIONTIME=', time
         do k=1,nz_global
            z = zzall(k)*hbar 
            write(930,*) z, umy(1,k), vmy(1,k), wmy1(1,k), pmy(1,k), &
                            sqrt(uf2my(1,k)), sqrt(vf2my(1,k)), sqrt(wf2my(1,k))
         end do
         close(930)
      end if

    end subroutine output_xy_average

   subroutine output_surface(time)
      use param, only:  xl, yl
      use grid, only : eta 
      use navier, only : u, v, w, pp
      use utils
      use MPI
      implicit none

      real(wp), intent(IN) :: time
      integer :: ierror, i, j
      real(wp) :: x, y
      real(wp), dimension(nx_global, ny_global) :: etaall, uall, vall, wall, ppall
      real(wp), dimension(nx_global, ny_global) :: tmp

      !if (istop) then
      if (isbot) then
         tmp = 0
         tmp(xst(1):xend(1), xst(2):xend(2)) = eta(:,:)
         call MPI_Reduce(tmp, etaall, nx_global*ny_global, real_type, &
                            MPI_SUM, 0, MPI_COMM_2D_COL, ierror)

         tmp = 0
         !tmp(xst(1):xend(1), xst(2):xend(2)) = u(:,:,xsz(3))
         tmp(xst(1):xend(1), xst(2):xend(2)) = u(:,:,1)
         call MPI_Reduce(tmp, uall, nx_global*ny_global, real_type, &
                            MPI_SUM, 0, MPI_COMM_2D_COL, ierror)

         tmp = 0
         !tmp(xst(1):xend(1), xst(2):xend(2)) = v(:,:,xsz(3))
         tmp(xst(1):xend(1), xst(2):xend(2)) = v(:,:,1)
         call MPI_Reduce(tmp, vall, nx_global*ny_global, real_type, &
                            MPI_SUM, 0, MPI_COMM_2D_COL, ierror)

         tmp = 0
         !tmp(xst(1):xend(1), xst(2):xend(2)) = w(:,:,xsz(3)-1)
         tmp(xst(1):xend(1), xst(2):xend(2)) = w(:,:,1)
         call MPI_Reduce(tmp, wall, nx_global*ny_global, real_type, &
                            MPI_SUM, 0, MPI_COMM_2D_COL, ierror)

         tmp = 0
         tmp(xst(1):xend(1), xst(2):xend(2)) = pp(:,:,1)
         call MPI_Reduce(tmp, ppall, nx_global*ny_global, real_type, &
                            MPI_SUM, 0, MPI_COMM_2D_COL, ierror)

         if (myid1 == 0) then
            write(93,*) "VARIABLES = x, y, eta, u, v, w, pp"
            write(93,"(A,F20.8,A,I5,A,I5,A,F20.8)") 'ZONE T="',time,'" I=',nx_global, &
                                           ' J=',ny_global,' SOLUTIONTIME=', time
            do j=1, ny_global
               do i=1, nx_global
                  x = (i-1)*xl/nx_global
                  y = (j-1)*yl/ny_global
                  write(93,*) x, y, etaall(i,j), uall(i,j), vall(i,j), wall(i,j), ppall(i,j)
               end do
            end do
         end if

      end if

   end subroutine output_surface
   
   subroutine output_plane
      use grid, only : zz
      use decomp, only : nx_global, ny_global, wp
      use navier
      use utils

      implicit none
      
      integer k
      real(wp), dimension(xsz(3)) :: um, vm, wm  
      real(wp), dimension(xsz(3)) :: uq, vq, wq  
      real(wp), dimension(xsz(1),xsz(2),xsz(3)) :: ufl, vfl, wfl  


      type(HDFObj) :: writer

      call IO_Init

      do k=1, xsz(3)
         um(k)=sum_plane_xy(u(:,:,k), mpi_comm_2d_col)/nx_global/ny_global
         vm(k)=sum_plane_xy(v(:,:,k), mpi_comm_2d_col)/nx_global/ny_global
         wm(k)=sum_plane_xy(w(:,:,k), mpi_comm_2d_col)/nx_global/ny_global
      end do 

      do k=1, xsz(3)
         ufl(:,:,k)=(u(:,:,k)-um(k))**2
         vfl(:,:,k)=(v(:,:,k)-vm(k))**2
         wfl(:,:,k)=(w(:,:,k)-wm(k))**2
      end do 

      do k=1, xsz(3)
         uq(k)=sum_plane_xy(ufl(:,:,k), mpi_comm_2d_col)/nx_global/ny_global
         vq(k)=sum_plane_xy(vfl(:,:,k), mpi_comm_2d_col)/nx_global/ny_global
         wq(k)=sum_plane_xy(wfl(:,:,k), mpi_comm_2d_col)/nx_global/ny_global
      end do 

      call writerOpen("plane.h5", MPI_COMM_2D_CART, writer)
      call write1D_z(writer, "zz", zz(1:xsz(3)), 0)
      call write1D_z(writer, "um", um(1:xsz(3)), 0)
      call write1D_z(writer, "vm", vm(1:xsz(3)), 0)
      call write1D_z(writer, "wm", wm(1:xsz(3)), 0)
      call write1D_z(writer, "uq", uq(1:xsz(3))**0.5_wp, 0)
      call write1D_z(writer, "vq", vq(1:xsz(3))**0.5_wp, 0)
      call write1D_z(writer, "wq", wq(1:xsz(3))**0.5_wp, 0)
      call writerClose(writer)

      call IO_Finalize

   end subroutine output_plane

   subroutine output_eta
      use grid, only : eta
      use decomp, only : nx_global, ny_global, wp
      use param, only : xl, yl
      use navier
      use utils
      use spectral

      implicit none

      real(wp), dimension(nx_global, ny_global) :: eta_all
      integer :: i, j

      call gather_2d_xy(eta, eta_all)

      if(myid==0)then 
         open(10, file="eta.dat")
         write(10,*)'variables=x , y, z ' 
         write(10,*) 'ZONE I=',nx_global,' J=',ny_global,' K=',1,' DATAPACKING=POINT'
         do j=1, ny_global 
            do i=1, nx_global 
               write(10,*) (i-nx_global/2-1)*xl/nx_global, (j-1)*yl/ny_global, eta_all(i,j) 
            end do
         end do
         close(10)
      end if 

   end subroutine 
end module

