module io_hos

  use hos_param
  use spectral_hos
  use utils_hos
  use hdf_io_hos

  implicit none
  
  integer,public :: nxmax,nymax

  public :: save_hos, read_hos
  public :: outsurf

  interface save_hos
     module procedure save_hos_new, save_hos_hdf5_1, save_hos_hdf5_2
  end interface save_hos

  interface read_hos
     module procedure read_hos_new, read_hos_hdf5_1, read_hos_hdf5_2
  end interface read_hos

  interface outsurf
     module procedure outsurf_gravity,outsurf_cap
  end interface outsurf

contains

  subroutine save_hos_new(eta,vps,pa)

    implicit none
    
    real(wp), intent(in), dimension(nxhos,nyhos/ncpu_hos) :: eta,vps,pa
    integer fileid
    
    fileid = 170000 + myid
    open (fileid)
    write (fileid,*) nxhos,nyhos,ncpu_hos,time_hos,fr2
    write (fileid,*) eta, vps, pa
    close(fileid)
  end subroutine save_hos_new

  subroutine read_hos_new(eta,vps,pa)

    implicit none

    real(wp), intent(inout), dimension(nxhos,nyhos/ncpu_hos) :: eta,vps,pa
    integer idum
    real(wp) dum
    integer fileid

    fileid = 160000 + myid
    open (fileid)
    read (fileid,*) idum,idum,idum,time_hos,dum
    read (fileid,*) eta, vps, pa
    close(fileid)
  end subroutine read_hos_new

  !--------------------------------------------------------------

  subroutine save_hos_hdf5_1(time, ioutd, ioutc, eta, vps, pa)
    
    implicit none
    
    integer,  intent(IN) :: ioutd, ioutc
    real(wp), intent(IN) :: time
    real(wp), intent(IN), dimension(:,:) :: eta,vps,pa

    type(HDFObj_hos) :: writer

    call IO_Init_hos
    
    if (myid == 0) then
       open(16000, file="restart_param_hos.dat")
       write(16000,*) time, ioutd, ioutc
       close(16000)
    end if
    
    call writerOpen_hos("restart_hos.h5", mpi_comm_world, writer)
    call write2D_xy_hos(writer, "eta_hos", eta)
    call write2D_xy_hos(writer, "vps_hos", vps)
    call write2D_xy_hos(writer, "pa_hos", pa)
    call writerClose_hos(writer)
    
    call IO_Finalize_hos
    
  end subroutine save_hos_hdf5_1
    
  !--------------------------------------------------------------

  subroutine save_hos_hdf5_2(time, ioutd, ioutc, eta, vps, pa, fparam, fdata)

    implicit none

    character(len=*), intent(IN) :: fparam, fdata

    integer,  intent(IN) :: ioutd, ioutc
    real(wp), intent(IN) :: time
    real(wp), intent(IN), dimension(:,:) :: eta,vps,pa

    type(HDFObj_hos) :: writer

    call IO_Init_hos

    if (myid == 0) then
       open(16000, file=fparam)
       write(16000,*) time, ioutd, ioutc
       close(16000)
    end if

    call writerOpen_hos(fdata, mpi_comm_world, writer)
    call write2D_xy_hos(writer, "eta_hos", eta)
    call write2D_xy_hos(writer, "vps_hos", vps)
    call write2D_xy_hos(writer, "pa_hos", pa)
    call writerClose_hos(writer)

    call IO_Finalize_hos

  end subroutine save_hos_hdf5_2

  !--------------------------------------------------------------

  subroutine read_hos_hdf5_1(time, ioutd, ioutc, eta, vps, pa)

    implicit none

    integer,  intent(out) :: ioutd, ioutc
    real(wp), intent(out) :: time
    real(wp), intent(out), dimension(:,:) :: eta,vps,pa

    type(HDFObj_hos) :: reader

    call IO_Init_hos

    open(16000, file="restart_param_hos.dat")
    read(16000,*) time, ioutd, ioutc
    close(16000)
    
    call readerOpen_hos("restart_hos.h5", mpi_comm_world, reader)
    call read2D_xy_hos(reader, "eta_hos", eta)
    call read2D_xy_hos(reader, "vps_hos", vps)
    call read2D_xy_hos(reader, "pa_hos", pa)
    call readerClose_hos(reader)

    call IO_Finalize_hos

  end subroutine read_hos_hdf5_1

  !--------------------------------------------------------------

  subroutine read_hos_hdf5_2(time, ioutd, ioutc, eta, vps, pa, fparam, fdata)

    implicit none

    character(len=*), intent(IN) :: fparam, fdata

    integer,  intent(out) :: ioutd, ioutc
    real(wp), intent(out) :: time
    real(wp), intent(out), dimension(:,:) :: eta,vps,pa

    type(HDFObj_hos) :: reader

    call IO_Init_hos

    open(16000, file=fparam)
    read(16000,*) time, ioutd, ioutc
    close(16000)

    call readerOpen_hos(fdata, mpi_comm_world, reader)
    call read2D_xy_hos(reader, "eta_hos", eta)
    call read2D_xy_hos(reader, "vps_hos", vps)
    call read2D_xy_hos(reader, "pa_hos", pa)
    call readerClose_hos(reader)

    call IO_Finalize_hos

  end subroutine read_hos_hdf5_2

  !--------------------------------------------------------------

  subroutine outsurf_cap(eta,pa,ioutd,ioutc,time)

    implicit none
    real(wp), intent(in) :: time
    real(wp), intent(in), dimension(:,:) :: eta, pa
    integer, intent(inout) :: ioutd,ioutc

    real(wp) x,y
    integer i,j
    real(wp), allocatable, dimension(:,:) :: eo, po
    real(wp), dimension(nxhos / 2) :: skx, spkx
    real(wp), dimension(nyhos / 2) :: sky

    ioutc = ioutc + 1
    if (ioutc == noutc_hos) then
       ioutc = 0

       !export surface elevation
       allocate(eo(nxhos,nyhos))
       allocate(po(nxhos,nyhos))

       call alltoone(eta,eo)
       call alltoone(pa,po)

       call spec_x(eo,skx)
       call spec_x(po,spkx)


       if (myid == 0) then
          open(93, file='fort.93.dat', action='write', status='replace')
          write(93,*) " VARIABLES = x,y,eta,pa_cap"
          write(93,*) ' ZONE T="',time,'" I=', nxhos, ' J=',nyhos,' F=POINT'

          

          do j = 1, nyhos
             y = (j - 1) * dy_hos
             do i = 1, nxhos
                x = (i - 1) * dx_hos
                write(93,'(25e12.4)') x,y,eo(i,j),po(i,j)
             end do
          end do
          close(93)

          !export energy spectrum

          open(32, file='fort.32.dat', action='write', status='replace')
          write(32,*) " VARIABLES = kx, Skx, Spkx"
          write(32,*) ' ZONE T="',time,'" I=', nxhos/2
          do i = 1, nxhos / 2
             write(32,'(25e12.4)') i * pex_hos, skx(i), spkx(i)
          end do
          close(32)
       end if
       deallocate(eo,po)
    end if

    
  end subroutine outsurf_cap

  !--------------------------------------------------------------

  subroutine outsurf_gravity(ioutd,ioutc,time)

    implicit none
    real(wp), intent(in) :: time
    integer, intent(inout) :: ioutd,ioutc

    real(wp) x,y
    integer i,j
    real(wp), allocatable, dimension(:,:) :: eo, vo
    real(wp), dimension(nxhos / 2) :: skx
    real(wp), dimension(nyhos / 2) :: sky
    real(wp), dimension(nxhos,nyhos) :: skxy
    real(wp), dimension(nxhos/2,-nyhos/2:nyhos/2) :: fskxy
    real(wp) :: vol, flux, pe, ke, ener, rms, hs, skew, kurt

    character(len=128) :: filen
    real(wp) :: sigma, period, time2

    if (ist_hos == 1) then
      sigma = (nswavex * pex_hos * (1/fr2 + (nswavex*pex_hos)**2/bond/fr2))**0.5
    else
      sigma = (nswavex * pex_hos / fr2)**0.5
    end if
    period = twopi / sigma
    time2 = time * period / dt_hos

    ioutc = ioutc + 1
    if (ioutc == noutc_hos) then
       ioutc = 0

       !export surface elevation
       allocate(eo(nxhos,nyhos))
       call alltoone(eta_hos,eo)
       call spec_x(eo,skx)
       call spec_xy(eo,skxy)
       call fullspec_xy(eo,fskxy)
       if (myid == 0) then
          write(filen, '(a,i0.6,a)') 'fort_93_', nint(time2),'.dat'
          open(93, file=filen, action='write', status='replace')
          write(93,*) " VARIABLES = x,y,eta"
          write(93,*) ' ZONE T="',time,'" I=', nxhos, ' J=',nyhos,' F=POINT'
          do j = 1, nyhos
             y = (j - 1) * dy_hos
             do i = 1, nxhos
                x = (i - 1) * dx_hos
                write(93,'(25e12.4)') x,y,eo(i,j)
             end do
          end do
          close(93)

          !export kx spectrum
          write(filen, '(a,i0.6,a)') 'fort_32_', nint(time2),'.dat'
          open(32, file=filen, action='write', status='replace')
          write(32,*) " VARIABLES = kx, Skx"
          write(32,"(A,F16.8,A,I5)") ' ZONE T="',time,'" I=', nxhos/2
          do i = 1, nxhos / 2
             write(32,'(25e12.4)') i * pex_hos, skx(i)
          end do
          close(32)

          !export 2d spectrum
          write(filen, '(a,i0.6,a)') 'fort_33_', nint(time2),'.dat'
          open(33, file=filen, action='write', status='replace')
          write(33,*) " VARIABLES = kx,ky,skxy"
          write(33,"(A,F16.8,A,I5,A,I5,A)") ' ZONE T="',time,'" I=', nxhos/2, ' J=',nyhos/2,' F=POINT'
          do j = 1, nyhos / 2
             y = (j - 1) * pey_hos
             do i = 1, nxhos / 2
                x = (i - 1) * pex_hos
                write(33,'(25e12.4)') x,y,skxy(i,j)
             end do
          end do
          close(33)

          !export 2d spectrum
          write(filen, '(a,i0.6,a)') 'fort_34_', nint(time2),'.dat'
          open(34, file=filen, action='write', status='replace')
          write(34,*) " VARIABLES = kx,ky,skxy"
          write(34,"(A,F16.8,A,I5,A,I5,A)") ' ZONE T="',time,'" I=', nxhos / 2, ' J=',nyhos + 1,' F=POINT'
          do j = -nyhos/2, nyhos / 2
             y = j * pey_hos
             do i = 1, nxhos / 2
                x = i * pex_hos
                write(34,'(25e12.4)') x,y,fskxy(i,j)
             end do
          end do
          close(34)
       end if
       deallocate(eo)

       call get_stats(vol,flux,pe,ke,ener,eta_hos,vps_hos,feta)
       call skwave(eta_hos,rms,hs,skew,kurt)

       if (myid == 1) then
          write(55,*) time, vol, flux, pe, ke, ener, rms,hs,skew, kurt
       end if
       
    end if

    ioutd = ioutd + 1
    if (ioutd == noutd_hos) then
       ioutd = 0

       !export surface elevation
       allocate(eo(nxhos,nyhos))
       allocate(vo(nxhos,nyhos))
       call alltoone(eta_hos,eo)
       call alltoone(vps_hos,vo)
       if (myid == 0) then
          open(340, file='fort.340.dat', action='write', status='replace')
          write(340,'(25e12.4)') eo, vo
          close(340)
       end if
       deallocate(eo,vo)
    end if
  end subroutine outsurf_gravity

end module io_hos
