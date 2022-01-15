module io_hos

  use hos_param
  use decomp, only : wp, mpi_comm_2d_col
  use spectral_hos
  use hdf_io_hos

  implicit none
  
  integer,public :: nxmax,nymax
  public :: save_hos, read_hos
  public :: outsurf

  interface save_hos
     module procedure save_hos_hdf5_1, save_hos_hdf5_2
  end interface save_hos

  interface read_hos
     module procedure read_hos_hdf5_1, read_hos_hdf5_2
  end interface read_hos

contains

  !--------------------------------------------------------------

  subroutine save_hos_hdf5_1(time, ioutd, ioutc, eta, vps, pa)

    implicit none

    integer,  intent(IN) :: ioutd, ioutc
    real(wp), intent(IN) :: time
    real(wp), intent(IN), dimension(:,:) :: eta,vps,pa

    type(HDFObj_hos) :: writer

    call IO_Init_hos

    if (myid_hos == 0) then
       open(16000, file="restart_param_hos.dat")
       write(16000,*) time, ioutd, ioutc
       close(16000)
    end if

    call writerOpen_hos("restart_hos.h5", mpi_comm_2d_col, writer)
    call write2D_xy_hos(writer, "eta_hos", eta)
    call write2D_xy_hos(writer, "vps_hos", vps)
    call write2D_xy_hos(writer, "pa_hos", pa)
    call writerClose_hos(writer)

    call IO_Finalize_hos

  end subroutine save_hos_hdf5_1

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

    call readerOpen_hos("restart_hos.h5", mpi_comm_2d_col, reader)
    call read2D_xy_hos(reader, "eta_hos", eta)
    call read2D_xy_hos(reader, "vps_hos", vps)
    call read2D_xy_hos(reader, "pa_hos", pa)
    call readerClose_hos(reader)

    call IO_Finalize_hos

  end subroutine read_hos_hdf5_1

  !--------------------------------------------------------------

  subroutine save_hos_hdf5_2(time, ioutd, ioutc, eta, vps, pa, fparam, fdata)

    implicit none

    character(len=*), intent(IN) :: fparam, fdata

    integer,  intent(IN) :: ioutd, ioutc
    real(wp), intent(IN) :: time
    real(wp), intent(IN), dimension(:,:) :: eta,vps,pa

    type(HDFObj_hos) :: writer

    call IO_Init_hos

    if (myid_hos == 0) then
       open(16000, file=fparam)
       write(16000,*) time, ioutd, ioutc
       close(16000)
    end if

    call writerOpen_hos(fdata, mpi_comm_2d_col, writer)
    call write2D_xy_hos(writer, "eta_hos", eta)
    call write2D_xy_hos(writer, "vps_hos", vps)
    call write2D_xy_hos(writer, "pa_hos", pa)
    call writerClose_hos(writer)

    call IO_Finalize_hos

  end subroutine save_hos_hdf5_2

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

    call readerOpen_hos(fdata, mpi_comm_2d_col, reader)
    call read2D_xy_hos(reader, "eta_hos", eta)
    call read2D_xy_hos(reader, "vps_hos", vps)
    call read2D_xy_hos(reader, "pa_hos", pa)
    call readerClose_hos(reader)

    call IO_Finalize_hos

  end subroutine read_hos_hdf5_2

  !--------------------------------------------------------------

  subroutine outsurf(eta,time)

    use hos
    use spectral_hos
    implicit none
    real(wp), intent(in) :: time
    real(wp), intent(in), dimension(:,:) :: eta

    real(wp) x,y
    integer i,j
    real(wp), allocatable, dimension(:,:) :: eo
    real(wp), dimension(nxhos / 2) :: skx
    real(wp), dimension(nyhos / 2) :: sky
    !print *, 'plyudebug, nxhos,nyhos:', nxhos, nyhos 
    allocate(eo(nxhos,nyhos))
    !print *,'plyudebug_nl, size(eta):',size(eta,1),size(eta,2)
    !print *,'plyudebug_nl, size(eo)',size(eo,1),size(eo,2)
    call alltoone(eta,eo)
    if (myid_hos == 0) then
       write(932,*) " VARIABLES = x,y,eta"
       write(932,*) ' ZONE T="',time,'" I=', nxhos, ' J=',nyhos,' F=POINT'
       do j = 1, nyhos
          y = (j - 1) * dy_hos
          do i = 1, nxhos
             x = (i - 1) * dx_hos
             write(932,'(25e12.4)') x,y,eo(i,j)
          end do
       end do

       !export energy spectrum
       call spec_x(eo,skx)
       write(32,*) " VARIABLES = kx, Skx"
       write(32,*) ' ZONE T="',time,'" I=', nxhos/2
       do i = 1, nxhos / 2
          write(32,'(25e12.4)') i * pex_hos, skx(i)
       end do
       
    end if
    deallocate(eo)
   end subroutine outsurf

end module io_hos
