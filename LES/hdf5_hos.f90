!-----------------------------------------------------------
! MODULE: hdf_io
!
!> This module provides parallel HDF5 IO facilities for HOS
!
!> @author Anqing Xuan
!> @modified by William Xuanting Hao
!-----------------------------------------------------------
module hdf_io_hos

  use hos_param
  use MPI
  use hdf5

  implicit none

  private        ! Make everything private unless declared public

  integer :: error

  type, public :: HDFObj_hos
     character(len=256) :: filename
     integer(hid_t) :: file_id
  end type HDFObj_hos

  public :: IO_init_hos, IO_finalize_hos
  public :: writerOpen_hos, writerClose_hos
  public :: readerOpen_hos, readerClose_hos
  public :: write2D_xy_hos, read2D_xy_hos

contains
  subroutine IO_init_hos
    implicit none

    call h5open_f(error)

  end subroutine IO_init_hos

  subroutine IO_finalize_hos
    implicit none

    call h5close_f(error)

  end subroutine IO_finalize_hos

  !--------------------------------------------------------------

  subroutine writerOpen_hos(filename, comm, writerObj)
    implicit none

    integer, intent(IN) :: comm
    character(len=*), intent(IN) :: filename
    type(HDFObj_hos), intent(OUT) :: writerObj

    integer(hid_t) :: plist_id, file_id
    integer :: info

    info = MPI_INFO_NULL
    ! create a property list for the newly opened file
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)

    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)

    call h5pclose_f(plist_id, error)

    writerObj%filename = trim(filename)
    writerObj%file_id = file_id

  end subroutine writerOpen_hos

  !--------------------------------------------------------------

  subroutine writerClose_hos(writerObj)
    implicit none

    type(HDFObj_hos) :: writerObj

    integer :: error

    call h5fclose_f(writerObj%file_id, error)

  end subroutine writerClose_hos

  !--------------------------------------------------------------

  subroutine readerOpen_hos(filename, comm, readerObj)
    implicit none

    integer, intent(IN) :: comm
    character(len=*), intent(IN) :: filename
    type(HDFObj_hos), intent(OUT) :: readerObj

    integer(hid_t) :: plist_id, file_id
    integer :: info

    info = MPI_INFO_NULL
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)

    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)

    call h5pclose_f(plist_id, error)

    readerObj%filename = trim(filename)
    readerObj%file_id = file_id

  end subroutine readerOpen_hos

  !--------------------------------------------------------------

  subroutine readerClose_hos(writerObj)
    implicit none

    type(HDFObj_hos) :: writerObj

    integer :: error

    call h5fclose_f(writerObj%file_id, error)

  end subroutine readerClose_hos

  !--------------------------------------------------------------

  subroutine write2D_xy_hos(writerObj, varname, data)
    implicit none
    
    type(HDFObj_hos), intent(IN) :: writerObj
    character(len=*), intent(IN) :: varname
    real(wp), dimension(nxhos,nyhos/ncpu_hos), intent(IN) :: data
    
    integer(hid_t) :: dset_id, filespace, memspace
    integer(hsize_t), dimension(2) :: count, offset, dims
    
    dims(1) = nxhos
    dims(2) = nyhos
    
    ! create dataset
    call h5screate_simple_f(2, dims, filespace, error)
    call h5dcreate_f(writerObj%file_id, trim(varname), H5T_NATIVE_DOUBLE, filespace, &
         dset_id, error)
    call h5sclose_f(filespace, error)
    
    ! create memory space
    count(1) = nxhos
    count(2) = nyhos/ncpu_hos
    offset(1) = 0
    offset(2) = myid_hos * nyhos/ncpu_hos
    
    call h5screate_simple_f(2, count, memspace, error)
    
    ! set hyperslab
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
    
    ! write data independently
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
         file_space_id=filespace, mem_space_id=memspace)
    
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
    
  end subroutine write2D_xy_hos

  !--------------------------------------------------------------

  subroutine read2D_xy_hos(writerObj, varname, data)
    implicit none

    type(HDFObj_hos), intent(IN) :: writerObj
    character(len=*), intent(IN) :: varname
    real(wp), dimension(nxhos,nyhos/ncpu_hos), intent(out) :: data

    integer(hid_t) :: dset_id, filespace, memspace
    integer(hsize_t), dimension(2) :: count, offset, dims

    dims(1) = nxhos
    dims(2) = nyhos

    ! open dataset
    call h5dopen_f(writerObj%file_id, trim(varname), dset_id, error)

    ! create memory space
    count(1) = nxhos
    count(2) = nyhos/ncpu_hos
    offset(1) = 0
    offset(2) = myid_hos * nyhos/ncpu_hos 

    call h5screate_simple_f(2, count, memspace, error)

    ! set hyperslab
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

    ! write data independently
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
         file_space_id=filespace, mem_space_id=memspace)

    call h5sclose_f(filespace, error)
    call h5dclose_f(dset_id, error)

  end subroutine read2D_xy_hos

end module hdf_io_hos
