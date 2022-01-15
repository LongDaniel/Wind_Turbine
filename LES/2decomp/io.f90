!-----------------------------------------------------------
! MODULE: hdf_io
!
!> This module provides parallel HDF5 IO facilities for based on
!! 2D decomposition.
!
!> @author Anqing Xuan
!-----------------------------------------------------------
module hdf_io

  use decomp
  use MPI
  use hdf5

  implicit none

  private        ! Make everything private unless declared public

  integer :: error

  type, public :: HDFObj
     character(len=256) :: filename
     integer(hid_t) :: file_id
  end type HDFObj

  public :: IO_init, IO_finalize
  public :: writerOpen, writerClose
  public :: readerOpen, readerClose
  public :: write1D_z, read1D_z
  public :: write2D_xy, read2D_xy
  public :: write3D_xyz, read3D_xyz

contains
   subroutine IO_init
      implicit none

      call h5open_f(error)

   end subroutine

   subroutine IO_finalize
      implicit none

      call h5close_f(error)

   end subroutine

   subroutine writerOpen(filename, comm, writerObj)
      implicit none

      integer, intent(IN) :: comm
      character(len=*), intent(IN) :: filename
      type(HDFObj), intent(OUT) :: writerObj

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

    end subroutine writerOpen

    subroutine writerClose(writerObj)
       implicit none

       type(HDFObj) :: writerObj

       integer :: error

       call h5fclose_f(writerObj%file_id, error)

    end subroutine

    subroutine write2D_xy(writerObj, varname, data, id2)
       implicit none

       type(HDFObj), intent(IN) :: writerObj
       character(len=*), intent(IN) :: varname
       real(wp), dimension(xsz(1),xsz(2)), intent(IN) :: data
       integer , intent(IN) :: id2

       integer(hid_t) :: dset_id, filespace, memspace
       integer(hsize_t), dimension(2) :: count, offset, dims

       dims(1) = xsz(1); dims(2) = ysz(1)

       ! create dataset
       call h5screate_simple_f(2, dims, filespace, error)
       call h5dcreate_f(writerObj%file_id, trim(varname), H5T_NATIVE_DOUBLE, filespace, &
                        dset_id, error)
       call h5sclose_f(filespace, error)

       ! create memory space
       if (myid2 == id2) then
          count(1) = xsz(1)
          count(2) = xsz(2)
          offset(1) = 0
          offset(2) = xst(2)-1
       else
          count(1) = 0
          count(2) = 0
          offset(1) = 0
          offset(2) = 0
       end if
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

    end subroutine write2D_xy

    subroutine write3D_xyz(writerObj, varname, data)
       implicit none

       type(HDFObj), intent(IN) :: writerObj
       character(len=*), intent(IN) :: varname
       real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(IN) :: data

       integer(hid_t) :: dset_id, filespace, memspace
       integer(hsize_t), dimension(3) :: count, offset, dims

       dims(1) = xsz(1); dims(2) = ysz(1); dims(3) = nz_global

       ! create dataset
       call h5screate_simple_f(3, dims, filespace, error)
       call h5dcreate_f(writerObj%file_id, trim(varname), H5T_NATIVE_DOUBLE, filespace, &
                        dset_id, error)
       call h5sclose_f(filespace, error)

       ! create memory space
       count(1) = xsz(1)
       count(2) = xsz(2)
       count(3) = xsz(3)
       offset(1) = 0
       offset(2) = xst(2)-1
       offset(3) = xst(3)-1
       call h5screate_simple_f(3, count, memspace, error)

       ! set hyperslab
       call h5dget_space_f(dset_id, filespace, error)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
       
       ! write data independently
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                       file_space_id=filespace, mem_space_id=memspace)
       
       call h5sclose_f(filespace, error)
       call h5sclose_f(memspace, error)
       call h5dclose_f(dset_id, error)

    end subroutine write3D_xyz

    subroutine write1D_z(writerObj, varname, data, id1)
       implicit none

       integer, intent(IN) :: id1
       type(HDFObj), intent(IN) :: writerObj
       character(len=*), intent(IN) :: varname
       real(wp), dimension(xsz(3)), intent(IN) :: data

       integer(hid_t) :: dset_id, filespace, memspace
       integer(hsize_t), dimension(1) :: count, offset, dims

       dims(1) = nz_global

       ! create dataset
       call h5screate_simple_f(1, dims, filespace, error)
       call h5dcreate_f(writerObj%file_id, trim(varname), H5T_NATIVE_DOUBLE, filespace, &
                        dset_id, error)
       call h5sclose_f(filespace, error)

       ! create memory space
       if (myid1 == id1) then
          count(1) = xsz(3)
          offset(1) = xst(3)-1
       else
          count(1) = 0
          offset(1) = 0
       end if
       call h5screate_simple_f(1, count, memspace, error)

       ! set hyperslab
       call h5dget_space_f(dset_id, filespace, error)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

       ! write data independently
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                       file_space_id=filespace, mem_space_id=memspace)

       call h5sclose_f(filespace, error)
       call h5sclose_f(memspace, error)
       call h5dclose_f(dset_id, error)

    end subroutine write1D_z

   subroutine readerOpen(filename, comm, readerObj)
      implicit none

      integer, intent(IN) :: comm
      character(len=*), intent(IN) :: filename
      type(HDFObj), intent(OUT) :: readerObj

      integer(hid_t) :: plist_id, file_id
      integer :: info

      info = MPI_INFO_NULL
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      if (error .ne. 0) print *, "h5pcreate_f error"
      call h5pset_fapl_mpio_f(plist_id, comm, info, error)
      if (error .ne. 0) print *, "h5pset_fapl_mpio_f error"

      call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
      if (error .ne. 0) print *, "h5fopen_f error"

      call h5pclose_f(plist_id, error)
      if (error .ne. 0) print *, "h5pclose_f error"

      readerObj%filename = trim(filename)
      readerObj%file_id = file_id

    end subroutine readerOpen

    subroutine readerClose(writerObj)
       implicit none

       type(HDFObj) :: writerObj

       integer :: error

       call h5fclose_f(writerObj%file_id, error)

    end subroutine readerClose

    subroutine read2D_xy(writerObj, varname, data, id2)
       implicit none

       type(HDFObj), intent(IN) :: writerObj
       character(len=*), intent(IN) :: varname
       real(wp), dimension(xsz(1),xsz(2)), intent(OUT) :: data
       integer , intent(IN) :: id2

       integer(hid_t) :: dset_id, filespace, memspace
       integer(hsize_t), dimension(2) :: count, offset, dims

       dims(1) = xsz(1); dims(2) = ysz(1)

       ! open dataset
       call h5dopen_f(writerObj%file_id, trim(varname), dset_id, error)

       ! create memory space
       if (id2 == -1 .or. myid2 == id2) then
          count(1) = xsz(1)
          count(2) = xsz(2)
          offset(1) = 0
          offset(2) = xst(2)-1
       else
          count(1) = 0
          count(2) = 0
          offset(1) = 0
          offset(2) = 0
       end if
       call h5screate_simple_f(2, count, memspace, error)

       ! set hyperslab
       call h5dget_space_f(dset_id, filespace, error)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

       ! read data independently
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                      file_space_id=filespace, mem_space_id=memspace)

       call h5sclose_f(filespace, error)
       call h5dclose_f(dset_id, error)

    end subroutine read2D_xy

    subroutine read3D_xyz(writerObj, varname, data)
       implicit none

       type(HDFObj), intent(IN) :: writerObj
       character(len=*), intent(IN) :: varname
       real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: data

       integer(hid_t) :: dset_id, filespace, memspace
       integer(hsize_t), dimension(3) :: count, offset, dims

       dims(1) = xsz(1); dims(2) = ysz(1); dims(3) = nz_global

       ! create dataset
       call h5dopen_f(writerObj%file_id, trim(varname), dset_id, error)
       if (error .ne. 0) print *, "error in h5dopen_f"

       ! create memory space
       count(1) = xsz(1)
       count(2) = xsz(2)
       count(3) = xsz(3)
       offset(1) = 0
       offset(2) = xst(2)-1
       offset(3) = xst(3)-1
       call h5screate_simple_f(3, count, memspace, error)
       if (error .ne. 0) print *, "error in h5screate_simple_f"

       ! set hyperslab
       call h5dget_space_f(dset_id, filespace, error)
       if (error .ne. 0) print *, "error in h5dget_space_f"
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
       if (error .ne. 0) print *, "error in h5sselect_hyperslab_f"

       ! read data independently
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                      file_space_id=filespace, mem_space_id=memspace)
       if (error .ne. 0) print *, "error in h5dread_f"

       call h5sclose_f(filespace, error)
       if (error .ne. 0) print *, "error in h5sclose_f"
       call h5dclose_f(dset_id, error)
       if (error .ne. 0) print *, "error in h5dclose_f"

    end subroutine read3D_xyz

    subroutine read1D_z(writerObj, varname, data, id1)
       implicit none

       integer, intent(IN) :: id1
       type(HDFObj), intent(IN) :: writerObj
       character(len=*), intent(IN) :: varname
       real(wp), dimension(xsz(3)), intent(OUT) :: data

       integer(hid_t) :: dset_id, filespace, memspace
       integer(hsize_t), dimension(1) :: count, offset, dims

       dims(1) = nz_global

       ! create dataset
       call h5dopen_f(writerObj%file_id, trim(varname), dset_id, error)

       ! create memory space
       if (id1 == -1 .or. myid1 == id1) then
          count(1) = xsz(3)
          offset(1) = xst(3)-1
       else
          count(1) = 0
          offset(1) = 0
       end if
       call h5screate_simple_f(1, count, memspace, error)

       ! set hyperslab
       call h5dget_space_f(dset_id, filespace, error)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

       ! read data independently
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                      file_space_id=filespace, mem_space_id=memspace)

       call h5sclose_f(filespace, error)
       call h5dclose_f(dset_id, error)

    end subroutine read1D_z

end module hdf_io
