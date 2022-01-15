module utils
   implicit none

   public :: sumall
   public :: max_abs_diff
   public :: check_nan
   public :: sum_plane_xy

   interface check_nan
      module procedure check_nan_sp2, check_nan_dp2
      module procedure check_nan_sp3, check_nan_dp3
   end interface

   interface sumall
      module procedure sumall_2d_scalar
      module procedure sumall_2d_1d
   end interface sumall

   interface get_average_y
      module procedure get_average_y_2d, get_average_y_3d
   end interface get_average_y

   interface get_average_x
      module procedure get_average_x_3d
   end interface get_average_x

contains

   function sum_plane_xy(f, mpicomm) result(s)
      use MPI
      use decomp, only : wp, real_type

      implicit none

      integer,  intent(IN) :: mpicomm
      real(wp), dimension(:,:), intent(IN) :: f
      real(wp) :: s

      integer  :: ierror
      real(wp) :: slocal

      slocal = sum(f(:,:))
      call MPI_Allreduce(slocal, s, 1, real_type, MPI_SUM, mpicomm, ierror)

   end function sum_plane_xy

   subroutine gather_1d_z(in, out, root2)
      use MPI
      use decomp

      implicit none

      integer, optional :: root2
      real(wp), dimension(xsz(3)), intent(IN) :: in
      real(wp), dimension(nz_global), intent(OUT) :: out

      integer :: root, ierror
      real(wp), dimension(nz_global) :: tmp

      if (present(root2)) then
         root = root2
      else
         root = -1
      end if

      tmp = 0
      tmp(xst(3):xend(3)) = in(:)
      if (root < 0 .or. root >= nproc2) then
         call MPI_Allreduce(tmp, out, nz_global, real_type, &
                            MPI_SUM, MPI_COMM_2D_ROW, ierror)
      else
         call MPI_Reduce(tmp, out, nz_global, real_type, &
                         MPI_SUM, root, MPI_COMM_2D_ROW, ierror)
      end if

   end subroutine gather_1d_z

   subroutine gather_2d_xy(in, out, root1)
      use MPI
      use decomp

      implicit none

      integer, optional :: root1
      real(wp), dimension(xsz(1),xsz(2)), intent(IN) :: in
      real(wp), dimension(nx_global,ny_global), intent(OUT) :: out

      integer :: root, ierror
      real(wp), dimension(nx_global,ny_global) :: tmpall

      if (present(root1)) then
         root = root1
      else
         root = -1
      end if

      tmpall = 0
      tmpall(:,xst(2):xend(2)) = in(:,:)
      if (root < 0 .or. root >= nproc1) then
         call MPI_Allreduce(tmpall, out, nx_global*ny_global, real_type, &
                            MPI_SUM, MPI_COMM_2D_COL, ierror)
      else
         call MPI_Reduce(tmpall, out, nx_global*ny_global, real_type, &
                            MPI_SUM, root, MPI_COMM_2D_COL, ierror)
      end if

   end subroutine gather_2d_xy

   subroutine gather_2d_xz(in, out, root2)
      use MPI
      use decomp

      implicit none

      integer, optional :: root2
      real(wp), dimension(xsz(1),xsz(3)), intent(IN) :: in
      real(wp), dimension(nx_global,nz_global), intent(OUT) :: out

      integer :: root, ierror
      real(wp), dimension(nx_global,nz_global) :: tmpall

      if (present(root2)) then
         root = root2
      else
         root = -1
      end if

      tmpall = 0
      tmpall(:,xst(3):xend(3)) = in(:,:)
      if (root < 0 .or. root >= nproc2) then
         call MPI_Allreduce(tmpall, out, nx_global*nz_global, real_type, &
                            MPI_SUM, MPI_COMM_2D_ROW, ierror)
      else
         call MPI_Reduce(tmpall, out, nx_global*nz_global, real_type, &
                            MPI_SUM, root, MPI_COMM_2D_ROW, ierror)
      end if

   end subroutine gather_2d_xz

   subroutine gather_3d_xyz(in, out, root0)
      use MPI
      use decomp

      implicit none

      integer, optional :: root0
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(IN) :: in
      real(wp), dimension(nx_global,ny_global,nz_global), intent(OUT) :: out

      integer :: i, ierror, root, coord(2)
      integer, dimension(-1:nproc1-1) :: xend2_all
      integer, dimension(-1:nproc2-1) :: xend3_all
      real(wp), dimension(:,:,:), allocatable :: tmpall

      if (present(root0)) then
         root = root0
      else
         root = -1
      end if

      if (root < 0 .or. root >= nproc) then
         allocate(tmpall(nx_global, ny_global, nz_global))
         tmpall = 0
         tmpall(xst(1):xend(1),xst(2):xend(2),xst(3):xend(3)) = in(:,:,:)
         call MPI_Allreduce(tmpall, out, nx_global*ny_global*nz_global, real_type, &
                            MPI_SUM, MPI_COMM_2D_COL, ierror)
      else
         xend2_all(0:) = ydist
         xend2_all(-1) = 0
         do i=1, nproc1-1
            xend2_all(i) = xend2_all(i-1)+ydist(i)
         end do
         xend3_all(0:) = zdist
         xend3_all(-1) = 0
         do i=1, nproc2-1
            xend3_all(i) = xend3_all(i-1)+zdist(i)
         end do

         if (myid /= root) then
            call MPI_Send(in, xsz(1)*xsz(2)*xsz(3), real_type, & 
                          root, myid, MPI_COMM_2D_CART, ierror)
         else
            do i=0, nproc-1
               call MPI_Cart_coords(MPI_COMM_2D_CART, i, 2, coord, ierror)
               allocate(tmpall(xsz(1),ydist(coord(2)),zdist(coord(1))))
               if (i/=root) then
                  call MPI_Recv(tmpall, xsz(1)*xsz(2)*xsz(3), real_type, &
                                i, i, MPI_COMM_2D_CART, MPI_STATUS_IGNORE, ierror)
               else
                  tmpall = in
               end if
               out(1:nx_global, &
                   xend2_all(coord(2)-1)+1:xend2_all(coord(2)), &
                   xend3_all(coord(1)-1)+1:xend3_all(coord(1))) = tmpall
               deallocate(tmpall)
            end do
         end if
         ! call MPI_Reduce(tmpall, out, nx_global*ny_global*nz_global, real_type, &
         !                    MPI_SUM, root, MPI_COMM_2D_COL, ierror)
      end if

   end subroutine gather_3d_xyz

   function sumall_2d_scalar(f, mpicomm) result(s)
      use MPI
      use decomp, only : wp, real_type

      implicit none

      integer,  intent(IN) :: mpicomm
      real(wp), dimension(:,:), intent(IN) :: f
      real(wp) :: s

      integer  :: ierror
      real(wp) :: slocal

      slocal = sum(f(:,:))
      call MPI_Allreduce(slocal, s, 1, real_type, MPI_SUM, mpicomm, ierror)

   end function sumall_2d_scalar

   function sumall_2d_1d(f, mpicomm, dims) result(s)
      use MPI
      use decomp

      implicit none

      integer, intent(IN) :: mpicomm, dims
      real(wp), dimension(:,:), intent(IN) :: f
      real(wp), dimension(:), allocatable :: s

      integer :: ierror
      real(wp), dimension(:), allocatable :: slocal

      if (dims == 1) then
         allocate(slocal(ny_global), s(ny_global))
         slocal = 0
         slocal(xst(2):xend(2)) = sum(f(:,:), dims)
         call MPI_Allreduce(slocal, s, ny_global, real_type, MPI_SUM, mpicomm, ierror)
      else if (dims == 2) then
         allocate(slocal(nx_global), s(nx_global))
         slocal = sum(f(:,:), dims)
         call MPI_Allreduce(slocal, s, nx_global, real_type, MPI_SUM, mpicomm, ierror)
      else
         return
      end if

   end function sumall_2d_1d

   subroutine get_average_x_3d(f, fm)
      use decomp
      use MPI
      implicit none

      real(wp), dimension(:,:,:), intent(IN)  :: f
      real(wp), dimension(:,:),   intent(OUT) :: fm

      integer :: k, ierror
      real(wp), dimension(:,:), allocatable :: tmpall

      allocate(tmpall(ny_global, nz_global))
      tmpall = 0
      do k=xst(3), xend(3)
         tmpall(1:ny_global, k) = sumall(f(:,:,k-xst(3)+1), MPI_COMM_2D_COL, 1)/nx_global
      end do
      call MPI_Allreduce(tmpall, fm, ny_global*nz_global, real_type, &
                         MPI_SUM, MPI_COMM_2D_ROW, ierror)

   end subroutine get_average_x_3d

   subroutine get_average_y_3d(f, fm)
      use decomp
      use MPI
      implicit none

      real(wp), dimension(:,:,:), intent(IN)  :: f
      real(wp), dimension(:,:),   intent(OUT) :: fm

      integer :: k, ierror
      real(wp), dimension(:,:), allocatable :: tmpall

      allocate(tmpall(nx_global, nz_global))
      tmpall = 0
      do k=xst(3), xend(3)
         tmpall(1:nx_global, k) = sumall(f(:,:,k-xst(3)+1), MPI_COMM_2D_COL, 2)/ny_global
      end do
      call MPI_Allreduce(tmpall, fm, nx_global*nz_global, real_type, &
                         MPI_SUM, MPI_COMM_2D_ROW, ierror)

   end subroutine get_average_y_3d

   subroutine get_average_y_2d(f, fm)
      use decomp
      use MPI
      implicit none

      real(wp), dimension(:,:), intent(IN)  :: f
      real(wp), dimension(:),   intent(OUT) :: fm

      ! real(wp), dimension(:), allocatable :: tmpall

      fm(1:nx_global) = sumall(f, MPI_COMM_2D_COL, 2)/ny_global

   end subroutine get_average_y_2d

   subroutine get_average_xy_2d(f, fm)
      use decomp
      use MPI
      implicit none

      real(wp), dimension(:,:), intent(IN)  :: f
      real(wp), intent(OUT) :: fm

      fm = sumall(f, MPI_COMM_2D_COL)/nx_global/ny_global

   end subroutine get_average_xy_2d

   subroutine get_average_xy_3d(f, fm)
      use decomp
      use MPI
      implicit none

      real(wp), dimension(:,:,:), intent(IN)  :: f
      real(wp), dimension(:),     intent(OUT) :: fm

      integer :: k
      real(wp), dimension(xsz(3)) :: tmp

      do k=1, xsz(3)
         tmp(k) = sumall(f(:,:,k), MPI_COMM_2D_COL)/nx_global/ny_global
      end do
      
      fm(1:xsz(3))=tmp(1:xsz(3))
      ! call gather_1d_z(tmp, fm)

   end subroutine get_average_xy_3d

   function max_abs_diff(a, b, mpicomm) result(maxdiff)
      use decomp, only : wp, real_type
      use MPI
      implicit none

      integer, intent(IN), optional :: mpicomm
      real(wp), dimension(:,:,:), contiguous, intent(IN) :: a, b
      real(wp) :: maxdiff

      integer  :: i,j,k
      real(wp) :: tmp

      maxdiff = 0
      do k=1, size(a,3)
         do j=1, size(a,2)
            do i=1, size(a,1)
               tmp = abs(a(i,j,k)-b(i,j,k))
               if (tmp > maxdiff) maxdiff = tmp
            end do
         end do
      end do

      if (present(mpicomm)) then
         call MPI_Allreduce(maxdiff, tmp, 1, real_type, MPI_MAX, mpicomm, i)
         maxdiff = tmp
      end if

   end function max_abs_diff

   function check_nan_sp2(f) result(existnan)
      use iso_fortran_env, only : wp=>REAL32
      implicit none

      real(wp), dimension(:,:), intent(IN) :: f
      logical :: existnan

      integer :: i, j

      existnan = .false.
      do j=1, size(f,2)
         do i=1, size(f,1)
            existnan = ((f(i,j)/=f(i,j)) .or. existnan)
         end do
      end do

   end function check_nan_sp2

   function check_nan_dp2(f) result(existnan)
      use iso_fortran_env, only : wp=>REAL64
      implicit none

      real(wp), dimension(:,:), intent(IN) :: f
      logical :: existnan

      integer :: i, j

      existnan = .false.
      do j=1, size(f,2)
         do i=1, size(f,1)
            existnan = ((f(i,j)/=f(i,j)) .or. existnan)
         end do
      end do

   end function check_nan_dp2

   function check_nan_sp3(f) result(existnan)
      use iso_fortran_env, only : wp=>REAL32
      implicit none

      real(wp), dimension(:,:,:), intent(IN) :: f
      logical :: existnan

      integer :: i, j, k

      existnan = .false.
      do k=1, size(f,3)
         do j=1, size(f,2)
            do i=1, size(f,1)
               existnan = ((f(i,j,k)/=f(i,j,k)) .or. existnan)
            end do
         end do
      end do

   end function check_nan_sp3

   function check_nan_dp3(f) result(existnan)
      use iso_fortran_env, only : wp=>REAL64
      implicit none

      real(wp), dimension(:,:,:), intent(IN) :: f
      logical :: existnan

      integer :: i, j, k

      existnan = .false.
      do k=1, size(f,3)
         do j=1, size(f,2)
            do i=1, size(f,1)
               existnan = ((f(i,j,k)/=f(i,j,k)) .or. existnan)
            end do
         end do
      end do

   end function check_nan_dp3

   subroutine error_abort(errorcode, msg)
      use MPI
      implicit none

      integer, intent(IN) :: errorcode
      character(len=*), intent(IN), optional :: msg

      integer :: ierror
      
      if (present(msg)) then
         write(*,'(A6,I2.2,A,A)') 'ERROR(',errorcode,'): ',msg
      end if
      call MPI_Abort(MPI_COMM_WORLD, errorcode, ierror)

   end subroutine error_abort

end module utils
