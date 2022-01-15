  !---------------------------------------------------------
  !> @brief Update ghost points out-of-place
  !
  !> @param[in]  in  the input array, should only include interior points
  !> @param[out] out the output array of the expanded size
  !> @param[in]  level the number of ghost grids to be expanded
  !---------------------------------------------------------
  subroutine update_ghost_outplace3(in, out, level)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    real(wp), dimension(:,:,1:), intent(IN) :: in    
    real(wp), dimension(:,:,1-level:), intent(OUT) :: out

    real(wp), allocatable, dimension(:,:,:) :: buf_s, buf_r
    integer :: n1, n2, n3

    integer :: k, ierror

    n1 = size(in, 1); n2 = size(in, 2); n3 = xsz(3)
    if (level > n3) call decomp_abort(10, 'Ghost cell exceeds local size')
    allocate(buf_s(n1,n2,level*2))
    allocate(buf_r(n1,n2,level*2))

    do k=1,level
       buf_s(:,:,k) = in(:,:,k)
       buf_s(:,:,k+level) = in(:,:,n3-k+1)
    end do

    buf_r = 0
    call MPI_NEIGHBOR_ALLTOALL(buf_s, n1*n2*level, real_type, &
                               buf_r, n1*n2*level, real_type, &
                               MPI_COMM_2D_ROW, ierror)

    out(:,:,1:n3) = in(:,:,1:n3)
    do k=1,level
       out(:,:,1-k) = buf_r(:,:,k)
       out(:,:,n3+k) = buf_r(:,:,k+level)
    end do

    deallocate(buf_s, buf_r)
  end subroutine update_ghost_outplace3

  !---------------------------------------------------------
  !> @brief Update ghost points in-place
  !
  !> @param[in,out] in the array to be expanded, the array should be of 
  !!   the expanded size
  !> @param[in]   level the number of ghost points
  !---------------------------------------------------------
  subroutine update_ghost_inplace3(in, level)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    real(wp), dimension(:,:,1-level:), intent(INOUT) :: in    

    real(wp), allocatable, dimension(:,:,:) :: buf_s, buf_r
    integer :: n1, n2, n3

    integer :: k, ierror

    n1 = size(in, 1); n2 = size(in, 2); n3 = xsz(3)
    if (level > n3) call decomp_abort(10, 'Ghost cell exceeds local size')
    allocate(buf_s(n1,n2,level*2))
    allocate(buf_r(n1,n2,level*2))

    do k=1,level
       buf_s(:,:,k) = in(:,:,k)
       buf_s(:,:,k+level) = in(:,:,n3-k+1)
    end do

    buf_r = 0
    call MPI_NEIGHBOR_ALLTOALL(buf_s, n1*n2*level, real_type, &
                               buf_r, n1*n2*level, real_type, &
                               MPI_COMM_2D_ROW, ierror)

    if (istop) then
       do k=1, level
          in(:,:,1-k) = buf_r(:,:,k)
       end do
    else if (isbot) then
       do k=1, level
          in(:,:,n3+k) = buf_r(:,:,k+level)
       end do
    else
       do k=1,level
          in(:,:,1-k) = buf_r(:,:,k)
          in(:,:,n3+k) = buf_r(:,:,k+level)
       end do
    end if

    deallocate(buf_s, buf_r)
  end subroutine update_ghost_inplace3

  subroutine update_ghost_outplace1(in, out, level)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    real(wp), dimension(1:), intent(IN) :: in    
    real(wp), dimension(1-level:), intent(OUT) :: out

    real(wp), allocatable, dimension(:) :: buf_s, buf_r
    integer :: n3

    integer :: k, ierror

    n3 = xsz(3)
    if (level > n3) call decomp_abort(10, 'Ghost cell exceeds local size')
    allocate(buf_s(level*2))
    allocate(buf_r(level*2))

    do k=1,level
       buf_s(k) = in(k)
       buf_s(k+level) = in(n3-k+1)
    end do

    buf_r = 0
    call MPI_NEIGHBOR_ALLTOALL(buf_s, level, real_type, &
                               buf_r, level, real_type, &
                               MPI_COMM_2D_ROW, ierror)

    out(1:n3) = in(1:n3)
    do k=1,level
       out(1-k) = buf_r(k)
       out(n3+k) = buf_r(k+level)
    end do

    deallocate(buf_s, buf_r)
  end subroutine update_ghost_outplace1

  subroutine update_ghost_inplace1(in, level)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    real(wp), dimension(1-level:), intent(INOUT) :: in    

    real(wp), allocatable, dimension(:) :: buf_s, buf_r
    integer :: n3

    integer :: k, ierror

    n3 = xsz(3)
    if (level > n3) call decomp_abort(10, 'Ghost cell exceeds local size')
    allocate(buf_s(level*2))
    allocate(buf_r(level*2))

    do k=1,level
       buf_s(k) = in(k)
       buf_s(k+level) = in(n3-k+1)
    end do

    buf_r = 0
    call MPI_NEIGHBOR_ALLTOALL(buf_s, level, real_type, &
                               buf_r, level, real_type, &
                               MPI_COMM_2D_ROW, ierror)

    if (istop) then
       do k=1, level
          in(1-k) = buf_r(k)
       end do
    else if (isbot) then
       do k=1, level
          in(n3+k) = buf_r(k+level)
       end do
    else
       do k=1,level
          in(1-k) = buf_r(k)
          in(n3+k) = buf_r(k+level)
       end do
    end if
    
    deallocate(buf_s, buf_r)
   
  end subroutine update_ghost_inplace1

