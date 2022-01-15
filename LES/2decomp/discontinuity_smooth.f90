module discontinuity_smooth
  !> Added by plyu.
  !! Each row of limit1d is a two-element pair of lower and upper indices of
  !! smoothing range. 

  use decomp, only : wp, myid
  use fll_mods_m !, only : fll_mkdir, fll_mk, fll_mv
   
!  type limit_type
!    integer :: lim(2)
!  end type
!
!  type arr_limit_1d
!    integer :: n !> number of ranges to be smoothed
!    integer, dimension(:,:), allocatable :: lim !> size=(n,2)
!  end type  
!  
!  type arr_limit_2d
!    integer :: n !> number of 1D lines to be smoothed
!    integer, dimension(:), allocatable :: line_index !> which line to be smoothed
!    type(arr_limit_1d), dimension(:), allocatable :: linelim !> detail of one 1D line smoothing
!  end type
!
!  type arr_limit_3d
!    integer :: n 
!    integer, dimension(:), allocatable :: plane_index
!    type(arr_limit_2d), dimension(:), allocatable :: planelim 
!  end type

  type(dnode), pointer, public, save :: tbn_lim_x, tbn_lim_y

  public :: cubic_smooth_1d, fdm_in_smooth_1d, recover_in_smooth_1d
  public :: cubic_smooth_3d, postproc_smooth_3d
  public :: discontinuity_smooth_init, discontinuity_smooth_end
  public :: analyze_limit1d_from_1d
  public :: analyze_limit3d_from_3d

  integer, public, save :: n_buf_ds
  integer, public, save :: n_pdfx, n_pdfy, i_div_ustar, i_grad_p, lim_dir
  real, parameter :: onetwelfth = 0.083333333333333333
  real, parameter :: twothird = 0.666666666666666667

contains 
  
  subroutine discontinuity_smooth_init
    implicit none
    !integer :: n

    n_buf_ds = 4
    i_div_ustar = 0
    i_grad_p = 0
    n_pdfx = 0
    n_pdfy = 0
    lim_dir = 0
  end subroutine discontinuity_smooth_init

  subroutine discontinuity_smooth_end()
    implicit none
  end subroutine discontinuity_smooth_end
  
  !> analyze limit1d from 1d field.
  subroutine analyze_limit1d_from_1d (flag1d, lim1d, nsegs)
    implicit none

    real(wp), dimension(:), intent(in) :: flag1d
    type(dnode), pointer, intent(out) :: lim1d
    integer, intent(out) :: nsegs

    integer :: n1, i, x1, x2
    logical :: i_add

    type(dnode), pointer :: ptmp, ptmp2
    type(func_data_set) :: fpar
    logical :: ok

    n1 = size(flag1d, 1)
    nsegs = 0
    x1 = -1; x2 = -1
    i_add = .false.

    lim1d => fll_mkdir('limits_1d', fpar)
    
    do i = 4, n1-2
      
      if ( flag1d(i-1)<0.5 .and. flag1d(i)>0.5) then
        x1 = i - n_buf_ds
      elseif (flag1d(i-1)>0.5 .and. flag1d(i)<0.5) then
        x2 = i-1 + n_buf_ds
        i_add = .true.
      !elseif (flag1d(i-1) .and. flag1d(i) .and. i.eq.n1) then
      !  x2 = i
      !  i_add = .true.
      endif

      if (i_add) then
        nsegs = nsegs + 1
        if (x1 < 3) x1 = 3
        if (x2 > (n1-2)) x2 = n1 - 2 
        ptmp => fll_mk('lim1d_unit', 'I', 1_lint, 2_lint, fpar)
        ptmp%i1(:) = (/x1, x2/) 
        ok = fll_mv(ptmp, lim1d, fpar)
        i_add = .false.
      endif
    enddo

    !> merge those regions that are too close to each other.
    if (nsegs>=2) then
      ptmp => lim1d%pchild
      ptmp2 => ptmp%pnext
      do i = 1, nsegs - 1
        if (ptmp%i1(2) >= (ptmp2%i1(1)-2)) then
          !> Merge infor of two nodes into the first node
          ptmp%i1(2) = ptmp2%i1(2)
          
          !> remove the second node
          call fll_rm(ptmp2, fpar)

          !> update ptmp and ptmp2 for next loop
          ptmp2 => ptmp%pnext
        else
          ptmp => ptmp2
          ptmp2 => ptmp%pnext
        endif
      enddo
      nsegs = lim1d%ndim
    endif
  end subroutine analyze_limit1d_from_1d

  subroutine analyze_limit3d_from_3d (flag3d, lim3d)
    use decomp, only : xst, yst
    implicit none

    real(wp), dimension(:,:,:), intent(in) :: flag3d
    type(dnode), pointer, intent(out) :: lim3d

    type(dnode), pointer :: p_2d_indices, lim1d, ptmp 
    type(func_data_set) :: fpar
    logical :: ok

    integer :: n1, n2, n3, j, k, nsegs
    character(len=64) :: fn_lim3d
    !logical, dimension(:) :: flag1d

    n1 = size(flag3d, 1); n2 = size(flag3d, 2); n3 = size(flag3d, 3)
    
    !> clear the old lim3d, if it exists
    if (associated(lim3d)) then
      call fll_rm(lim3d, fpar)
    endif
    !> create the new list
    lim3d => fll_mkdir('lim3d', fpar)
    
    !> size(2d_index) = (1,3), contains (iy, iz, and n_x_segs)
    p_2d_indices => fll_mkdir('2d_indices', fpar) 
    ok = fll_mv(p_2d_indices, lim3d, fpar)

    do k = 1, n3
      do j = 1, n2
        call analyze_limit1d_from_1d (flag3d(:,j,k), lim1d, nsegs)
        if (nsegs > 0) then
          ptmp => fll_mk('2d_index', 'I', 3_lint, 1_lint, fpar)
          ptmp%i1(:) = (/j,k,nsegs/)
          ok = fll_mv (ptmp, p_2d_indices, fpar)

          ok = fll_mv (lim1d, lim3d, fpar)
        endif          
      enddo
    enddo

    !> for debug
    !if (n_pdfy == 39 .and. lim3d%pchild%ndim>0 ) then
    !  if (lim_dir .eq. 1) print *, 'xst', myid, lim3d%pchild%ndim, xst(1:3)
    !  if (lim_dir .eq. 2) print *, 'yst', myid, lim3d%pchild%ndim, yst(1:3)
    !  write(fn_lim3d, '(a, i0.1, a, i0.4, a)') 'lim3d_',lim_dir, '_', myid, '.fll'
    !  ok = fll_write(lim3d, fn_lim3d, 1000, 'A', fpar)
    !endif

  end subroutine analyze_limit3d_from_3d

  subroutine cubic_smooth_1d (fin, lim1d, idx_j, idx_k)
    implicit none
 
    real(wp), dimension(:), intent(inout) :: fin
    type(dnode), pointer, intent(in) :: lim1d
    integer :: idx_j, idx_k
    
    type(dnode), pointer :: lim1d_unit
 
    integer :: i, j, n1, i1, i2, n_t !, n2
    real(wp) :: xtmp1(4), ftmp1(4), a(4), xtmp2(10)
    real(wp), allocatable, dimension(:) :: xtmp3
    n1 = lim1d%ndim
    lim1d_unit => lim1d%pchild
    do i = 1, n1
      i1 = lim1d_unit%i1(1); i2 = lim1d_unit%i1(2)
      lim1d_unit => lim1d_unit%pnext

      xtmp1(1:2) = dble((/-2, -1/)) 
      xtmp1(3:4) = dble((/i2-i1+1, i2-i1+2/))
      ftmp1(1:2) = fin((i1-2):(i1-1))
      ftmp1(3:4) = fin((i2+1):(i2+2))
      call solve_cubic_interp_1d (xtmp1, ftmp1, a)
      
      !print *, 'myid=', myid, ', j=', idx_j, ', k=',idx_k, &
      !  ', i=', i, ', i1=', i1,&
      !  ', i2=', i2, ', ftmp=', ftmp1, ', cubic_coeff=', a

      n_t = i2-i1+1
      if (n_t<=10) then
        do j = 1, n_t
          xtmp2(j) = dble(j-1)
        enddo
        call func_cubic (xtmp2(1:n_t), fin(i1:i2), a)
        !print *, 'smoothed:', fin(i1:i2)
      else
        allocate(xtmp3(n_t))
        do j = 1, n_t
          xtmp3(j) = dble(j-1)
        enddo
        call func_cubic (xtmp3, fin(i1:i2), a)
        !print *, 'smoothed:', fin(i1:i2)
        deallocate(xtmp3)
      endif
    enddo
  end subroutine cubic_smooth_1d

  !> In-place cubic smoothing for 3D data
  subroutine cubic_smooth_3d (input, lim3d)
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT),  contiguous :: input
    type(dnode), pointer, intent(in) :: lim3d

    integer  :: i, j, k
    type(dnode), pointer :: p_2d_indices, p_2d_index, lim1d

    p_2d_indices => lim3d%pchild
    lim1d => p_2d_indices%pnext
    p_2d_index => p_2d_indices%pchild
    do i = 1, p_2d_indices%ndim
      j = p_2d_index%i1(1)
      k = p_2d_index%i1(2)
      !> plyunote: both input and flag3d consist of no ghost cell in vertical
      !!           direction, namely their size is usually
      !!           (xsz(1),xsz(2),xsz(3))
      call cubic_smooth_1d(input(:,j,k), lim1d, j, k)
      p_2d_index => p_2d_index%pnext
      lim1d => lim1d%pnext
    enddo
  end subroutine cubic_smooth_3d
  
  !> In the smoothed region, use 2nd order finite difference
  !! method to estimate the derivate dy/dx
  subroutine fdm_in_smooth_1d (fin, fout, lim1d, dx, order)
    implicit none
 
    real(wp), dimension(:), intent(in) :: fin
    real(wp), dimension(:), intent(inout) :: fout
    type(dnode), pointer, intent(in) :: lim1d
    real(wp), intent(in) :: dx
    integer, intent(in) :: order
    
    type(dnode), pointer :: lim1d_unit
 
    integer :: i, j, n1, i1, i2, n_t !, n2
    real(wp) :: dydx(10)
    real(wp), allocatable, dimension(:) :: dydx2
    n1 = lim1d%ndim
    lim1d_unit => lim1d%pchild
    do i = 1, n1
      !> for dydx, the range need to be corrected is 1 grid wider than y(x) on
      !! each side
      i1 = lim1d_unit%i1(1) - 1; i2 = lim1d_unit%i1(2) + 1 
      lim1d_unit => lim1d_unit%pnext

      n_t = i2-i1+1
      if (n_t<=10) then
        if (order .eq. 2) then
          do j = 1, n_t
            dydx(j) = (fin(i1+j) - fin(i1+j-2)) / 2.0_wp
          enddo
        else if (order .eq. 4) then
          do j = 1, n_t
            dydx(j) = - onetwelfth*fin(i1+j+1) + twothird*fin(i1+j) &
              - twothird*fin(i1+j-2) + onetwelfth*fin(i1+j-3)
          enddo
        endif
        fout(i1:i2) = dydx(1:n_t) / dx
      else
        allocate(dydx2(n_t))
        if (order .eq. 2) then
          do j = 1, n_t
            dydx2(j) = (fin(i1+j) - fin(i1+j-2)) / 2.0_wp
          enddo
        else if (order .eq. 4) then
          do j = 1, n_t
            dydx2(j) = - onetwelfth*fin(i1+j+1) + twothird*fin(i1+j) &
              - twothird*fin(i1+j-2) + onetwelfth*fin(i1+j-3)
          enddo
        endif
        fout(i1:i2) = dydx2(1:n_t) / dx
        deallocate(dydx2)
      endif
    enddo
  end subroutine fdm_in_smooth_1d

  subroutine recover_in_smooth_1d (fin, fout, lim1d)
    implicit none
 
    real(wp), dimension(:), intent(in) :: fin
    real(wp), dimension(:), intent(inout) :: fout
    type(dnode), pointer, intent(in) :: lim1d
    
    type(dnode), pointer :: lim1d_unit
 
    integer :: i, j, n1, i1, i2, n_t !, n2
    real(wp) :: dydx(10)
    real(wp), allocatable, dimension(:) :: dydx2
    n1 = lim1d%ndim
    lim1d_unit => lim1d%pchild
    do i = 1, n1
      i1 = lim1d_unit%i1(1); i2 = lim1d_unit%i1(2)
      fout(i1:i2) = fin(i1:i2)
      lim1d_unit => lim1d_unit%pnext
    enddo
  end subroutine recover_in_smooth_1d

  subroutine postproc_smooth_3d (input, output, lim3d, mode, dx)
    implicit none

    real(wp), dimension(:,:,:), intent(in), contiguous :: input
    real(wp), dimension(:,:,:), intent(inout), contiguous :: output
    type(dnode), pointer, intent(in) :: lim3d
    integer :: mode
    real(wp), optional :: dx

    integer  :: i, j, k
    type(dnode), pointer :: p_2d_indices, p_2d_index, lim1d

    p_2d_indices => lim3d%pchild
    lim1d => p_2d_indices%pnext
    p_2d_index => p_2d_indices%pchild
    
    if (mode .eq. 0) then
      !> mode 0: in the smoothed range, recover output to input 
      do i = 1, p_2d_indices%ndim
        j = p_2d_index%i1(1)
        k = p_2d_index%i1(2)
        call recover_in_smooth_1d(input(:,j,k), output(:,j,k), lim1d)
        p_2d_index => p_2d_index%pnext
        lim1d => lim1d%pnext
      enddo
    elseif (mode .eq. 1) then
      !> mode 1: in the smoothed range, replace output with dy/dx of input
      if (.not. present(dx)) then
        print *, 'The optional parameter dx is needed in mode 1 of ', &
          'postproc_smooth_3d'
        call exit(1)
      endif
      do i = 1, p_2d_indices%ndim
        j = p_2d_index%i1(1)
        k = p_2d_index%i1(2)
        !> plyunote: both input and flag3d consist of no ghost cell in vertical
        !!           direction.
        call fdm_in_smooth_1d(input(:,j,k), output(:,j,k), lim1d, dx, 4)
        p_2d_index => p_2d_index%pnext
        lim1d => lim1d%pnext
      enddo
    endif

  end subroutine postproc_smooth_3d
  
  subroutine func_cubic (xin, fout, a)
    implicit none
    real(wp), dimension(:), intent(in) :: xin, a
    real(wp), dimension(:), intent(out) :: fout
 
    fout = a(1) + a(2) * xin + a(3) * xin**2 + a(4) * xin**3
  end subroutine func_cubic
  
  !> Do the cubic fit of fin = a1 + a2*xin + a3*xin**2 + a4*xin**3
  subroutine solve_cubic_interp_1d (xin, fin, a)
    implicit none
 
    real(wp), intent(in) :: xin(4), fin(4)
    real(wp), intent(out) :: a(4)
 
    real(wp) :: al_A(4,4), al_A_inv(4,4)
    integer :: i, j
 
    do i = 1, 4
      al_A(i, 1) = 1.0_wp
      do j = 2, 4
        al_A(i, j) = xin(i)**(j-1)
      enddo
    enddo
    al_A_inv = matinv4(al_A)
    a = matmul(al_A_inv, fin)
  end subroutine solve_cubic_interp_1d
 
 !> Copied from http://fortranwiki.org/fortran/show/Matrix+inversion
 pure function matinv4(A) result(B)
   !! Performs a direct calculation of the inverse of a 4×4 matrix.
   real(wp), intent(in) :: A(4,4)   !! Matrix
   real(wp)             :: B(4,4)   !! Inverse matrix
   real(wp)             :: detinv
 
   ! Calculate the inverse determinant of the matrix
   detinv = &
     1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
      - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
      + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
      - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))
 
   ! Calculate the inverse of the matrix
   B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
   B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
   B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
   B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
   B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
   B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
   B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
   B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
   B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
   B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
   B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
   B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
   B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
   B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
   B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
   B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
 end function
 
end module discontinuity_smooth
