  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routine to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! X-pencil real arrays
  subroutine alloc_x(var, expand, opt_global)

    implicit none

    real(wp), allocatable, dimension(:,:,:) :: var
    integer, intent(IN), optional :: expand
    logical, intent(IN), optional :: opt_global

    logical :: global
    integer :: level , errorcode

    if (present(expand)) then
       level = expand
    else
       level = 0
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    if (global) then
       allocate(var(xst(1):xend(1),  &
                    xst(2):xend(2),  &
                    xst(3)-level:xend(3)+level), &
                    stat=errorcode)
    else
       allocate(var(1:xsz(1),1:xsz(2),1-level:xsz(3)+level), stat=errorcode)
    end if
    
    if (errorcode /= 0) then
       call decomp_abort(8, 'Out of memory creating new arrays')
    end if

  end subroutine alloc_x

  ! Y-pencil real arrays
  subroutine alloc_y(var, expand, opt_global)

    implicit none

    real(wp), allocatable, dimension(:,:,:) :: var
    integer, intent(IN), optional :: expand
    logical, intent(IN), optional :: opt_global

    logical :: global
    integer :: level, errorcode

    if (present(expand)) then
       level = expand
    else
       level = 0
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    if (global) then
       allocate(var(yst(1):yend(1),  &
                    yst(2):yend(2),  & 
                    yst(3)-level:yend(3)+level), &
                    stat=errorcode)
    else
       allocate(var(1:ysz(1),1:ysz(2),1-level:ysz(3)-level), stat=errorcode)
    end if
    
    if (errorcode /= 0) then
       call decomp_abort(8, 'Out of memory creating new arrays')
    end if

  end subroutine alloc_y

