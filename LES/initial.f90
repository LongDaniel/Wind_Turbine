subroutine add_turbulence_velocity(pu, pv, pw, u, v, w, expand)
    use decomp
    use param
    use grid
    !use navier
    use spectral
#ifdef __INTEL_COMPILER
    use ifport
#endif

    implicit none

    integer :: i,j,k
    integer, intent(IN) :: expand
    real(wp) :: amg
    real(wp), dimension(xsz(1),xsz(2),1-expand:xsz(3)+expand) :: pu, pv, pw
    real(wp), dimension(xsz(1),xsz(2),1-expand:xsz(3)+expand) :: u, v, w
!print *, myid, 0.31
    do k=1,xsz(3)
      do j=1,xsz(2)
          do i=1,xsz(1)
              call random_number(pu(i,j,k))
              call random_number(pv(i,j,k))
              call random_number(pw(i,j,k))
           enddo 
        enddo 
    enddo
!print *, myid, 0.32
    call cutoff(pu(:,:,1:xsz(3)))
    call cutoff(pv(:,:,1:xsz(3)))
    call cutoff(pw(:,:,1:xsz(3)))
!print *, myid, 0.33
    do k=1, xsz(3)
       if(k==1) amg=0
       call ampran(zz(k)-zl,amg)
       pu(:,:,k)=amg*pu(:,:,k)
       pv(:,:,k)=amg*pv(:,:,k)
       pw(:,:,k)=amg*pw(:,:,k)
    end do 
!print *, myid, 0.34
    call update_ghost(pu,expand)
    call update_ghost(pv,expand)
    call update_ghost(pw,expand)
!print *, myid, 0.35
    call soleon(pu, pv, pw, u, v, w, expand)
!print *, myid, 0.36
    call rescale 

end subroutine add_turbulence_velocity

!checked
subroutine add_log_velocity(u,expand)
    use decomp
    use param
    !use navier 
    use grid, only : zz

    implicit none

    integer, intent(IN) :: expand
    integer :: k, nz2
    real(wp) :: zlbot, utop, z, zs, su, u1, u2

    real(wp), dimension(xsz(1),xsz(2),1-expand:xsz(3)+expand), intent(INOUT) :: u

    nz2=xsz(3)

    zlbot=hbar*resbot/(resbot+restop)

    if(restop.le.1.e-6) then
       utop=usbot*(2.5_wp*log(hbar/z0))
    else
       utop=usbot*((2.5_wp*log(hbar*resbot**2/(resbot+restop))+5.0_wp) &
           +(2.5_wp*log(hbar*restop**2/(resbot+restop))+5.0_wp)*restop &
           /resbot)
    endif

    do k=1,nz2
       z=zz(k)*hbar

       if(z.le.zlbot) then
          zs=z/zlsbot
          su=0
          u1=usbot*zs
          u2=usbot*(2.5_wp*log(z/z0))
          if(u1.lt.u2) su=u1
          if(u1.ge.u2) su=u2
          if(u2.lt.0) su=u1
       endif

       if(z.gt.zlbot) then
          !> plyunote: for cases restop<1.e-6, zlstop is 0, will cause error
          if(abs(zlstop).le.1.0e-6) then
            zlstop = 2.0e-6_wp
          endif

          zs=(hbar-z)/zlstop
          if(abs(zs).lt.1.e-6) then
             su=utop
          else
             u1=ustop*zs
             u2=ustop*(2.5_wp*log(z/z0))
             if(u1.lt.u2) su=utop-u1
             if(u1.ge.u2) su=utop-u2
             if(u2.lt.0) su=utop-u1
          endif
       endif
       
       if (isbot.and.k==1)then 
          u(:,:,k)=u(:,:,k)
       else
          u(:,:,k)=u(:,:,k)+su
       end if 

    enddo
end subroutine add_log_velocity

!checked
subroutine soleon(pu, pv, pw, u, v, w, expand)
    use decomp
    use param
    use grid
    use spectral 
    use solver_common
    !use navier
    
    implicit none
   
    integer :: nz2, k
    integer, intent(IN) :: expand
    real(wp), dimension(xsz(1), xsz(2),1-expand:xsz(3)+expand), intent(IN) :: pu, pv, pw
    real(wp), dimension(xsz(1), xsz(2),1-expand:xsz(3)+expand), intent(OUT) :: u, v, w 
    real(wp), dimension(xsz(1), xsz(2),1-expand:xsz(3)+expand) :: puy, pvx, pwx, pwy
    real(wp), dimension(xsz(1), xsz(2),1-expand:xsz(3)+expand) :: tmp1, tmp2, tmp3 
   
    nz2=xsz(3)

    ! compute dpu/dpsi
    call pdfy_x(pu(:,:,1:xsz(3)), puy(:,:,1:xsz(3)), pey)
    call update_ghost(puy, expand)
    ! compute dpv/dxi
    call pdfx(pv(:,:,1:xsz(3)), pvx(:,:,1:xsz(3)), pex)
    call update_ghost(pvx, expand)
    ! compute dpw/dxi
    call pdfx(pw(:,:,1:xsz(3)), pwx(:,:,1:xsz(3)), pex)
    call update_ghost(pwx, expand)
    ! compute dpw/dpsi
    call pdfy_x(pw(:,:,1:xsz(3)), pwy(:,:,1:xsz(3)), pey)
    call update_ghost(pwy, expand)

    ! compute dpu/dzeta
    call calc_uzeta(pu, tmp1(:,:,1:), expand)
    if(isbot)then 
        tmp1(:,:,2)=(pu(:,:,3)-pu(:,:,1))/(2.0_wp*dz(1)+dz(2))
    end if 
    if(istop)then 
        tmp1(:,:,nz2-1)=(pu(:,:,nz2)-pu(:,:,nz2-2))/(2.0_wp*dz(nz2-1)+dz(nz2-2))
    end if 
    call update_ghost(tmp1,expand)

    ! compute dpv/dzeta
    call calc_uzeta(pv, tmp2(:,:,1:), expand)
    if(isbot)then 
        tmp2(:,:,2)=(pv(:,:,3)-pv(:,:,1))/(2.0_wp*dz(1)+dz(2))
    end if 
    if(istop)then 
        tmp2(:,:,nz2-1)=(pv(:,:,nz2)-pv(:,:,nz2-2))/(2.0_wp*dz(nz2-1)+dz(nz2-2))
    end if 
    call update_ghost(tmp2,expand)

    ! compute dpw/dzeta
    call calc_uzeta(pw, tmp3(:,:,1:), expand)
    if(isbot)then 
        tmp3(:,:,2)=(pw(:,:,3)-pw(:,:,1))/(2.0_wp*dz(1)+dz(2))
    end if 
    if(istop)then 
        tmp3(:,:,nz2-1)=(pw(:,:,nz2)-pw(:,:,nz2-2))/(2.0_wp*dz(nz2-1)+dz(nz2-2))
    end if 
    call update_ghost(tmp3,expand)

    call dealiasxy(tmp1(:,:,1:nz2))
    call dealiasxy(tmp2(:,:,1:nz2))
    call dealiasxy(tmp3(:,:,1:nz2))



    do k=1, nz2
       ! u=dpw/dy-dpv/dz
       u(:,:,k)=pwy(:,:,k)+zetay(:,:,k)*tmp3(:,:,k)-her(:,:)*tmp2(:,:,k)
       ! u=dpu/dz-dpw/dx
       v(:,:,k)=her(:,:)*tmp1(:,:,k)-pwx(:,:,k)-zetax(:,:,k)*tmp3(:,:,k)
       ! u=dpv/dx-dpu/dy
       w(:,:,k)=pvx(:,:,k)+zetax(:,:,k)*tmp2(:,:,k)-puy(:,:,k)-zetay(:,:,k)*tmp1(:,:,k)
    end do 

    ! call dealiasxy(u(:,:,1:xsz(3)))
    ! call dealiasxy(v(:,:,1:xsz(3)))
    ! call dealiasxy(w(:,:,1:xsz(3)))

    call update_ghost(w,expand)
    
    do k=1,nz2
       w(:,:,k)=0.5_wp*(w(:,:,k)+w(:,:,k+1))
    end do 

    if(isbot)then 
        u(:,:,1)=0
        v(:,:,1)=0
        w(:,:,1)=0
    end if 

    if(istop)then 
        u(:,:,nz2)=0
        v(:,:,nz2)=0
        w(:,:,nz2)=0
        w(:,:,nz2-1)=0 
    end if 

end subroutine soleon 

!checked
subroutine rescale
    use decomp
    use param
    use navier 
    use mpi

    implicit none

    integer :: i, j, k, ierr
    real(wp) :: rmsmax_all, rmsmax, rmsmax_all2, rms, scal

    rmsmax_all=0
    rmsmax=0
    do k=1,xsz(3)
       
       rms=0
       do i=1,xsz(1)
          do j=1,xsz(2)
             rms=rms+u(i,j,k)**2+v(i,j,k)**2+0.25_wp*(w(i,j,k)+w(i,j,k-1))**2
          enddo
       enddo
       
       call mpi_allreduce(rms,rmsmax_all,1,mpi_double_precision, &
            mpi_sum,mpi_comm_2d_col,ierr)

       if(rmsmax.lt.rmsmax_all) then
          rmsmax=rmsmax_all
       endif
 
    enddo

    !why calculate rmsmax_all2
    call mpi_allreduce(rmsmax,rmsmax_all2,1,mpi_double_precision, &
         mpi_max,mpi_comm_world,ierr)

    rmsmax_all2=(rmsmax_all2/nx_global/ny_global/3.0_wp)**0.5_wp

    !> plyunote: avoid divided by zero
    scal=arm/(rmsmax_all2+1.0d-20)

    do k=1,xsz(3)
        u(:,:,k)=u(:,:,k)*scal
        v(:,:,k)=v(:,:,k)*scal
        w(:,:,k)=w(:,:,k)*scal
    enddo

end subroutine rescale

subroutine cutoff(f)
    ! can only be applied when parallel is in y direction 
  use MPI
  use decomp
  use param
  use fft_hos
  use spectral
  use utils

  implicit none
  
  real(wp), dimension(xsz(1), xsz(2), xsz(3)) :: f
  
  real(wp), dimension(nx_global,ny_global,nz_global) :: fall, tmpf
  real(wp), dimension(ny_global,nx_global,nz_global) :: ft1
  real(wp), dimension(nz_global,nx_global,ny_global) :: ft2
  
  integer*8 :: plan1, plan2
  
  real(wp) :: amx, amy, amz
  
  integer :: i,j,k,root0, ierror
  integer :: l, ll, m, mm, n, nn
  integer :: modx, mody, modz
  
  modx=18; mody=18; modz=7
  
  root0 = 0; !print *, myid, 0.321; call mpi_barrier(mpi_comm_world, ierror)
  call gather_3d_xyz(f,fall,root0)
  !print *, myid, 0.322; call mpi_barrier(mpi_comm_world, ierror)
  if (myid == root0) then
     ! 3Dfft
     do k = 1, nz_global
        call fft_for_x_hos(fall(:,:,k),nx_global,ny_global)       
        ft1(:,:,k) = transpose(fall(:,:,k))
        call fft_for_x_hos(ft1(:,:,k),ny_global,nx_global)
     end do
     
     do i = 1, nx_global
        !       do j = 1, ny_global
        !          do k = 1, nz_global
        !             ft2(k,i,j) = ft1(j,i,k)
        !          end do
        !       end do
        ft2(:,i,:) = transpose(ft1(:,i,:))
     end do
     
     do j = 1, ny_global
        call fft_for_x_hos(ft2(:,:,j),nz_global,nx_global)
     end do
     
     ! Smooth out small-scale structures
     do l=1,nx_global
        ll=(l-1)/2
        amx=exp(-(float(ll)/float(modx))**2)
        do m=1,ny_global
           mm=(m-1)/2
           amy=exp(-(float(mm)/float(mody))**2)
           do n=1,nz_global
              nn=(n-1)/2
              amz=exp(-(float(nn)/float(modz))**2)
              ft2(n,l,m)=amx*amy*amz*ft2(n,l,m)
              if(l.le.2.or.m.le.2.or.n.le.2) ft2(n,l,m)=0
           enddo
        enddo
     enddo
     
     ! inverse 3Dfft
     do j = 1, ny_global
        call fft_bac_x_hos(ft2(:,:,j),nz_global,nx_global)
     end do
     
     do i = 1, nx_global
        ft1(:,i,:) = transpose(ft2(:,i,:))
     end do
     
     do k = 1, nz_global
        call fft_bac_x_hos(ft1(:,:,k),ny_global,nx_global)
        fall(:,:,k) = transpose(ft1(:,:,k))
        call fft_bac_x_hos(fall(:,:,k),nx_global,ny_global)
     end do
  end if
  
  if (myid /= root0) fall=0.0
  !print *, myid, 0.323; call mpi_barrier(mpi_comm_world, ierror)
  
  call MPI_ALLREDUCE(fall,tmpf,nx_global*ny_global*nz_global,mpi_double_precision,mpi_sum, &
       mpi_comm_world,ierror)
  !print *, myid, 0.324; call mpi_barrier(mpi_comm_world, ierror)
    
  f = tmpf(xst(1):xend(1),xst(2):xend(2),xst(3):xend(3))

end subroutine cutoff
!checked
subroutine ampran(z,amg)
    use decomp, only : wp
    ! amplitude of random velocity field as function of depth (z)
  
    implicit none

    real(wp), parameter :: thick =0.1_wp
    real(wp), intent(IN) :: z
    real(wp), intent(OUT) :: amg
    ! data thick/0.1/
      
    amg=1.0_wp-exp(-(z/thick)**2)-exp(-((z+1.0_wp)/thick)**2)-0.005_wp

    if(amg.lt.0) amg=0

end subroutine ampran
