module smooth

  use hos_param
  use fft_hos
  use spectral_hos
  implicit none

  public :: filter_lp, filter_hp

  interface filter_lp
     module procedure filter_lp1, filter_lp2, filter_lp3, filter_lp4, filter_lp3_aux
  end interface filter_lp

contains

  subroutine filter_hp(eta,vps,fac)

    implicit none
    real(wp), intent(in) :: fac
    real(wp), intent(inout), dimension(nxhos,nyhos/ncpu_hos) :: eta,vps
    real(wp) nxpat,nypat
    integer i,j

    i = 0
    j = 0

    call fft_for_x_hos(eta)
    call fft_for_x_hos(vps)

    nxpat = nxhos * fac
    nypat = nyhos * fac

    do i = 1, nxhos
       do j = 1, nyhos/ncpu_hos
          if (i > nxpat) then
             eta(i,j) = 0
             vps(i,j) = 0
          end if
       end do
    end do

    call transpose_2d(eta,bufb_hos,nxhos,nyhos/ncpu_hos)
    call fft_for_x_hos(bufb_hos)
    do j = 1, nyhos
       do i = 1, nxhos/ncpu_hos
          if (j > nypat) then
             bufb_hos(j,i) = 0
          end if
       end do
    end do
    call fft_bac_x_hos(bufb_hos)
    call transpose_2d(bufb_hos,eta,nyhos,nxhos/ncpu_hos)
    call fft_bac_x_hos(eta)

    call transpose_2d(vps,bufb_hos,nxhos,nyhos/ncpu_hos)
    call fft_for_x_hos(bufb_hos)
    do j = 1, nyhos
       do i = 1, nxhos/ncpu_hos
          if (j > nypat) then
             bufb_hos(j,i) = 0
          end if
       end do
    end do
    call fft_bac_x_hos(bufb_hos)
    call transpose_2d(bufb_hos,vps,nyhos,nxhos/ncpu_hos)
    call fft_bac_x_hos(vps)

  end subroutine filter_hp

  !-----------------------------------------------------------

  subroutine filter_lp1(eta,vps,fac)

    implicit none
    real(wp), intent(in) :: fac
    real(wp), intent(inout), dimension(nxhos,nyhos/ncpu_hos) :: eta,vps
    real(wp) nxpat,nypat
    integer i,j

    i = 0
    j = 0

    call fft_for_x_hos(eta)
    call fft_for_x_hos(vps)

    nxpat = nxhos * fac
    nypat = nyhos * fac

    do i = 1, nxhos
       do j = 1, nyhos/ncpu_hos
          if (i > nxpat) then
             eta(i,j) = 0
             vps(i,j) = 0
          end if
       end do
    end do

    call transpose_2d(eta,bufb_hos,nxhos,nyhos/ncpu_hos)
    call fft_for_x_hos(bufb_hos)
    do j = 1, nyhos
       do i = 1, nxhos/ncpu_hos
          if (j > nypat) then
             bufb_hos(j,i) = 0             
          end if
       end do
    end do
    call fft_bac_x_hos(bufb_hos)
    call transpose_2d(bufb_hos,eta,nyhos,nxhos/ncpu_hos)
    call fft_bac_x_hos(eta)

    call transpose_2d(vps,bufb_hos,nxhos,nyhos/ncpu_hos)
    call fft_for_x_hos(bufb_hos)
    do j = 1, nyhos
       do i = 1, nxhos/ncpu_hos
          if (j > nypat) then
             bufb_hos(j,i) = 0
          end if
       end do
    end do
    call fft_bac_x_hos(bufb_hos)
    call transpose_2d(bufb_hos,vps,nyhos,nxhos/ncpu_hos)
    call fft_bac_x_hos(vps)

  end subroutine filter_lp1

  !-----------------------------------------------------------

  !-----------------------------------------------------------

  subroutine filter_lp2(eta,vps,wvn,ratio)

    implicit none
    real(wp), intent(in) :: ratio
    real(wp), intent(inout), dimension(nxhos,nyhos/ncpu_hos) :: eta,vps
    real(wp), intent(in), dimension(nxhos,nyhos/ncpu_hos,npw) :: wvn
    real(wp)  fact_fun, ax, ay, wvmax, wvf, delta, dwvn, factor
    integer modex,modey
    integer i,j, l, m

    call fft_for_xy_hos(eta)
    call fft_for_xy_hos(vps)

    modex=((nxhos+1)/2)/3*2
    modey=((nyhos+1)/2)/3*2
    ax=pex_hos*modex
    ay=pey_hos*modey
    wvmax=min(ax,ay)
    wvf=ratio*wvmax
    delta=2.*(wvmax-wvf)

    do m=1,nyhos/ncpu_hos
       do l=1,nxhos
          if(wvn(l,m,1).ge.1.e-6) then
             dwvn=wvmax-wvn(l,m,1)
             if(dwvn.le.0.) dwvn=0.
             if(dwvn.gt.delta) then
                factor=1.
             else
                factor=462.*(dwvn/delta)**6-1980.*(dwvn/delta)**7+3465.*(dwvn/delta)**8&
                     -3080.*(dwvn/delta)**9+1386.*(dwvn/delta)**10-252.*(dwvn/delta)**11
             endif
             eta(l,m)=eta(l,m)*factor
             vps(l,m)=vps(l,m)*factor
          end if
       end do
    end do

    call fft_bac_xy_hos(eta)
    call fft_bac_xy_hos(vps)

  end subroutine filter_lp2

  !---------------------------------------------------------------

  subroutine filter_lp3(eta,vps,beta1,beta2)
     
    implicit none

    real(wp), intent(in) :: beta1,beta2
    real(wp), intent(inout), dimension(:,:) :: eta,vps

    real(wp) :: maxa, tmpa, con, tmpexp, tmpsum
    integer l,m,kxp,kyp,kx,ky
    real(wp), allocatable, dimension(:,:) :: etaall
    real(wp) :: kp, filter
    
    allocate(etaall(nxhos,nyhos))
    
    call fft_for_xy_hos(eta)
    call fft_for_xy_hos(vps)
    
    call alltoone(eta,etaall)

    kp = 0
    tmpsum = 0
    if (myid == 0) then
       maxa = -1.0
       do l = 1, nxhos - 1, 2
          kx = (l - 1) / 2
          do m = 1, nyhos - 1, 2
             ky = (m - 1) / 2
             tmpa = etaall(l,m)**2 + etaall(l+1,m)**2 + etaall(l,m+1)**2 +  etaall(l+1,m+1)**2
             con = 4.0
             if(l.eq.0.or.m.eq.0) con=2.
             if(l.eq.0.and.m.eq.0) con=1.
             tmpa = tmpa * con

             kp = kp + tmpa**4 * sqrt((kxp*pex_hos)**2 + (kyp*pey_hos)**2)
             tmpsum = tmpsum + tmpa**4
          enddo
       enddo
       kp = kp / tmpsum
    endif

    call mpi_bcast(kp,1,mpi_double_precision,0,mpi_comm_world,ierr)

    do l = 1, nxhos
       kx = (l - 1) / 2
       do m = 1, nyhos / ncpu_hos
          ky = (m + myid * nyhos / ncpu_hos - 1) / 2
          tmpexp = abs(sqrt((kx*pex_hos)**2+(ky*pey_hos)**2) / beta1 / kp)**beta2
          filter = exp(-tmpexp)
          eta(l,m) = eta(l,m) * filter
          vps(l,m) = vps(l,m) * filter
          if (kx .eq. 0 .and. ky .eq. 0) then
             eta(l,m) = 0
             vps(l,m) = 0
          endif
       enddo
    enddo
    deallocate(etaall)
    call fft_bac_xy_hos(eta)
    call fft_bac_xy_hos(vps)
    
  end subroutine filter_lp3

  subroutine filter_lp3_aux(eta,vps,beta1,beta2,sds)

    implicit none

    real(wp), intent(in) :: beta1,beta2
    real(wp), intent(inout), dimension(:,:) :: eta,vps
    real(wp), intent(out), dimension(:,:) :: sds

    real(wp) :: maxa, tmpa, con, tmpexp, tmpsum
    integer l,m,kxp,kyp,kx,ky
    real(wp), allocatable, dimension(:,:) :: etaall
    real(wp) :: kp, filter

    sds = 0.0

    allocate(etaall(nxhos,nyhos))

    call fft_for_xy_hos(eta)
    call fft_for_xy_hos(vps)

    call alltoone(eta,etaall)

    kp = 0
    tmpsum = 0
    if (myid == 0) then
       maxa = -1.0
       do l = 1, nxhos - 1, 2
          kx = (l - 1) / 2
          do m = 1, nyhos - 1, 2
             ky = (m - 1) / 2
             tmpa = etaall(l,m)**2 + etaall(l+1,m)**2 + etaall(l,m+1)**2 +  etaall(l+1,m+1)**2
             con = 4.0
             if(l.eq.0.or.m.eq.0) con=2.
             if(l.eq.0.and.m.eq.0) con=1.
             tmpa = tmpa * con

             kp = kp + tmpa**4 * sqrt((kxp*pex_hos)**2 + (kyp*pey_hos)**2)
             tmpsum = tmpsum + tmpa**4
          enddo
       enddo
       kp = kp / tmpsum
    endif

    call mpi_bcast(kp,1,mpi_double_precision,0,mpi_comm_world,ierr)

    do l = 1, nxhos
       kx = (l - 1) / 2
       do m = 1, nyhos / ncpu_hos
          ky = (m + myid * nyhos / ncpu_hos - 1) / 2
          tmpexp = abs(sqrt((kx*pex_hos)**2+(ky*pey_hos)**2) / beta1 / kp)**beta2
          filter = exp(-tmpexp)
          if (tmpexp >= 20.0) filter = 0
          eta(l,m) = eta(l,m) * filter
          vps(l,m) = vps(l,m) * filter
          if (kx .eq. 0 .and. ky .eq. 0) then
             eta(l,m) = 0
             vps(l,m) = 0
          endif
       enddo
    enddo

    if (myid == 0) then
       do l = 3, nxhos - 1, 2
          kx = (l - 1) / 2
          do m = 3, nyhos - 1, 2
             ky = (m - 1) / 2
             tmpa = etaall(l,m)**2 + etaall(l+1,m)**2 + etaall(l,m+1)**2 +  etaall(l+1,m+1)**2
             con = 4.0
             if(l.eq.0.or.m.eq.0) con=2.
             if(l.eq.0.and.m.eq.0) con=1.
             tmpa = tmpa * con
             tmpexp = abs(sqrt((kx*pex_hos)**2+(ky*pey_hos)**2) / beta1 / kp)**beta2
             filter = exp(-tmpexp)
             if (tmpexp >= 20.0) filter = 0 
             sds(kx,ky) = 1 - filter 
          enddo
       enddo
    end if

    deallocate(etaall)
    call fft_bac_xy_hos(eta)
    call fft_bac_xy_hos(vps)

  end subroutine filter_lp3_aux


  subroutine filter_lp4(eta,vps,ex,ey)

    implicit none

    real(wp), intent(in), dimension(:,:) :: ex,ey
    real(wp), intent(inout), dimension(:,:) :: eta,vps

    real(wp), allocatable, dimension(:,:) :: eo, vo, exo, eyo
    real(wp) :: elimit, exmax, eymax
    integer i,j,nsmooth
    logical ibreak, ismooth
    

    allocate(eo(nxhos,nyhos),vo(nxhos,nyhos))
    allocate(exo(nxhos,nyhos),eyo(nxhos,nyhos))

    elimit = 0.7

    call alltoone(eta,eo)
    call alltoone(vps,vo)
    call alltoone(ex,exo)
    call alltoone(ey,eyo)
    
    if (myid == 0) then
       nsmooth = 0
       do 
!          exo = abs(exo)
!          eyo = abs(eyo)
          exmax = maxval(abs(exo))
          eymax = maxval(abs(eyo))
          nsmooth = nsmooth+1

          ismooth = .false.
          do j = 1, nyhos
             do i = 1, nxhos
                if (exo(i,j) > elimit .or. eyo(i,j) > elimit) then
                   ismooth = .true.
                   call smooth2_sub(nxhos,nyhos,eo,i,j)
                   call smooth2_sub(nxhos,nyhos,vo,i,j)
                end if
             end do
          end do
          
          if (ismooth) then
             if (nsmooth .gt. 50) then
                ibreak = .true.
                exit
             end if
          else
             exit
          end if
          
          call pdfx_hos(eo,exo,pex_hos,nxhos,nyhos)
          call pdfy_hos(eo,eyo,pey_hos,nxhos,nyhos)
       end do
    end if

    call onetoall(eo,eta)
    call onetoall(vo,vps)

    deallocate(eo, vo, exo, eyo)

  end subroutine filter_lp4
  
  subroutine smooth2_sub(nxmod,nymod,f,i,j)


    !     by di yang, 08/2007
    
    !     9-points local average

    implicit none
    
    integer i,j,ks,i1,j1,ip,jp,im,jm,idex,jdex
    integer nxmod,nymod
    
    real(wp) f(nxmod,nymod), tmp(nxmod,nymod)

    ks = 32
    do j1=j-ks/2,j+ks/2
       jp=j1
       if(jp.gt.nymod) jp=jp-nymod
       if(jp.lt.1) jp=jp+nymod
       do i1=i-ks/2,i+ks/2
          ip=i1
          if(ip.gt.nxmod) ip=ip-nxmod
          if(ip.lt.1) ip=ip+nxmod
          tmp(ip,jp)=0.
          do jm=-1,1
             do im=-1,1
                idex=ip+im
                jdex=jp+jm
                if(idex.gt.nxmod) idex=idex-nxmod
                if(idex.lt.1) idex=idex+nxmod
                if(jdex.gt.nymod) jdex=jdex-nymod
                if(jdex.lt.1) jdex=jdex+nymod
                tmp(ip,jp)=tmp(ip,jp)+f(idex,jdex)
             enddo
          enddo
          f(ip,jp)=tmp(ip,jp)/9.
       enddo
    enddo

    return
  end subroutine smooth2_sub



end module smooth
