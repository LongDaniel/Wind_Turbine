module utils_hos

  use hos_param, only : wp, twopi, nxhos, nyhos, ncpu_hos, fr2, pex_hos, pey_hos
  use fft_hos
  use MPI

  implicit none

  private

  public :: skwave, kk2comp, fullspec_xy, get_moment, get_stats
  public :: get_edf, get_pdf, gram_charlier, get_peak, get_width
  public :: hilbert, get_wave_coherent
  public :: fit_jonswap

  interface get_peak
     module procedure get_peak_wave_energy, get_peak_wave_action
  end interface
    
contains

  subroutine gram_charlier(x, p, skew, kurt)

    implicit none

    real(wp), intent(in) :: x, skew, kurt
    real(wp), intent(out) :: p
    
    real(wp) h3, h4, k3, k4

    h3 = x**3-3*x
    h4 = x**4 - 6*x**2+3

    k3 = skew
    k4 = kurt - 3
    
    p = (1+k3*h3/6+k4*h4/24)*exp(-0.5*x**2) / sqrt(twopi)

  end subroutine gram_charlier

  !-------------------------------------------------------------------------

  subroutine get_edf(edf, etaall, rms)

    implicit none
    
    real(wp), intent(inout), dimension(:,:) :: edf
    real(wp), intent(in) :: rms
    real(wp), intent(in), dimension(:,:) :: etaall
    
    real(wp), allocatable,  dimension(:,:) :: tmp
    real(wp) x, amp
    integer i,j, it

    allocate(tmp(nxhos, nyhos))

    tmp = etaall / rms
    edf(:,2) = 0
    do j = 1, nyhos
       do i = 1, nxhos
          do it = 1, nxhos
             x = edf(it,1)
             if (tmp(i,j) > x) edf(it,2) = edf(it,2) + 1
          end do
       end do
    end do
    edf(:,2) = edf(:,2) / nxhos / nyhos 

    deallocate(tmp)
  end subroutine get_edf

  !-------------------------------------------------------------------------

  subroutine get_pdf(pdf, etaall, rms)

    implicit none

    real(wp), intent(inout), dimension(:,:) :: pdf
    real(wp), intent(in) :: rms
    real(wp), intent(in), dimension(:,:) :: etaall

    real(wp), allocatable,  dimension(:,:) :: tmp
    real(wp) x, amp, dx
    integer i,j, it

    allocate(tmp(nxhos, nyhos))

    tmp = etaall / rms
    dx = (pdf(nxhos,1) - pdf(1,1)) / (nxhos - 1.0)
    pdf(:,2) = 0
    do j = 1, nyhos
       do i = 1, nxhos
          do it = 1, nxhos
             x = pdf(it,1)
             if (tmp(i,j) > x - 0.5 * dx .and. tmp(i,j) < x + 0.5 * dx) then
                pdf(it,2) = pdf(it,2) + 1
             end if
          end do
       end do
    end do
    pdf(:,2) = pdf(:,2) / nxhos / nyhos / dx

    deallocate(tmp)
  end subroutine get_pdf

  !-------------------------------------------------------------------------

  subroutine get_width(sf, nmax, dwf, width)

    implicit none

    real(wp), intent(out) :: width
    integer, intent(in) :: nmax
    real(wp), intent(in) :: dwf
    real(wp), intent(in), dimension(nmax) :: sf

    integer i, imax
    real(wp) sf0, f1, f2

    imax = maxval(maxloc(sf))
    sf0 = sf(imax)

    i = 1
    do while (sf(i) <= 0.5*sf0)
       i = i + 1
    end do
    f1 = i * dwf

    i = imax
    do while (sf(i) >= 0.5 * sf0)
       i = i + 1
    end do
    f2 = i * dwf

    width = (f2 - f1) / (imax * dwf)
    
  end subroutine get_width

  !-------------------------------------------------------------------------

  subroutine skwave(eta, rms, hs, skew, kurt)

    implicit none

    real(wp), intent(out) :: rms, hs, skew, kurt
    real(wp), intent(in), dimension(:,:) :: eta
    
    real(wp) :: sum1, sum2, sum3
    real(wp) :: tmp1, tmp2, tmp3

    integer i,j

    tmp1 = 0.0_wp
    tmp2 = 0.0_wp
    tmp3 = 0.0_wp
!    tmp1 = sum(eta*eta)
!    tmp2 = sum(eta*eta*eta)
!    tmp3 = sum(eta*eta*eta*eta)
    do i = 1, nxhos
       do j = 1, nyhos / ncpu_hos
          tmp1 = tmp1 + eta(i,j) ** 2
          tmp2 = tmp2 + eta(i,j) ** 3
          tmp3 = tmp3 + eta(i,j) ** 4
       end do
    end do
    
    call mpi_allreduce(tmp1,sum1,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(tmp2,sum2,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(tmp3,sum3,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    rms = sqrt(sum1 / nxhos / nyhos)
    hs = 4 * rms
    skew = (sum2 / nxhos / nyhos) / rms**3
    kurt = (sum3 / nxhos / nyhos) / rms**4
    
  end subroutine skwave

  !-------------------------------------------------------------------------
 
  subroutine kk2comp(a,acomp)

    implicit none
    
    real(wp), intent(in), dimension(:,:) :: a
    real(wp), intent(out), dimension(:,:) :: acomp

    integer kx, ky, l, m

    do kx = 0, nxhos / 2
       l = kx*2+1
       do ky = 0, nyhos / 2
          m = ky*2+1
          acomp(l,m) = a(l,m) + a(l+1,m+1)
          acomp(l,m+1) = -a(l,m+1) - a(l+1,m)
          acomp(l+1,m) = sqrt(acomp(l,m)**2+acomp(l,m+1)**2)
          acomp(l+1,m+1) = atan2(acomp(l,m+1),acomp(l,m))
       end do
    end do
    
  end subroutine kk2comp

  !-------------------------------------------------------------------------

  subroutine spec2phy(a,b,c,d,a1,a2,alpha,beta)
    
    !     by xuanting hao,
    !     04/23/2014   first edition
    
    !     cc/4, -sc/4      a   b
    !     -cs/4,ss/4       c   d
    
    !      4a  =  cc = a1*cos(alpha) + a2*cos(beta)
    !     -4b  =  sc = -a1*sin(alpha) - a2*sin(beta)
    !     -4c  =  cs = -a1*sin(alpha) + a2*sin(beta)
    !      4d  =  ss = -a1*cos(alpha) + a2*cos(beta)
    !     a1^2 + a2^2 = 8 * (a^2+b^2+c^2+d^2)
    
    implicit none
    
    real(wp), intent(in) :: a,b,c,d
    real(wp), intent(out) :: a1,a2,alpha,beta
    
    a1 = 2*sqrt((b + c)**2 + (a - d)**2)
    a2 = 2*sqrt((b - c)**2 + (a + d)**2)
    
    alpha = sign(acos((a - d) / a1),b + c)
    beta = sign(acos((a + d) / a2),b - c)   
    
  end subroutine spec2phy
  
  !-------------------------------------------------------------------------

  subroutine get_moment(moment,skxy,n_mom)

    implicit none

    real(wp), intent(out) :: moment
    integer, intent(in) :: n_mom
    real(wp), intent(in), dimension(nxhos/2,-nyhos/2:nyhos/2) :: skxy
    
    integer i,j
    real(wp) wkx,wky,omega

    moment = 0.0_wp
    
    do j = -nyhos/2, nyhos/2
       wky = j * pey_hos
       do i = 1, nxhos/2
          wkx = i * pex_hos
          omega = sqrt(sqrt(wkx**2+wky**2) / fr2)
          moment = moment + omega**n_mom * skxy(i,j) * pex_hos * pey_hos
       end do
    end do
  end subroutine get_moment

  !-------------------------------------------------------------------------
  
  subroutine fullspec_xy(f,skxy)
    
    implicit none
    
    real(wp), intent(in), dimension(:,:) :: f
    real(wp), intent(out), dimension(nxhos/2,-nyhos/2:nyhos/2) :: skxy
    
    real(wp), allocatable, dimension(:,:) :: tmp
    real(wp) a1, a2, alpha, beta, wkx,wky
    integer i,j,l,m,kx,ky, root
    
    allocate(tmp(nxhos,nyhos))
    root = 0

    tmp = f
    !     call fft_for_xy_hos(tmp,nxhos,nyhos)
    call fft_for_xy_hos(f,tmp,nxhos,nyhos)
    
    skxy = 0.0_wp
    
    do l = 3, nxhos - 1, 2
       kx = (l - 1) / 2
       wkx = kx * pex_hos
       do m = 1, nyhos - 1, 2
          ky = (m - 1) / 2
          wky = ky * pey_hos
          if (abs(kx)+abs(ky) /= 0) then
             if (kx /= 0) then
                call spec2phy(tmp(l,m),tmp(l+1,m),tmp(l,m+1),tmp(l+1,m+1),a1,a2,alpha,beta)
                skxy(kx,ky) = a1**2 / 2 / pex_hos / pey_hos
                skxy(kx,-ky) = a2**2 / 2 / pex_hos / pey_hos
                if (ky == 0) then
                   skxy(kx,ky) = (a1**2 + a2**2) / 8 / pex_hos / pey_hos
                endif
             endif
          endif
       enddo
    enddo
    
    j = nxhos/2 * (nyhos+1)
    call mpi_bcast(skxy, j , mpi_double_precision, root, mpi_comm_world,ierr)

    deallocate(tmp)
    
  end subroutine fullspec_xy

  !-------------------------------------------------------------------------

  subroutine get_peak_wave_energy(skxy, kpx, kpy, kp, tp, omgp, cp, charv)

    implicit none

    real(wp), intent(in), dimension(nxhos/2,-nyhos/2:nyhos/2) :: skxy
    real(wp), intent(out) :: kpx, kpy, kp, tp, omgp, cp, charv

    integer i,j    
    real(wp) wkx, wky, tmp, tscale, cg, rek

    kpx = 0.0
    kpy = 0.0
    charv = 0
    tscale = sum(abs(skxy)) / nxhos / nyhos / 2

    do i = 1, nxhos / 2
       wkx = pex_hos * i
       do j = -nyhos / 2, nyhos / 2
          wky = pey_hos * j
          rek = sqrt(wkx**2+wky**2)
          cg = 0.5 * sqrt(1.0 / fr2 / rek)
          charv = charv + skxy(i,j) * cg
          kpx = kpx + (skxy(i,j)/tscale)**4 * wkx
          kpy = kpy + (skxy(i,j)/tscale)**4 * wky
          tmp = tmp + (skxy(i,j)/tscale)**4
       end do
    end do

    charv = charv / sum(skxy)
    kpx = kpx / tmp
    kpy = kpy / tmp
    kp = sqrt(kpx**2 + kpy**2)
    omgp = sqrt(kp / fr2)
    cp = omgp / kp
    tp = twopi / omgp

  end subroutine get_peak_wave_energy

  !-------------------------------------------------------------------------

  subroutine get_peak_wave_action(skxy, kpx, kpy, kp, tp, omgp, cp)

    implicit none

    real(wp), intent(in), dimension(nxhos/2,-nyhos/2:nyhos/2) :: skxy
    real(wp), intent(out) :: kpx, kpy, kp, tp, omgp, cp

    integer i,j
    real(wp) wkx, wky, tmp, tscale, cg, rek, omega
    real(wp), allocatable, dimension(:,:) :: nkxy

    kpx = 0.0
    kpy = 0.0

    allocate(nkxy(nxhos/2,-nyhos/2:nyhos/2))
    do i = 1, nxhos / 2
       wkx = pex_hos * i
       do j = -nyhos / 2, nyhos / 2
          wky = pey_hos * j
          rek = sqrt(wkx**2+wky**2)
          omega = sqrt(rek / fr2)
          nkxy(i,j) = skxy(i,j) / omega
       end do
    end do

    tscale = sum(abs(nkxy)) / nxhos / nyhos / 2

    do i = 1, nxhos / 2
       wkx = pex_hos * i
       do j = -nyhos / 2, nyhos / 2
          wky = pey_hos * j
          rek = sqrt(wkx**2+wky**2)
          omega = sqrt(rek / fr2)
          cg = 0.5 * sqrt(1.0 / fr2 / rek)
          kpx = kpx + (nkxy(i,j)/tscale)**4 * wkx
          kpy = kpy + (nkxy(i,j)/tscale)**4 * wky
          tmp = tmp + (nkxy(i,j)/tscale)**4
       end do
    end do
    
    kpx = kpx / tmp
    kpy = kpy / tmp
    kp = sqrt(kpx**2 + kpy**2)
    omgp = sqrt(kp / fr2)
    cp = omgp / kp
    tp = twopi / omgp

    deallocate(nkxy)

  end subroutine get_peak_wave_action

  !-------------------------------------------------------------------------
  
  subroutine get_stats(vol,flux,pe,ke,ener,eta,vps,feta)
    
    implicit none
    
    real(wp), intent(in), dimension(:,:) :: eta, vps
    real(wp), intent(in), dimension(:,:,:) :: feta
    real(wp), intent(out) :: vol, flux, pe, ke, ener
    
    integer i,j,root
    real(wp), allocatable, dimension(:,:) :: t1, t2
    
    allocate(t1(nxhos,nyhos/ncpu_hos))
    allocate(t2(nxhos,nyhos/ncpu_hos))
    
    root = 0
    
    ! volume
    t1 = eta
    call fft_for_xy_hos(t1)
    vol = (twopi/pex_hos) * (twopi/pey_hos) * t1(1,1)
    
    ! volume flux 
    t2 = feta(:,:,1)
    call fft_for_xy_hos(t2)
    flux = (twopi/pex_hos) * (twopi/pey_hos) * t2(1,1)
    
    ! potential and kinetic energy
    t1 = eta * eta
    t2 = vps * feta(:,:,1)
    call fft_for_xy_hos(t1)
    call fft_for_xy_hos(t2)
    pe = 0.5*(twopi/pex_hos) * (twopi/pey_hos) * t1(1,1) / fr2
    ke = 0.5*(twopi/pex_hos) * (twopi/pey_hos) * t2(1,1)
    ener = pe + ke
    
    call mpi_bcast(vol, 1, mpi_double_precision, root, mpi_comm_world,ierr)
    call mpi_bcast(flux, 1, mpi_double_precision, root, mpi_comm_world,ierr)
    call mpi_bcast(pe, 1, mpi_double_precision, root, mpi_comm_world,ierr)
    call mpi_bcast(ke, 1, mpi_double_precision, root, mpi_comm_world,ierr)
    call mpi_bcast(ener, 1, mpi_double_precision, root, mpi_comm_world,ierr)
    
    deallocate(t1,t2)
    
  end subroutine get_stats

  !-------------------------------------------------------------------------

  subroutine fit_jonswap(u10_fit,fetch_fit,sf,nfhos,dfq,clen,ctime)

    implicit none

    integer, intent(in) :: nfhos
    real(wp), intent(in), dimension(:) :: sf
    real(wp), intent(in) :: dfq, clen, ctime
    real(wp), intent(out) :: u10_fit,fetch_fit

    real(wp), allocatable, dimension(:) :: sj

    integer i, j, nu,nf
    real(wp) umin,umax,fmin,fmax
    parameter (nu=4096, nf=4096)
    real(wp) jcost,jmin, myu, myf

    allocate(sj(nfhos))

    umin = 3.0
    umax = 15.0
    fmin = 3000
    fmax = 100000

    jmin = 5e8

!    print *, myid, nfhos

    do i = 1, nu
       myu = umin + (i-1) * (umax - umin) / (nu - 1.0)
       do j = 1, nf
          myf = fmin + (j - 1) * (fmax - fmin) / (nf - 1.0)
          call jonswap_generator(myu, myf, dfq/ctime, nfhos, sj)
          jcost = sum(abs(sf(1:nfhos)*clen**2*ctime-sj(1:nfhos)))
          if (jcost <= jmin) then
             jmin = jcost
             u10_fit = myu
             fetch_fit = myf
          end if         
       end do
    end do

    call jonswap_generator(u10_fit, fetch_fit, dfq/ctime, nfhos, sj)

    if (myid == 0) then
!       do i = 1, nfhos
!          write(295,'(25e12.4)') i*dfq/ctime, sf(i)*clen**2*ctime, sj(i)
!       end do
    end if
 
    deallocate(sj)
    
  end subroutine fit_jonswap

  !-------------------------------------------------------------------------
  
  subroutine jonswap_generator(u10, fetch, dfq, nfhos, sj)

    implicit none

    integer, intent(in) :: nfhos
    real(wp), intent(in) :: u10, fetch, dfq
    real(wp), intent(out), dimension(:) :: sj

    integer i
    real(wp) alpha, omgp, gamma, sigma, omega

    alpha = 0.076 * (u10**2/fetch/9.8)**0.22
    omgp = 22 * (9.8**2/u10/fetch) **(1.0/3)
    gamma = 3.3
    do i = 1, nfhos
       omega = twopi * i * dfq
       if ( omega .lt. omgp) then
          sigma = 0.07
       else
          sigma = 0.09
       end if
       sj(i) = twopi * alpha * 9.8**2 / omega**5 * exp(-1.25*(omega/omgp)**(-4.))&
            * 3.3**(exp(-(omega-omgp)**2/2/sigma**2/omgp**2))
    end do

  end subroutine jonswap_generator

  !-------------------------------------------------------------------------

  subroutine get_wave_coherent(f,f_wave,eta,n,dtheta)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in), dimension(:) :: f, eta 
    real(wp), intent(out), dimension(:) :: f_wave
    real(wp), intent(out) :: dtheta

    integer i, j
    real(wp), allocatable, dimension(:) :: heta, deta, dheta, test
    real(wp) sum1, sum2, amp, tmp, alp, f_ave

    allocate(heta(n))
    allocate(deta(n))
    allocate(dheta(n))
    allocate(test(n))

    alp = 0.00001
    call hilbert(eta,heta,n)

    do i = 1, n - 1
       deta(i) = (eta(i+1)-eta(i)) 
       dheta(i) = (heta(i+1)-heta(i)) 
    end do
    deta(n) = deta(n-1)
    dheta(n) = dheta(n-1)

    f_ave = sum(f(1:n))/n
    do j = 1, n
       sum1 = 0
       sum2 = 0
       do i = 1, n
          amp = sqrt(deta(i)**2 + dheta(i)**2)
          tmp = sqrt((eta(i)-eta(j))**2 + (heta(i)-heta(j))**2)
          sum1 = sum1 + f(i) * amp * delta(tmp,alp)
          sum2 = sum2 + amp * delta(tmp,alp) 
       end do

       f_wave(j) = sum1 / sum2 - f_ave
    end do

    call fft_for_hos(f_wave,test,n)
    test(n/4+1:n) = 0
    call fft_bac_hos(test,f_wave,n)

    sum1 = 0
    sum2 = 0
    do i = 1, n
       sum1 = sum1 + f_wave(i) * heta(i)
       sum2 = sum2 + f_wave(i) * eta(i)
    end do

    dtheta = atan2(sum1,sum2)
    
    deallocate(heta,deta,dheta,test)

  end subroutine get_wave_coherent

  !-------------------------------------------------------------------------

  subroutine hilbert(f,hf,n)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in), dimension(:) :: f
    real(wp), intent(out), dimension(:) :: hf

    integer l
    real(wp), allocatable, dimension(:) :: tmp

    allocate(tmp(n))

    call fft_for_hos(f,tmp,n)
    
    hf(1) = 0
    hf(2) = 0
    do l = 3, n - 1, 2
       hf(l) = tmp(l+1)
       hf(l+1) = -tmp(l)
    end do
    
    tmp = hf(1:n)

    call fft_bac_hos(tmp,hf,n)

    deallocate(tmp)

  end subroutine hilbert

  !-------------------------------------------------------------------------

  subroutine get_inst(f,hf,iamp,iphs,n)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in), dimension(:) :: f, hf
    real(wp), intent(out), dimension(:) :: iamp, iphs

    integer i

    do i = 1, n
       iamp(i) = sqrt(f(i)**2 + hf(i)**2)
       iphs(i) = atan(hf(i), f(i))
    end do    

  end subroutine get_inst

  !-------------------------------------------------------------------------

  real(wp) function delta(x,alpha)
    
    implicit none

    real(wp), intent(in) :: x, alpha

    real(wp) tmp
    if (alpha <= 0) then
       print *,"Wrong alpha value!"
       stop
    end if

    tmp = -x**2/alpha**2
    
    if (tmp <= -600) then 
       delta = 0
    else 
       delta = (1/alpha/sqrt(0.5*twopi)) * exp(tmp)
    end if

  end function delta

end module utils_hos
 
