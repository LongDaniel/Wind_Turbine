module src_stat

  use hos_param, only : wp, pi, twopi, nxhos, nyhos, ncpu_hos, fr2, pex_hos, pey_hos
  use fft_hos
  use spectral_hos
  use MPI

  implicit none

  private

  real(wp), public :: fql, fqh, thl, thh, dfq, dth
  integer, public :: maxnfq, minnfq, maxnth, minnth, nfhos, nthos
  public :: src_init,fkk2omg, snl_lavr, snl_webb, twod2oned, check_snl

contains

  subroutine src_init

    implicit none
    real(wp) tmpk, tmpomg

!    nfhos = 128
!    nthos = 128
    read(13,*) nfhos, nthos

    if (mod(nthos,ncpu_hos) /= 0) then
       nthos = (nthos/ncpu_hos) * ncpu_hos
       print *, "Nthos changed to:", nthos
    end if

    tmpk = sqrt((pex_hos*nxhos/2)**2 + (pey_hos*nyhos/2)**2)
    tmpomg = sqrt(tmpk/fr2)
    dfq = tmpomg / twopi / nfhos
    dth = twopi / nthos
   
    fql = 3.5*dfq
    fqh = nfhos*dfq / 1.1
    thl = -pi
    thh = pi

    maxnfq = floor(fqh / dfq) + 1
    minnfq = ceiling(fql / dfq) + 1
!    maxnfq = nfhos
!    minnfq = 1
    maxnth = nthos
    minnth = 1


    if (myid .eq. 0) then
       print *,"dfq=", dfq
       print *,"dth=", dth
       
       print *,"fql =",fql
       print *,"fqh =",fqh
       print *,"thl =",thl
       print *,"thh =",thh
       print *,"maxnfq = ",maxnfq
       print *,"minnfq = ",minnfq
       print *,"maxnth = ",maxnth
       print *,"minnth = ",minnth
    endif

  end subroutine src_init

  subroutine fkk2omg(skxy, fk)

    implicit none

    real(wp), intent(in), dimension(nxhos/2,-nyhos/2:nyhos/2) :: skxy
    real(wp), intent(out), dimension(nfhos,nthos) :: fk

    real(wp) freq, theta, kx, ky, tmpk(2), omg, sumene
    integer i,j, nj, ni, tmpi, tmpj, ui, li, uj, lj 

    fk = 0
!    do nj = 1, nthos
!       theta = (nj - 1) * dth
!       do ni = 1, nfhos
!          freq = (ni - 1) * dfq         
!          sumene = 0
!          do j = -nyhos / 2, nyhos / 2
!             ky = j * pey_hos
!             do i = 1, nxhos / 2
!                kx = i * pex_hos
!                tmpk(1) = kx
!                tmpk(2) = ky
!                call kk2omgthe(tmpk,omg,theta,fr2)
!                freq = omg / twopi
                
!                if (freq < (ni+0.5)*dfq .and. freq >= (ni-0.5)*dfq &
!                   .and. theta < (nj+0.5)*dth - pi .and. theta >= (nj-0.5)*dth - pi) then
!                   sumene = sumene + skxy(i,j) * pex_hos * pey_hos 
!                end if
!             end do
!          end do
!          fk(ni,nj) = sumene / dfq / dth
!       end do
!    end do


    do j = -nyhos / 3, nyhos / 3
       ky = j * pey_hos
       do i = 1, nxhos / 3
          kx = i * pex_hos
          tmpk(1) = kx
          tmpk(2) = ky
          call kk2omgthe(tmpk,omg,theta,fr2)
          freq = omg / twopi
          tmpi = int(freq / dfq)
          tmpj = int((theta+pi)/dth)
          if (myid == 0) then
             write (385,*) freq, theta
          end if
          ui = tmpi + 2
          li = tmpi - 2
          uj = tmpj + 2
          lj = tmpj - 2

          if (ui > nfhos) ui = nfhos
          if (li < 1) li = 1
          if (uj > nthos) uj = nthos
          if (lj < 1) lj = 1

          do nj = lj, uj
             do ni = li, ui
                if (freq < (ni+0.5)*dfq .and. freq >= (ni-0.5)*dfq &
                     .and. theta < (nj+0.5)*dth - pi .and. theta >= (nj-0.5)*dth - pi) then
                   fk(ni,nj) = fk(ni,nj) + skxy(i,j) * pex_hos * pey_hos / dfq / dth
                end if
             end do
          end do
       end do
    end do
    
  end subroutine fkk2omg

  subroutine omgthe2kk(k,omg,theta,fr2)
    !     by xuanting hao
    
    !     this subroutine converts omg-theta coordinate to kx-ky coordinate
    
    !     03/16/2014  first edition
    !     04/03/2014  use module
    !     07/01/2014  delete module
    
    !      use m_const
    implicit none
    real(wp),intent(in)::fr2,omg,theta
    real(wp) k(2)
    
    k(1) = fr2 * omg**2 * cos(theta)
    k(2) = fr2 * omg**2 * sin(theta)
    
  end subroutine omgthe2kk
  
!--------------------------------------------

  subroutine kk2omgthe(k,omg,theta,fr2)
    !     by xuanting hao
    
    !     this subroutine converts kx-ky coordinate to omg-theta coordinate
    !     range of theta:[-pi,pi)
    
    !     03/16/2014  first edition
    !     03/20/2014  use another method for theta
    !     04/02/2014  use module
    !     07/01/2014  delete module
    
    !      use m_const
    implicit none
    real(wp),intent(in)::k(2),fr2
    real(wp) omg,theta,pi
    
    pi = acos(-1.0)
    omg = sqrt(sqrt(k(1)**2 + k(2)**2) / fr2)
    
    theta = sign(acos(k(1) / sqrt(k(1)**2 + k(2)**2)),k(2))
    
    do while (theta .lt. -pi)
       theta = theta + 2*pi
    enddo
    
    do while (theta .ge. pi)
       theta = theta - 2*pi
    enddo
    
  end subroutine kk2omgthe



!-------------------------------------------------------------------------                                   
  subroutine regulate(theta)
    !     by xuanting hao     
    !     this subroutine regulates theta to the domain [-pi,pi)         
    !     05/28/2014  first edition                                                                              
    implicit none
    
    real(wp), intent(inout) ::  theta
    real(wp) eps
    
    eps = 1.0e-15
    
    do while (theta .lt. -pi)
       theta = theta + 2*pi
    enddo
    
    do while (theta .ge. pi)
       theta = theta - 2*pi
    enddo
    
    if (abs(theta-pi) .lt. eps) then
       theta = -pi
    endif
    
  end subroutine regulate
  
  subroutine checkresonant(k,k1,k2,k3,omg,omg1,omg2,omg3,ite)
    !     by xuanting hao                                              
    !     this subroutine checks the resonant condition for the              
    !     nonlinear wave interactions                              
    !     06/04/2014   first edition               
    !     06/17/2014   second edition                                
    implicit none
    
    real(wp),dimension(2),intent(in)::k,k1,k2,k3
    real(wp),intent(in)::omg,omg1,omg2,omg3
    integer,intent(in)::ite
    
    real(wp) err1,err2,tmp(2),eps, test(2)
    
    eps = 1.0e-5
    tmp = k + k1 - k2 - k3
    test = abs(k) + abs(k1) + abs(k2) + abs(k3)
    err1 = sqrt(tmp(1)**2+tmp(2)**2)
    
    err2 = omg + omg1 - omg2 - omg3
    
    if (err1 > eps * sqrt(test(1)**2+test(2)**2)) then
       print *,"ite = ",ite
       print *,"err k resonant! err1=",err1
    endif
    if (err2 > eps * max(omg,omg1,omg2,omg3)) then
       print *,"ite = ",ite
       print *,"err omega resonant! err2=",err2
    endif
    
  end subroutine checkresonant

!------------------------------------------------------------------------                                    
  subroutine check_snl(fk,snl,act_tot,mom_x_tot,mom_y_tot,ene_tot,ent_tot)

    implicit none
    
    real(wp), intent(out) :: act_tot,mom_x_tot,mom_y_tot,ene_tot,ent_tot
    real(wp), intent(in), dimension(:,:) :: fk, snl
    real(wp), allocatable, dimension(:,:) :: fkall, snlall
    
    integer i,j
    real(wp) freq, theta, wvn, temp, eps

    allocate(fkall(nfhos,nthos))
    allocate(snlall(nfhos,nthos))
    call alltoone(snl,snlall,nfhos,nthos)
    
    temp = sum(abs(fkall)) / nfhos / nthos

    ! wave action
    act_tot = 0
    mom_x_tot = 0
    mom_y_tot = 0
    ene_tot = 0
    ent_tot = 0
    eps = 1.0e-6

    do j = 1, nthos
       theta = j * dth - pi
       do i = 1, nfhos
          freq = i * dfq
          wvn = fr2 * (twopi * freq)**2
          act_tot = act_tot + (snlall(i,j) / twopi / freq) * dfq * dth
          ene_tot = ene_tot + snlall(i,j) * dfq * dth
          if (fkall(i,j) >= eps * temp) then
             ent_tot = ent_tot + (snlall(i,j) / fkall(i,j)) * dfq * dth
          end if
          mom_x_tot = mom_x_tot + (snlall(i,j)*wvn*cos(theta)/twopi/freq) * dfq * dth
          mom_y_tot = mom_y_tot + (snlall(i,j)*wvn*sin(theta)/twopi/freq) * dfq * dth
       end do
    end do

    deallocate(fkall,snlall)

  end subroutine check_snl
!------------------------------------------------------------------------ 
  subroutine snl_lavr(fk,snl)

    !     by xuanting hao                                                                     
    
    !     this subroutine calculates snl using lavrenov's method                              
    !     time complexity o(nfmod^2*ntmax/ncpu*n^2)                                           
    !     03/16/2014 the first edition                                                        
    !     03/28/2014 the second edition introduce module                                      
    !     05/20/2014 complete the parallel version                                            
    !     07/01/2014 remove the use of module                                                 

    implicit none
    
    real(wp), dimension(:,:),intent(in)::fk
    real(wp), dimension(:,:)::snl

    real(wp), allocatable, dimension(:,:)::fkall
    integer i,j,n1,n2,isign,rej
    integer m1  !integration over omg1                                                  
    integer m2  !integration over theta1                                                
    integer m3  !integration over omg2                                                  
    real(wp) omg,theta,omg1,theta1,omg2,theta2,omg3,theta3
    real(wp) omga,ka,epsa
    real(wp),dimension(2)::k,k1,k2,k3,kav
    real(wp) maxomg2,minomg2
    !     note: f(omg,the)=f(fq,the) / twopi                                                  
    real(wp) f,f1,f2,f3
    real(wp) hfun,coef,den,beta,a,domg
    real(wp) sum1,sum2,sum3
    
    n1 = 20
    n2 = 20
    domg = twopi * dfq
    snl = 0.0

    allocate(fkall(nfhos,nthos))
    call alltoone(fk,fkall,nfhos,nthos)

    do i = nfhos, 2, -1 
       fkall(i,:) = fkall(i-1,:)
    end do
    fkall(1,:) = 0
    
    do i = minnfq, maxnfq
       do j = 1, nthos / ncpu_hos
          rej = j + myid * nthos / ncpu_hos
          omg = twopi * (i - 1) *dfq
          theta = (rej - 1) * dth - pi
          f = fkall(i,rej) / twopi
          call omgthe2kk(k,omg,theta,fr2)
          
          !     integration over omg1                 
          do m1 = minnfq, maxnfq
             omg1 = twopi * (m1 - 1) * dfq
             call int_the1(sum1,omg,theta,k,omg1,f,fkall)
             snl(i,j) = snl(i,j) + sum1 * domg * twopi
          enddo
       enddo
    enddo
    
    call alltoone(snl,fkall,nfhos,nthos)

    do i = 1, nfhos - 1
       fkall(i,:) = fkall(i+1,:)
    end do
    
    call onetoall(fkall,snl,nfhos,nthos)

    deallocate(fkall)
    
  end subroutine snl_lavr
!-------------------------------------------------------------------------   

  subroutine snl_webb(fk,snl)

    !     by xuanting hao
    
    !     this subroutine calculates snl using webb's method
    
    !     06/16/2014   the first edition
    !     07/01/2014   remove the use of module

    implicit none

    real(wp), dimension(:,:),intent(in)::fk
    real(wp), dimension(:,:)::snl

    real(wp), allocatable, dimension(:,:)::fkall
    integer i,j,isign,rej
    integer m1,m2

    real(wp) omg1,theta1,omg2,theta2,omg3,theta3,omg4,theta4
    real(wp),dimension(2)::k1,k2,k3,k4
    real(wp) f1,f2,f3,f4
    real(wp) n1,n2,n3,n4
    real(wp) t13,maxomg2,minomg2
    !     note: f(omg,the)=f(fq,the) / twopi
    real(wp) hfun,coef,den,beta,a,domg

    domg = twopi * dfq
    snl = 0.0

    allocate(fkall(nfhos,nthos))
    call alltoone(fk,fkall,nfhos,nthos)
    do i = nfhos, 2, -1
       fkall(i,:) = fkall(i-1,:)
    end do
    fkall(1,:) = 0

    do i = minnfq, maxnfq
       do j = 1, nthos / ncpu_hos
          rej = j + myid * nthos / ncpu_hos
print *, 5, 1, i, j
          omg1 = twopi * (i - 1) *dfq
          theta1 = (rej - 1) * dth - pi
          f1 = fkall(i,rej) / twopi
          call omgthe2kk(k1,omg1,theta1,fr2)
print *, 5, 2 
          n1 = f1 / 2 / omg1**4 / fr2**2

          do m1 = minnfq, maxnfq
             do m2 = minnth, maxnth
                omg3 = twopi * (m1 - 1) * dfq
                theta3 = (m2 - 1) * dth - pi
                call omgthe2kk(k3,omg3,theta3,fr2)
             !print *, 5, 3, m1, m2   
                f3 = fkall(m1,m2) / twopi
                n3 = f3 / 2 / omg3**4 / fr2**2
                
                if (abs(i - m1) + abs(rej - m2) .ne. 0) then
                   call calt13(t13,k1,k3,omg1,omg3,theta1,theta3,n1,n3,fkall,i,rej,m1,m2)
                   snl(i,j) = snl(i,j) + (4 * pi * omg1**4 * fr2**2) * t13 &
                        * (2 * omg3**3 * fr2**2) * domg * dth
                endif
             end do
          end do
       end do
    end do

    call alltoone(snl,fkall,nfhos,nthos)

    do i = 1, nfhos - 1
       fkall(i,:) = fkall(i+1,:)
    end do

    call onetoall(fkall,snl,nfhos,nthos)
    deallocate(fkall)

  end subroutine snl_webb
  !-------------------------------------------------------------------------

  subroutine twod2oned(f2d,f1d)
    
    implicit none

    real(wp), intent(in), dimension(:,:) :: f2d
    real(wp), intent(out), dimension(:) :: f1d

    integer i,j
    real(wp) theta
    real(wp), allocatable, dimension(:,:) :: f2dall

    allocate(f2dall(nfhos,nthos))
    call alltoone(f2d,f2dall,nfhos,nthos)
    
    f1d = 0.0
    do j = 1, nthos
       do i = 1, nfhos
          theta = (j - 1) * dth - pi
          if (theta <= thh .and. theta >= thl) then
             f1d(i) = f1d(i) + f2dall(i,j)
          end if
       end do
       f1d(i) = dth * f1d(i)
    end do
    
    deallocate(f2dall)
  end subroutine twod2oned

  !-------------------------------------------------------------------------
  integer function heavi(k1,k2,k3,k4)

    implicit none

    real(wp), dimension(2), intent(in) :: k1,k2,k3,k4
    real(wp) tmp1,tmp2,x
    tmp1 = dot_product(k1-k4,k1-k4)
    tmp2 = dot_product(k1-k3,k1-k3)
    x = tmp1 - tmp2
    if (x > 0) then
       heavi = 1
    else
       heavi = 0
    endif

  end function heavi

!-------------------------------------------------------------------------

  subroutine calt13(t13,k1,k3,omg1,omg3,theta1,theta3,n1,n3,fkall,i1,j1,i2,j2)
    
    !     by xuanting hao
    
    !     implicie none
    !     this subroutine calculates t13(k1,k3)
    
    !     06/23/2014   first edition use a different method to calculate ds and k2
    !     06/26/2014   calculate the special case
    
    implicit none
    
    real(wp), intent(out) :: t13
    real(wp), intent(in) :: k1(2),k3(2),omg1,omg3,theta1,theta3,n1,n3
    integer, intent(in) :: i1,j1,i2,j2     
    real(wp), intent(in), dimension(:,:) :: fkall
    
    real(wp), dimension(2) :: k2,k4,kc,kp 
    real(wp) omg2,omg4,theta2,theta4,thep
    real(wp) n2,n4,f2,f4,k2old(2)
    real(wp) coef,jacob,pi,twopi,ka,kb,mp,mq,mk2,ds
    integer i,j,isign,is1,js1,is2,js2
    real(wp) locus(-nfhos:nfhos,8),the2next,k2next(2)
    real(wp) eps, domg
    integer count
    
    eps = 1.0e-13
    domg = twopi * dfq
    kp = k1 - k3
    call kk2omgthe(kp,mp,thep,fr2)
    mp = sqrt(dot_product(kp,kp))
    kc = -kp
    mq = (k1(1)**2+k1(2)**2)**0.25 - (k3(1)**2+k3(2)**2)**0.25
    
    ka = 0.25 * (-mq + sqrt(2*mp - mq**2))**2
    if (mq .lt. 0) then
       kb = ((-mp - mq**2) / (2 * mq))**2
    else
       kb = ((mp - mq**2) / (2 * mq))**2
    end if
    
    t13 = 0.0
    count = 0
    is1 = nfhos / 4
    js1 = nthos / 6
    is2 = nfhos / 3
    js2 = nthos * 4 / 7
    !     special case
    if (i1 .eq. i2 .and. j1 .ne. j2) then
       ka = 0.5 * sqrt(dot_product(kp,kp))
       kb = sqrt((twopi * fqh)**2 * fr2)
       if (ka <= kb) then
          do isign = -1,+1,2
             do i = 1, nfhos
                mk2 = ka + (kb - ka) * (i - 1) / (nfhos - 1)
                omg2 = sqrt(mk2 / fr2)
                call solveang_loci_s(k2,theta2,ds,k1,k3,k4,omg1,omg2,omg3,omg4,theta1,&
                     theta3,theta4,mk2,mp,kp,thep,k2old,domg,isign)
                call checkresonant(k1,k2,k3,k4,omg1,omg2,omg3,omg4,1)
                if (isnan(k2(1))) then
                   print *,"incorrect k2 here!"
                   print *,"i-j",i1,i2,j1,j2
                   print *,"k1-k4",k1,k2,k3,k4
                   print *,"mk2",mk2
                   print *,"omg2",omg2
                   print *,"ka",ka
                   print *,"kb",kb
                   print *,"isign",isign
                   stop
                end if
                locus(isign * i,1:2) = k2
                locus(isign * i,3:4) = k4
                locus(isign * i,5) = omg2
                locus(isign * i,6) = omg4
                locus(isign * i,7) = theta2
                locus(isign * i,8) = theta4
             end do
          end do
          do i = -nfhos,nfhos
             if (i /= 0) then
                k2 = locus(i,1:2)
                k4 = locus(i,3:4)
                omg2 = locus(i,5)
                omg4 = locus(i,6)
                theta2 = locus(i,7)
                theta4 = locus(i,8)
                if (i == nfhos) then
                   k2next = locus(-nfhos,1:2)
                   the2next = locus(-nfhos,7)
                else
                   if (i == -1) then
                      k2next = locus(1,1:2)
                      the2next = locus(1,7)
                   else
                      k2next = locus(i + 1,1:2)
                      the2next = locus(i + 1,7)
                   end if
                end if
                ds = sqrt(dot_product(k2 - k2next,k2 - k2next))
                if (i1 == is1 .and. j1 == js1 .and. i2 == is2 .and. j2 == js2) then
                   write(58,*) i,k2(1),k2(2),k1(1),k1(2)
                   write(59,*) i,k4(1),k4(2),k3(1),k3(2)
                   write(52,*) i,ds
                end if
                call snlcoef(coef,k1,k2,k3,k4)
                call caljacob(jacob,k1,k2,k3,fr2,omg1,omg2,omg3)
                call f_int(f2,fkall,omg2,theta2)
                call f_int(f4,fkall,omg4,theta4)
                n2 = f2 / 2 / omg2**4 / fr2**2
                n4 = f4 / 2 / omg4**4 / fr2**2
                t13 = t13 + 2 * ds * heavi(k1,k2,k3,k4) * coef * jacob &
                     * (n1 * n3 * (n4 - n2) + n2 * n4 * (n3 - n1))
                if (isnan(t13)) then
                   print *,"nan t13 here!"
                end if
                count = count + 1
             end if
          end do
       end if
    else
       do isign = -1,+1,2
          do i = 1, nfhos
             mk2 = ka + (kb - ka) * (i - 1) / (nfhos - 1)
             omg2 = sqrt(mk2 / fr2)
             call solveang_loci(k2,theta2,ds,k1,k3,k4,omg1,omg2,omg3,omg4,theta1,theta3,&
                  theta4,mk2,mp,kp,thep,k2old,domg,isign)
             call checkresonant(k1,k2,k3,k4,omg1,omg2,omg3,omg4,1)                
             if (isnan(k2(1))) then
                print *,"incorrect k2 here!"
                print *,"i-j",i1,i2,j1,j2
                print *,"k1-k4",k1,k2,k3,k4
                print *,"mk2",mk2
                print *,"omg2",omg2
                print *,"ka",ka
                print *,"kb",kb
                print *,"isign",isign
                
             end if
             locus(isign * i,1:2) = k2
             locus(isign * i,3:4) = k4
             locus(isign * i,5) = omg2
             locus(isign * i,6) = omg4     
             locus(isign * i,7) = theta2
             locus(isign * i,8) = theta4          
          end do
       end do
       do i = -nfhos, nfhos
          if (i /= 0) then
             k2 = locus(i,1:2)
             k4 = locus(i,3:4)
             omg2 = locus(i,5)
             omg4 = locus(i,6)
             theta2 = locus(i,7)
             theta4 = locus(i,8)
             if (i == nfhos) then
                k2next = locus(-nfhos,1:2)
                the2next = locus(-nfhos,7)
             else
                if (i == -1) then
                   k2next = locus(1,1:2)
                   the2next = locus(1,7)
                else
                   k2next = locus(i + 1,1:2)
                   the2next = locus(i + 1,7)
                end if
             end if
             ds = sqrt(dot_product(k2 - k2next,k2 - k2next))
             if (i1 == is1 .and. j1 == js1 .and. i2 == is2 .and. j2 == js2) then
                write(55,*) i,k2(1),k2(2),k1(1),k1(2)
                write(54,*) i,k4(1),k4(2),k3(1),k3(2)
                write(57,*) i,ds
             endif
             call snlcoef(coef,k1,k2,k3,k4)
             call caljacob(jacob,k1,k2,k3,fr2,omg1,omg2,omg3)
             call f_int(f2,fkall,omg2,theta2)
             call f_int(f4,fkall,omg4,theta4)
             n2 = f2 / 2 / omg2**4 / fr2**2
             n4 = f4 / 2 / omg4**4 / fr2**2
             t13 = t13 + 2 * ds * heavi(k1,k2,k3,k4) * coef * jacob &
                  * (n1 * n3 * (n4 - n2) + n2 * n4 * (n3 - n1))
             if (isnan(t13)) then
                print *,"nan t13 here!"
             end if
             count = count + 1
          end if
       end do
    end if
    
  end subroutine calt13
  !-----------------------------------------------------------------------

  subroutine solveang_loci(k2,theta2,ds,k1,k3,k4,omg1,omg2,omg3,omg4,theta1,theta3,&
       theta4,mk2,mp,kp,thep,k2old,domg,isign)
    
    !     by xuanting hao
    
    !     this subroutine calculates one point on the locus
    !     refer to van vledder 2006 paper
    
    !     06/18/2014 first edition
    
    implicit none
    
    real(wp), intent(in) :: k1(2),k3(2),omg1,omg2,omg3,theta1,theta3
    integer, intent(in) :: isign
    real(wp), intent(in) :: mk2,mp,kp(2),thep,k2old(2), domg
    
    real(wp), intent(out) :: ds,k2(2),k4(2),omg4,theta2,theta4
    real(wp) eps,tmp
    
    omg4 = omg1 + omg2 - omg3
    tmp = (omg4**4 * fr2**2 - mk2**2 - mp**2) / 2.0 / mk2 / mp
    eps = 1.0e-10
    if (abs(tmp+1) <= eps) then
       theta2 = thep - isign * pi
    else if (abs(tmp-1) <= eps) then
       theta2 = thep
    else
       theta2 = thep + isign * acos(tmp)
    end if
    call regulate(theta2)
    call omgthe2kk(k2,omg2,theta2,fr2)
    k4 = k1 + k2 - k3
    call kk2omgthe(k4,omg4,theta4,fr2)
    call regulate(theta4)
    omg4 = omg1 + omg2 - omg3
    
  end subroutine solveang_loci
!-------------------------------------------------------------------------

!--------------------------------------------------------------------------

  subroutine solveang_loci_s(k2,theta2,ds,k1,k3,k4,omg1,omg2,omg3,omg4,theta1,theta3,&
       theta4,mk2,mp,kp,thep,k2old,domg,isign)
    
    !     by xuanting hao
    
    !     this subroutine calculates one point on a linear locus 
    !     refer to van vledder 2006 paper
    
    !     06/26/2014 first edition
    
    implicit none
    
    real(wp), intent(in) :: k1(2),k3(2),omg1,omg2,omg3,theta1,theta3
    integer, intent(in) :: isign
    real(wp), intent(in) :: mk2,mp,kp(2),thep,k2old(2), domg
    
    real(wp), intent(out) :: ds,k2(2),k4(2),omg4,theta2,theta4
    
    real(wp) eps,tmp, thetas
    
    omg4 = omg1 + omg2 - omg3
    eps = 1.0e-12
    thetas = 0.5 * (theta1 + theta3)
    tmp = 1 - mp**2 / 2 / mk2**2
    if (abs(tmp+1) <= eps) then
       theta2 = thetas - 0.5 * isign * pi
    else if (abs(tmp-1) <= eps) then
       theta2 = thetas
    else
       theta2 = thetas + 0.5 * isign * acos(tmp)
    endif
    call regulate(theta2)
    call omgthe2kk(k2,omg2,theta2,fr2)
    k4 = k1 + k2 - k3
    call kk2omgthe(k4,omg4,theta4,fr2)
    call regulate(theta4)
    omg4 = omg1 + omg2 - omg3  
    
  end subroutine solveang_loci_s

!-------------------------------------------------------------------------  
  subroutine int_the1(sum1,omg,theta,k,omg1,f,fkall)
    
    !     by xuanting hao                                                                     
    
    !     this subroutine integrates over theta1                                              
    !     time complexity o(n^2)                                                              
    !     05/28/2014    first edition                                                         
    !     07/01/2014    remove the use of module                                              
    
    
    implicit none
    
    real(wp), intent(in)::omg,theta,f,omg1
    real(wp), dimension(2),intent(in)::k
    real(wp), intent(in),dimension(:,:) :: fkall
    real(wp), intent(out) :: sum1
    
    real(wp) k1(2),theta1,f1
    real(wp) sum2, afun,lambda,xmin,xmax,x
    integer i,j,n,isign
    
    sum1 = 0.0
    lambda = omg/omg1
    afun = ((lambda + 1)**4 - 4 * (lambda**4 + 1)) / (8*lambda**2)
    n = nthos / 4
    do isign = -1,1,2
       !     [-1,1]                                                                              
       if (afun .lt. -1) then
          xmin = -1.0
          xmax = 1.0
          do i = 1, n
             x = (xmax + xmin) / 2 + (xmax - xmin) / 2 * cos((2*i - 1) * pi / 2 / n)
             theta1 = theta + isign * acos(x)
             call regulate(theta1)
             call f_int(f1,fkall,omg1,theta1)
             call omgthe2kk(k1,omg1,theta1,fr2)
             if (isnan(theta1)) then
                print *,"here!"
             endif
             call int_omg2(sum2,omg,theta,omg1,theta1,k,k1,f,f1,fkall)
             sum1 = sum1 + (pi / n) * sum2 * sqrt(abs(cos(theta-theta1)-afun)) / sqrt(abs(afun-x))
             if (isnan(sum1)) then
                print *,"here!"
             endif
          enddo
       else
          !     [-1,afun]                                                                           
          xmin = -1.0
          xmax = afun
          do i = 1, n
             x = (xmax + xmin) / 2 + (xmax - xmin) / 2 * cos((2*i - 1) * pi / 2 / n)
             theta1 = theta + isign * acos(x)
             call regulate(theta1)
             call f_int(f1,fkall,omg1,theta1)
             call omgthe2kk(k1,omg1,theta1,fr2)
             if (isnan(theta1)) then
                print *,"here!"
             endif
             call int_omg2(sum2,omg,theta,omg1,theta1,k,k1,f,f1,fkall)
             sum1 = sum1 + (pi / n) * sum2 * sqrt(abs(cos(theta-theta1)-afun)) / sqrt(1 - x)
          enddo
          !     [afun,1]                                                                            
          xmin = afun
          xmax = 1.0
          do i = 1, n
             x = (xmax + xmin) / 2 + (xmax - xmin) / 2 * cos((2*i - 1) * pi / 2 / n)
             theta1 = theta + isign * acos(x)
             call regulate(theta1)
             call f_int(f1,fkall,omg1,theta1)
             call omgthe2kk(k1,omg1,theta1,fr2)
             if (isnan(theta1)) then
                print *,"here!"
             endif
             call int_omg2(sum2,omg,theta,omg1,theta1,k,k1,f,f1,fkall)
             sum1 = sum1 + (pi / n) * sum2 * sqrt(abs(cos(theta-theta1)-afun)) / sqrt(1 + x)
          enddo
       endif
    enddo
    
  end subroutine int_the1
!-----------------------------------------------------------

  subroutine int_omg2(sum,omg,theta,omg1,theta1,k,k1,f,f1,fkall)
    
    !     by xuanting hao                                                                                        
    
    !     this subroutine integrates over omg2                                                                   
    !     time complexity o(n)                                                                                   
    !     05/27/2014   first edition                                                                             
    !     07/01/2014   remove the use of module      
    
    implicit none
    
    real(wp) :: omg,theta,omg1,theta1,f,f1
    real(wp), dimension(2),intent(in)::k,k1
    
    real(wp), intent(in), dimension(:,:)::fkall
    real(wp), intent(out) :: sum
    
    real(wp) f2,f3
    real(wp), dimension(2)::kav,k2,k3
    real(wp) omga,ka,epsa,tmpa
    real(wp) omg2,omg3,theta2,theta3
    integer i,j,n,isign
    real(wp) minomg2,maxomg2
    real(wp) hfun,den,coef
    integer ite,nite
    
    nite = 1
    n = nthos / 4
    omga = omg + omg1 
    kav = k + k1
    ka = sqrt(kav(1)**2 + kav(2)**2)
    epsa = 2 * ka / omga**2 / fr2
    minomg2 = 0.5 * omga * (1 - 0.5 * epsa)
    if (epsa .gt. 1.0) then
       maxomg2 = 0.5 * omga * (1 - sqrt(epsa - 1.0))
    else
       maxomg2 = 0.5 * omga
    endif
    
    sum = 0.0
    
    if (epsa .gt. 1) then
       do i = 1, n
          omg2 = (maxomg2 + minomg2) / 2 + (maxomg2 - minomg2) / 2 * cos((2*i - 1) * pi / 2 / n)
          omg3 = omga - omg2
          
          do isign = -1,1,2
             do ite = 1, nite
                call solveang(tmpa,theta2,theta3,omg,omg1,omg2,omg3,theta,theta1,k2,k3,ka,kav,isign)
                call checkresonant(k,k1,k2,k3,omg,omg1,omg2,omg3,ite)
             enddo
             call f_int(f2,fkall,omg2,theta2)
             call f_int(f3,fkall,omg3,theta3)
             call snlcoef(coef,k,k1,k2,k3)
             
             hfun = f2*f3*(f*omg1**4 + f1*omg**4) - f1*f*(f2*omg3**4 + f3*omg2**4)
             den = omg1 * omg2 * omg3 * sqrt(omga) * sqrt((ka/fr2 + omg3**2)**2 - omg2**4) * sqrt(omga - omg2 - maxomg2)
             sum = sum + (pi / n) * 4 * hfun * coef / den
             if (isnan(sum)) then
                call dealnan(k,k1,k2,k3,omg,omg1,omg2,omg3,theta,theta1,theta2,theta3)
             endif
          enddo
       enddo
       
    else 
       do i = 1, n
          omg2 = (omga * epsa / 4) * ((i - 1.0) / n)**2 + minomg2
          omg3 = omga - omg2
          do isign = -1,1,2
             do ite = 1, nite
                call solveang(tmpa,theta2,theta3,omg,omg1,omg2,omg3,theta,theta1,k2,k3,ka,kav,isign)
                call checkresonant(k,k1,k2,k3,omg,omg1,omg2,omg3,ite)
             enddo
             call f_int(f2,fkall,omg2,theta2)
             call f_int(f3,fkall,omg3,theta3)
             call snlcoef(coef,k,k1,k2,k3)
             hfun = f2*f3*(f*omg1**4 + f1*omg**4) - f1*f*(f2*omg3**4 + f3*omg2**4)
             den = omg1*omg2*omg3*sqrt(omga) * sqrt((ka/fr2 + omg3**2)**2 - omg2**4) &
                  * sqrt((omg2 - 0.5*omga)**2 - 0.25 * omga**2*(epsa - 1))
             sum = sum +  4 * hfun * coef / den * 2 * sqrt(maxomg2 - minomg2) / n
             if (isnan(sum)) then
                call dealnan(k,k1,k2,k3,omg,omg1,omg2,omg3,theta, theta1,theta2,theta3)
             endif
          enddo
       enddo
    endif
    
  end subroutine int_omg2
  !--------------------------------------------------------------------------- 
  
  !-------------------------------------------------------------------------                                   
  
  subroutine snlcoef(c,k1,k2,k3,k4)
    
    !     by xuanting hao                                                                                        
    
    !     this subroutine is a second edition to calculate the coefficient                                       
    !     c(k1,k2,k3,k4) for deep water wave                                                                     
    !     see lavrenov's paper                                                                                   
    
    !     04/07/2014   first edition!                                                                            
    !     05/26/2014   second edition increase robustness                                                        
        
    implicit none
    
    real(wp),dimension(2),intent(in)::k1,k2,k3,k4
    real(wp), intent(out) :: c

    real(wp) d
    real(wp) tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
    real(wp) km1,km2,km3,km4
    real(wp) omg1,omg2,omg3,omg4
    real(wp) den1,den2,den3
    real(wp) eps
    
    eps = 1.e-16

    km1 = sqrt(dot_product(k1,k1))
    km2 = sqrt(dot_product(k2,k2))
    km3 = sqrt(dot_product(k3,k3))
    km4 = sqrt(dot_product(k4,k4))
    omg1 = sqrt(km1 / fr2)
    omg2 = sqrt(km2 / fr2)
    omg3 = sqrt(km3 / fr2)
    omg4 = sqrt(km4 / fr2)
    
    den1 = sqrt(dot_product(k1+k2,k1+k2)) / fr2 - (omg1 + omg2)**2
    if (abs(den1) .gt. eps) then
       tmp1 = 2*(omg1 + omg2)**2 * (km1*km2 - dot_product(k1,k2)) * (km3*km4 - dot_product(k3,k4)) / den1
    else
       tmp1 = 0
    endif
    den2 = sqrt(dot_product(k1-k3,k1-k3)) / fr2 - (omg1 - omg3)**2
    if (abs(den2) .gt. eps) then
       tmp2 = 2*(omg1 - omg3)**2 * (km1*km3 + dot_product(k1,k3)) * (km2*km4 + dot_product(k2,k4)) / den2
    else
       tmp2 = 0
    endif
    den3 = sqrt(dot_product(k1-k4,k1-k4)) / fr2 - (omg1 - omg4)**2
    if (abs(den3) .gt. eps)  then
       tmp3 = 2*(omg1 - omg4)**2 * (km1*km4 + dot_product(k1,k4)) * (km2*km3 + dot_product(k2,k3)) / den3
    else
       tmp3 = 0
    endif
    tmp4 = 0.5*(dot_product(k1,k2) * dot_product(k3,k4)) + 0.5*(dot_product(k1,k3) * dot_product(k2,k4)) &
         + 0.5*(dot_product(k1,k4) * dot_product(k2,k3))
    tmp5 = -0.25*fr2**2*(dot_product(k1,k2) + dot_product(k3,k4)) * (omg1 + omg2)**4 &
         + 0.25*fr2**2*(dot_product(k1,k3) + dot_product(k2,k4)) * (omg1 - omg3)**4 &
         + 0.25*fr2**2*(dot_product(k1,k4) + dot_product(k2,k3)) * (omg1 - omg4)**4
    tmp6 = fr2**3 * (omg1 + omg2)**2 *(omg1 - omg3)**2 * (omg1 - omg4)**2 * (km1 + km2 + km3 + km4)
    tmp7 = 2.5 * km1 * km2 * km3 * km4
    d = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7
    c = 0.25 * pi * d**2 / fr2**4 / omg1 / omg2 / omg3 / omg4
    
  end subroutine snlcoef

!-------------------------------------------------------------------------        

!--------------------------------------------------------------------------                                  
  subroutine dealnan(k,k1,k2,k3,omg,omg1,omg2,omg3,theta,theta1,theta2,theta3)
    !     by xuanting hao                                                                                        
    !     this subroutine exports wave components when the code blows up.                                        
    implicit none
    
    real(wp),dimension(2),intent(in)::k,k1,k2,k3
    real(wp),intent(in)::omg,omg1,omg2,omg3,theta,theta1,theta2,theta3
    
    print *,"myid=",myid
    print *,"k-k3"
    print *,k,k1,k2,k3
    print *,"omg-omg3"
    print *,omg,omg1,omg2,omg3
    print *,"theta-theta3"
    print *,theta,theta1,theta2,theta3
    
    stop
    
  end subroutine dealnan
  !----------------------------------------------------------------------- 
  
  !-------------------------------------------------------- 
  
  subroutine solveang(tmpa,theta2,theta3,omg,omg1,omg2,omg3,theta,theta1,k2,k3,ka,kav,isign)
    
    !     by xuanting hao       
    !     this subroutine calculates tmpa,theta2 and theta3      
    !     so that the resonant condition "omg + omg1 = omg2 + omg3" is satisfied. 
    !     see laverenov's paper                                    
    !     06/17/2014   second edition    
    implicit none
    
    real(wp), intent(in)::omg,omg1,omg2,theta,theta1,ka
    real(wp), dimension(2),intent(in)::k2, kav
    integer, intent(in) :: isign
    
    real(wp), intent(out) :: tmpa,theta2,theta3,omg3
    real(wp), intent(out), dimension(2) :: k3
    
    real(wp) eps,num1,num2,num3

    eps = 1.0e-10
    num1 = (omg**2*cos(theta)+omg1**2*cos(theta1)) * fr2 / ka
    num2 = (ka**2 + omg2**4 * fr2**2 - omg3**4 * fr2**2) / (2*ka*omg2**2*fr2)
    num3 = (ka**2 + omg3**4*fr2**2 - omg2**4 * fr2**2) / (2*ka*omg3**2*fr2)
    if (abs(num1+1) <= eps) then
       tmpa = -pi
    else if (abs(num1-1) <= eps) then
       tmpa = 0
    else
       tmpa = sign(acos(num1),kav(2))
    end if
    if (abs(num2+1) < eps) then
       theta2 = tmpa + isign*(-pi)
    else if (abs(num2-1) < eps) then
       theta2 = tmpa
    else
       theta2 = tmpa + isign*acos(num2)
    end if
    if (isnan(theta2)) then
       print *,"num2=",num2
       print *,"omg-omg3",omg,omg1,omg2,omg3
       print *,"ka",ka
       print *,"kav",kav
    end if
    call regulate(theta2)
    
    call omgthe2kk(k2,omg2,theta2,fr2)
    k3 = kav - k2
    call kk2omgthe(k3,omg3,theta3,fr2)
    omg3 = omg + omg1 - omg2
    
  end subroutine solveang

!------------------------------------------------------------------------       
  
  subroutine caljacob(jacob,k1,k2,k3,fr2,omg1,omg2,omg3)
    !     by xuanting hao
    
    !     this subroutine calculates the jacobian, see van vledder's paper
    
    !     06/18/2014  first edition
    
    implicit none
    real(wp), intent(out) :: jacob
    real(wp), intent(in), dimension(2) :: k1,k2,k3
    real(wp), intent(in) :: omg1,omg2,omg3,fr2
    
    real(wp) cg1(2),cg2(2),k4(2),omg4,dcg(2)
    
    k4 = k1 + k2 - k3
    omg4 = omg1 + omg2 - omg3
    cg1 = k2 / fr2**2 / 2.0 / omg2**3
    cg2 = k4 / fr2**2 / 2.0 / omg4**3
    dcg = cg1 - cg2
    jacob = 1 / sqrt(dcg(1)**2 + dcg(2)**2)
    
  end subroutine caljacob

!-------------------------------------------------------------------------    
  subroutine f_int(f,fk,omg,theta)
    !     by xuanting hao                                                        
    !     this subroutine can calculate f when the point   
    !     is not exactly on the grid.  
    !     theta = [-pi,pi)                
    !     03/17/2014  first edition       
    !     04/01/2014  second edition  
    !     05/18/2014  complete the parallel version
    !     06/06/2014  remove the usage of module  
    
    implicit none
    
    !     output---f=f(omg,the)   
    real(wp), intent(out) :: f
    !     input---fk = fk(fq,the)      
    real(wp), intent(in), dimension(:,:) :: fk
    
    real(wp) :: omg,theta,fq
    real(wp) :: x1,y1,x2,y2
    real(wp) :: w1,w2,w3,w4
    integer i1,j1,i2,j2,i,j
    real(wp) :: ftmp
    real(wp), allocatable, dimension(:,:) :: extfk

    fq = omg / 2 / pi
    
    allocate(extfk(size(fk,1)+2,size(fk,2)+2))
    extfk(1:nfhos,1:nthos) = fk
    extfk(1:nfhos,nthos + 1) = fk(:,1)
        
    !     calculate the indexes                                                                                  
    i1 = ceiling(fq / dfq) + 1
    i2 = floor(fq / dfq) + 1
    if (i1 .eq. i2) then
       i2 = nint(fq / dfq) + 1
       i1 = i2 + 1
    endif
    j1 = ceiling(theta / dth + nthos / 2) + 1
    j2 = floor(theta / dth + nthos / 2) + 1
    if (j1 .eq. j2) then
       j2 = nint(theta / dth + nthos / 2) + 1
       j1 = j2 + 1
    endif

!    i1 = 1
!    do while (fq >= i1 * dfq) 
!       i1 = i1 + 1
!    end do
!    i2 = i1 - 1
!    j1 = 1
!    do while (theta >= (j1 - 1) * dth - pi)
!       j1 = j1 + 1
!    end do
!    j2 = j1 - 1
    
    if (theta > thl .and. theta < thh .and. i1 <= maxnfq .and. i2 >= minnfq) then
       if (j1 .gt. nthos + 1.or. j2 .lt. 1) then
          print *,"j incorrect!"
          print *,"fq,theta,j1,j2",fq,theta,j1,j2
          print *,"111"
       endif
       x1 = (i1 - 1) * dfq
       x2 = (i2 - 1) * dfq
       y1 = (j1 - 1) * dth - pi
       y2 = (j2 - 1) * dth - pi
       !     evaluate weighting factors using 2d linear interpolation 
       call intp2d(w1,w2,w3,w4,x1,y1,x2,y2,fq,theta)
       f = (w1 * extfk(i1,j1) + w2 * extfk(i1,j2) &
       + w3 * extfk(i2,j2) + w4 * extfk(i2,j1)) / twopi
    else
       f = 0.0
       !     extrapolate the spectrum assuming a f^-4 tail                           
       if (fq .gt. (maxnfq-1)*dfq) then
          if (j1 .gt. nthos + 1 .or. j2 .lt. 1) then
             print *,"j incorrect!"
             print *,"fq,theta,j1,j2",fq,theta,j1,j2
             print *,"222"
          end if
          
          ftmp = extfk(maxnfq,j2) * ((j1-1)*dth-pi - theta) / dth &
               + extfk(maxnfq,j1) * (theta - (j2-1)*dth + pi) / dth
          f = ftmp * ((maxnfq-1)*dfq / fq)**4 / twopi
       end if
    endif

    deallocate(extfk)
  end subroutine f_int

!------------------------------------------------------------------------  

  subroutine intp2d(w1,w2,w3,w4,x1,y1,x2,y2,x,y)
    !     by xuanting hao  
    !     this subroutine determines the weighting factor using the bilinear interpolation   
    !     02/16/2014 first edition                                                      
    !     02/22/2014 improve robustness                                       
    
    !                                                                                                            
    !                  x2,y1 w4        x1,y1  w1                                    
    !                       ---------------                                           
    !                       |             |                                        
    !                       |   x,y       |                                        
    !                       |             |                                              
    !                       |             |                                            
    !                       ---------------                                          
    !               x2,y2  w3        x1,y2  w2                                         
    !                                                                                     
    !   2d interpolation formula:                       
    !   f(x,y) = w1 * f(x1,y1) + w2 * f(x2,y2) + w3 * f(x3,y3) + w4 * f(x4,y4)            
    !       where:                                                          
    !       w1 = (x - x2) * (y - y2) / (x2 - x1) /(y2 - y1)                
    !       w2 = (x1 - x) * (y - y2) / (x2 - x1) /(y2 - y1)                 
    !       w3 = (x1 - x) * (y1 - y) / (x2 - x1) /(y2 - y1)                    
    !       w4 = (x - x2) * (y1 - y) / (x2 - x1) /(y2 - y1)                       
    implicit none
    
    real(wp) :: x1,x2,y1,y2,x,y,w1,w2,w3,w4
    if (x1 .ne. x2 .and. y1 .ne. y2) then
       w1 = (x - x2) * (y - y2) / (x2 - x1) /(y2 - y1)
       w2 = (x1 - x) * (y - y2) / (x2 - x1) /(y2 - y1)
       w3 = (x1 - x) * (y1 - y) / (x2 - x1) /(y2 - y1)
       w4 = (x - x2) * (y1 - y) / (x2 - x1) /(y2 - y1)
    endif
    
    if (x1 .eq. x2 .and. y1 .ne. y2) then
       w1 = (y - y2) / (y1 - y2) / 2
       w2 = (y1 - y) / (y1 - y2) / 2
       w3 = w2
       w4 = w1
    endif
    
    if (x1 .ne. x2 .and. y1 .eq. y2) then
       w1 = (x - x2) / (x1 - x2) / 2
       w2 = w1
       w3 = (x1 - x) / (x1 - x2) / 2
       w4 = w3
    endif
    
    if (x1 .eq. x2 .and. y1 .eq. y2) then
       w1 = 0.25
       w2 = 0.25
       w3 = 0.25
       w4 = 0.25
    endif

  end subroutine intp2d
  
!-----------------------------------------------------------------------  
  
end module src_stat

