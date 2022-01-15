module post_proc

  use hos
  use hdf_io_hos
  use io_hos
  use src_stat

  public :: post_hos_resio, post_wind_wave

contains

  subroutine post_hos_resio

    implicit none

    integer it,icon
    real(wp) :: sigma, period
    real(wp) :: tmp
    real(wp), allocatable, dimension(:,:) :: eo, vo
    real(wp), allocatable, dimension(:) :: sk, sf, sk1, sk2
    real(wp) :: vol, flux, pe, ke, ener, rms, hs, skew, kurt
    real(wp), allocatable, dimension(:) :: skx, skx_ave, snl1d
    real(wp), allocatable, dimension(:,:) :: fskxy, fs_ave
    real(wp), allocatable, dimension(:,:) :: sds
    real(wp) beta1, beta2
    
    real(wp) :: t_ave, fetch_ave, dwk, omega, sw
    real(wp) :: tp,kp,kpx, kpy,omgp, cp, rek,cg, charv,t2f, old_time, x,y,  m0, m1, m2, nu
    integer i,j,root, i_ave, nmax, icount, valueRSS
    character (len=100) :: fileid, fparam, fdata
    real(wp) :: clen, ctime, cspeed, ratio

    icon = 16
    root = 0
    call jonswap_init
    allocate(eo(nxhos,nyhos))
    allocate(vo(nxhos,nyhos))
    
    allocate(fskxy(nxhos/2,-nyhos/2:nyhos/2),fs_ave(nxhos/2,-nyhos/2:nyhos/2))
    allocate(sds(nxhos,nyhos))
    allocate(skx(nxhos/2),skx_ave(nxhos/2))

    beta1 = 16.0
    beta2 = 30.0
    ratio = 0.6

    if (ist_hos == 1) then
       sigma = (nswavex * pex_hos * (1/fr2 + (nswavex*pex_hos)**2/bond/fr2))**0.5
    else
       sigma = (nswavex * pex_hos / fr2)**0.5
    end if
    
    period = twopi / sigma
    
    if (dt_hos <= 0) then
       dt_hos = period / ntp
    end if
    
    nmax = nxhos/2-1
    dwk = pex_hos
    allocate(sk(nmax), sf(nfhos), sk1(nmax), sk2(nmax))
    icount = 0
    ioutd_hos = 0
    time_hos = 0

!    write(fileid,*) it
!    fileid = trim(adjustl(fileid))
    
!    fparam = trim("restart_param_hos.dat"//fileid)
!    fdata = trim("restart_hos.h5"//fileid)
    
    call read_hos(time_hos, ioutd_hos, ioutc_hos, eta_hos,vps_hos,pa_hos) !,fparam, fdata)

    call wavenum(wvn)
    call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
    call zeta(eta_hos,zp)
    call boundvp(vps_hos,r_hos,zp)
    call wsurf(w_hos,r_hos,zp,wvn)
    call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
    call righ(eta_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos,w_hos,feta,fvps,pa_hos,pa_st,1)

    call alltoone(eta_hos,eo)
    call alltoone(vps_hos,vo)

    call fullspec_xy(eo,fskxy)
    call get_moment(m0,fskxy,0)
    call get_moment(m1,fskxy,1)
    call get_moment(m2,fskxy,2)

    call get_peak(fskxy, kpx, kpy, kp, tp, omgp, cp, charv)

    call spec_1df(eo, sf, nfhos, dfq)

    if (myid == 0) then
       open(667,FILE="hos_eta.dat")
       open(668,FILE="hos_E_kx_ky.dat")
       open(669,FILE="hos_E_f.dat")
       
       write(667,*) nxhos, nyhos
       write(667,*) time_hos
       do j = 1, nyhos
          y = (j - 1) * dy_hos
          do i = 1, nxhos
             x = (i - 1) * dx_hos
             write(667,'(25e12.4)') x,y,eo(i,j)
          end do
       end do
       
       write(668,'(25e12.4)') ratio * pex_hos * nxhos / 2.0, ratio * pey_hos * nyhos / 2.0       
       write(668,*) nxhos/2,nyhos/2
       do j = -nyhos/2, nyhos / 2
          y = j * pey_hos
          do i = 1, nxhos / 2
             x = i * pex_hos
             write(668,'(25e12.4)') x,y,fskxy(i,j)
          end do
       end do

       write(669,*) nfhos
       do i = 1, nfhos
          write(669,'(25e12.4)') i*dfq, sf(i)
       end do
       close(667)
       close(668)
       close(669)
    end if
    
  end subroutine post_hos_resio

  subroutine post_wind_wave
    
    implicit none
    
    integer it,icon
    real(wp) :: sigma, period
    real(wp) :: tmp
    real(wp), allocatable, dimension(:,:) :: eo, vo
    real(wp), allocatable, dimension(:) :: sk, sf, sk1, sk2
    real(wp) :: vol, flux, pe, ke, ener, rms, hs, skew, kurt
    real(wp), allocatable, dimension(:) :: skx, skx_ave, snl1d
    real(wp), allocatable, dimension(:,:) :: fskxy, fs_ave, fk_fit, snl, snlall
    real(wp), allocatable, dimension(:,:) :: sds
    real(wp) beta1, beta2
    
    real(wp) :: t_ave, fetch_ave, dwk, omega, sw
    real(wp) :: tp,kpe, kpn,kpx, kpy,omgp, cp, rek,cg, charv
    real(wp) :: t2f, old_time, x,y,  m0, m1, m2, nu
    integer i,j,root, i_ave, nmax, icount, valueRSS, n_ave
    character (len=100) :: fileid, fparam, fdata
    real(wp) :: clen, ctime, cspeed, u10_fit, fetch_fit, width, bfi, ursell
    real(wp) act_tot,mom_x_tot,mom_y_tot,ene_tot,ent_tot
    
    !  real(wp) :: eta0(nxhos),myeta(nxhos)                                                    
    !  real(wp) :: err                                                                         
    
    icon = 16
    root = 0
    call jonswap_init
    allocate(eo(nxhos,nyhos))
    allocate(vo(nxhos,nyhos))
    
    allocate(fskxy(nxhos/2,-nyhos/2:nyhos/2),fs_ave(nxhos/2,-nyhos/2:nyhos/2))
    allocate(sds(nxhos,nyhos))
    allocate(skx(nxhos/2),skx_ave(nxhos/2))

    read(14,*) clen, cspeed, ctime
    
    if (myid == 0) then
       write(55,*) " VARIABLES = t, fetch, ener, rms,hs,skew, kurt,tpe,kpe,kpn,omgp, cp,charv"
       write(55,"(A,F16.8,A,I5)") ' ZONE T="',time_hos,'" I=', ntime_hos
       write(55,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
       write(55,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'

       write(56,*) " VARIABLES = t,  act_tot,mom_x_tot,mom_y_tot,ene_tot,ent_tot, bfi"
       write(56,"(A,F16.8,A,I5)") ' ZONE T="',time_hos,'" I=', ntime_hos
       write(56,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
       write(56,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'

    end if
    
    t2f = 9500/clen
    old_time = 550.0
    i_ave = 0
    t_ave = 0
    fetch_ave = 0
    skx_ave = 0
    fs_ave = 0
    
    beta1 = 7.0
    beta2 = 20.0
    
    if (myid == 0) then
       !export kx-ky spectrum                                                                                               
       write(34,*) " VARIABLES = kx,ky,skxy"
       write(34,"(A,F16.8,A,I5,A,I5,A)") ' ZONE T="',t_ave,'" I=', nxhos/2, ' J=',nyhos+1,' F=POINT'
       write(34,"(A,F16.8,A)") ' AUXDATA t="',t_ave,'"'
       write(34,"(A,F16.8,A)") ' AUXDATA fetch="',fetch_ave,'"'
       write(34,"(A,F16.8,A)") ' AUXDATA tp="',tp,'"'
       write(34,"(A,F16.8,A)") ' AUXDATA kpe="',kpe,'"'
       write(34,"(A,F16.8,A)") ' AUXDATA kpn="',kpn,'"'
       write(34,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
       write(34,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'
       
       do j = -nyhos/3, nyhos / 3
          y = j * pey_hos
          do i = 1, nxhos / 3
             x = i * pex_hos
             write(34,'(25e12.4)') x,y,fs_ave(i,j)
          end do
       end do

       write(35,*) " VARIABLES = freq,theta,snl"
       write(35,"(A,F16.8,A,I5,A,I5,A)") ' ZONE T="',time_hos,'" I=', nfhos, ' J=',nthos,' F=POINT'
       write(35,"(A,F16.8,A)") ' AUXDATA t="',time_hos,'"'
       write(35,"(A,F16.8,A)") ' AUXDATA kpe="',kpe,'"'
       write(35,"(A,F16.8,A)") ' AUXDATA kpn="',kpn,'"'
       write(35,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
       write(35,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'
       
       do j = 1, nthos
          do i = 1, nfhos
             write(35,'(25e12.4)') i*dfq,j*dth,0.0
          end do
       end do

    end if
    
    if (ist_hos == 1) then
       sigma = (nswavex * pex_hos * (1/fr2 + (nswavex*pex_hos)**2/bond/fr2))**0.5
    else
       sigma = (nswavex * pex_hos / fr2)**0.5
    end if
    
    period = twopi / sigma
    
    if (dt_hos <= 0) then
       dt_hos = period / ntp
    end if
    
    if (myid == 0) then
       print *,"dt_hos=",dt_hos
       
       write(37,*) " VARIABLES = f, sf, omg, somg, k, sk, nomg, nk"
       write(37,"(A,F16.8,A,I5)") ' ZONE T="',time_hos,'" I=', nfhos
       write(37,"(A,F16.8,A)") ' AUXDATA t="',t_ave,'"'
       write(37,"(A,F16.8,A)") ' AUXDATA fetch="',fetch_ave,'"'
       write(37,"(A,F16.8,A)") ' AUXDATA tp="',tp,'"'
       write(37,"(A,F16.8,A)") ' AUXDATA kpe="',kpe,'"'
       write(37,"(A,F16.8,A)") ' AUXDATA kpn="',kpn,'"'
       write(37,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
       write(37,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'
       
       do i = 1, nfhos
          omega = i*dfq * twopi
          sw = alpha / fr2**2 / omega**5 * exp(-1.25*(omega/omega0)**(-4.0)) &
               * gamma**(exp(-(omega-omega0)**2/2/sigma**2/omega0**2)) * twopi
          
          write(37,*) i*dfq, sw, omega, sw/twopi, omega**2*fr2, sw/twopi/fr2/2.0/omega, &
               sw/twopi / omega,  sw/twopi/fr2/2.0/omega**2
       end do
    end if
    
    !  do it = 501, 900                                                                                                       
    
    nmax = nxhos/2-1
    dwk = pex_hos
    allocate(sk(nmax), sf(nfhos), sk1(nmax), sk2(nmax))
    allocate(fk_fit(nfhos,nthos),snlall(nfhos,nthos))
    allocate(snl(nfhos,nthos/ncpu_hos), snl1d(nfhos))
    icount = 0
    ioutd_hos = 0
    time_hos = 0

    n_ave = 5
    do it = 601, 600 + ntime_hos
       !     i_ave = i_ave + 1                                                                                                   
       
       write(fileid,*) it
       fileid = trim(adjustl(fileid))
       
       fparam = trim("restart_param_hos.dat"//fileid)
       fdata = trim("restart_hos.h5"//fileid)
       
       
       call read_hos(time_hos, ioutd_hos, ioutc_hos, eta_hos,vps_hos,pa_hos,fparam, fdata)
       
       !     time_hos = time_hos + dt_hos                                                                                        
       
       !     ioutd_hos = ioutd_hos + 1                                                                                           
       
       !     if (ioutd_hos == noutd_hos) then                                                                                    
       !        if (myid == 0) print *,"i_ave=", i_ave                                                                           
       i_ave = i_ave + 1
       ioutd_hos = 0
       !        read(340,*) eo,vo                                                                                                
       !        if (myid == 0) then                                                                                              
       !           print *,"it=",it," reread successful!"                                                                        
       !        end if                                                                                                           
       icount = icount + 1
       !        eta_hos = eo(:,myid*nyhos/ncpu_hos+1:(myid+1)*nyhos/ncpu_hos)                                                    
       !        vps_hos = vo(:,myid*nyhos/ncpu_hos+1:(myid+1)*nyhos/ncpu_hos)                                                    
       
       if (myid == 0) then
          print *,"t=",time_hos," reread successful!"
       end if
       
       !        call filter_lp(eta_hos,vps_hos,beta1,beta2,sds)                                                                  
       ! Get other quantities                                                                                            
       call wavenum(wvn)
       call derivh(eta_hos,vps_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
       call zeta(eta_hos,zp)
       call boundvp(vps_hos,r_hos,zp)
       call wsurf(w_hos,r_hos,zp,wvn)
       call uvsurf(u_hos,v_hos,w_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos)
       call righ(eta_hos,ex_hos,ey_hos,vpsx_hos,vpsy_hos,w_hos,feta,fvps,pa_hos,pa_st,1)
       
       call alltoone(eta_hos,eo)
       call alltoone(vps_hos,vo)
       
       call spec_1dk(eo, sk, nmax, dwk)      
       call spec_1df(eo, sf, nfhos, dfq)
       
       call spec_x(eo,skx)
       call fullspec_xy(eo,fskxy)
       call get_moment(m0,fskxy,0)
       call get_moment(m1,fskxy,1)
       call get_moment(m2,fskxy,2)
       nu = sqrt(m0*m2/m1**2-1)
       
       call get_peak(fskxy, kpx, kpy, kpe, tp, omgp, cp, charv)
       call get_peak(fskxy, kpx, kpy, kpn, tp, omgp, cp)
 
       if (mod(it,100) == 1) then
          call fit_jonswap(u10_fit,fetch_fit,sf,nfhos,dfq,clen,ctime)
          call jonswap_dir(fk_fit,nfhos,nthos,u10_fit,fetch_fit,dfq, dth,ctime)
          call snl_webb(fk_fit,snl)
          call twod2oned(snl,snl1d)
          call check_snl(fk_fit,snl,act_tot,mom_x_tot,mom_y_tot,ene_tot,ent_tot)
       end if

       call alltoone(snl,snlall,nfhos,nthos)

!       if (myid == 0) print *,u10_fit, fetch_fit

       !        if (time_hos >= 300.0) then                                                                                      
       t2f = t2f + charv * (time_hos - old_time)
       old_time = time_hos
       !        end if                                                                                                           
       
       t_ave = t_ave + time_hos / n_ave
       fetch_ave = fetch_ave + t2f / n_ave
       skx_ave = skx_ave + skx / n_ave
       fs_ave = fs_ave + fskxy / n_ave
       
       call get_stats(vol,flux,pe,ke,ener,eta_hos,vps_hos,feta)
       call skwave(eta_hos,rms,hs,skew,kurt)
       
       call get_width(sf, nfhos, dfq, width)
       bfi = sqrt(kpe**2*rms**2)*sqrt(2.0)/width

       if (myid == 0) then
          
          write(55,'(25e12.4)') time_hos, t2f, ener, rms, hs,skew, kurt,tp,kpe,kpn,omgp, cp, charv

          write(56,'(25e12.4)') time_hos, act_tot,mom_x_tot,mom_y_tot,ene_tot,ent_tot, bfi
          
          write(37,*) " VARIABLES = f, sf, omg, somg, k, sk, nomg, nk"
          write(37,"(A,F16.8,A,I5)") ' ZONE T="',time_hos,'" I=', nfhos
          write(37,"(A,F16.8,A)") ' AUXDATA t="',time_hos,'"'
          write(37,"(A,F16.8,A)") ' AUXDATA fetch="',t2f,'"'
          write(37,"(A,F16.8,A)") ' AUXDATA tp="',tp,'"'
          write(37,"(A,F16.8,A)") ' AUXDATA kpe="',kpe,'"'
          write(37,"(A,F16.8,A)") ' AUXDATA kpn="',kpn,'"'
          write(37,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
          write(37,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'
          do i = 1, nfhos
             omega = twopi *i * dfq
             write(37,*) i*dfq, sf(i), omega, sf(i)/twopi, omega**2*fr2, sf(i)/twopi/fr2/2.0/omega, &
                  sf(i)/twopi / omega,  sf(i)/twopi/fr2/2.0/omega**2
          end do
        
          !export snl spectrum        
          if (mod(it,100) == 1) then
             write(38,*) " VARIABLES = f, snl1d"
             write(38,"(A,F16.8,A,I5)") ' ZONE T="',time_hos,'" I=', nfhos
             write(38,"(A,F16.8,A)") ' AUXDATA t="',time_hos,'"'
             write(38,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
             write(38,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'
             do i = 1, nfhos
                write(38,*) i*dfq, snl1d(i)
             end do


             write(35,*) " VARIABLES = freq,theta,snl"
             write(35,"(A,F16.8,A,I5,A,I5,A)") ' ZONE T="',time_hos,'" I=', nfhos, ' J=',nthos,' F=POINT'
             write(35,"(A,F16.8,A)") ' AUXDATA t="',time_hos,'"'
             write(35,"(A,F16.8,A)") ' AUXDATA kpe="',kpe,'"'
             write(35,"(A,F16.8,A)") ' AUXDATA kpn="',kpn,'"'
             write(35,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
             write(35,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'
             write(35,*) "VARSHARELIST=([1-2]=1)"
             write(35,'(25e12.4)') snlall
          end if
       end if

       

       if (myid == -1) then
          write(38,*) " VARIABLES = kx,ky,sds"
          write(38,"(A,F16.8,A,I5,A,I5,A)") ' ZONE T="',t_ave,'" I=', nxhos/2, ' J=',nyhos/2,' F=POINT'
          do j = 1, nyhos / 2
             y = j * pey_hos
             do i = 1, nxhos / 2
                x = i * pex_hos
                write(38,'(25e12.4)') x, y, sds(i,j)
             end do
          end do
       end if
       if (i_ave == n_ave) then
          if (myid  == 0) then
             !export kx spectrum                                                  
             write(32,*) " VARIABLES = kx, nkx"
             write(32,"(A,F16.8,A,I5)") ' ZONE T="',t_ave,'" I=', nxhos/2
             write(32,"(A,F16.8,A)") ' AUXDATA t="',t_ave,'"'
             write(32,"(A,F16.8,A)") ' AUXDATA fetch="',fetch_ave,'"'
             write(32,"(A,F16.8,A)") ' AUXDATA tp="',tp,'"'
             write(32,"(A,F16.8,A)") ' AUXDATA kpe="',kpe,'"'
             write(32,"(A,F16.8,A)") ' AUXDATA kpn="',kpn,'"'
             write(32,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
             write(32,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'

                         
             do i = 1, nxhos / 2
                write(32,'(25e12.4)') i * pex_hos, skx_ave(i) / sqrt(i*pex_hos/fr2)
             end do
             
             !export kx-ky spectrum                                                                                      
             write(34,*) " VARIABLES = kx,ky,skxy"
             write(34,"(A,F16.8,A,I5,A,I5,A)") ' ZONE T="',t_ave,'" I=', nxhos/2, ' J=',nyhos+1,' F=POINT'
             write(34,"(A,F16.8,A)") ' AUXDATA t="',t_ave,'"'
             write(34,"(A,F16.8,A)") ' AUXDATA fetch="',fetch_ave,'"'
             write(34,"(A,F16.8,A)") ' AUXDATA tp="',tp,'"'
             write(34,"(A,F16.8,A)") ' AUXDATA kpe="',kpe,'"'
             write(34,"(A,F16.8,A)") ' AUXDATA kpn="',kpn,'"'
             write(34,"(A,F16.8,A)") ' AUXDATA clen="',clen,'"'
             write(34,"(A,F16.8,A)") ' AUXDATA ctime="',ctime,'"'
             write(34,*) "VARSHARELIST=([1-2]=1)"
             write(34,'(25e12.4)') fs_ave(1:nxhos/3,-nyhos/3:nyhos/3)
          end if
          i_ave = 0
          t_ave = 0
          fetch_ave = 0
          skx_ave = 0
          fs_ave = 0
       end if
       !       call system_mem_usage(valueRSS)
       !        write(36700+myid,*)"it=",it,"valueRSS=",valueRSS                                                                 
    end do
    
    deallocate(eo,vo)
    deallocate(sk,sf,sk1,sk2)
    deallocate(fskxy, fs_ave)
    deallocate(sds)
    deallocate(skx,skx_ave)
    
  end subroutine post_wind_wave
      
end module post_proc
