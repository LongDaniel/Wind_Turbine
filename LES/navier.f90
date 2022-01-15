module navier
   use decomp
   use param
   use solver_common, only : tmp_x1, tmp_x2, tmp_x3
   use solver_common, only : tmp_x4, tmp_x5, tmp_x6, tmp_x7
   use solver_common, only : tmp_y1
   use solver_common, only : les_filter, calc_uzeta, calc_wzeta
   use hos_param, only: npw, eta_hos, u_hos, v_hos, w_hos
   use discontinuity_smooth, only : tbn_lim_x, tbn_lim_y

   implicit none

   private

   integer :: level

   real(wp), allocatable, dimension(:,:,:), public :: u, v, w, pp
   real(wp), allocatable, dimension(:,:,:), public :: hu, hv, hw
   real(wp), allocatable, dimension(:,:,:), public :: fu, fv, fw
   real(wp), allocatable, dimension(:,:),   public :: pf0

   !> plyunote: add wtforce in y/z direction for AL model
   real(wp), allocatable, dimension(:,:,:), public :: wtforce
   real(wp), allocatable, dimension(:,:,:), public :: wtforce_y
   real(wp), allocatable, dimension(:,:,:), public :: wtforce_z
   real(wp), public :: fturbinex, fturbiney
   !> plyunote: although the flag can be as simple as a logical type, 
   !!           I choose to use real type so that I can reuse the existing
   !!           functions. Another way to reuse the real-type function for
   !!           logical type might be generic programming:
   !!           http://fortranwiki.org/fortran/show/Generic+programming
   real(wp), allocatable, dimension(:,:,:), public, save :: flag_de
      
   ! HOS variables
   ! real(wp), allocatable, dimension(:,:),   public :: etab
   ! real(wp), allocatable, dimension(:,:),   public :: vps
   
   ! real(wp), allocatable, dimension(:,:),   public :: pas
   ! real(wp), allocatable, dimension(:,:),   public :: pas0

   ! real(wp), allocatable, dimension(:,:),   public :: ebx, eby
   
   ! real(wp), allocatable, dimension(:,:),   public :: vpsx, vpsy

   ! real(wp), allocatable, dimension(:,:,:),   public :: wvn
   ! real(wp), allocatable, dimension(:,:,:),   public :: zp
   ! real(wp), allocatable, dimension(:,:,:),   public :: r

   ! real(wp), allocatable, dimension(:,:),   public :: ubs, vbs, wbs  
   
   !end
   
   
   real(wp), allocatable, dimension(:,:),   public :: uzfs
   real(wp), allocatable, dimension(:,:),   public :: vzfs
   real(wp), allocatable, dimension(:,:),   public :: wzfs

   real(wp), allocatable, dimension(:,:),   public :: ut, vt, wt
   real(wp), allocatable, dimension(:,:),   public :: ub, vb, wb

   real(wp), allocatable, dimension(:,:,:),   public :: pu,pv,pw

   
   real(wp), allocatable, dimension(:,:,:) :: u_xi
   real(wp), allocatable, dimension(:,:,:) :: u_psi
   real(wp), allocatable, dimension(:,:,:), public, save :: u_zeta
   real(wp), allocatable, dimension(:,:,:) :: ux, uy, uz_u
   real(wp), allocatable, dimension(:,:,:) :: vx, vy, vz_u
   real(wp), allocatable, dimension(:,:,:) :: wx, wy, wz_w
   
   real(wp), allocatable, dimension(:,:,:) :: lap_u_al, lap_u_ao, lap_u_ar
   real(wp), allocatable, dimension(:,:,:) :: lap_w_al, lap_w_ao, lap_w_ar
   real(wp), allocatable, dimension(:,:,:) :: lap_p_al, lap_p_ao, lap_p_ar

   
   ! LES variables
   real(wp), dimension(:,:,:), allocatable, public :: uf, vf, wf
   real(wp), dimension(:,:,:), allocatable, public :: nut, nutw
   real(wp), dimension(:,:,:), allocatable, public, save :: s11, s12, s13, s22, s23, s33, ss
   real(wp), dimension(:,:,:), allocatable, public :: s11w, s12w, s13w, s22w, s23w, s33w
   real(wp), dimension(:,:,:), allocatable, public :: t11, t12, t13, t22, t23, t33
   real(wp), dimension(:,:,:), allocatable, public :: t11w, t12w, t13w, t22w, t23w, t33w
   real(wp), allocatable, dimension(:,:),   public :: tauwx, tauwy
   
   real(wp), allocatable, dimension(:,:,:), public :: cc

   !navier.f
   public :: navier_init
   public :: get_max_divu
   public :: bottom2, bottom_hos_les
   !> added by plyu
   public :: gradu, gradw, calc_divu_center
   
   !navier_fu.f
   public :: calc_F
   !> added by plyu
   public :: wall_model_v3

   !navier_pressure.f
   public :: init_p_poisson, press_g, correction_us, bc_lnr, volume_lnr2


   !navier_les.f: LES subroutines
   public :: les_init
   public :: div_tau, div_tau_w
   !> added by plyu
   public :: les_filter, get_strain, get_nut, get_SGS_stress
   public :: remove_plane_mean

   
contains

   !========================================================
   !> @brief Allocate and initialize variables in module navier.
   !========================================================


  !input: level_
  !subroutine common variable: level
  subroutine navier_init(level_)
      implicit none

      integer, intent(IN) :: level_

      level = level_

      allocate(u(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(v(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(hu(xsz(1),xsz(2),xsz(3)))
      allocate(hv(xsz(1),xsz(2),xsz(3)))
      allocate(hw(xsz(1),xsz(2),xsz(3)))
      allocate(fu(xsz(1),xsz(2),xsz(3)))
      allocate(fv(xsz(1),xsz(2),xsz(3)))
      allocate(fw(xsz(1),xsz(2),xsz(3)))

      allocate(pp(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(pf0(xsz(1),xsz(2)))

      ! if(myid==2)print*,'xsz(3)',size(u,3) 
      !HOS variables
      ! allocate(pas(nxs,nys/np1))
      ! allocate(pas0(nxs,nys/np1))

      ! allocate(etab(nxs,nys/np1))
      ! allocate(vps(nxs,nys/np1))

      ! allocate(ebx(nxs,nys/np1))
      ! allocate(eby(nxs,nys/np1))

      ! allocate(vpsx(nxs,nys/np1))
      ! allocate(vpsy(nxs,nys/np1))      

      ! allocate(wvn(nxs,nys/np1,npw))      
      ! allocate(zp(nxs,nys/np1,npw))
      ! allocate(r(nxs,nys/np1,npw))

      ! allocate(ubs(nxs,nys/np1))
      ! allocate(vbs(nxs,nys/np1))
      ! allocate(wbs(nxs,nys/np1))      
      !end HOS
      
      allocate(uzfs(xsz(1),xsz(2)))
      allocate(vzfs(xsz(1),xsz(2)))
      allocate(wzfs(xsz(1),xsz(2)))
      
      allocate(ut(xsz(1),xsz(2)))
      allocate(vt(xsz(1),xsz(2)))
      allocate(wt(xsz(1),xsz(2)))
      allocate(ub(xsz(1),xsz(2)))
      allocate(vb(xsz(1),xsz(2)))
      allocate(wb(xsz(1),xsz(2)))

      ! plyunote: modified. conclude from line navier_fu.f90:276
      ! size of wtforce should be the same with fu  
      !allocate(wtforce(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(wtforce(xsz(1),xsz(2),xsz(3)))
      allocate(wtforce_y(xsz(1),xsz(2),xsz(3)))
      allocate(wtforce_z(xsz(1),xsz(2),xsz(3)))
      allocate(flag_de(xsz(1),xsz(2),xsz(3)))
      wtforce = 0.0_wp; wtforce_y = 0.0_wp; wtforce_z = 0.0_wp
      fturbinex = 0.0_wp; fturbiney = 0.0_wp
      flag_de(:,:,:) = 0.0_wp 
      
      allocate(pu(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(pv(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(pw(xsz(1),xsz(2),1-level:xsz(3)+level))     
      
      allocate(u_xi(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(u_psi(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(u_zeta(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(ux(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(uy(xsz(1),xsz(2),1-level:xsz(3)+level))
      !allocate(uz_w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(uz_u(xsz(1),xsz(2),1-level:xsz(3)+level))      
      allocate(vx(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(vy(xsz(1),xsz(2),1-level:xsz(3)+level))
      !allocate(vz_w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(vz_u(xsz(1),xsz(2),1-level:xsz(3)+level))      
      allocate(wx(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(wy(xsz(1),xsz(2),1-level:xsz(3)+level))
      !allocate(wz_u(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(wz_w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(tmp_x1(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(tmp_x2(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(tmp_x3(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(tmp_x4(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(tmp_x5(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(tmp_x6(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(tmp_x7(xsz(1),xsz(2),1-level:xsz(3)+level))

      allocate(tmp_y1(ysz(1), ysz(2), 1-level:ysz(3)+level))
      ! allocate(tmp_y2(ysz(1), ysz(2), 1-l:ysz(3)+l))

      allocate(tauwx(xsz(1),xsz(2)))
      allocate(tauwy(xsz(1),xsz(2))) 
      
      allocate(cc(xsz(1),xsz(2),8))

      allocate(lap_u_ao(ysz(1),ysz(2),ysz(3)))
      allocate(lap_u_al(ysz(1),ysz(2),ysz(3)))
      allocate(lap_u_ar(ysz(1),ysz(2),ysz(3)))
      allocate(lap_w_ao(ysz(1),ysz(2),ysz(3)))
      allocate(lap_w_al(ysz(1),ysz(2),ysz(3)))
      allocate(lap_w_ar(ysz(1),ysz(2),ysz(3)))
      allocate(lap_p_ao(ysz(1),ysz(2),ysz(3)))
      allocate(lap_p_al(ysz(1),ysz(2),ysz(3)))
      allocate(lap_p_ar(ysz(1),ysz(2),ysz(3)))


    end subroutine navier_init

#include "navier_les.f90"
#include "navier_fu.f90"
#include "navier_pressure.f90"
    
    
    !checked
    !input: u
    !output: ux, uy, uz_u
    subroutine gradu(u, ux, uy, uz_u)
      use grid, only : dz, zetax, zetay, her
      use param, only: pex, pey
      use spectral
      implicit none

      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN)  :: u
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(OUT) :: ux, uy, uz_u

      integer :: k
      real(wp), allocatable, dimension(:,:,:) :: tmp

      allocate(tmp(xsz(1),xsz(2),1-level:xsz(3)+level))
      !--------------------
      ! Calculate du/dzeta
      !--------------------
      call calc_uzeta(u, tmp(:,:,1:xsz(3)), level)
      if (istop)then
         tmp(:,:,xsz(3)-1) = (u(:,:,xsz(3)+1) - u(:,:,xsz(3)-2))/(2*dz(xsz(3)-2))
         tmp(:,:,xsz(3)) = (u(:,:,xsz(3)+1) - u(:,:,xsz(3)-1))/dz(xsz(3)-2)
      endif
      
      call update_ghost(tmp, level)

      !-----------------
      ! Calculate du/dx
      !-----------------
      call pdfx(u(:,:,1:xsz(3)), ux(:,:,1:xsz(3)), pex)
      call update_ghost(ux, level)
      if (istop) then
         call pdfx(u(:,:,xsz(3)+1), ux(:,:,xsz(3)+1), pex)
      end if
      ux = ux + zetax(:,:,:)*tmp(:,:,:)

      !-----------------
      ! Calculate du/dy
      !-----------------
      call pdfy_x(u(:,:,1:xsz(3)), uy(:,:,1:xsz(3)), pey)
      call update_ghost(uy, level)
      if (istop) then
         call pdfy_x(u(:,:,xsz(3)+1), uy(:,:,xsz(3)+1), pey)
      end if
      uy = uy + zetay(:,:,:)*tmp(:,:,:)

      !-----------------
      ! Calculate du/dz at u,v
      !-----------------
      do k=1-level, xsz(3)+level
         uz_u(:,:,k) = her(:,:)*tmp(:,:,k)
      end do     

      deallocate(tmp)
   end subroutine gradu

   !checked
   !input: w
   !output: wx, wy, wz_w
   subroutine gradw(w, wx, wy, wz_w)
     use grid, only :   zetax_w, zetay_w, her
     use param, only: pex, pey
     use spectral
     implicit none

      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN)  :: w
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(OUT) :: wx, wy, wz_w

      integer :: k
      real(wp), allocatable, dimension(:,:,:) :: tmp

      allocate(tmp(xsz(1),xsz(2),1-level:xsz(3)+level))
      !--------------------
      ! Calculate dw/dzeta
      !--------------------
      call calc_wzeta(w, tmp(:,:,1:xsz(3)), level)
      call update_ghost(tmp, level)
      
      !-----------------
      ! Calculate dw/dx
      !-----------------
      call pdfx(w(:,:,1:xsz(3)), wx(:,:,1:xsz(3)), pex)
      call update_ghost(wx, level)
      if (isbot) then
         wx(:,:,1:) = wx(:,:,1:) + zetax_w(:,:,1:)*tmp(:,:,1:)
      else
         wx(:,:,:) = wx(:,:,:) + zetax_w(:,:,:)*tmp(:,:,:)
      end if

      !-----------------
      ! Calculate dw/dy
      !-----------------
      call pdfy_x(w(:,:,1:xsz(3)), wy(:,:,1:xsz(3)), pey)
      call update_ghost(wy, level)
      if (isbot) then
         wy(:,:,1:) = wy(:,:,1:) + zetay_w(:,:,1:)*tmp(:,:,1:)
      else
         wy(:,:,:) = wy(:,:,:) + zetay_w(:,:,:)*tmp(:,:,:)
      end if

      !-----------------
      ! Calculate dw/dz at w
      !-----------------
      
      do k=1-level, xsz(3)+level
         wz_w(:,:,k) = her(:,:)*tmp(:,:,k)
      end do 

      deallocate(tmp)
   end subroutine gradw

   
   !unchecked
   subroutine calc_divu(divu)
     use grid, only : dz, dzw, zetax, zetay, her
     use param, only: pex, pey
     use spectral
     
     implicit none
     
     real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: divu
     
     integer :: k
     
     !---------------
     ! compute du/dx
     !---------------
     call pdfx(u(:,:,1:xsz(3)),u_xi(:,:,1:xsz(3)),pex)
     call calc_uzeta(u, u_zeta(:,:,1:), level)
     if (isbot) then
        u_zeta(:,:,1)=0
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-2)=(u(:,:,xsz(3)-1)-u(:,:,xsz(3)-3))/2/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)-1)=(4*u(:,:,xsz(3))-3*u(:,:,xsz(3)-1)-u(:,:,xsz(3)-2))/3/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)) = (9*u(:,:,xsz(3)+1)+8*u(:,:,xsz(3))  &
                             -18*u(:,:,xsz(3)-1)+u(:,:,xsz(3)-2)) &
                             /12/dz(xsz(3)-2)
      end if
      u_zeta(:,:,1:xsz(3)) = u_zeta(:,:,1:xsz(3))*zetax(:,:,1:xsz(3))
      divu = u_xi(:,:,1:xsz(3)) + u_zeta(:,:,1:xsz(3))

      !---------------
      ! compute dv/dy
      !---------------
      call pdfy_x(v(:,:,1:xsz(3)),u_psi(:,:,1:xsz(3)),pey)
      call calc_uzeta(v, u_zeta(:,:,1:), level)
      if (isbot) then
         u_zeta(:,:,1)=0
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-2)=(v(:,:,xsz(3)-1)-v(:,:,xsz(3)-3))/2/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)-1)=(4*v(:,:,xsz(3))-3*v(:,:,xsz(3)-1)-v(:,:,xsz(3)-2))/3/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)) = (9*v(:,:,xsz(3)+1)+8*v(:,:,xsz(3))  &
                             -18*v(:,:,xsz(3)-1)+v(:,:,xsz(3)-2)) &
                             /12/dz(xsz(3)-2)
      end if
      u_zeta(:,:,1:xsz(3)) = u_zeta(:,:,1:xsz(3))*zetay(:,:,1:xsz(3))
      divu = divu + u_psi(:,:,1:xsz(3)) + u_zeta(:,:,1:xsz(3))

      !---------------
      ! compute dw/dz
      !---------------
      do k=1, xsz(3)
         u_zeta(:,:,k) = (w(:,:,k)-w(:,:,k-1))/dzw(k-1)
      end do
      if (isbot) then
         u_zeta(:,:,1) = 2*w(:,:,1)/dz(1)
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-1)=(w(:,:,xsz(3)-1)-w(:,:,xsz(3)-2))/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)) = (2*w(:,:,xsz(3))+3*w(:,:,xsz(3)-1)-6*w(:,:,xsz(3)-2)+w(:,:,xsz(3)-3))/6/dz(xsz(3)-2)
      end if
      do k=1,xsz(3)
         u_zeta(:,:,k) = u_zeta(:,:,k)*her(:,:)
      end do
      divu = divu + u_zeta(:,:,1:xsz(3))

   end subroutine calc_divu

   !> added by plyu. Calculate divergence when u, v, w, are co-located at the same place.
   subroutine calc_divu_center(divu, u_, v_, w_)
     use grid, only : dz, dzw, zetax, zetay, her
     use param, only: pex, pey
     use spectral
     
     implicit none
     
     real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: divu
     real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(in) :: u_, v_, w_
     
     integer :: k
     
     !---------------
     ! compute du/dx
     !---------------
     call pdfx(u_(:,:,1:xsz(3)),u_xi(:,:,1:xsz(3)),pex)
     call calc_uzeta(u_, u_zeta(:,:,1:), level)
     if (isbot) then
        u_zeta(:,:,1)=0
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-2)=(u_(:,:,xsz(3)-1)-u_(:,:,xsz(3)-3))/2/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)-1)=(4*u_(:,:,xsz(3))-3*u_(:,:,xsz(3)-1)-u_(:,:,xsz(3)-2))/3/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)) = (9*u_(:,:,xsz(3)+1)+8*u_(:,:,xsz(3))  &
                             -18*u_(:,:,xsz(3)-1)+u_(:,:,xsz(3)-2)) &
                             /12/dz(xsz(3)-2)
      end if
      u_zeta(:,:,1:xsz(3)) = u_zeta(:,:,1:xsz(3))*zetax(:,:,1:xsz(3))
      divu = u_xi(:,:,1:xsz(3)) + u_zeta(:,:,1:xsz(3))

      !---------------
      ! compute dv/dy
      !---------------
      call pdfy_x(v_(:,:,1:xsz(3)),u_psi(:,:,1:xsz(3)),pey)
      call calc_uzeta(v_, u_zeta(:,:,1:), level)
      if (isbot) then
         u_zeta(:,:,1)=0
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-2)=(v_(:,:,xsz(3)-1)-v_(:,:,xsz(3)-3))/2/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)-1)=(4*v_(:,:,xsz(3))-3*v_(:,:,xsz(3)-1)-v_(:,:,xsz(3)-2))/3/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)) = (9*v_(:,:,xsz(3)+1)+8*v_(:,:,xsz(3))  &
                             -18*v_(:,:,xsz(3)-1)+v_(:,:,xsz(3)-2)) &
                             /12/dz(xsz(3)-2)
      end if
      u_zeta(:,:,1:xsz(3)) = u_zeta(:,:,1:xsz(3))*zetay(:,:,1:xsz(3))
      divu = divu + u_psi(:,:,1:xsz(3)) + u_zeta(:,:,1:xsz(3))

      !---------------
      ! compute dw/dz (w_in is defined at u, v location)
      !---------------
      call calc_uzeta(w_, u_zeta(:,:,1:xsz(3)), level)
      if (istop)then
         u_zeta(:,:,xsz(3)-1) = (w_(:,:,xsz(3)+1) - w_(:,:,xsz(3)-2))/(2*dz(xsz(3)-2))
         u_zeta(:,:,xsz(3)) = (w_(:,:,xsz(3)+1) - w_(:,:,xsz(3)-1))/dz(xsz(3)-2)
      endif
      do k = 1, xsz(3)
        u_zeta(:,:,k) = u_zeta(:,:,k) * her(:,:)
      enddo
      call update_ghost(u_zeta, level)

      divu = divu + u_zeta(:,:,1:xsz(3))

   end subroutine calc_divu_center

   !unchecked
   function get_max_divu(maxdiv_loc) result(maxdiv)
      use MPI
      use spectral, only : dealiasxy
      implicit none

      integer, dimension(3), intent(OUT), optional :: maxdiv_loc
      real(wp) :: maxdiv

      integer  :: max_ijk(3), ierror
      real(wp) :: tmp

      call calc_divu(tmp_x1(:,:,1:xsz(3)))
      call dealiasxy(tmp_x1(:,:,1:xsz(3)))
      tmp_x1 = abs(tmp_x1)
      max_ijk = maxloc(tmp_x1(:,:,1:xsz(3)))
      tmp = tmp_x1(max_ijk(1),max_ijk(2),max_ijk(3))
      call MPI_Allreduce(tmp, maxdiv, 1, real_type, MPI_MAX, MPI_COMM_2D_CART, ierror)

      if (present(maxdiv_loc)) then
         if (maxdiv == tmp) then
            max_ijk(1:3) = max_ijk(1:3)+xst(1:3)-1
         else
            max_ijk(1:3) = 0
         end if
         call MPI_Allreduce(max_ijk, maxdiv_loc, 3, MPI_INTEGER, MPI_SUM,  &
                            MPI_COMM_2D_CART, ierror)
      end if

    end function get_max_divu

  
   !checked
   !---------------------------------------------------------------------
   !input: time, timewavy, iwavy, nx, dx, tcoef, hk, ha, homeg, pex, pey
   !output: hh, ht, hx, hy, hxy, hxx, hyy, eta0(set zero),eta
   !output: u, v, w, ub, vb, wb
   !---------------------------------------------------------------------
   subroutine bottom2(time)
     
     use MPI
     use spectral
     use grid, only: hh, ht, hx, hy, hxy, hxx, hyy, eta0, eta
     use param
     use fft
     
     implicit none
     integer :: i,j 
     real(wp), intent(IN) :: time  
     real(wp) :: x, fsin, fcos, fexp, time0
     real(wp), dimension(ysz(1),ysz(2)) :: tmpy1, tmpy2
     
     if(time .lt. timewavy) then 
        time0 = 0
     else
        time0 = time - timewavy
     end if
     
     hh=0
     ht=0
     ub=0
     vb=0
     wb=-ht
     
     ! iwavy=1: solid wavy wall 
     if (iwavy==1)then
        do i=1,xsz(1)
           do j=1,xsz(2)
              x=(i-nx/2-1)*dx
              fsin=sin(hk*x)
              fcos=cos(hk*x)
              fexp=exp(-tcoef*time0**2)
              
              ! hh(i,j)=-ha*fsin*(1.0_wp-fexp)
              ! ht(i,j)=-ha*fsin*2.0_wp*tcoef*time0*fexp
              hh(i,j)=-ha*fsin
              ht(i,j)=0
              
              ub(i,j)=0
              vb(i,j)=0
              wb(i,j)=-ht(i,j)
           enddo
        enddo
        do i=1, xsz(1)
           write(102, *)i, hh(i,1)
        end do 
        goto 999 
   endif
   
   ! iwavy=3: water wave surface
   if(iwavy.eq.3) then
      
      do i=1,xsz(1)
         do j=1,xsz(2)
            x=(i-nx/2-1)*dx
            fsin=sin(hk*x-homeg*time0)
            fcos=cos(hk*x-homeg*time0)
            fexp=exp(-tcoef*time0**2)
            
            hh(i,j)=-ha*fsin*(1.0_wp-fexp)
            ht(i,j)=ha*homeg*fcos*(1.0_wp-fexp)-ha*fsin*2.0_wp*tcoef*time0*fexp
            
            ub(i,j)=ha*homeg*fsin*(1.0_wp-fexp)
            vb(i,j)=0
            wb(i,j)=-ht(i,j)
         enddo
      enddo
      !if(myid==0) then
      !  print *, 'time = ', time, ', time0 = ', time0
      !  print *, 'hk = ', hk, ', homeg = ', homeg
      !  print *, 'x(1:5)=', (1-nx/2-1)*dx, (2-nx/2-1)*dx, (3-nx/2-1)*dx,&
      !    (4-nx/2-1)*dx, (5-nx/2-1)*dx
      !  print *, 'hh(1:5)=', hh(1:5, 1)
      !endif      
      goto 999
   end if
   
999 continue
   
   call pdfx(hh, hx, pex)
   call pdfy_x(hh, hy, pey)
   call pdfy_x(hx, hxy, pey)

   ! compute hh_xx   
   call fft_r2c_x(hh,hxx)
   call pdfxx_(hxx,hxx,1,pex)
   call fft_c2r_x(hxx)

   ! compute hh_yy
   call transpose_xy(hh,tmpy1)
   !tmpy1 saves the eta in spectral domain in y-pencil
   call fft_r2c_y(tmpy1)
   call pdfyy_(tmpy1,tmpy2,1,pey)
   call fft_c2r_y(tmpy2)
   call transpose_yx(tmpy2,hyy)
   
   eta0=0
   eta(:,:)=hh(:,:)+0
 ! ------------------------------------------------------
 ! > give the bottom boundary condition
 ! ------------------------------------------------------
   if(isbot)then 
      u(:,:,1)=ub(:,:)
      v(:,:,1)=vb(:,:)
      w(:,:,1)=wb(:,:)
   end if 
   
 end subroutine bottom2

 ! checked
 !---------------------------------------------------------------------
 !input: etas, ubs, vbs, wbs, nxs, nys, np1
 !input: timewavy, time, tcoef, pex, pey
 !output: eta, hh, ht, hx, hy, hxy, hxx, hyy
 !output: u, v, w, ub, vb, wb
 !---------------------------------------------------------------------
 subroutine bottom_hos_les(time)
   use MPI
   use spectral
   use param
   use grid, only: eta, hh, ht, hx, hy, hxy, hxx, hyy
   use fft
   
   implicit none
   
   real(wp), intent(IN) :: time  
   real(wp) ::  fexp, time0
   real(wp), dimension(ysz(1),ysz(2)) :: tmpy1, tmpy2
   real(wp), dimension(xsz(1),xsz(2)) :: tmp1,tmp2
   
   integer :: i, j

   call filtering(nxs,nys,np1,eta_hos,eta)
   call filtering(nxs,nys,np1,u_hos,ub)
   call filtering(nxs,nys,np1,v_hos,vb)
   call filtering(nxs,nys,np1,w_hos,wb)
   
   !print *, "nxs, nys, np1=", nxs, nys, np1
   !print *, "eta_hos=", eta_hos(1:5,1) , ", eta=", eta(1:5,1)
   !do i = 1, size(eta_hos, 1)
   !  do j = 1, size(eta_hos, 2)
   !    if (isnan(eta_hos(i,j))) print *, "nan in eta_hos_in_bottom_hos_les:", i, j
   !  enddo
   !enddo
   !print *, "all eta_hos=", eta_hos
   
   call pdfx(eta, tmp1, pex)
   call pdfy_x(eta, tmp2, pey)
   
   if(time.lt.timewavy)then 
      time0=0
   else
      time0=time-timewavy
   endif

   fexp=exp(-tcoef*time0**2)

   hh = -eta * (1.0_wp-fexp)
   ht = -(wb-ub*tmp1-vb*tmp2)*(1.0_wp-fexp)
   hx = -tmp1 *(1.0_wp-fexp)
   hy = -tmp2 *(1.0_wp-fexp)
   
   call dealiasxy(ub)
   call dealiasxy(vb)
   call dealiasxy(wb)
   call dealiasxy(hh)
   call dealiasxy(ht)
   
   u(:,:,1)=ub*(1.0_wp-fexp)
   v(:,:,1)=vb*(1.0_wp-fexp)
   w(:,:,1)=wb*(1.0_wp-fexp)
   ub = u(:,:,1)
   vb = v(:,:,1)
   wb = w(:,:,1)

   call pdfy_x(hx, hxy, pey)

   ! compute hh_xx   
   call fft_r2c_x(hh,hxx)
   call pdfxx_(hxx,hxx,1,pex)
   call fft_c2r_x(hxx)

   ! compute hh_yy
   call transpose_xy(hh,tmpy1)
   !tmpy1 saves the eta in spectral domain in y-pencil
   call fft_r2c_y(tmpy1)
   call pdfyy_(tmpy1,tmpy2,1,pey)
   call fft_c2r_y(tmpy2)
   call transpose_yx(tmpy2,hyy)
      
   ! if(myid==0)print*, 'hh', hh(1:8,1)
 end subroutine bottom_hos_les

 subroutine remove_plane_mean(egg)
   use mpi
   implicit none

   real(wp), dimension(:,:,:), intent(inout) :: egg
   integer :: k
   
   !> for local averaging
   real(wp) :: vm_tmp, vmm

   !> for global averaging
   real(wp), allocatable, dimension(:) :: avg_local, avg_global
   integer :: ierror

   !!> local averaging
   !vmm = 0.0_wp
   !do k = 1, xsz(3)
   !  vm_tmp = sum(egg(:,:,k)) / xsz(1) / xsz(2)
   !  egg(:,:,k) = egg(:,:,k) - vm_tmp
   !  vmm = vmm + vm_tmp
   !enddo
   !vmm = vmm / xsz(3)
   !if (myid.eq.0) print *, 'CPU0 removed mean flow: ', vmm

   !> global averaging
   if (.not. allocated(avg_local)) then
     allocate(avg_local(nz_global), avg_global(nz_global))
   endif
   avg_local(:) = 0.0_wp; avg_global(:) = 0.0_wp
   do k = 1, xsz(3)   
     avg_local(xst(3)+k-1) = avg_local(xst(3)+k-1) + sum(egg(:,:,k))
   enddo
   !if (myid.eq.0) print *, 'local avg: ', avg_local
   
   !> plyunote: about mpi_comm_2d_col, see collect_grid in turbine_model.f90
   call MPI_Allreduce(avg_local, avg_global, nz_global, mpi_double_precision, &
     MPI_SUM, MPI_COMM_2D_COL, ierror)
   avg_local(:) = avg_global(:) / nx_global / ny_global
   !if (myid.eq.0) print *, 'global avg 1: ', avg_local
   avg_global(:) = 0.0_wp
   call MPI_Allreduce(avg_local, avg_global, nz_global, mpi_double_precision, &
     MPI_SUM, MPI_COMM_2D_ROW, ierror)
   !if (myid.eq.0) print *, 'global avg 2: ', avg_global

   do k = 1, xsz(3)
     egg(:,:,k) = egg(:,:,k) - avg_global(xst(3)+k-1)
   enddo
   
   call update_ghost(egg, level)
   if (myid.eq.0) print *, 'plyunote: remove_plane_mean is called.'

 end subroutine remove_plane_mean
   
end module navier
