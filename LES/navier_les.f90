!whole .f checked
!input & output checked
subroutine les_init(level_)
      implicit none

      integer :: level_

      level = level_

      allocate(nut(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(nutw(xsz(1),xsz(2),1-level:xsz(3)+level))

      allocate(s11(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s12(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s13(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s22(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s23(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s33(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t11(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t12(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t13(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t22(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t23(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t33(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(ss(xsz(1),xsz(2),1:xsz(3)))

      allocate(s11w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s12w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s13w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s22w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s23w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(s33w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t11w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t12w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t13w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t22w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t23w(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(t33w(xsz(1),xsz(2),1-level:xsz(3)+level))

      allocate(uf(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(vf(xsz(1),xsz(2),1-level:xsz(3)+level))
      allocate(wf(xsz(1),xsz(2),1-level:xsz(3)+level))

   end subroutine les_init

   !checked
   !input: deltax, deltay, u, v, w, uf, vf, wf
   !input: hbar, cs0, resbot, aplus, z0, ivandriest, icsc
   !input: zz
   !input: nz_global
   !output: s11~s33, ss
   !output: nut
   subroutine get_nut(u, v, w, uf, vf, wf, deltax, deltay)
      use grid, only : zz
      use spectral
      use MPI
      use param, only: hbar, z0, cs0, icsc, aplus, resbot, ivandriest
      implicit none

      real(wp), intent(IN) :: deltax, deltay
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: u, v, w
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: uf, vf, wf

      integer  :: i, j, k
      real(wp) :: dzm, delta, deltat
      real(wp), allocatable, dimension(:) :: mlr, mmr, mlp, mmp, csp
      real(wp), allocatable, dimension(:,:,:) :: m11, m12, m13, m22, m23, m33
      real(wp), allocatable, dimension(:,:,:) :: l11, l12, l13, l22, l23, l33

      allocate(mlr(xsz(3)), mmr(xsz(3)), mlp(xsz(3)), mmp(xsz(3)), csp(xsz(3)))
      allocate(m11(xsz(1),xsz(2),xsz(3)), m12(xsz(1),xsz(2),xsz(3)))
      allocate(m13(xsz(1),xsz(2),xsz(3)), m22(xsz(1),xsz(2),xsz(3)))
      allocate(m23(xsz(1),xsz(2),xsz(3)), m33(xsz(1),xsz(2),xsz(3)))
      allocate(l11(xsz(1),xsz(2),xsz(3)), l12(xsz(1),xsz(2),xsz(3)))
      allocate(l13(xsz(1),xsz(2),xsz(3)), l22(xsz(1),xsz(2),xsz(3)))
      allocate(l23(xsz(1),xsz(2),xsz(3)), l33(xsz(1),xsz(2),xsz(3)))

      dzm = hbar/real(nz_global-2, wp)
      delta = (deltax*deltay*dzm)**(1.0_wp/3)
      deltat = (4*deltax*deltay*dzm)**(1.0_wp/3)

      call get_strain_mag(s11(:,:,1:), s12(:,:,1:), s13(:,:,1:), s22(:,:,1:),  & 
                          s23(:,:,1:), s33(:,:,1:), ss(:,:,1:))

      call get_Lij(u, v, w, uf, vf, wf, &
                   l11, l12, l13, l22, l23, l33, &
                   mfilt, deltax, deltay)

      call get_Mij(s11(:,:,1:), s12(:,:,1:), s13(:,:,1:), s22(:,:,1:), &
                   s23(:,:,1:), s33(:,:,1:), ss(:,:,1:), &
                   m11(:,:,1:), m12(:,:,1:), m13(:,:,1:),&
                   m22(:,:,1:), m23(:,:,1:), m33(:,:,1:), &
                   mfilt, deltax, deltay, delta, deltat)

      mlp = 0; mmp = 0
      do k=1, xsz(3)
         do j=1, xsz(2)
            do i=1, xsz(1)
               mlp(k) = mlp(k)+m11(i,j,k)*l11(i,j,k)+m22(i,j,k)*l22(i,j,k)+m33(i,j,k)*l33(i,j,k) &
                   +2*m12(i,j,k)*l12(i,j,k)+2*m13(i,j,k)*l13(i,j,k)+2*m23(i,j,k)*l23(i,j,k)
               mmp(k) = mmp(k)+m11(i,j,k)*m11(i,j,k)+m22(i,j,k)*m22(i,j,k)+m33(i,j,k)*m33(i,j,k) &
                   +2*m12(i,j,k)*m12(i,j,k)+2*m13(i,j,k)*m13(i,j,k)+2*m23(i,j,k)*m23(i,j,k)
            end do
         end do
      end do
      call MPI_Allreduce(mlp, mlr, xsz(3), real_type, MPI_SUM, MPI_COMM_2D_COL, i)
      call MPI_Allreduce(mmp, mmr, xsz(3), real_type, MPI_SUM, MPI_COMM_2D_COL, i)

      csp = mlr/mmr

      if (ivandriest .eq. 1)then
         do k=1, xsz(3)
            csp(k) = (cs0 * (1.0_wp-exp(-zz(k)*resbot/0.5/aplus)))**2.0_wp
         enddo
      elseif(ivandriest .eq. 2)then
         do k=1, xsz(3)
            csp(k) = ((cs0 * delta)**(-3.0_wp) + (0.4 * (zz(k) * hbar + z0)) &
                 ** (-3.0_wp)) ** (-2./3.) /delta **2.0_wp
         enddo
      else
         continue
      endif

      if(icsc .eq. 1)then
         do k =1, xsz(3)
            if (csp(k) .lt. 0)then
               csp(k) = ((cs0 * delta) ** (-3.0_wp) + (0.4_wp*(zz(k) * hbar +z0))&
                    ** (-3.0_wp)) ** (-2.0_wp/3.0_wp) / delta**2.0_wp
            endif
         enddo
      endif

      do k = 1, xsz(3)
         do j =1, xsz(2)
            do i =1, xsz(1)
               nut(i,j,k) = csp(k)* delta**2 * ss(i,j,k)
            end do
         end do
      end do

      call dealiasxy(nut(:,:,1:xsz(3)))

      call update_ghost(nut, level)
      
      !boundary condition
      if (istop) then
         nut(:,:,xsz(3)+1) = nut(:,:,xsz(3)-1)
      end if
      
      if (isbot) then
         nut(:,:,1) = 0
      end if
      
      ! print*, 'nut', nut(1,1,1:xsz(3))
      deallocate(mlr, mmr, mlp, mmp, csp)
      deallocate(m11, m12, m13, m22, m23, m33)
      deallocate(l11, l12, l13, l22, l23, l33)

   end subroutine get_nut

   !checked
   !input: deltax, deltay, u, v, w, uf, vf, wf
   !input: hbar, cs0, resbot, aplus, z0, ivandriest, icsc
   !input: zw
   !input: nz_global
   !output: s11w~s33w, ss
   !output: nutw
   subroutine get_nutw(u, v, w, uf, vf, wf, deltax, deltay)
      use grid, only : zw
      use spectral
      use MPI
      use param, only: hbar, z0, cs0, icsc, aplus, resbot, ivandriest
      implicit none

      real(wp), intent(IN) :: deltax, deltay
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: u, v, w
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: uf, vf, wf

      integer  :: i, j, k
      real(wp) :: dzm, delta, deltat

      real(wp), allocatable, dimension(:) :: mlr, mmr, mlp, mmp, csp
      real(wp), allocatable, dimension(:,:,:) :: m11, m12, m13, m22, m23, m33
      real(wp), allocatable, dimension(:,:,:) :: l11, l12, l13, l22, l23, l33

      allocate(mlr(xsz(3)), mmr(xsz(3)), mlp(xsz(3)), mmp(xsz(3)), csp(xsz(3)))
      allocate(m11(xsz(1),xsz(2),xsz(3)), m12(xsz(1),xsz(2),xsz(3)))
      allocate(m13(xsz(1),xsz(2),xsz(3)), m22(xsz(1),xsz(2),xsz(3)))
      allocate(m23(xsz(1),xsz(2),xsz(3)), m33(xsz(1),xsz(2),xsz(3)))
      allocate(l11(xsz(1),xsz(2),xsz(3)), l12(xsz(1),xsz(2),xsz(3)))
      allocate(l13(xsz(1),xsz(2),xsz(3)), l22(xsz(1),xsz(2),xsz(3)))
      allocate(l23(xsz(1),xsz(2),xsz(3)), l33(xsz(1),xsz(2),xsz(3)))


      dzm = hbar/real(nz_global-2, wp)
      delta = (deltax*deltay*dzm)**(1.0_wp/3)
      deltat = (4*deltax*deltay*dzm)**(1.0_wp/3)

      call get_strain_mag(s11w(:,:,1:xsz(3)), s12w(:,:,1:xsz(3)), s13w(:,:,1:xsz(3)),&
           s22w(:,:,1:xsz(3)), s23w(:,:,1:xsz(3)), s33w(:,:,1:xsz(3)), ss)

      call get_Lij_w(u, v, w, uf, vf, wf, &
                     l11, l12, l13, l22, l23, l33, &
                     mfilt, deltax, deltay)

      call get_Mij(s11w(:,:,1:), s12w(:,:,1:), s13w(:,:,1:), & 
                   s22w(:,:,1:), s23w(:,:,1:), s33w(:,:,1:), ss, &
                   m11(:,:,1:), m12(:,:,1:), m13(:,:,1:), & 
                   m22(:,:,1:), m23(:,:,1:), m33(:,:,1:), &
                   mfilt, deltax, deltay, delta, deltat)

      ! if (myid==0)print*, m13(1,1,1:xsz(3))
      ! if (myid==1)print*, m13(1,1,1:xsz(3))

      mlp = 0; mmp = 0
      do k=1, xsz(3)
         do j=1, xsz(2)
            do i=1, xsz(1)
               mlp(k) = mlp(k)+m11(i,j,k)*l11(i,j,k)+m22(i,j,k)*l22(i,j,k)+m33(i,j,k)*l33(i,j,k) &
                   +2*m12(i,j,k)*l12(i,j,k)+2*m13(i,j,k)*l13(i,j,k)+2*m23(i,j,k)*l23(i,j,k)
               mmp(k) = mmp(k)+m11(i,j,k)*m11(i,j,k)+m22(i,j,k)*m22(i,j,k)+m33(i,j,k)*m33(i,j,k) &
                   +2*m12(i,j,k)*m12(i,j,k)+2*m13(i,j,k)*m13(i,j,k)+2*m23(i,j,k)*m23(i,j,k)
            end do
         end do
      end do
      call MPI_Allreduce(mlp, mlr, xsz(3), real_type, MPI_SUM, MPI_COMM_2D_COL, i)
      call MPI_Allreduce(mmp, mmr, xsz(3), real_type, MPI_SUM, MPI_COMM_2D_COL, i)

      ! print*, mlr(1:xsz(3))

      csp = mlr/mmr

      if (ivandriest .eq. 1)then
         do k=1, xsz(3)
            csp(k) = (cs0 * (1.0_wp-exp(-zw(k)*resbot/0.5/aplus)))**2.0_wp
         enddo
      elseif(ivandriest .eq. 2)then
         do k=1, xsz(3)
            csp(k) = ((cs0 * delta)**(-3.0_wp) + (0.4 * (zw(k) * hbar + z0)) &
                 ** (-3.0_wp)) ** (-2.0_wp/3.0_wp) /delta **2.0_wp
         enddo
      else
         continue
      endif

      if(icsc .eq. 1)then
         do k =1, xsz(3)
            if (csp(k) .lt. 0)then
               csp(k) = ((cs0 * delta) ** (-3.0_wp) + (0.4_wp*(zw(k) * hbar +z0))&
                    ** (-3.0_wp)) ** (-2.0_wp/3.0_wp) / delta**2.0_wp
            endif
         enddo
      endif

      ! print*, 'csp', csp(1:xsz(3))

      csp = 0.16

      do k = 1, xsz(3)
         do j =1, xsz(2)
            do i =1, xsz(1)
               nutw(i,j,k) = csp(k)* delta**2 * ss(i,j,k)
            end do
         end do
      end do

      call dealiasxy(nutw(:,:,1:xsz(3)))
      call update_ghost(nutw, level)

      !boundary condition
      if (istop) then
         nutw(:,:,xsz(3)) = nutw(:,:,xsz(3)-2) 
      end if

      deallocate(mlr, mmr, mlp, mmp, csp)
      deallocate(m11, m12, m13, m22, m23, m33)
      deallocate(l11, l12, l13, l22, l23, l33)

      ! print*, nutw(1,1,1:xsz(3))
   end subroutine get_nutw

   !checked
   subroutine get_Mij(s11, s12, s13, s22, s23, s33, ss, &
                      m11, m12, m13, m22, m23, m33, &
                      mfilt, deltax, deltay, delta, deltat)
      use spectral
      implicit none

      integer,  intent(IN) :: mfilt
      real(wp), intent(IN) :: deltax, deltay, delta, deltat
      real(wp), dimension(xsz(1),xsz(2),1:xsz(3)), intent(IN) :: s11, s12, s13, s22, s23, s33, ss
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: m11, m12, m13, m22, m23, m33

      real(wp), allocatable, dimension(:,:,:) :: ss11f, ss12f, ss13f, ss22f, ss23f, ss33f, ssf

      integer :: i, j, k

      allocate(ss11f(xsz(1),xsz(2),xsz(3)))
      allocate(ss12f(xsz(1),xsz(2),xsz(3)))
      allocate(ss13f(xsz(1),xsz(2),xsz(3)))
      allocate(ss22f(xsz(1),xsz(2),xsz(3)))
      allocate(ss23f(xsz(1),xsz(2),xsz(3)))
      allocate(ss33f(xsz(1),xsz(2),xsz(3)))
      allocate(ssf(xsz(1),xsz(2),xsz(3)))

      ss11f(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) * s11(:,:,1:xsz(3))
      ss12f(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) * s12(:,:,1:xsz(3))
      ss13f(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) * s13(:,:,1:xsz(3))
      ss22f(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) * s22(:,:,1:xsz(3))
      ss23f(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) * s23(:,:,1:xsz(3))
      ss33f(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) * s33(:,:,1:xsz(3))

      call dealiasxy(ss11f(:,:,1:))
      call dealiasxy(ss12f(:,:,1:))
      call dealiasxy(ss13f(:,:,1:))
      call dealiasxy(ss22f(:,:,1:))
      call dealiasxy(ss23f(:,:,1:))
      call dealiasxy(ss33f(:,:,1:))

      call les_filter(ss11f, mfilt, deltax, deltay)
      call les_filter(ss12f, mfilt, deltax, deltay)
      call les_filter(ss13f, mfilt, deltax, deltay)
      call les_filter(ss22f, mfilt, deltax, deltay)
      call les_filter(ss23f, mfilt, deltax, deltay)
      call les_filter(ss33f, mfilt, deltax, deltay)

      m11(:,:,1:xsz(3))=s11(:,:,1:xsz(3))
      m12(:,:,1:xsz(3))=s12(:,:,1:xsz(3))
      m13(:,:,1:xsz(3))=s13(:,:,1:xsz(3))
      m22(:,:,1:xsz(3))=s22(:,:,1:xsz(3))
      m23(:,:,1:xsz(3))=s23(:,:,1:xsz(3))
      m33(:,:,1:xsz(3))=s33(:,:,1:xsz(3))

      call les_filter(m11, mfilt, deltax, deltay)
      call les_filter(m12, mfilt, deltax, deltay)
      call les_filter(m13, mfilt, deltax, deltay)
      call les_filter(m22, mfilt, deltax, deltay)
      call les_filter(m23, mfilt, deltax, deltay)
      call les_filter(m33, mfilt, deltax, deltay)

      call get_strain_mag(m11(:,:,1:), m12(:,:,1:), m13(:,:,1:), m22(:,:,1:), & 
                          m23(:,:,1:), m33(:,:,1:), ssf)

      do k=1, xsz(3)
         do j=1, xsz(2)
            do i=1, xsz(1)
               m11(i,j,k) = 2*delta**2 * ss11f(i,j,k) - 2*deltat**2 * ssf(i,j,k) * m11(i,j,k)
               m12(i,j,k) = 2*delta**2 * ss12f(i,j,k) - 2*deltat**2 * ssf(i,j,k) * m12(i,j,k)
               m13(i,j,k) = 2*delta**2 * ss13f(i,j,k) - 2*deltat**2 * ssf(i,j,k) * m13(i,j,k)
               m22(i,j,k) = 2*delta**2 * ss22f(i,j,k) - 2*deltat**2 * ssf(i,j,k) * m22(i,j,k)
               m23(i,j,k) = 2*delta**2 * ss23f(i,j,k) - 2*deltat**2 * ssf(i,j,k) * m23(i,j,k)
               m33(i,j,k) = 2*delta**2 * ss33f(i,j,k) - 2*deltat**2 * ssf(i,j,k) * m33(i,j,k)
            end do
         end do
      end do

      call dealiasxy(m11)
      call dealiasxy(m12)
      call dealiasxy(m13)
      call dealiasxy(m22)
      call dealiasxy(m23)
      call dealiasxy(m33)

      deallocate(ss11f, ss12f, ss13f, ss22f, ss23f, ss33f, ssf)

   end subroutine get_Mij

   !checked
   subroutine get_Lij(u, v, w, uf, vf, wf, &
                      l11, l12, l13, l22, l23, l33, &
                      mfilt, deltax, deltay)
      use spectral
      implicit none

      integer,  intent(IN) :: mfilt
      real(wp), intent(IN) :: deltax, deltay
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: u, v, w
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: uf, vf, wf
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: l11, l12, l13, l22, l23, l33

      !real(wp), dimension(xsz(1),xsz(2),xsz(3)) :: u13, u23, u33
      real(wp), allocatable, dimension(:,:,:) :: w0
      integer :: k

      allocate(w0(xsz(1),xsz(2),xsz(3)))

      do k=1, xsz(3)
         w0(:,:,k) = (w(:,:,k) + w(:,:,k-1)) / 2
      end do
      if (isbot) w0(:,:,1) = 0
      if (istop) w0(:,:,xsz(3)) = w(:,:,xsz(3)-1)
      l11(:,:,1:xsz(3)) = u(:,:,1:xsz(3))**2
      l12(:,:,1:xsz(3)) = u(:,:,1:xsz(3)) * v(:,:,1:xsz(3))
      l22(:,:,1:xsz(3)) = v(:,:,1:xsz(3))**2
      l13(:,:,1:xsz(3)) = u(:,:,1:xsz(3)) * w0(:,:,1:xsz(3))
      l23(:,:,1:xsz(3)) = v(:,:,1:xsz(3)) * w0(:,:,1:xsz(3))
      l33(:,:,1:xsz(3)) = w0(:,:,1:xsz(3))**2
      call les_filter(l11, mfilt, deltax, deltay)
      call les_filter(l12, mfilt, deltax, deltay)
      call les_filter(l13, mfilt, deltax, deltay)
      call les_filter(l22, mfilt, deltax, deltay)
      call les_filter(l23, mfilt, deltax, deltay)
      call les_filter(l33, mfilt, deltax, deltay)

      do k=1, xsz(3)
         w0(:,:,k) = (wf(:,:,k) + wf(:,:,k-1)) / 2
      end do
      if (isbot) w0(:,:,1) = 0
      if (istop) w0(:,:,xsz(3)) = wf(:,:,xsz(3)-1)
      l11(:,:,1:xsz(3)) = l11(:,:,1:xsz(3)) - uf(:,:,1:xsz(3))**2
      l12(:,:,1:xsz(3)) = l12(:,:,1:xsz(3)) - uf(:,:,1:xsz(3))*vf(:,:,1:xsz(3))
      l22(:,:,1:xsz(3)) = l22(:,:,1:xsz(3)) - vf(:,:,1:xsz(3))**2
      l13(:,:,1:xsz(3)) = l13(:,:,1:xsz(3)) - uf(:,:,1:xsz(3))*w0(:,:,1:xsz(3))
      l23(:,:,1:xsz(3)) = l23(:,:,1:xsz(3)) - vf(:,:,1:xsz(3))*w0(:,:,1:xsz(3))
      l33(:,:,1:xsz(3)) = l33(:,:,1:xsz(3)) - w0(:,:,1:xsz(3))**2

      call dealiasxy(l11)
      call dealiasxy(l12)
      call dealiasxy(l13)
      call dealiasxy(l22)
      call dealiasxy(l23)
      call dealiasxy(l33)

      ! print*, 'L23',L23(1,1,2:xsz(3))

      deallocate(w0)
   end subroutine get_Lij

!checked   
   subroutine get_Lij_w(u, v, w, uf, vf, wf, &
                        l11, l12, l13, l22, l23, l33, &
                        mfilt, deltax, deltay)
      use spectral
      implicit none

      integer,  intent(IN) :: mfilt
      real(wp), intent(IN) :: deltax, deltay
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: u, v, w
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: uf, vf, wf
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: l11, l12, l13, l22, l23, l33

      real(wp), allocatable, dimension(:,:,:) :: u0, v0
      integer :: k

      allocate(u0(xsz(1),xsz(2),xsz(3)))
      allocate(v0(xsz(1),xsz(2),xsz(3)))

      do k=1, xsz(3)
         u0(:,:,k) = (u(:,:,k) + u(:,:,k+1)) / 2.0_wp
         v0(:,:,k) = (v(:,:,k) + v(:,:,k+1)) / 2.0_wp
      end do
      if (isbot) then
         u0(:,:,1) = u(:,:,1)
         v0(:,:,1) = v(:,:,1)
      end if
      if (istop) then
         u0(:,:,xsz(3)-1) = u(:,:,xsz(3))
         v0(:,:,xsz(3)-1) = v(:,:,xsz(3))
      end if
      l11(:,:,1:xsz(3)) = u0(:,:,1:xsz(3))**2
      l12(:,:,1:xsz(3)) = u0(:,:,1:xsz(3))*v0(:,:,1:xsz(3))
      l13(:,:,1:xsz(3)) = u0(:,:,1:xsz(3))*w(:,:,1:xsz(3))
      l22(:,:,1:xsz(3)) = v0(:,:,1:xsz(3))**2
      l23(:,:,1:xsz(3)) = v0(:,:,1:xsz(3))*w(:,:,1:xsz(3))
      l33(:,:,1:xsz(3)) = w(:,:,1:xsz(3))**2
      call les_filter(l11, mfilt, deltax, deltay)
      call les_filter(l12, mfilt, deltax, deltay)
      call les_filter(l13, mfilt, deltax, deltay)
      call les_filter(l22, mfilt, deltax, deltay)
      call les_filter(l23, mfilt, deltax, deltay)
      call les_filter(l33, mfilt, deltax, deltay)

      do k=1, xsz(3)
         u0(:,:,k) = (uf(:,:,k) + uf(:,:,k+1)) / 2.0_wp
         v0(:,:,k) = (vf(:,:,k) + vf(:,:,k+1)) / 2.0_wp
      end do
      if (isbot) then
         u0(:,:,1) = uf(:,:,1)
         v0(:,:,1) = vf(:,:,1)
      end if
      if (istop) then
         u0(:,:,xsz(3)-1) = uf(:,:,xsz(3))
         v0(:,:,xsz(3)-1) = vf(:,:,xsz(3))
      end if
      l11(:,:,1:xsz(3)) = l11(:,:,1:xsz(3)) - u0(:,:,1:xsz(3))**2
      l12(:,:,1:xsz(3)) = l12(:,:,1:xsz(3)) - u0(:,:,1:xsz(3))*v0(:,:,1:xsz(3))
      l13(:,:,1:xsz(3)) = l13(:,:,1:xsz(3)) - u0(:,:,1:xsz(3))*wf(:,:,1:xsz(3))
      l22(:,:,1:xsz(3)) = l22(:,:,1:xsz(3)) - v0(:,:,1:xsz(3))**2
      l23(:,:,1:xsz(3)) = l23(:,:,1:xsz(3)) - v0(:,:,1:xsz(3))*wf(:,:,1:xsz(3))
      l33(:,:,1:xsz(3)) = l33(:,:,1:xsz(3)) - wf(:,:,1:xsz(3))**2

      call dealiasxy(l11)
      call dealiasxy(l12)
      call dealiasxy(l13)
      call dealiasxy(l22)
      call dealiasxy(l23)
      call dealiasxy(l33)

      ! print*, 'l13', l13(1,1,1:xsz(3))

      deallocate(u0,v0)

   end subroutine get_Lij_w

   !checked
   !input:u, v, w
   !output: s11~s33, s11w~s33w
   !output: u_xi, u_psi, u_zeta, tmp_x5, tmp_x6 ALL ZERO!
   subroutine get_strain(u, v, w)
      use grid, only : dz,  zetax, zetay, zetax_w, zetay_w, her
      use spectral
      use param, only: pex, pey

      implicit none

      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: u, v, w
      integer :: k

      s11 = 0; s12 = 0; s13 = 0; s22 = 0; s23 = 0; s33 = 0
      s11w = 0; s12w = 0; s13w = 0; s22w = 0; s23w = 0; s33w = 0

      ! compute du/dxi
      call pdfx(u(:,:,1:xsz(3)), u_xi(:,:,1:xsz(3)), pex)
      call update_ghost(u_xi,level)
      if (istop) then
         call pdfx(u(:,:,xsz(3)+1),u_xi(:,:,xsz(3)+1),pex)
      end if

      ! compute du/dpsi
      call pdfy_x(u(:,:,1:xsz(3)), u_psi(:,:,1:xsz(3)), pey)
      call update_ghost(u_psi,level)
      if (istop) then
         call pdfy_x(u(:,:,xsz(3)+1),u_psi(:,:,xsz(3)+1),pey)
      end if

      ! compute du/dzeta
      call calc_uzeta(u, u_zeta(:,:,1:), level)
      call update_ghost(u_zeta,level)
      if (isbot) then
         u_zeta(:,:,2)=(-u(:,:,4) + 4*u(:,:,3) - 3*u(:,:,2))&
              /(2*dz(2))
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-1)=(u(:,:,xsz(3)+1)-u(:,:,xsz(3)-2))/(2*dz(xsz(3)-2))
         u_zeta(:,:,xsz(3))=(u(:,:,xsz(3)+1)-u(:,:,xsz(3)-1))/dz(xsz(3)-2)
      end if

      !> added by plyu, to debug the abnormal NaN when calling s11 outside navier
      !if (isnan(u_xi(1,1,1))) stop 'u_xi is NaN'
      !if (isnan(zetax(1,1,1))) stop 'zetax is NaN'
      !if (isnan(u_zeta(1,1,1))) stop 'u_zeta is NaN'
      !if (isnan(s11(1,1,1))) stop 's11 is NaN'
      
      s11 = u_xi + zetax*u_zeta
      s12 = u_psi + zetay*u_zeta
      do k=1, xsz(3)
         s13(:,:,k) = her(:,:)*u_zeta(:,:,k)
      end do
      ! compute du/dzeta at w points
      do k=1, xsz(3)
         u_zeta(:,:,k) = (u(:,:,k+1) - u(:,:,k))/dz(k)
      end do

      !if(isbot) then
      !   print*, "u_xi=",u_xi(1:5,1,1),",zetax=",zetax(1:5,1,1),&
      !     ",u_zeta=",u_zeta(1:5,1,1), ",s11=", s11(1:5,1,1)
      ! endif

      if (istop) then
         u_zeta(:,:,xsz(3)-1) = (u(:,:,xsz(3)+1) - u(:,:,xsz(3)-1))/dz(xsz(3)-2)
      end if
      ! compute du/dx, du/dy, du/dz at w points
      do k=1, xsz(3)
         tmp_x5(:,:,k) = (u_xi(:,:,k) + u_xi(:,:,k+1))/2
         tmp_x6(:,:,k) = (u_psi(:,:,k) + u_psi(:,:,k+1))/2
      end do

      if(istop)then
         tmp_x5(:,:,xsz(3)-1) = u_xi(:,:,xsz(3))
         tmp_x6(:,:,xsz(3)-1) = u_psi(:,:,xsz(3))
      endif

      do k=1, xsz(3)
         s11w(:,:,k) = tmp_x5(:,:,k) + zetax_w(:,:,k)*u_zeta(:,:,k)
         s12w(:,:,k) = tmp_x6(:,:,k) + zetay_w(:,:,k)*u_zeta(:,:,k)
         s13w(:,:,k) = her(:,:) * u_zeta(:,:,k)
      end do

      u_xi=0
      u_psi=0
      u_zeta=0
      tmp_x5=0
      tmp_x6=0
      
!!$      if (istop) then
!!$         s11w(:,:,xsz(3)-1) = s11(:,:,xsz(3))
!!$         s12w(:,:,xsz(3)-1) = s12(:,:,xsz(3))
!!$      end if

      ! compute dv/dxi
      call pdfx(v(:,:,1:xsz(3)), u_xi(:,:,1:xsz(3)), pex)
      call update_ghost(u_xi,level)
      if (istop) then
         call pdfx(v(:,:,xsz(3)+1),u_xi(:,:,xsz(3)+1),pex)
      end if

      ! compute dv/dpsi
      call pdfy_x(v(:,:,1:xsz(3)), u_psi(:,:,1:xsz(3)), pey)
      call update_ghost(u_psi,level)
      if (istop) then
         call pdfy_x(v(:,:,xsz(3)+1),u_psi(:,:,xsz(3)+1),pey)
      end if

      ! compute dv/dzeta
      call calc_uzeta(v, u_zeta(:,:,1:), level)
      call update_ghost(u_zeta,level)
      
      if (isbot) then
         u_zeta(:,:,2)=(-v(:,:,4) + 4*v(:,:,3) - 3*v(:,:,2))&
              /(2*dz(2))
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-1)=(v(:,:,xsz(3)+1)-v(:,:,xsz(3)-2))/(2*dz(xsz(3)-2))
         u_zeta(:,:,xsz(3))=(v(:,:,xsz(3)+1)-v(:,:,xsz(3)-1))/dz(xsz(3)-2)
      end if

      s12 = s12 + u_xi + zetax*u_zeta
      s22 = u_psi + zetay*u_zeta
      do k=1, xsz(3)
         s23(:,:,k) = her(:,:)*u_zeta(:,:,k)
      end do
      ! compute dv/dzeta at w points
      do k=1, xsz(3)
         u_zeta(:,:,k) = (v(:,:,k+1) - v(:,:,k))/dz(k)
      end do
      if (istop) then
         u_zeta(:,:,xsz(3)-1) = (v(:,:,xsz(3)+1) - v(:,:,xsz(3)-1))/dz(xsz(3)-2)
      end if
      ! compute dv/dx, dv/dy, dv/dz at w points
      do k=1, xsz(3)
         tmp_x5(:,:,k) = (u_xi(:,:,k) + u_xi(:,:,k+1))/2
         tmp_x6(:,:,k) = (u_psi(:,:,k) + u_psi(:,:,k+1))/2
      end do
      
      if(istop)then
         tmp_x5(:,:,xsz(3)-1) = u_xi(:,:,xsz(3))
         tmp_x6(:,:,xsz(3)-1) = u_psi(:,:,xsz(3))
      endif
      
      do k=1, xsz(3)
         s12w(:,:,k) = s12w(:,:,k) + tmp_x5(:,:,k) + zetax_w(:,:,k)*u_zeta(:,:,k)
         s22w(:,:,k) = tmp_x6(:,:,k) + zetay_w(:,:,k)*u_zeta(:,:,k)
         s23w(:,:,k) = her(:,:) * u_zeta(:,:,k)
      end do
      
!!$      if (istop) then
!!$         s12w(:,:,xsz(3)-1) = s12(:,:,xsz(3))
!!$         s22w(:,:,xsz(3)-1) = s22(:,:,xsz(3))
!!$      end if

      u_xi=0
      u_psi=0
      u_zeta=0
      tmp_x5=0
      tmp_x6=0
      
      ! compute dw/dxi
      call pdfx(w(:,:,1:xsz(3)), u_xi(:,:,1:xsz(3)), pex)
      call update_ghost(u_xi, level)

      ! compute dw/dpsi
      call pdfy_x(w(:,:,1:xsz(3)), u_psi(:,:,1:xsz(3)), pey)
      call update_ghost(u_psi, level)

      ! compute dw/dzeta at w points
!!$      call calc_wzeta(w, u_zeta(:,:,1:), level)

      do k=1, xsz(3)
         u_zeta(:,:,k) = (w(:,:,k+1)-w(:,:,k-1))/(dz(k+1)/2.+dz(k)+dz(k-1)/2.)
      end do
      
      call update_ghost(u_zeta,level)
      
      if(isbot)then
         u_zeta(:,:,2) = (w(:,:,3)-w(:,:,1))/(dz(3)/2.+dz(2)+dz(1))
      endif
      
      if(istop)then
         u_zeta(:,:,xsz(3)-1) = (w(:,:,xsz(3))-w(:,:,xsz(3)-2))/(2.*dz(xsz(3)-2))
         u_zeta(:,:,xsz(3)-2) = (w(:,:,xsz(3)-1)-w(:,:,xsz(3)-3))&
              / (dz(xsz(3)-1) + dz(xsz(3)-2) + dz(xsz(3)-3)/2.)
      endif

      s13w = s13w + u_xi + zetax_w*u_zeta
      s23w = s23w + u_psi + zetay_w*u_zeta
      do k=1, xsz(3)
         s33w(:,:,k) = her(:,:)*u_zeta(:,:,k)
      end do

      ! compute dw/dzeta at u points
      do k=1, xsz(3)
         u_zeta(:,:,k) = (w(:,:,k) - w(:,:,k-1))/(dz(k)/2+dz(k-1)/2)
      end do
      if (isbot) then
         u_zeta(:,:,2) = (w(:,:,2) - w(:,:,1))/(dz(2)/2+dz(1))
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-1) = (w(:,:,xsz(3)-1)-w(:,:,xsz(3)-2))/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)) = (w(:,:,xsz(3))-w(:,:,xsz(3)-2))/2/dz(xsz(3)-2)
      end if

      do k=1, xsz(3)
         tmp_x5(:,:,k) = (u_xi(:,:,k-1) + u_xi(:,:,k)) /2
         tmp_x6(:,:,k) = (u_psi(:,:,k-1)+ u_psi(:,:,k)) /2
      end do

      if(istop)then
         tmp_x5(:,:,xsz(3)) = u_xi(:,:,xsz(3)-1)
         tmp_x6(:,:,xsz(3)) = u_psi(:,:,xsz(3)-1)
      endif

      do k=1, xsz(3)
         s13(:,:,k) = s13(:,:,k) + tmp_x5(:,:,k) + zetax(:,:,k)*u_zeta(:,:,k)
         s23(:,:,k) = s23(:,:,k) + tmp_x6(:,:,k) + zetay(:,:,k)*u_zeta(:,:,k)
         s33(:,:,k) = her(:,:)*u_zeta(:,:,k)
      end do

      u_xi=0
      u_psi=0
      u_zeta=0
      tmp_x5=0
      tmp_x6=0
      
      s12 = s12/2; s12w = s12w/2
      s13 = s13/2; s13w = s13w/2
      s23 = s23/2; s23w = s23w/2
      
      call dealiasxy(s11(:,:,1:xsz(3)))
      call dealiasxy(s12(:,:,1:xsz(3)))
      call dealiasxy(s13(:,:,1:xsz(3)))
      call dealiasxy(s22(:,:,1:xsz(3)))
      call dealiasxy(s23(:,:,1:xsz(3)))
      call dealiasxy(s33(:,:,1:xsz(3)))
      call dealiasxy(s11w(:,:,1:xsz(3)))
      call dealiasxy(s12w(:,:,1:xsz(3)))
      call dealiasxy(s13w(:,:,1:xsz(3)))
      call dealiasxy(s22w(:,:,1:xsz(3)))
      call dealiasxy(s23w(:,:,1:xsz(3)))
      call dealiasxy(s33w(:,:,1:xsz(3)))

      ! print*, 's23', s23(1,1,1:xsz(3))

      ! apply boundary condition

      if(istop)then 
         s11(:,:,xsz(3)+1)=s11(:,:,xsz(3)-1)
         s12(:,:,xsz(3)+1)=s12(:,:,xsz(3)-1)
         s13(:,:,xsz(3)+1)=-s13(:,:,xsz(3)-1)
         s22(:,:,xsz(3)+1)=s22(:,:,xsz(3)-1)
         s23(:,:,xsz(3)+1)=-s23(:,:,xsz(3)-1)
         s33(:,:,xsz(3)+1)=s33(:,:,xsz(3)-1)

         s11w(:,:,xsz(3))=s11w(:,:,xsz(3)-2)
         s12w(:,:,xsz(3))=s12w(:,:,xsz(3)-2)
         s13w(:,:,xsz(3))=-s13w(:,:,xsz(3)-2)
         s22w(:,:,xsz(3))=s22w(:,:,xsz(3)-2)
         s23w(:,:,xsz(3))=-s23w(:,:,xsz(3)-2)
         s33w(:,:,xsz(3))=s33w(:,:,xsz(3)-2)  
      end if 


   end subroutine get_strain
!checked
   subroutine get_strain_mag(s11, s12, s13, s22, s23, s33, ss)
      use spectral
      implicit none

      real(wp), dimension(:,:,:), intent(IN)  :: s11, s12, s13, s22, s23, s33
      real(wp), dimension(:,:,:), intent(OUT) :: ss
      
      ss(:,:,1:xsz(3)) = s11(:,:,1:xsz(3))**2
      ss(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) + s22(:,:,1:xsz(3))**2
      ss(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) + s33(:,:,1:xsz(3))**2
      ss(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) + 2 * s12(:,:,1:xsz(3))**2
      ss(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) + 2 * s13(:,:,1:xsz(3))**2
      ss(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) + 2 * s23(:,:,1:xsz(3))**2
      ss(:,:,1:xsz(3)) = ss(:,:,1:xsz(3)) * 2
      
      ss = sqrt(ss)
      call dealiasxy(ss(:,:,1:xsz(3)))

   end subroutine get_strain_mag

   !checked
   !input: s11~s33, s11w~s33w, nut, nutw
   !output: t11~t33, t11w~t33w
   subroutine get_SGS_stress
      use spectral
      implicit none

      t11(:,:,1:xsz(3))  = -2 * nut(:,:,1:xsz(3))  * s11(:,:,1:xsz(3))
      t12(:,:,1:xsz(3))  = -2 * nut(:,:,1:xsz(3))  * s12(:,:,1:xsz(3))
      t13(:,:,1:xsz(3))  = -2 * nut(:,:,1:xsz(3))  * s13(:,:,1:xsz(3))
      t22(:,:,1:xsz(3))  = -2 * nut(:,:,1:xsz(3))  * s22(:,:,1:xsz(3))
      t23(:,:,1:xsz(3))  = -2 * nut(:,:,1:xsz(3))  * s23(:,:,1:xsz(3))
      t33(:,:,1:xsz(3))  = -2 * nut(:,:,1:xsz(3))  * s33(:,:,1:xsz(3))
      t11w(:,:,1:xsz(3)) = -2 * nutw(:,:,1:xsz(3)) * s11w(:,:,1:xsz(3))
      t12w(:,:,1:xsz(3)) = -2 * nutw(:,:,1:xsz(3)) * s12w(:,:,1:xsz(3))
      t13w(:,:,1:xsz(3)) = -2 * nutw(:,:,1:xsz(3)) * s13w(:,:,1:xsz(3))
      t22w(:,:,1:xsz(3)) = -2 * nutw(:,:,1:xsz(3)) * s22w(:,:,1:xsz(3))
      t23w(:,:,1:xsz(3)) = -2 * nutw(:,:,1:xsz(3)) * s23w(:,:,1:xsz(3))
      t33w(:,:,1:xsz(3)) = -2 * nutw(:,:,1:xsz(3)) * s33w(:,:,1:xsz(3))

      call dealiasxy(t11(:,:,1:xsz(3)))
      call dealiasxy(t11w(:,:,1:xsz(3)))
      call dealiasxy(t12(:,:,1:xsz(3)))
      call dealiasxy(t12w(:,:,1:xsz(3)))
      call dealiasxy(t13(:,:,1:xsz(3))) 
      call dealiasxy(t13w(:,:,1:xsz(3)))
      call dealiasxy(t22(:,:,1:xsz(3)))
      call dealiasxy(t22w(:,:,1:xsz(3)))
      call dealiasxy(t23(:,:,1:xsz(3)))
      call dealiasxy(t23w(:,:,1:xsz(3)))
      call dealiasxy(t33(:,:,1:xsz(3)))
      call dealiasxy(t33w(:,:,1:xsz(3)))
      call update_ghost(t11, level)
      call update_ghost(t11w, level)
      call update_ghost(t12, level)
      call update_ghost(t12w, level)
      call update_ghost(t13, level)
      call update_ghost(t13w, level)
      call update_ghost(t22, level)
      call update_ghost(t22w, level)
      call update_ghost(t23, level) 
      call update_ghost(t23w, level)
      call update_ghost(t33, level)
      call update_ghost(t33w, level)

      ! print*,'t11',  t11(1,1,2:xsz(3)-1)
      ! print*,'t12',  t12(1,1,2:xsz(3)-1)
      ! print*,'t13',  t13w(1,1,2:xsz(3)-1)

    end subroutine get_SGS_stress
    
   !checked
   subroutine div_tau(t1, t2, t3w, divtau)
      use grid, only : zetax, zetay, her, dz
      use spectral
      use param, only: pex, pey
      implicit none

      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: t1, t2
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) ::t3w
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: divtau

      integer :: k

      tmp_x1 = 0.
      ! compute d(tau11)/dx
      call pdfx(t1(:,:,1:xsz(3)), tmp_x1(:,:,1:xsz(3)), pex)
      divtau(:,:,1:xsz(3)) = tmp_x1(:,:,1:xsz(3))

      tmp_x1 = 0.
      call calc_uzeta(t1, tmp_x1(:,:,1:xsz(3)), level)
      if (isbot) then
         tmp_x1(:,:,1) = 0
         tmp_x1(:,:,2) = (t1(:,:,3) - t1(:,:,2))/dz(2)
      end if
      if (istop) then
         tmp_x1(:,:,xsz(3)-1) = (t1(:,:,xsz(3)-1) - t1(:,:,xsz(3)-2))/dz(xsz(3)-2)/2.0_wp
      end if
      
      divtau(:,:,1:xsz(3)) = divtau(:,:,1:xsz(3)) + zetax(:,:,1:xsz(3))*tmp_x1(:,:,1:xsz(3))

      tmp_x1=0.
      ! compute d(tau12)/dy
      call pdfy_x(t2(:,:,1:xsz(3)), tmp_x1(:,:,1:xsz(3)), pey)
      divtau(:,:,1:xsz(3)) = divtau(:,:,1:xsz(3)) + tmp_x1(:,:,1:xsz(3))
      call calc_uzeta(t2, tmp_x1(:,:,1:xsz(3)), level)
      if (isbot) then
         tmp_x1(:,:,1) = 0
         tmp_x1(:,:,2) = (t2(:,:,3) - t2(:,:,2))/dz(2)
      end if
      if (istop) then
         tmp_x1(:,:,xsz(3)-1) = (t2(:,:,xsz(3)-1) - t2(:,:,xsz(3)-2))/dz(xsz(3)-2)/2.0_wp
      end if

      divtau(:,:,1:xsz(3)) = divtau(:,:,1:xsz(3)) + zetay(:,:,1:xsz(3))*tmp_x1(:,:,1:xsz(3))

      tmp_x1=0
      ! compute d(tau13)/dz
      do k=1, xsz(3)
         tmp_x1(:,:,k) = (t3w(:,:,k)-t3w(:,:,k-1))*2.0_wp/(dz(k)+dz(k-1))
      end do
      if (isbot) then
         tmp_x1(:,:,1) = 0
         tmp_x1(:,:,2) = (t3w(:,:,2) - t3w(:,:,1))/dz(2)    
      end if
      if (istop) then
         tmp_x1(:,:,xsz(3)-1) = (t3w(:,:,xsz(3)-1)-t3w(:,:,xsz(3)-2)) / dz(xsz(3)-2)
      end if
      
      do k=1, xsz(3)
         divtau(:,:,k) = divtau(:,:,k) + her(:,:) * tmp_x1(:,:,k)
      enddo

   end subroutine div_tau

   !checked
   subroutine div_tau_w(t1, t2, t3, t1w, t2w,  divtau)
      use grid, only : zetax_w, zetay_w, pex, pey, her, dz
      use spectral
      implicit none

      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: t1, t2, t3
      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: t1w, t2w
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: divtau

      integer :: k

      tmp_x1 = 0.
      ! compute d(tau31)/dx
      call pdfx(t1w(:,:,1:xsz(3)), tmp_x1(:,:,1:xsz(3)), pex)
      divtau(:,:,1:xsz(3)) = tmp_x1(:,:,1:xsz(3))

      tmp_x1 = 0.
      do k=1,xsz(3)
         tmp_x1(:,:,k) = (t1(:,:,k+1)-t1(:,:,k))/dz(k)
      enddo
      
      divtau(:,:,1:xsz(3)) = divtau(:,:,1:xsz(3)) + zetax_w(:,:,1:xsz(3))*tmp_x1(:,:,1:xsz(3))

      tmp_x1=0.
      ! compute d(tau32)/dy
      call pdfy_x(t2w(:,:,1:xsz(3)), tmp_x1(:,:,1:xsz(3)), pey)
      divtau(:,:,1:xsz(3)) = divtau(:,:,1:xsz(3)) + tmp_x1(:,:,1:xsz(3))

      tmp_x1=0.
      do k=1, xsz(3)
         tmp_x1(:,:,k) = (t2(:,:,k+1)-t2(:,:,k))/dz(k)
      end do      

      divtau(:,:,1:xsz(3)) = divtau(:,:,1:xsz(3)) + zetay_w(:,:,1:xsz(3))*tmp_x1(:,:,1:xsz(3))

      tmp_x1=0.
      ! compute d(tau33)/dz
      do k=1, xsz(3)
         tmp_x1(:,:,k) = her(:,:) *(t3(:,:,k+1)-t3(:,:,k)) / dz(k)
      end do
      
      divtau(:,:,1:xsz(3)) = divtau(:,:,1:xsz(3)) + tmp_x1(:,:,1:xsz(3))

   end subroutine div_tau_w
 
