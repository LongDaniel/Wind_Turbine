!checked
subroutine init_p_poisson()
      use grid, only : pex, pey, hbar, dz
      use poisson

      implicit none

      integer  :: k, m, l, nz, modex, modey
      real(wp) :: c, ax, ay

      nz = ysz(3)

      do m = 1,ysz(2)
         ax = pex * ((yst(2)+m-2)/2)
         do l = 1,ysz(1)
            ay = pey * ((l-1)/2)
            c = -ax**2-ay**2

            modex=(yst(2)+m-2)/2
            modey=(l-1)/2

            ! coefficients for pressure 
            do k=1, xsz(3)
               lap_p_al(l,m,k)=1/dz(k-1)
               lap_p_ao(l,m,k)=c*hbar**2*(dz(k)+dz(k-1))/2-(1/dz(k-1)+1/dz(k))
               lap_p_ar(l,m,k)=1/dz(k)
            end do

            if (isbot) then
               lap_p_al(l,m,2)=0
               lap_p_ao(l,m,2)=c*hbar**2*dz(2)-1/dz(2)
               lap_p_ar(l,m,2)=1/dz(2)
   
               lap_p_al(l,m,3)=1/dz(2)
               lap_p_ao(l,m,3)=c*hbar**2*dz(2)-2/dz(2)
               lap_p_ar(l,m,3)=1/dz(2)

             if(modex.eq.0.and.modey.eq.0)then
                 lap_p_al(1,1,2)=0
                 lap_p_ao(1,1,2)=1.0_wp
                 lap_p_ar(1,1,2)=0
             end if
            end if

            if (istop) then
               lap_p_al(l,m,nz-2)=1/dz(nz-2)
               lap_p_ao(l,m,nz-2)=c*hbar**2*dz(nz-2)-2/dz(nz-2)
               lap_p_ar(l,m,nz-2)=1/dz(nz-2)
               
               lap_p_al(l,m,nz-1)=1/dz(nz-2)
               lap_p_ao(l,m,nz-1)=c*hbar**2*dz(nz-2)-1/dz(nz-2)
               lap_p_ar(l,m,nz-1)=0
            end if

         end do
      end do
      ! perform LU factorization
      if (istop.and.isbot) then
         call tridiag_lu(lap_p_ao(:,:,2:nz-1),lap_p_al(:,:,2:nz-1),lap_p_ar(:,:,2:nz-1), &
                         nz-2,ysz(1)*ysz(2),MPI_COMM_2D_ROW)
      else if (istop.and.(.not.isbot)) then
         call tridiag_lu(lap_p_ao(:,:,1:nz-1),lap_p_al(:,:,1:nz-1),lap_p_ar(:,:,1:nz-1), &
                         nz-1,ysz(1)*ysz(2),MPI_COMM_2D_ROW)
      else if (isbot.and.(.not.istop)) then
         call tridiag_lu(lap_p_ao(:,:,2:nz),lap_p_al(:,:,2:nz),lap_p_ar(:,:,2:nz), &
                         nz-1,ysz(1)*ysz(2),MPI_COMM_2D_ROW)
      else
         call tridiag_lu(lap_p_ao,lap_p_al,lap_p_ar,nz,ysz(1)*ysz(2),MPI_COMM_2D_ROW)
      end if

    end subroutine init_p_poisson

    !checked
    subroutine press_g(time)
      use MPI
      use spectral
      use utils
      use param, only: itmax, timewavy
      use discontinuity_smooth, only: i_div_ustar !> for debug

      implicit none

      !real(wp), dimension(xsz(1), xsz(2)), intent(IN) :: ut, vt, wt
      ! real(wp), dimension(xsz(1), xsz(2)), intent(OUT) :: pf0
      
      real(wp), intent(IN) :: time
      integer  :: it, nz, ierr
      real(wp) :: err1, err2, ome
      real(wp) :: er1, er2

      nz = xsz(3)
      if (myid == 0) print *,'Solving p, eps=', erlim

      ! rhs: Div(u)
      i_div_ustar = 1 
      call source_s_cnn(tmp_x5(1,1,1))
      i_div_ustar = 0
      ! print*, tmp_x5(1,1,1:xsz(3))

      !-----------------
      ! Iteration start
      !-----------------
      do it = 1, itmax
      ! do it = 1, 1

         call poisson_du_cnn(tmp_x5(1,1,1), tmp_x6(1,1,1))

         err1=0
         if (isbot.and.istop) then 
            err1 = max_abs_diff(pp(:,:,2:nz-1), tmp_x6(:,:,2:nz-1))
         else if (isbot.and.(.not.istop)) then 
            err1 = max_abs_diff(pp(:,:,2:nz), tmp_x6(:,:,2:nz))
         else if (istop.and.(.not.isbot)) then
            err1 = max_abs_diff(pp(:,:,1:nz-1), tmp_x6(:,:,1:nz-1))
         else
            err1 = max_abs_diff(pp(:,:,1:nz), tmp_x6(:,:,1:nz))
         end if
         call mpi_allreduce(abs(err1),er1,1,mpi_double_precision, &
              mpi_max,mpi_comm_world,ierr)

         call update_ghost(tmp_x6, level)
         pp=tmp_x6

         !IF FLAT BOTTOM, NO NEED FOR ITERATION
         if (time.le.timewavy) goto 999
         
         call poisson_du_cnn(tmp_x5(1,1,1), tmp_x7(1,1,1))

         err2=0
         if (isbot.and.istop) then 
            err2 = max_abs_diff(tmp_x6(:,:,2:nz-1), tmp_x7(:,:,2:nz-1))
         else if (isbot.and.(.not.istop)) then 
            err2 = max_abs_diff(tmp_x6(:,:,2:nz), tmp_x7(:,:,2:nz))
         else if (istop.and.(.not.isbot)) then
            err2 = max_abs_diff(tmp_x6(:,:,1:nz-1), tmp_x7(:,:,1:nz-1))
         else
            err2 = max_abs_diff(tmp_x6(:,:,1:nz), tmp_x7(:,:,1:nz))
         end if
         call mpi_allreduce(abs(err2),er2,1,mpi_double_precision, &
              mpi_max,mpi_comm_world,ierr)

         if(xsz(2).eq.0) then
            ome=0
            er2=0
            goto 888
         endif

         if(abs(er2)+abs(er1).eq.0) then
            print*, 'error at |er2+er1|'
            stop
         endif

         ome=abs(er2)/(abs(er2)+abs(er1))
         
 888     continue

         if (isbot.and.istop) then 
            pp(:,:,2:nz-1)=tmp_x7(:,:,2:nz-1)-ome*(tmp_x7(:,:,2:nz-1)-tmp_x6(:,:,2:nz-1))
         else if (isbot.and.(.not.istop)) then 
            pp(:,:,2:nz)=tmp_x7(:,:,2:nz)-ome*(tmp_x7(:,:,2:nz)-tmp_x6(:,:,2:nz))
         else if (istop.and.(.not.isbot)) then 
            pp(:,:,1:nz-1)=tmp_x7(:,:,1:nz-1)-ome*(tmp_x7(:,:,1:nz-1)-tmp_x6(:,:,1:nz-1))
         else 
            pp(:,:,1:nz)=tmp_x7(:,:,1:nz)-ome*(tmp_x7(:,:,1:nz)-tmp_x6(:,:,1:nz))
         end if 

         call update_ghost(pp, level)

         if(myid.eq.0) write(*,*) it,er2,erlim
         if(abs(er2).lt.erlim) goto 999

      end do 

 999  continue

      return

   end subroutine press_g
    


   !checked
   !input: ut, vt, wt
   !output: divu
   subroutine source_s_cnn(divu)
     use grid,only : dz, her, zetax, zetay
     use param, only: pex, pey
      use spectral
      !use navier, only: ub,vb,wb
      implicit none

      !real(wp), intent(IN) :: dt
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: divu

      integer :: k

      ! ---------------
      ! calculate du/dx
      ! ---------------
      call pdfx(u(:,:,1:xsz(3)),u_xi(:,:,1:xsz(3)),pex)
      call update_ghost(u_xi,level)
      call calc_uzeta(u, u_zeta(:,:,1:), level)
      call update_ghost(u_zeta,level)
      if (isbot) then
         u_zeta(:,:,2)=(u(:,:,2)-ub(:,:)) /dz(1)
         u_zeta(:,:,3)=(u(:,:,4)-u(:,:,2)) /dz(2)/2.0_wp         
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-2)=(u(:,:,xsz(3)-1)-u(:,:,xsz(3)-3))/2.0_wp/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)-1)=(u(:,:,xsz(3)-1)-u(:,:,xsz(3)-2))/2.0_wp/dz(xsz(3)-2)
      end if
      
      do k=1, xsz(3)
         divu(:,:,k)=u_xi(:,:,k)+zetax(:,:,k)*u_zeta(:,:,k)
      end do
      
      ! ---------------
      ! calculate dv/dy
      ! ---------------
      call pdfy_x(v(:,:,1:xsz(3)),u_psi(:,:,1:xsz(3)),pey)
      call update_ghost(u_psi,level)
      call calc_uzeta(v, u_zeta(:,:,1:), level)
      call update_ghost(u_zeta,level)
      if (isbot) then
         u_zeta(:,:,2)=(v(:,:,2)-vb(:,:)) /dz(1)
         u_zeta(:,:,3)=(v(:,:,4)-v(:,:,2)) /dz(2)/2.0_wp
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-2)=(v(:,:,xsz(3)-1)-v(:,:,xsz(3)-3))/2.0_wp/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)-1)=(v(:,:,xsz(3)-1)-v(:,:,xsz(3)-2))/2.0_wp/dz(xsz(3)-2)
      end if

      do k=1, xsz(3)
         divu(:,:,k)=divu(:,:,k)+u_psi(:,:,k)+zetay(:,:,k)*u_zeta(:,:,k)
      end do

      ! ---------------
      ! calcualte dw/dz
      ! ---------------
      do k=1, xsz(3)
         u_zeta(:,:,k)=2.0_wp*(w(:,:,k)-w(:,:,k-1))/(dz(k-1)+dz(k))
      end do 
      call update_ghost(u_zeta,level)
      if (isbot) then 
         u_zeta(:,:,2)=(w(:,:,2)-wb(:,:))/dz(2)
         u_zeta(:,:,3)=(w(:,:,3)-w(:,:,2))/dz(2)
      end if 
      if (istop) then
         u_zeta(:,:,xsz(3)-2)=(w(:,:,xsz(3)-2)-w(:,:,xsz(3)-3))/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)-1)=(w(:,:,xsz(3)-1)-w(:,:,xsz(3)-2))/dz(xsz(3)-2)
      end if

      do k=1, xsz(3)
         divu(:,:,k)=divu(:,:,k)+her(:,:)*u_zeta(:,:,k)
      end do
      
      divu=divu/dt

      call dealiasxy(divu(:,:,1:xsz(3)))

      ! print*, divu(1,1,2:xsz(3)-1)

   end subroutine source_s_cnn



 !checked  
   subroutine poisson_du_cnn(divu,sigma)
      use MPI
      use fft
      use poisson
      use spectral
      use param,only: hbar
      use grid, only: dz

      implicit none

      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(IN) :: divu
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: sigma

      integer :: k, l, m, nz  , ierr
      real(wp) :: sa, sai

      nz=ysz(3)

      call source_e_cnn(pp,sigma)

      ! print*,'sigma', sigma(1,1,2:xsz(3)-1)

      sigma=-sigma+divu
      
      call dealiasxy(sigma)
      
      !> modified by plyu
      if (idsp .eq. 1 .or. idsp .eq. 2) then
        call fft_r2c_xy_out_s32 (sigma(:,:,1:nz), tmp_y1(:,:,1:nz))
      else
        call fft_r2c_xy(sigma(:,:,1:nz),tmp_y1(:,:,1:nz))
      endif
 
      tmp_y1 = tmp_y1/nx_global/ny_global

      !> Following lines are moved to 2decomp/fft_fftw3.f90:fft_r2c_xy_out_s32
      !!> added by plyu, do a smoothing to remove high frequency part when a notable
      !!!   discontinuity was introduced into the fluid field by turbine force
      !if (idsp .eq. 1) then
      !  !> dimension of temp_y1 is (ysz(1), ysz(2), 1-level:ysz(3)+level)
      !  do m = 1, ysz(2)/2
      !    tmp_y1(:,(2*m-1):(2*m),:) = co_s32_x((yst(2)+1)/2+m-1) * tmp_y1(:,(2*m-1):(2*m),:)
      !  enddo
      !  do l = 1, ysz(1)/2
      !    tmp_y1((2*l-1):(2*l),:,:) = co_s32_y(l) * tmp_y1((2*l-1):(2*l),:,:)
      !  enddo
      !endif

      ! compute rhs 
      m=1; l=nz
      if (isbot) then
          tmp_y1(:,:,2)=tmp_y1(:,:,2)*hbar**2*dz(2)
          tmp_y1(:,:,3)=tmp_y1(:,:,3)*hbar**2*dz(2)
          m=4
      end if 
      if (istop) then
          tmp_y1(:,:,nz-2)=tmp_y1(:,:,nz-2)*hbar**2*dz(nz-2)
          tmp_y1(:,:,nz-1)=tmp_y1(:,:,nz-1)*hbar**2*dz(nz-2)
          l=nz-3
      end if 
      do k=m, l
          tmp_y1(:,:,k)=tmp_y1(:,:,k)*hbar**2*(dz(k-1)+dz(k))/2.0_wp 
      end do

      ! compute the average
      if (myid1==0)then 
          sai=0.0_wp
          if (isbot.and.istop) then
              do k=2, nz-1
                 sai=sai+tmp_y1(1,1,k)
              end do 
          else if (isbot.and.(.not.istop)) then
              do k=2, nz
                 sai=sai+tmp_y1(1,1,k)
              end do 
          else if (istop.and.(.not.isbot)) then 
              do k=1, nz-1
                 sai=sai+tmp_y1(1,1,k)
              end do 
          else
              do k=1, nz
                 sai=sai+tmp_y1(1,1,k)
              end do 
          end if
         
          sa=0.0_wp
          call MPI_Allreduce(sai, sa, 1, real_type, MPI_SUM, MPI_COMM_2D_ROW, ierr)
          sa=sa/(nz_global-2)
    
          if (isbot) then
              tmp_y1(1,1,2)=0
              tmp_y1(1,1,3)=tmp_y1(1,1,3)-sa
              m=4
          end if 
          if (istop) then
              tmp_y1(1,1,nz-2)=tmp_y1(1,1,nz-2)-sa
              tmp_y1(1,1,nz-1)=tmp_y1(1,1,nz-1)-sa
              l=nz-3
          end if 
          do k=m, l
              tmp_y1(1,1,k)=tmp_y1(1,1,k)-sa
          end do
      end if 


      if (istop.and.isbot) then
         call tridiag_lu_solve(lap_p_ao(:,:,2:nz-1),lap_p_al(:,:,2:nz-1),lap_p_ar(:,:,2:nz-1), &
                         tmp_y1(:,:,2:nz-1),nz-2,ysz(1)*ysz(2),MPI_COMM_2D_ROW)
      else if (istop.and.(.not.isbot)) then
         call tridiag_lu_solve(lap_p_ao(:,:,1:nz-1),lap_p_al(:,:,1:nz-1),lap_p_ar(:,:,1:nz-1), &
                         tmp_y1(:,:,1:nz-1),nz-1,ysz(1)*ysz(2),MPI_COMM_2D_ROW)
      else if (isbot.and.(.not.istop)) then
         call tridiag_lu_solve(lap_p_ao(:,:,2:nz),lap_p_al(:,:,2:nz),lap_p_ar(:,:,2:nz), &
                         tmp_y1(:,:,2:nz),nz-1,ysz(1)*ysz(2),MPI_COMM_2D_ROW)
      else
         call tridiag_lu_solve(lap_p_ao(:,:,1:nz),lap_p_al(:,:,1:nz),lap_p_ar(:,:,1:nz), &
                         tmp_y1(:,:,1:nz),nz,ysz(1)*ysz(2),MPI_COMM_2D_ROW)
      end if

      ! dealias and inverse transform
      m = ((xsz(1)/2)*2/3)*2+1
      l = ((ysz(1)/2)*2/3)*2+1
      ! tmp_y1(l:,:,:) = 0

      !> added by plyu: remove high frequency part after solving the pressure
      !! poisson equation
      if (idsp .eq. 2) then
        do l = 1, ysz(1)/2
          tmp_y1((2*l-1):(2*l),:,:) = co_s32_y(l) * tmp_y1((2*l-1):(2*l),:,:)
        enddo
      endif

      call fft_c2r_y(tmp_y1(:,:,1:nz))
      call transpose_yx(tmp_y1(:,:,1:nz), sigma(:,:,1:nz))
      ! sigma(m:,:,:) = 0
      
      !> added by plyu: remove high frequency part after solving the pressure
      !! poisson equation
      if (idsp .eq. 2) then
        do l = 1, xsz(1)/2
          sigma((2*l-1):(2*l),:,:) = co_s32_x(l) * sigma((2*l-1):(2*l),:,:)
        enddo
      endif
      
      call fft_c2r_x(sigma(:,:,1:nz))

      tmp_y1=0
      
   end subroutine poisson_du_cnn

!checked
   subroutine source_e_cnn(p,sigma)
     use grid, only : dz, zetax, zetay, her
     use param, only: pex, pey, hbar
      use spectral

      implicit none

      real(wp), dimension(xsz(1),xsz(2),1-level:xsz(3)+level), intent(IN) :: p
      real(wp), dimension(xsz(1),xsz(2),xsz(3)), intent(OUT) :: sigma
      real(wp), dimension(xsz(1),xsz(2)) :: pt, ptx, pty

      integer :: k

      call pdfx(p(:,:,1:xsz(3)),u_xi(:,:,1:xsz(3)),pex)
      call update_ghost(u_xi,level)
      if (istop) then
         call pdfx(p(:,:,xsz(3)+1),u_xi(:,:,xsz(3)+1),pex)
      end if

      call pdfy_x(p(:,:,1:xsz(3)),u_psi(:,:,1:xsz(3)),pey)
      call update_ghost(u_psi,level)
      if (istop) then
         call pdfy_x(p(:,:,xsz(3)+1),u_psi(:,:,xsz(3)+1),pey)
      end if

      ! compute dp/dzeta and impose pressure boundary conditions
      call calc_uzeta(p, u_zeta(:,:,1:), level)
      call update_ghost(u_zeta,level)
      if (isbot) then
         !u_zeta(:,:,1)=(-p(:,:,3)+9.0_wp*p(:,:,2)-8.0_wp*p(:,:,1))/3.0_wp/dz(2) !for test
         !u_zeta(:,:,2)=(-p(:,:,4)+4.0_wp*p(:,:,3)-3.0_wp*p(:,:,2))/2.0_wp/dz(2)
         u_zeta(:,:,2)=(p(:,:,3)-p(:,:,2))/dz(2)
         u_zeta(:,:,3)=(p(:,:,4)-p(:,:,2))/2.0_wp/dz(2)
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-2)=(p(:,:,xsz(3)-1)-p(:,:,xsz(3)-3))&
              /2.0_wp/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3)-1)=(3.0_wp*p(:,:,xsz(3)-1)-4.0_wp*p(:,:,xsz(3)-2) &
                                +p(:,:,xsz(3)-3))/2.0_wp/dz(xsz(3)-2)
         u_zeta(:,:,xsz(3))=(2.0_wp*p(:,:,xsz(3)-1)-3.0_wp*p(:,:,xsz(3)-2) &
                                +p(:,:,xsz(3)-3))/dz(xsz(3)-2)
      end if

      ! compute dp/dx and store in array tmp_x1
      tmp_x3(:,:,:) = zetax(:,:,:)*u_zeta(:,:,:)
      tmp_x1(:,:,:) = u_xi(:,:,:)+tmp_x3(:,:,:)
      ! compute d(dp/dx)/dzeta and store in array tmp_x2
      call calc_uzeta(tmp_x1, tmp_x2(:,:,1:), level)
      if (isbot) then
         !tmp_x2(:,:,2)=(tmp_x1(:,:,3)+3.0_wp*tmp_x1(:,:,2)) / 3.0_wp/dz(2)
         tmp_x2(:,:,2)=tmp_x1(:,:,2) / dz(1)
      end if
      if (istop) then
         pt(:,:)=(35.0_wp*p(:,:,xsz(3)-1)-35.0_wp*p(:,:,xsz(3)-2)+21.0_wp*p(:,:,xsz(3)-3) &
                 -5.0_wp*p(:,:,xsz(3)-4)) / 16.0_wp

         call pdfx(pt(:,:),ptx(:,:),pex)
         ptx(:,:)=ptx(:,:)+tmp_x3(:,:,xsz(3))
         tmp_x2(:,:,xsz(3)-1)= (4.0_wp*ptx(:,:)-3.0_wp*tmp_x1(:,:,xsz(3)-1) &
                              -tmp_x1(:,:,xsz(3)-2)) / 3.0_wp/dz(xsz(3)-2)
      end if
      call dealiasxy(tmp_x2(:,:,1:xsz(3)))
      ! compute zeta_x * d(dp/dx)/dzeta
      tmp_x2(:,:,1:xsz(3)) = tmp_x2(:,:,1:xsz(3))*zetax(:,:,1:xsz(3))

      sigma(:,:,1:xsz(3)) = tmp_x2(:,:,1:xsz(3))

      ! compute d(zetax*u_zeta)/dxi 
      call pdfx(tmp_x3(:,:,1:xsz(3)), tmp_x2(:,:,1:xsz(3)), pex)
      sigma(:,:,1:xsz(3)) = sigma(:,:,1:xsz(3)) + tmp_x2(:,:,1:xsz(3))

      ! compute dp/dy and store in array tmp_x1
      tmp_x3(:,:,:) = zetay(:,:,:)*u_zeta(:,:,:)
      tmp_x1(:,:,:) = u_psi(:,:,:)+tmp_x3(:,:,:)
      ! compute d(dp/dy)/dzeta and store in array tmp_x2
      call calc_uzeta(tmp_x1, tmp_x2(:,:,1:), level)
      if (isbot) then
         tmp_x2(:,:,2)=tmp_x1(:,:,2)/dz(1) 
      end if
      if (istop) then
         call pdfy_x(pt(:,:),pty(:,:),pey)
         pty(:,:)=pty(:,:)+tmp_x3(:,:,xsz(3))
         tmp_x2(:,:,xsz(3)-1)= (4.0_wp*pty(:,:)-3.0_wp*tmp_x1(:,:,xsz(3)-1) &
                              -tmp_x1(:,:,xsz(3)-2)) / 3.0_wp/dz(xsz(3)-2)
      end if
      call dealiasxy(tmp_x2(:,:,1:xsz(3)))

      ! compute zeta_y * d(dp/dy)/dzeta
      tmp_x2(:,:,1:xsz(3)) = tmp_x2(:,:,1:xsz(3))*zetay(:,:,1:xsz(3))
      sigma(:,:,1:xsz(3)) = sigma(:,:,1:xsz(3)) + tmp_x2(:,:,1:xsz(3))

      ! compute d(zetay*u_zeta)/dpsi
      call pdfy_x(tmp_x3(:,:,1:xsz(3)), tmp_x2(:,:,1:xsz(3)), pey)
      sigma(:,:,1:xsz(3)) = sigma(:,:,1:xsz(3)) + tmp_x2(:,:,1:xsz(3))

      ! compute d(dp/zeta)/dzeta and store in tmp_x1
      do k=1,xsz(3)
         tmp_x1(:,:,k) = (p(:,:,k+1)-p(:,:,k))/dz(k)-(p(:,:,k)-p(:,:,k-1))/dz(k-1)
         tmp_x1(:,:,k) = tmp_x1(:,:,k)*2.0_wp/(dz(k)+dz(k-1))
      end do
      if (isbot) then
         tmp_x1(:,:,2) = (p(:,:,3)-p(:,:,2))/dz(2)**2
         tmp_x1(:,:,3) = (p(:,:,4)-2.*p(:,:,3)+p(:,:,2))/dz(2)**2
      end if
      if (istop) then
         ! tmp_x1(:,:,xsz(3)-1) = (p(:,:,xsz(3)-2)-p(:,:,xsz(3)-1))/dz(xsz(3)-2)**2
         tmp_x1(:,:,xsz(3)-2) = (p(:,:,xsz(3)-1)-2*p(:,:,xsz(3)-2)+p(:,:,xsz(3)-3))/dz(xsz(3)-2)**2         
         tmp_x1(:,:,xsz(3)-1) = -(p(:,:,xsz(3)-1)-p(:,:,xsz(3)-2))/dz(xsz(3)-2)**2
      end if
      do k=1,xsz(3)
         tmp_x2(:,:,k) = tmp_x1(:,:,k)*her(:,:)
      end do
      tmp_x1 = tmp_x1/hbar**2

      ! print*, her(:,1)*tmp_x2(:,1,xsz(3)-1) -tmp_x1(:,1,xsz(3)-1)
      do k=1,xsz(3)
         sigma(:,:,k) = sigma(:,:,k)+her(:,:)*tmp_x2(:,:,k) -tmp_x1(:,:,k)
      end do
      ! print*, sigma(:,1,xsz(3)-1)
      call dealiasxy(sigma(:,:,1:xsz(3)))

      tmp_x1=0
      tmp_x2=0
      tmp_x3=0

      
   end subroutine source_e_cnn
   
   !checked
   !input: pp
   !output: u, v, w
   subroutine correction_us
     use grid, only : dz, her, zetax, zetay
     use param, only: pex, pey, dt
     use spectral
     use discontinuity_smooth, only : i_grad_p
     
     implicit none
     
     
      integer :: k, nz

      nz=xsz(3)

      if (isbot.and.istop)then 
         pp(:,:,1)=(3*pp(:,:,2)-pp(:,:,3))/2
         pp(:,:,nz)=pp(:,:,nz-1)
      else if (isbot.and.(.not.istop))then 
         pp(:,:,1)=(3*pp(:,:,2)-pp(:,:,3))/2         
      else if (istop.and.(.not.isbot))then 
         pp(:,:,nz)=pp(:,:,nz-1)
      end if

      !i_grad_p = 1
      call pdfx(pp(:,:,1:xsz(3)),u_xi(:,:,1:xsz(3)),pex)
      !i_grad_p = 0
      call update_ghost(u_xi,level)

      i_grad_p = 1
      call pdfy_x(pp(:,:,1:xsz(3)),u_psi(:,:,1:xsz(3)),pey)
      i_grad_p = 0
      call update_ghost(u_psi,level)

      ! compute dp/dzeta
      call calc_uzeta(pp, u_zeta(:,:,1:), level)
      call update_ghost(u_zeta,level)
      if (isbot) then
         u_zeta(:,:,1)=0
         u_zeta(:,:,2)=(pp(:,:,3)-pp(:,:,2))/dz(2)
      end if
      if (istop) then
         u_zeta(:,:,nz-1)=(pp(:,:,nz-1)-pp(:,:,nz-2))/2.0_wp/dz(nz-2)
         u_zeta(:,:,nz)=0
      end if

      ! compute dp/dx and correct u
      tmp_x3(:,:,:) = zetax(:,:,:)*u_zeta(:,:,:)
      u(:,:,:) = u(:,:,:) - (u_xi(:,:,:)+tmp_x3(:,:,:))*dt

      ! compute dp/dy and correct v
      tmp_x3(:,:,:) = zetay(:,:,:)*u_zeta(:,:,:)
      v(:,:,:) = v(:,:,:) - (u_psi(:,:,:)+tmp_x3(:,:,:))*dt

      ! compute dp/dz and correct w
      do k=1, xsz(3)
         u_zeta(:,:,k) = (pp(:,:,k+1)-pp(:,:,k))/dz(k)
      end do

      do k=1, xsz(3)
         u_zeta(:,:,k) = u_zeta(:,:,k) * her(:,:)
      end do
      
      if (istop)then
          w(:,:,1:nz-2) = w(:,:,1:nz-2) - dt*u_zeta(:,:,1:nz-2)
          w(:,:,nz-1)=0
      else 
          w(:,:,1:xsz(3)) = w(:,:,1:xsz(3)) - dt*u_zeta(:,:,1:xsz(3))
      end if 
       
      call dealiasxy(u(:,:,1:xsz(3)))
      call dealiasxy(v(:,:,1:xsz(3)))
      call dealiasxy(w(:,:,1:xsz(3)))

      if(istop)then 
         call dealiasxy(u(:,:,xsz(3)+1))
         call dealiasxy(v(:,:,xsz(3)+1))
         call dealiasxy(w(:,:,xsz(3)+1))
      end if 

      call update_ghost(u, level)
      call update_ghost(v, level)
      call update_ghost(w, level)

      ! print*, u(1,1,2:xsz(3)+1)
   end subroutine correction_us

   !checked
   !input: uzfs, vzfs
   !input: ub, vb, wb
   !output: u, v, w
   subroutine bc_lnr
      use spectral
      use grid,only : her, dz, zz, hyr, eyr, hxr, exr
      !use navier, only: ub, vb, wb
      use param, only: pex, pey
      implicit none

      real(wp), dimension(xsz(1),xsz(2)) :: ux, vy
      
      integer :: nz

      nz=xsz(3)

      if (isbot)then 
          u(:,:,1)=ub(:,:)
          v(:,:,1)=vb(:,:)
          w(:,:,1)=wb(:,:)
      end if

      if (istop) then 
         call pdfx(u(:,:,nz),ux,pex)
         call pdfy_x(v(:,:,nz),vy,pey)

         ux(:,:)=ux(:,:)+(hxr(:,:)-zz(nz)*exr(:,:))*uzfs(:,:)
         vy(:,:)=vy(:,:)+(hyr(:,:)-zz(nz)*eyr(:,:))*vzfs(:,:)

         u(:,:,nz+1)=(12*dz(nz-2)*uzfs(:,:)-u(:,:,nz-2)+18*u(:,:,nz-1)-8*u(:,:,nz))/9.
         v(:,:,nz+1)=(12*dz(nz-2)*vzfs(:,:)-v(:,:,nz-2)+18*v(:,:,nz-1)-8*v(:,:,nz))/9.
         w(:,:,nz-1)=0
         w(:,:,nz)=(6*dz(nz-2)*(-ux(:,:)-vy(:,:))/her(:,:)-w(:,:,nz-3)+6*w(:,:,nz-2)-3*w(:,:,nz-1))/2.

         ! if(myid==0)print*,'bc_lnr',(hxr(1,1)-zz(nz)*exr(1,1))
      end if 

   end subroutine bc_lnr

   ! function max_abs_diff(a, b) result(maxdiff)
   !    implicit none

   !    real(wp), dimension(:,:,:), contiguous, intent(IN) :: a, b
   !    real(wp) :: maxdiff

   !    integer  :: i,j,k
   !    real(wp) :: tmp

   !    maxdiff = 0
   !    do k=1, size(a,3)
   !       do j=1, size(a,2)
   !          do i=1, size(a,1)
   !             tmp = abs(a(i,j,k)-b(i,j,k))
   !             if (tmp > maxdiff) maxdiff = tmp
   !          end do
   !       end do
   !    end do

   ! end function max_abs_diff

   !input: u, v, w
   subroutine volume_lnr2(time)
     use grid,only : dz, her, zetax, zetay
     use param, only: pex, pey, nx, ny, nz
      use spectral
      use mpi

      implicit none

      integer :: i, j, k, ierr
      integer :: imax, jmax, kmax

      real(wp) :: div_all, divmax
      real(wp), dimension(xsz(1),xsz(2),xsz(3)) :: divu

      real(wp), intent(IN) :: time
      
      !> added by plyu
      real(wp) :: norm1_div, norm2_div, norm1_div_all, norm2_div_all
      norm1_div = 0.0_wp; norm2_div = 0.0_wp
      norm1_div_all = 0.0_wp; norm2_div_all = 0.0_wp  

      div_all=0; divmax=0; imax=0; jmax=0; kmax=0

      ! ---------------
      ! calculate du/dx
      ! ---------------
      call pdfx(u(:,:,1:xsz(3)),u_xi(:,:,1:xsz(3)),pex)
      call update_ghost(u_xi,level)
      call calc_uzeta(u, u_zeta(:,:,1:), level)
      call update_ghost(u_zeta,level)
      if (isbot) then
         u_zeta(:,:,2)=(u(:,:,2)-u(:,:,1))/dz(1)
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-1)=(u(:,:,xsz(3)-1)-u(:,:,xsz(3)-2))/2.0_wp/dz(xsz(3)-2)
      end if
      
      do k=1, xsz(3)
         divu(:,:,k)=u_xi(:,:,k)+zetax(:,:,k)*u_zeta(:,:,k)
      end do
      
      ! ---------------
      ! calculate dv/dy
      ! ---------------
      call pdfy_x(v(:,:,1:xsz(3)),u_psi(:,:,1:xsz(3)),pey)
      call update_ghost(u_psi,level)
      call calc_uzeta(v, u_zeta(:,:,1:), level)
      call update_ghost(u_zeta,level)
      if (isbot) then
         u_zeta(:,:,2)=(v(:,:,2)-v(:,:,1))/dz(1)
      end if
      if (istop) then
         u_zeta(:,:,xsz(3)-1)=(v(:,:,xsz(3)-1)-v(:,:,xsz(3)-2))/2.0_wp/dz(xsz(3)-2)
      end if

      do k=1, xsz(3)
         divu(:,:,k)=divu(:,:,k)+u_psi(:,:,k)+zetay(:,:,k)*u_zeta(:,:,k)
      end do

      ! ---------------
      ! calcualte dw/dz
      ! ---------------
      do k=1, xsz(3)
         u_zeta(:,:,k)=2.0_wp*(w(:,:,k)-w(:,:,k-1))/(dz(k-1)+dz(k))
      end do 
      call update_ghost(u_zeta,level)
      if (isbot) then 
         u_zeta(:,:,2)=(w(:,:,2)-w(:,:,1))/dz(2)
      end if 
      if (istop) then
         u_zeta(:,:,xsz(3)-1)=(w(:,:,xsz(3)-1)-w(:,:,xsz(3)-2))/dz(xsz(3)-2)
      end if

      do k=1, xsz(3)
         divu(:,:,k)=divu(:,:,k)+her(:,:)*u_zeta(:,:,k)
      end do
      
      call dealiasxy(divu(:,:,1:xsz(3)))

      ! if (istop)then 
      !     call dealiasxy(divu(:,:,xsz(3)+1))
      ! end if 

      if (istop.and.isbot)then 
          do k=2, xsz(3)-1
             do j=1, xsz(2)
                do i=1, xsz(1)
                   if (abs(divu(i, j, k)).gt.divmax)then 
                       divmax=abs(divu(i,j,k))
                       imax=i
                       jmax=j
                       kmax=k
                    end if 
                end do 
             end do 
          end do 

          !> added by plyu
          norm1_div = sum(abs(divu(:,:,2:(xsz(3)-1))))
          norm2_div = sum(divu(:,:,2:(xsz(3)-1))**2)

      else if (istop.and.(.not.isbot))then 
          do k=1, xsz(3)-1
             do j=1, xsz(2)
                do i=1, xsz(1)
                   if (abs(divu(i, j, k)).gt.divmax)then 
                       divmax=abs(divu(i,j,k))
                       imax=i
                       jmax=j
                       kmax=k
                   end if 
                end do 
             end do 
          end do 
          
          !> added by plyu
          norm1_div = sum(abs(divu(:,:,1:(xsz(3)-1))))
          norm2_div = sum(divu(:,:,1:(xsz(3)-1))**2)
      else if (isbot.and.(.not.istop))then 
          do k=2, xsz(3)
             do j=1, xsz(2)
                do i=1, xsz(1)
                   if (abs(divu(i, j, k)).gt.divmax)then 
                       divmax=abs(divu(i,j,k))
                       imax=i
                       jmax=j
                       kmax=k
                    end if
                end do 
             end do 
          end do 
          
          !> added by plyu
          norm1_div = sum(abs(divu(:,:,2:xsz(3))))
          norm2_div = sum(divu(:,:,2:xsz(3))**2)
       else 
          do k=1, xsz(3)
             do j=1, xsz(2)
                do i=1, xsz(1)
                   if (abs(divu(i, j, k)).gt.divmax)then 
                       divmax=abs(divu(i,j,k))
                       imax=i
                       jmax=j
                       kmax=k
                    end if
                end do 
             end do 
          end do 
          
          !> added by plyu
          norm1_div = sum(abs(divu(:,:,1:xsz(3))))
          norm2_div = sum(divu(:,:,1:xsz(3))**2)
        end if
      
      call mpi_allreduce(divmax,div_all,1,mpi_double_precision, &
           mpi_max,mpi_comm_world,ierr)

      call mpi_allreduce(norm1_div, norm1_div_all, 1, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)
      call mpi_allreduce(norm2_div, norm2_div_all, 1, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)
      norm1_div_all = norm1_div_all / nx / ny / (nz-2)
      norm2_div_all = sqrt(norm2_div_all / nx / ny / (nz-2))

      !> modifided by plyu
      if(myid.eq.0) then
         write(*, '(a,f10.5,3e11.3,3i5)') "Divergence_info: ",time, &
           norm1_div_all, norm2_div_all, div_all, imax,jmax,kmax
      endif

 end subroutine volume_lnr2
