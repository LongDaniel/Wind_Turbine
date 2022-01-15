module wavecontrol
   use decomp
   use fft
   use param, only:  pex, pey, fr2, dt
   implicit none

   real(wp), public, save :: at1, aa, ak

contains

   function getp0(kx,pex,ky,pey,fr2,etahat,what,thramp,thrph) result(p0)
!     BY ANQING
!     This function return the pressure needed to attenuate
!     standing wave. If the wave is within the threshold or
!     the wave phase is not suitable for wave removal, the
!     pressure return will be 0.

      implicit none

      integer,  intent(IN) :: kx,ky
      real(wp), intent(IN) :: pex,pey,etahat,what,thramp,thrph,fr2
      real(wp) :: p0

      real(wp) :: a,s2,sigma,beta

      p0=0

      sigma = sqrt( ( kx*pex )**2 + ( ky*pey )**2 )
      sigma = sqrt( sigma/fr2 )

      if (kx==0 .or. ky==0) then
         beta = 2
      else
         beta = 4
      endif

      a = sqrt(etahat**2+(what/sigma)**2)
      if (beta*a > thramp) then
         ! compute cos(alpha)
         s2 = abs(what/sigma/a)
                 
         if (s2 > thrph) then
         ! add pressure
            p0 = sign(a,what)/sigma/fr2
         endif
      endif

   end function

!==========================================================================
   SUBROUTINE KILLSTANDINGWAVE_ALL(kx,ky,at,W,PA)
!C--BY ANQING-------------------------------------------------CC
!     DELTA PRESSURE
!C------------------------------------------------------------CC
      use MPI
      use utils
      use grid, only: eta
      !use param, only: pex, pey, fr2, dt

      implicit none
   !C--INPUT
      integer,  intent(IN) :: kx, ky
      !REAL(wp), intent(IN) :: at, FR2, DELTA, PEX, PEY
      REAL(wp), intent(IN) :: at

      !REAL(wp), dimension(xsz(1),xsz(2)), intent(IN) :: ETA, W
      REAL(wp), dimension(xsz(1),xsz(2)), intent(IN) ::  W

   !C--OUTPUT
      REAL(wp), dimension(xsz(1),xsz(2)), intent(INOUT) :: PA
   !C--INTERNAL
      integer :: i, j, nswavex, nswavey, kx1, kx2, ky1, ky2
      real(wp), dimension(xsz(1),xsz(2)) :: tmpx
      real(wp), dimension(ysz(1),ysz(2)) :: tmpy
      real(wp), dimension(nx_global, ny_global) :: PATMP, etaall, wsall
      ! The amplititude threshold for removing wave and the phase from
      ! which additional pressure is added
      real(wp) :: thramp, thrph

      thramp=0.010_wp*at; thrph=7/real(8,wp)

      tmpx = eta
      call fft_r2c_xy(tmpx, tmpy)
      tmpy = tmpy/nx_global/ny_global
      call transpose_yx(tmpy, tmpx)
      call gather_2d_xy(tmpx, etaall)

      tmpx = w
      call fft_r2c_xy(tmpx, tmpy)
      tmpy = tmpy/nx_global/ny_global
      call transpose_yx(tmpy, tmpx)
      call gather_2d_xy(tmpx, wsall)

      PATMP = 0

      ! Determine the wavenumbers to suppress
      kx1 = 0
      kx2 = nint(kx*1.8)
      ky1 = 0
      ky2 = ky+8

      DO NSWAVEX = kx1, kx2
         DO NSWAVEY = ky1, ky2
            if (nswavex==0 .and. nswavey==0) cycle
            if (nswavex==kx .and. nswavey==ky) cycle

            I = 2*NSWAVEX+1
            J = 2*NSWAVEY+1
            PATMP(I,J)=GETP0(NSWAVEX,PEX,NSWAVEY,PEY,FR2, &
                       ETAALL(I,J),WSALL(I,J),THRAMP,THRPH)/DT
            IF (istop .and. MYID1==0 .AND. PATMP(I,J)/=0) THEN
               PRINT *,'Remove standing wave at',I,J
            ENDIF

            I = 2*NSWAVEX+1
            J = 2*NSWAVEY+2
            PATMP(I,J)=GETP0(NSWAVEX,PEX,NSWAVEY,PEY,FR2, &
                       ETAALL(I,J),WSALL(I,J),THRAMP,THRPH)/DT
            IF (istop .and. MYID1==0 .AND. PATMP(I,J)/=0) THEN
               PRINT *,'Remove standing wave at',I,J
            ENDIF

            I = 2*NSWAVEX+2
            J = 2*NSWAVEY+1
            PATMP(I,J)=GETP0(NSWAVEX,PEX,NSWAVEY,PEY,FR2, &
                       ETAALL(I,J),WSALL(I,J),THRAMP,THRPH)/DT
            IF (istop .and. MYID1==0 .AND. PATMP(I,J)/=0) THEN
               PRINT *,'Remove standing wave at',I,J
            ENDIF

            I = 2*NSWAVEX+2
            J = 2*NSWAVEY+2
            PATMP(I,J)=GETP0(NSWAVEX,PEX,NSWAVEY,PEY,FR2, &
                       ETAALL(I,J),WSALL(I,J),THRAMP,THRPH)/DT
            IF (istop .and. MYID1==0 .AND. PATMP(I,J)/=0) THEN
               PRINT *,'Remove standing wave at',I,J
            ENDIF

            IF (NSWAVEX == 0) THEN
               PATMP(2*NSWAVEX+2,2*NSWAVEY+1)=0
               PATMP(2*NSWAVEX+2,2*NSWAVEY+2)=0
            ENDIF
            IF (NSWAVEY == 0) THEN
               PATMP(2*NSWAVEX+1,2*NSWAVEY+2)=0
               PATMP(2*NSWAVEX+2,2*NSWAVEY+2)=0
            ENDIF

         ENDDO
      ENDDO

      do j = 1,ysz(2)
         tmpy(:,j) = patmp(yst(2)+j-1,:)
      end do
      call fft_c2r_xy(tmpy, tmpx)

      pa = pa+tmpx

   END SUBROUTINE KILLSTANDINGWAVE_ALL

!==========================================================================
   SUBROUTINE KILLSTANDINGWAVE_ONE(NSWAVEX,NSWAVEY,eta, W,PA)
!C--BY XIN GUO------------------------------------------------CC
!     DELTA PRESSURE
!C------------------------------------------------------------CC
      use MPI
      use utils
      !use grid, only: eta
      !use param, only: pex, pey, fr2, dt

      implicit none
   !C--INPUT
      INTEGER,  intent(IN) :: NSWAVEX, NSWAVEY
      !REAL(wp), intent(IN) :: FR2, DELTA, PEX, PEY
      REAL(wp), dimension(xsz(1),xsz(2)), intent(IN) :: ETA, W
      !REAL(wp), dimension(xsz(1),xsz(2)), intent(IN) ::  W

   !C--OUTPUT
      REAL(wp), dimension(xsz(1),xsz(2)), intent(OUT) :: PA(xsz(1),xsz(2))
!C--INTERNAL
      INTEGER :: J
      real(wp), dimension(xsz(1),xsz(2)) :: tmpx
      real(wp), dimension(ysz(1),ysz(2)) :: tmpy
      REAL(wp), dimension(nx_global, ny_global) :: PATMP, etaall, wsall
      REAL(wp) :: SIGMA
      REAL(wp) :: S1, S2

      tmpx = eta
      call fft_r2c_xy(tmpx, tmpy)
      tmpy = tmpy/nx_global/ny_global
      call transpose_yx(tmpy, tmpx)
      call gather_2d_xy(tmpx, etaall)

      tmpx = w
      call fft_r2c_xy(tmpx, tmpy)
      tmpy = tmpy/nx_global/ny_global
      call transpose_yx(tmpy, tmpx)
      call gather_2d_xy(tmpx, wsall)

      PATMP = 0

      IF ( NSWAVEY /= 0 ) THEN
      
         SIGMA = sqrt( ( NSWAVEX * PEX )**2 + ( NSWAVEY * PEY )**2 )
         SIGMA = sqrt( SIGMA / FR2 )
      
         S1 = sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
                   + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) / SIGMA )**2 )
      
         S2 = sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+2)**2 &
                   + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+2) / SIGMA )**2 )
      
         IF ( S1 /= 0 ) THEN
            S1 = ETAALL(2*NSWAVEX+1,2*NSWAVEY+1) &
            / sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
            + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) &
            / SIGMA )**2 )
         
            IF ( abs(S1) < sqrt(3.0_wp) / 2 ) then
               S1 = ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) / SIGMA ) &
               / sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
               + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) &
               / SIGMA )**2 )
               IF ( S1 < -0.5_wp ) THEN
                  S1 = - sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
                  + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) &
                  / SIGMA )**2 )
               ELSE IF ( S1 > 0.5_wp ) THEN
                  S1 = sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
                  + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) &
                  / SIGMA )**2 )
               ELSE
                  S1 = 0
               END IF
            ELSE
               S1 = 0
            END IF
         ELSE
            S1 = 0
         END IF
      
         IF ( S2 /= 0 ) THEN
            S2 = ETAALL(2*NSWAVEX+1,2*NSWAVEY+2) &
            / sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+2)**2 &
            + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+2) &
            / SIGMA )**2 )
         
            IF ( abs(S2) < sqrt(3.0_wp)/2 ) then
               S2 = ( WSALL(2*NSWAVEX+1,2*NSWAVEY+2) &
               / SIGMA ) &
               / sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+2)**2 &
               + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+2) &
               / SIGMA )**2 )
               IF ( S2 < -0.5_wp ) THEN
                  S2 = - sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+2)**2 &
                  + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+2) &
                  / SIGMA )**2 )
               ELSE IF ( S2 > 0.5_wp ) THEN
                  S2 = sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+2)**2 &
                  + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+2) &
                  / SIGMA )**2 )**0.5
               ELSE
                  S2 = 0
               END IF
            ELSE
               S2 = 0
            END IF
         ELSE
            S2 = 0
         END IF
      
         PATMP(2*NSWAVEX+1,2*NSWAVEY+1) = S1 / FR2 / SIGMA  &
         / 4 / DT
         PATMP(2*NSWAVEX+1,2*NSWAVEY+2) = S2 / FR2 / SIGMA  &
         / 4 / DT
         PATMP(2*NSWAVEX+2,2*NSWAVEY+1) = S2 / FR2 / SIGMA  &
         / 4 / DT
         PATMP(2*NSWAVEX+2,2*NSWAVEY+2) = - S1 / FR2 / SIGMA  &
         / 4 / DT
      
         IF ( NSWAVEX == 0 ) THEN
            PATMP(2*NSWAVEX+1,2*NSWAVEY+1) = 2 &
            * PATMP(2*NSWAVEX+1,2*NSWAVEY+1)
            PATMP(2*NSWAVEX+2,2*NSWAVEY+1) = 0
            PATMP(2*NSWAVEX+1,2*NSWAVEY+2) = 2 &
            * PATMP(2*NSWAVEX+1,2*NSWAVEY+2)
            PATMP(2*NSWAVEX+2,2*NSWAVEY+2) = 0
         END IF
      ELSE
         SIGMA = sqrt( ( NSWAVEX * PEX )**2 + ( NSWAVEY * PEY )**2 )
         SIGMA = sqrt( SIGMA / FR2 )
      
         S1 = sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
         + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) &
         / SIGMA )**2 )
      
         S2 = sqrt( ETAALL(2*NSWAVEX+2,2*NSWAVEY+1)**2 &
         + ( WSALL(2*NSWAVEX+2,2*NSWAVEY+1) &
         / SIGMA )**2 )
      
         IF ( S1 /= 0 ) THEN
            S1 = ETAALL(2*NSWAVEX+1,2*NSWAVEY+1) &
            / sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
            + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) &
            / SIGMA )**2 )
         
            IF ( abs(S1) < sqrt(3.0_wp)/2 ) then
               S1 = ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) / SIGMA ) &
               / sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
               + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) &
               / SIGMA )**2 )
               IF ( S1 < -0.5_wp ) THEN
                  S1 = - sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
                  + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) &
                  / SIGMA )**2 )
               ELSE IF ( S1 > 0.5_wp ) THEN
                  S1 = sqrt( ETAALL(2*NSWAVEX+1,2*NSWAVEY+1)**2 &
                  + ( WSALL(2*NSWAVEX+1,2*NSWAVEY+1) &
                  / SIGMA )**2 )
               ELSE
                  S1 = 0
               END IF
            ELSE
               S1 = 0
            END IF
         ELSE
            S1 = 0
         END IF
      
         IF ( S2 /= 0 ) THEN
            S2 = ETAALL(2*NSWAVEX+2,2*NSWAVEY+1) &
            / sqrt( ETAALL(2*NSWAVEX+2,2*NSWAVEY+1)**2 &
            + ( WSALL(2*NSWAVEX+2,2*NSWAVEY+1) &
            / SIGMA )**2 )
         
            IF ( abs(S2) < sqrt(3.0_wp)/2 ) then
               S2 = ( WSALL(2*NSWAVEX+2,2*NSWAVEY+1) / SIGMA ) &
               / sqrt( ETAALL(2*NSWAVEX+2,2*NSWAVEY+1)**2 &
               + ( WSALL(2*NSWAVEX+2,2*NSWAVEY+1) &
               / SIGMA )**2 )
               IF ( S2 < -0.5_wp ) THEN
                  S2 = - sqrt( ETAALL(2*NSWAVEX+2,2*NSWAVEY+1)**2 &
                  + ( WSALL(2*NSWAVEX+2,2*NSWAVEY+1) &
                  / SIGMA )**2 )
               ELSE IF ( S2 > 0.5_wp ) THEN
                  S2 = sqrt( ETAALL(2*NSWAVEX+2,2*NSWAVEY+1)**2 &
                  + ( WSALL(2*NSWAVEX+2,2*NSWAVEY+1) &
                  / SIGMA )**2 )
               ELSE
                  S2 = 0
               END IF
            ELSE
               S2 = 0
            END IF
         ELSE
            S2 = 0
         END IF
      
         PATMP(2*NSWAVEX+1,2*NSWAVEY+1) = S1 / FR2 / SIGMA / 4 &
         / DT * 2
         PATMP(2*NSWAVEX+1,2*NSWAVEY+2) = 0
         PATMP(2*NSWAVEX+2,2*NSWAVEY+1) = S2 / FR2 / SIGMA / 4 &
         / DT * 2
         PATMP(2*NSWAVEX+2,2*NSWAVEY+2) = 0
      END IF

      do j = 1,ysz(2)
         tmpy(:,j) = patmp(yst(2)+j-1,:)
      end do
      call fft_c2r_xy(tmpy, tmpx)

      pa = tmpx

   END SUBROUTINE KILLSTANDINGWAVE_ONE

!==========================================================================
   SUBROUTINE MAINTAIN3(ETA,AMP,AT,SIGMA,AK,PA)
!C--BY XIN GUO------------------------------------------------CC
!     MAINTAIN WAVE WITH ENERGY
!C------------------------------------------------------------CC
      use MPI
      use constants, only : TWOPI
      !use param, only: pex, pey, dt

      IMPLICIT NONE

      !REAL(wp) :: AMP, AK, AT, SIGMA, PEX, PEY, DT
      
      REAL(wp) :: AMP, AK, AT, SIGMA      
      REAL(wp), dimension(xsz(1),xsz(2)), intent(IN)    :: ETA
      REAL(wp), dimension(xsz(1),xsz(2)), intent(INOUT) :: PA
   !C--INTERNAL
      INTEGER :: J, NSWAVEX
      REAL(wp) :: TMP, TMPALL, XL, YL
      REAL(wp) :: TMP1, TMP2, TMP1ALL, TMP2ALL
      REAL(wp), dimension(xsz(1),xsz(2)) :: etak, PPA

      XL = TWOPI / PEX
      YL = TWOPI / PEY

      NSWAVEX = NINT(AK/PEX)
      call fft_r2c_x(eta, etak)
      etak = etak/nx_global

      TMP = 0
      TMP1 = 0
      TMP2 = 0
      DO J = 1, xsz(2)
         TMP = TMP + sqrt( ETAk(2*NSWAVEX+1,J)**2 + ETAk(2*NSWAVEX+2,J)**2 )
         TMP1 = TMP1 + ETAk(2*NSWAVEX+1,J)
         TMP2 = TMP2 + ETAk(2*NSWAVEX+2,J)
      END DO
   ! C
      CALL MPI_ALLREDUCE(TMP,TMPALL,1,real_type,MPI_SUM, &
         MPI_COMM_2D_COL,j)
      CALL MPI_ALLREDUCE(TMP1,TMP1ALL,1,real_type,MPI_SUM, &
         MPI_COMM_2D_COL,j)
      CALL MPI_ALLREDUCE(TMP2,TMP2ALL,1,real_type,MPI_SUM, &
         MPI_COMM_2D_COL,j)
   !
      TMPALL = TMPALL / ny_global
      TMP1ALL = TMP1ALL / ny_global
      TMP2ALL = TMP2ALL / ny_global

      AMP = AMP / TMPALL / AT / SIGMA / DT * 2 / XL / YL / 2

   !C--TEST
   !      WRITE(*,*) AMP, TMPALL, AT, SIGMA, DT, XL , YL
   !      STOP
   !C--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   !C---STANDING WAVE
      ppa = 0
      ppa(2*nswavex+1,:) = -amp*tmp2all
      ppa(2*nswavex+2,:) = amp*tmp1all

      call fft_c2r_x(ppa)
      pa = pa+ppa

   END SUBROUTINE MAINTAIN3

!==========================================================================
   SUBROUTINE DIVTRAVSTAND(NX,NY,us,vs,ws,AS,AT,SI,CO, &
   ETATMP,ETTMP,ETAT,ETT)
!     BY XIN GUO
      use MPI
      use grid, only : eta, ehx, ehy
      !use param, only: pex, pey, fr2
    
      implicit none

      !C--INPUT
      INTEGER :: NX, NY
      !REAL(wp) :: PEX, PEY, FR2
      !REAL(wp), dimension(xsz(1),xsz(2)), intent(IN) :: ETA, us, vs, ws
      REAL(wp), dimension(xsz(1),xsz(2)), intent(IN) :: us, vs, ws

   !C--OUTPUT
      REAL(wp) :: AS, AT, SI, CO
      REAL(wp), dimension(xsz(1),xsz(2)) :: ETATMP, ETTMP, ETAT, ETT
   !C--INTERNAL
      INTEGER :: J
      real(wp), dimension(ysz(1),ysz(2)) :: tmpy
      REAL(wp) :: SIGMA
      REAL(wp) :: BR, BI

      sigma = sqrt( sqrt( ( nx * pex )**2 + ( ny * pey )**2 ) / Fr2 )

      etatmp = eta
      call fft_r2c_xy(etatmp, tmpy)
      tmpy = tmpy/nx_global/ny_global
      tmpy(1:2*ny,:) = 0; tmpy(2*ny+3:,:) = 0
      ! do j=1, ysz(2)
      !    if (yst(2)-1+j /= 2*nx+1 .and. yst(2)-1+j /= 2*nx+2) tmpy(:,j) = 0
      ! end do
      call fft_c2r_y(tmpy)
      call transpose_yx(tmpy, etatmp)
      etatmp(1:2*nx,:) = 0; etatmp(2*nx+3:,:) = 0

      ettmp = ws-us*ehx-vs*ehy
      call fft_r2c_xy(ettmp, tmpy)
      tmpy = tmpy/nx_global/ny_global
      tmpy(1:2*ny,:) = 0; tmpy(2*ny+3:,:) = 0
      ! do j=1, ysz(2)
      !    if (yst(2)-1+j /= 2*nx+1 .and. yst(2)-1+j /= 2*nx+2) tmpy(:,j) = 0
      ! end do
      call fft_c2r_y(tmpy)
      call transpose_yx(tmpy, ettmp)
      ettmp(1:2*nx,:) = 0; ettmp(2*nx+3:,:) = 0

      as = 2*sqrt( ( sigma*etatmp(2*nx+1,1) + ettmp(2*nx+2,1) )**2 &
               + ( sigma*etatmp(2*nx+2,1) - ettmp(2*nx+1,1) )**2 ) / abs(sigma)

      si =  2*( sigma*etatmp(2*nx+1,1) + ettmp(2*nx+2,1) )/sigma/as
      co = -2*( sigma*etatmp(2*nx+2,1) - ettmp(2*nx+1,1) )/sigma/as
   
      br = 2*etatmp(2*nx+1,1) - as / 2 * si
      bi = - 2*etatmp(2*nx+2,1) - as / 2 * co
   
      at = as / 2 + sqrt( br**2 + bi**2 )
   
      if ( abs(as-at*2) > 1d-12 ) then
         si = br / ( at - as / 2 )
         co = bi / ( at - as / 2 )
      else
         SIGMA = - SIGMA
         SI = sqrt( ETATMP(2*NX+1,1)**2 + ETATMP(2*NX+2,1)**2 )
         CO = ETATMP(2*NX+1,1) / SI
         SI = - ETATMP(2*NX+2,1) / SI
      end if
   
      ETATMP(2*NX+1,1) = ETATMP(2*NX+1,1) - AT * SI / 2
      ETATMP(2*NX+2,1) = ETATMP(2*NX+2,1) + AT * CO / 2
      ETTMP(2*NX+1,1) = ETTMP(2*NX+1,1) + AT * CO / 2 * SIGMA
      ETTMP(2*NX+2,1) = ETTMP(2*NX+2,1) + AT * SI / 2 * SIGMA
      etat = 0; ett = 0
      ETAT(2*NX+1,1) = AT * SI / 2
      ETAT(2*NX+2,1) = - AT * CO / 2
      ETT(2*NX+1,1) = -AT * CO / 2 * SIGMA
      ETT(2*NX+2,1) = -AT * SI / 2 * SIGMA
   
      DO J = 2, xsz(2)
         ETATMP(2*NX+1,J) = ETATMP(2*NX+1,1)
         ETATMP(2*NX+2,J) = ETATMP(2*NX+2,1)
         ETTMP(2*NX+1,J) = ETTMP(2*NX+1,1)
         ETTMP(2*NX+2,J) = ETTMP(2*NX+2,1)
         ETAT(2*NX+1,J) = ETAT(2*NX+1,1)
         ETAT(2*NX+2,J) = ETAT(2*NX+2,1)
         ETT(2*NX+1,J) = ETT(2*NX+1,1)
         ETT(2*NX+2,J) = ETT(2*NX+2,1)
      END DO

      call fft_c2r_x(etatmp)
      call fft_c2r_x(ettmp)
      call fft_c2r_x(etat)
      call fft_c2r_x(ett)

   END SUBROUTINE DIVTRAVSTAND

end module wavecontrol
