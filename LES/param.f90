module param
   use decomp, only : wp
   implicit none

   private

   real(wp), public :: pex,pey,zl,hbar,dx,dy,xl,yl
   integer, public :: istart,iscalar,iturbine,iwcontrol
   real(wp), public :: z0
   real(wp), public :: dt
   integer, public :: ntime
   integer, public :: itmax
   real(wp), public :: erlim
   real(wp), public :: Fr2   
   real(wp), public :: RWe
   integer, public :: noutd, noutc
   !real(wp), public :: aka
   real(wp), public :: HKA,CRAT
   integer, public :: NWAVE
   integer, public :: iwavy
   real(wp),public :: timewavy,tcoef

   integer, public :: nx,ny,nz
   integer, public :: nxs,nys
   !real(wp), public :: time
   real(wp), public :: resbot,restop
   integer, public :: np1,np2
   real(wp), public :: arm
   REAL(wp), public :: ERVFILT
   INTEGER, public :: NFILT,IFILT
   integer, public :: ipa,nth
   real(wp), public :: timep,tcp,timew,rdgl,timeturb
   real(wp), public :: clbeta,clgama
   !real(wp), public :: ak,aa
   real(wp), public :: zlsbot,zlstop,usbot,ustop
   real(wp), public :: re
   real(wp), public :: bforce

   real(wp), public :: HA,HK,VPHASE,HOMEG
   
   !wave control parameter
   real(wp), public :: aka
   integer, public :: nswave

   !-- LES parameter --
   integer, public :: mfilt
   integer, public :: ICSC, ILASD, IVANDRIEST,IWAVEBOT
   real(wp), public :: TLASD, TLASD0
   real(wp), public :: cs0, aplus,ZCS0

   
 end module param
 
module constants
   use decomp, only : wp

   implicit none

   real(wp), parameter :: PI=3.141592653589793238_wp
   real(wp), parameter :: TWOPI = PI*2.0_WP

end module constants
