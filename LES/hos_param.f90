module hos_param
  
  use MPI
  use decomp, only:wp
  
  implicit none
    
  !-- Constants --                                                                                  
  ! integer, parameter, public :: wp = kind(0.0d0)
  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  
  ! real(wp), parameter :: PI=acos(-1.0_wp)
  ! real(wp), parameter :: TWOPI = PI*2.0_WP
  
  !-- MPI parameters --                                                                             
  ! integer, public:: myid, numprocs, ierr
  integer, public:: myid_hos, ierr_hos
  
  !-- Numerical parameters --                                                                       
  integer, public :: nxhos, nyhos, ncpu_hos, npw
  real(wp), public :: pex_hos, pey_hos, dx_hos, dy_hos
  
  !-- Physical parameters --                                                                        
  ! real(wp), public :: fr2, g, weber
  real(wp), public :: fr2_hos, g, weber
  
  !-- Execution control --                                                                          
  integer, public :: istart_hos, ntime_hos, noutd_hos, noutc_hos, ntp
  integer, public :: ioutd_hos, ioutc_hos, ist_hos
  real(wp), public :: dt_hos, time_hos
  
  !-- Wave parameters--                                                                             
  real(wp), public :: akax, akay, aka_hos, ustar, u10, uss, gamma, fetch, phi
  integer, public :: nwave_hos, iswell, nswellx, nswelly, nswavex, nswavey
  real(wp), public :: alpha, omega_p, lambda_p, omega0
  
  real(wp), allocatable, dimension(:,:), public :: eta_hos, ex_hos, ey_hos
  real(wp), allocatable, dimension(:,:), public :: vps_hos, vpsx_hos, vpsy_hos
  real(wp), allocatable, dimension(:,:,:), public :: feta, fvps
  real(wp), allocatable, dimension(:,:), public :: u_hos, v_hos, w_hos, pa_hos, pa_st
  real(wp), allocatable, dimension(:,:), public :: pa0_hos
  ! real(wp), allocatable, dimension(:,:,:), public :: wvn, r_hos, zp
  real(wp), allocatable, dimension(:,:,:), public :: wvn_hos, r_hos, zp_hos
  
end module hos_param
