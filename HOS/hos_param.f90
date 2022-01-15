module hos_param
  
  use MPI
  
  implicit none
    
  !-- Constants --                                                                                  
  integer, parameter, public :: wp = kind(0.0d0)
  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  
  real(wp), parameter :: PI=acos(-1.0_wp)
  real(wp), parameter :: TWOPI = PI*2.0_WP
  
  !-- MPI parameters --                                                                             
  integer, public:: myid, numprocs, ierr
  
  !-- Numerical parameters --                                                                       
  integer, public :: nxhos, nyhos, ncpu_hos, npw
  real(wp), public :: pex_hos, pey_hos, dx_hos, dy_hos
  
  !-- Physical parameters --                                                                        
  real(wp), public :: fr2, g, bond
  
  !-- Execution control --                                                                          
  integer, public :: istart_hos, ntime_hos, noutd_hos, noutc_hos, ntp
  integer, public :: ioutd_hos, ioutc_hos, ist_hos, icap_hos
  real(wp), public :: dt_hos, time_hos
  
  !-- Wave parameters--                                                                             
  real(wp), public :: akax, akay, aka, ustar, u10, uss, gamma, fetch, phi
  integer, public :: nwave, iswell, nswellx, nswelly, nswavex, nswavey
  real(wp), public :: alpha, omega_p, lambda_p, omega0
  
  real(wp), allocatable, dimension(:,:), public :: eta_hos, ex_hos, ey_hos
  real(wp), allocatable, dimension(:,:), public :: vps_hos, vpsx_hos, vpsy_hos
  real(wp), allocatable, dimension(:,:,:), public :: feta, fvps
  real(wp), allocatable, dimension(:,:), public :: u_hos, v_hos, w_hos, pa_hos, pa_st
  real(wp), allocatable, dimension(:,:,:), public :: wvn, r_hos, zp
  
end module hos_param
