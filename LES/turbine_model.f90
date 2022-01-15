!!< actuator disk model is grabed from Di Yang's Wind_LES_HOS_TURBINE code.

module turbine_model
  !> don't use navier as module, to avoid circular dependency
  !use navier, only : wtforce
  
  !use io_hos, only : nx, ny
  use param, only : pex, pey, xl, yl, zl, dx, dy, dt, nx, ny, nz, &
    hbar, dt, noutc, bforce, usbot
  use constants, only : pi, twopi
  use decomp, only : myid, wp, xsz, ysz, zsz, xst, isbot, istop,&
    MPI_COMM_2D_CART, MPI_COMM_2D_ROW, MPI_COMM_2D_COL, &
    nx_global, ny_global, nz_global, nproc1, ids, idsd, idsp
  use hos_param, only : ncpu_hos, nxhos, nyhos
  use grid, only : eta, hh, zz, zw, hx, hy, dz, her
  use discontinuity_smooth
  use mpi
  implicit none
  
  !< integer, save :: myid, status

  !< **************************************
  !< Following are Actuator Disk Model of Di Yang
  !< Transplanted by plyu 
  !< begin: 12-17-2016 modified: 12-28-2016
  !< **************************************

  !< nxmod: nxmod is a deprecated variable, not used in other modules.
  !< so nxmod is renamed as nxmod0 in this module. nxmod0 = nxmod  = nx-1
  integer, save :: nxmod0, nymod0

  integer :: ierr_wt

  !< nxwt: number of rows of wind turbine
  !< nywt: number of colums of wind turbine
  !< rdisc: radius of wind turbine disc
  !< central height of wind turbine disc
  integer, public, save :: nxwt, nywt, iwt
  real(wp), public, save :: rdisc, hdisc

  !< i_thrust: method of thrust force calculation
  !         1, based on undisturbed upstream time-space-averaged velocity
  !         2, based on undisturbed upstream real-time space-averaged velocity
  !         3, based on disturbed rotor-plane time-space-averaged velocity
  !         4, based on disturbed rotor-plane real-time space-averaged velocity
  ! c_thrust: thrust coefficient based on undisturbed upstream velocity
  ! c_thrust_induced: thrust coefficient based on disturbed rotor-plane velocity
  integer, public, save :: i_thrust
  real(wp), public, save :: c_thrust, c_thrust_induced

  !< ndisc: disc number index (which disc the grid point belongs to)
  !< idisc: =1 if inside a turbine disc, =0 if outside any turbine disc
  !<integer, public, save :: ndisc(nx, ny, xsz(3))
  !<integer, public, save :: idisc(nx, ny, xsz(3))
  integer, allocatable, dimension(:,:,:), public, save :: ndisc, idisc
  real(wp), allocatable, dimension(:,:,:), public, save :: wtgamma0
  integer, allocatable, dimension(:,:,:), public, save :: ndisc_ref, idisc_ref
  real(wp), allocatable, dimension(:,:,:), public, save :: wtgamma0_ref

  !< wtgamma0: fraction of area overlap between mesh cell and turbine disc
  !< wtgamma0 is without domain decomposition
  !< real(wp), public, save :: wtgamma0(nx, ny, xsz(3))

  !< ndp: number of grid point for each disc
  !< wndp: ndp weighted by wtgamma0
  integer, public, save :: ndp, ndp_ref
  real(wp), allocatable, dimension(:), public, save :: wndp, wndp_ref

  !< output index for y-z, x-z, x-y planes
  !< icut: y-z plane at one diameter downstream of the turbine
  !< jcut: x-z plane at crossing the center of the turbine
  !< kcut: x-y plane crossing the center of the turbine hub (disc center)
  !< integer, public, save :: icut(nxwt), jcut(nywt), kcut
  integer, allocatable, dimension(:) :: icut, jcut
  integer, allocatable, dimension(:) :: icut_ref, jcut_ref
  integer :: kcut
  
  !< ud: disk-averaged time-averaged local velocity
  !< ud1: disk-averaged local velocity
  !< real(wp), public, save :: ud(nxwt*nywt)
  real(wp), allocatable, dimension(:) :: ud, ud1
  real(wp), allocatable, dimension(:) :: ud_ref, ud1_ref
  
  
  !< zzo: collect zz from all cpu
  !< zz1: temporary value of zzo
  !< cartx, carty, cartz: cartisian coordinate
  real(wp), allocatable, dimension(:), public, save :: zzo, zz1, zwo, zw1, cartx, carty
  !real(wp), allocatable, dimension(:,:,:), public, save :: cartz, cartzw

  real(wp), allocatable, dimension(:,:), public, save :: eo, hho, hxo, hyo


  public :: read_ad_dy_param, reread_admnr
  public :: fraction_v2, discloca_v3, veldisc_initial
  public :: veldisc_v2, wind_turbine_force 
  
  !< ***********************************
  !< Following are Actuator Line model from VWiS
  !< Tranplanted by plyu
  !< begin: 12-28-2016 modified:
  !< update history:
  !< 12-28-2016, project begin
  !< 01-03-2017, finished reading input part
  !< 01-12-2017, finished serial version
  !< 01-18-2017, finished and tested AL v1.0
  !< 01-23-2017, added Actuator Surface Model
  !< 01-24-2017, added Nacelle Model
  !< 02-07-2017, add different kind of delta function
  !< 02-07-2017, check Calc_F_lagr_ACL step-by-step with hand calculation
  !< 02-07-2017, add background force to whole domain to compensate turbine force
  !< 02-08-2017, add interpolation method to line2surface related subroutines  
  !< 02-09-2017, add output for line and surface location in tecplot format(vtk next?)
  !< 02-09-2017, fixed a bug: when use n_elmt instead of ibm(ibi)%n_elmt for
  !              convenience, forget assign value for n_elmt
  !< 02-10-2017, add tip loss correction of Shen 
  !< 02-11-2017, add undisturbed incident angle correction based on PhD-Schito-2011
  !< 02-14-2017, add mesh refinement for turbine region in z direction (grid.f90)
  !< 02-20-2017, add temporary solution for restarting a turbine case, valid
  !              only for fixed angular velocity case
  !< 02-22-2017, add Actuator Disk Model with Rotation, without Torque control
  !< 02-23-2017, add output for Power Coefficient
  !< 02-24-2017, validtion of Clipper turbine is successful, in power
  !              coefficient
  !< 02-27-2017, fixed a bug in Calc_F_eul. The height of a grid is not dz(k),
  !              it should be dz(k)*HBAR=dz(k)/her(i,j). After this fix, cases
  !              of high TSR get reasonable power coefficient. Thrust need
  !              improvement.
  !< 03-03-2017, add velocity restoration zone at upstream (add a force together
  !              with turbine force, before prediction step)
  !< 03-09-2017, add first version of post-process for hos-les-turbine case
  !< 03-10-2017, add another version of velocity restoration zone (after
  !              prediction step)
  !< 03-10-2017, add a backdoor for turbulent flow pre-generation
  !              ( set turbine_model=-100 to enable it)
  !< 03-15-2017, add a grid generation method, which allow refinement in inner
  !              layer (However when y+ < 20, it will blow up)
  !< 03-25-2017, improve post-process for turbine case, so that it can output
  !              both 1D line and 2D plane information, together with mean
  !              value.
  !< 03-26-2017, modify background force to compensate turbine drag to constant
  !              value
  !< 03-30-2017, fix some bugs in Di Yang's actuator disk model
  !< 03-31-2017, add restart ability to ADM-NR 
  !< 04-03-2017, add disk-average velocity output for disk-plane and 2D upstream
  !< 04-04-2017, add different ways to calculate thrust in ADM-NR; output_inlet
  !              is 3D mode now.
  !< 04-05-2017, Fix some bugs to make Nacelle Model get its fist run
  !< 04-06-2017, Add a adjust_CT process in both rotor model and nacelle model,
  !              to force CT to be prescribed value.
  !< 04-21-2017, Add more inlet import I-O interface. It can read results from 
  !              precursor simulation and synthetic turbulent time series now.
  !< 04-24-2017, Add consideration for extrapolation of CD and CL table
  !< 04-27-2017, Fix a vital bug in variable transfer of wtforce
  !< 05-03-2017, Modify Eulerian-Lagrangian convertion to simpler way (correct?)
  !< 05-11-2017, Angular velocity now changable. Fix bug in rotro_Rot.
  !< 06-14-2017, hos-les couple: fix bugs of grid sending from hos to les
  !< 06-27-2017, Add prescribe motion of nacelles and rotors.
  !< 07-04-2017, Add time-averaged U_ref to stablize Turbine angular speed 
  
  !< Todo: 
  !< 0. output hos info at cout
  !< 0. Extract time-average function to a separate one so that we can call it
  !     for many other variables, such as Uref, wind_direction
  !< 0. Add export for Shear force and Moment.
  !< 0. Turbine control algorithm: pitch, angular speed, yaw.
  !< 0. Inflow import from WRF or LiDAR. Consider Jie's adjoint model.
  !< 0. change adjust_CT to base it on a time-averaged result. Make the thrust
  !     vary smoothly.
  !     are distributed unevenly in post-process results
  !< 0. add CT and CP to output file
  !< 2. postprocess module
  !< 3. have a look at FEM code
  !< 4. turbine platform motion(only wind force)
  !< 5. unsteady AL model
  !< 6. turbine platform motion(wave and mooring force)
  !< 7. Nacelle Model with wall model


  !<   mainly refer to main.c and rotor_model.c in VWiS
  !< ***********************************
  
  !< **** variables read from control.dat
  !< rotor_model: 0, inactive; 1, AD with induction factor;
  !<   3, AL; 4, AD with thrust coefficient; 5, AS; 
  !<   6, AL with an additional line for computing U_ref
  !<   7, ADM-Rotation
  integer, public, save :: rotor_model, NumberOfTurbines, num_blade, num_foiltype
  real(wp) :: loc_refvel
  integer :: FixTipSpeedRatio, FixTurbineAngVel, URefDynamicAverage

  real(wp) :: halfwidth_dfunc
  integer :: Itpwidth
  integer :: deltafunc_U, deltafunc_F
  real(wp), save :: sqrt_pi

  !< corrections for rotor model
  integer :: correction3D, correction3D_CH, correction3D_DS
  real(wp) :: c1_CH, c2_CH, c3_CH, CD_0

  !< Tip loss correction
  integer :: Shen1_AL
  real(wp) :: Shen1_AL_tipcorrection, Shen1_AL_tipcorrectionratio_Fa
  
  !< Inlet velocity forcing
  !< InletRelaxation: 0: disabled.
  !<   1: nudge after prediction, in navier_fu.f90:Add_Velocity_Relaxation
  !<   2: nudge before prediction, in turbine_model.f90:Add_Inlet_Forcing
  !<   3: nudge after prediction, using ?
  !<   4: nudge after prediction, using imported inlet data
  !< IR_dataindex: index of inlet data source
  !< IR_scale: rescale the velocity
  integer, public, save :: InletRelaxation, IR_dataindex, IR_bw
  integer, public, save :: IR_ti_shift
  real(wp), public, save :: IR_scale 

  !< c_bforce is used to add an additional Background Force to keep the flow
  !  rate, normally it is set to be 4/3 which is thrust coefficient.
  !  if Umean_hub upstream is not your desired value, adjust c_bforce.
  real(wp), public, save :: c_bforce
 
  !< whether to adjust CT or CP to desired value
  integer :: adjust_CT_ACL, adjust_CP_ACL
  integer :: adjust_CT_nacelle

  !< whether to remove lateral ACL force
  integer :: remove_lateral_ACL_force
  
  !< whether to output info at certain wake locations
  integer :: output_wake_switch, output_wake_step
  real(wp) :: wake_rel_xmin, wake_rel_xmax, wake_rel_y
   
  !< nacelle_model: 
  !< 7: Calc_Nut_eul(Nut_lagr, lnut_eul), fx_t=0.0
  !< 1: force in simple IB formula.
  !< 4: ShearDesired, frictionfactor(l), Comput_desiredshear_nacelle1
  !< 5: ShearDesired, friction_factor,
  integer :: nacelle_model, rotate_nacelle
  integer :: NumberOfNacelle, NumNacellePerLoc, NumLoc

  !< FSI
  integer :: fsitype
  integer, dimension(6) :: fsi_dof
  real(wp), dimension(6) :: fsi_eqpos, fsi_amp, fsi_omega
  
  real(wp) :: dA_mean, dh_fsi, dh_uref, dh_mean
  real(wp), public, save :: first_step_time
  integer, public, save :: ti, ti_first

  !< **** variable read from Turbine.inp
  !< variable definition is grabbed from variable.h
  type FSInfo
    real(wp) :: nx_tb, ny_tb, nz_tb !< direction vector of rotor axis
    real(wp) :: nx_tb0, ny_tb0, nz_tb0 !< original direction vector of rotor axis
    real(wp), dimension(3) :: angvel
    
    !< ang_axia: real-time angle. angvel_axis: Omega
    !< ang0: initial angle
    real(wp) :: angvel_axis, ang_axis, ang0_yaw, ang0_rot, angvel0
    real(wp) :: x_c, y_c, z_c
    real(wp) :: r_rotor
    real(wp) :: Tipspeedratio, angvel_fixed

    !< desired value of thrust coefficient CT0, and power coeff CP0
    real(wp) :: CT0, CP0, area_CT, ratio_CT0

    real(wp) :: Torque_aero, Force_axis

    real(wp) :: x_c0, y_c0, z_c0

    integer :: rotate_alongaxis
    real(wp) :: xnacelle_upstreamend, ynacelle_upstreamend, znacelle_upstreamend

    integer :: rotor_model_type

    !< six degree location
    real(wp) :: x_co, y_co, z_co
    real(wp), dimension(6) :: disp, vel
    
    !< translation displacement and velocity
    !real(wp), dimension(6) :: S_new, S_old
    !real(wp), dimension(6) :: S_ang_n, S_ang_o 
  end type

  !< **** variable read from grid file in UCD format
  !< such as UrefData%03d, acldata%03d, acsdata%03d
  type IBMNodes
    !> number of ?
    integer :: nbnumber
    !> number of vertices and elements
    integer :: n_v, n_elmt
    !> node index of each triangle
    integer, allocatable, dimension(:) :: nv1, nv2, nv3
    !> normal direction
    real(wp), allocatable, dimension(:) :: nf_x, nf_y, nf_z
    real(wp), allocatable, dimension(:) :: ns_x, ns_y, ns_z
    real(wp), allocatable, dimension(:) :: nt_x, nt_y, nt_z
    !> coordinate of IBM surface nodes
    !> _i is initial value，
    !> _o is offseted value with x_c in Turbine.inp
    real(wp), allocatable, dimension(:) :: x_bp, y_bp, z_bp
    real(wp), allocatable, dimension(:) :: x_bp_i, y_bp_i, z_bp_i
    real(wp), allocatable, dimension(:) :: x_bp_o, y_bp_o, z_bp_o
    real(wp), allocatable, dimension(:) :: x_bp0, y_bp0, z_bp0

    !>center of the element
    real(wp), allocatable, dimension(:) :: cent_x, cent_y, cent_z 
    real(wp), allocatable, dimension(:,:) :: u, uold, urm1
    !> area of an element
    real(wp), allocatable, dimension(:) :: dA  
    !> force at the IB surface points (lagrange points)
    real(wp), allocatable, dimension(:) :: F_lagr_x, F_lagr_y, F_lagr_z
    real(wp), allocatable, dimension(:) :: U_lagr_x, U_lagr_y, U_lagr_z
    
    !> point index for quick search
    integer, allocatable, dimension(:) :: i_min, j_min, k_min
    integer, allocatable, dimension(:) :: i_max, j_max, k_max
    integer, allocatable, dimension(:) :: iclose, jclose, kclose
    !> for Interpolation Points in wall model
    integer, allocatable, dimension(:) :: iIP_min, jIP_min, kIP_min
    integer, allocatable, dimension(:) :: iIP_max, jIP_max, kIP_max
    
    !> relative incoming velocity for actuator model
    !> Urel: vector of the relative incoming velocity
    !> Uinduced: vector of the induced velocity
    !> circulation: circulation vector on the blade
    real(wp), allocatable, dimension(:,:)::Urel,Uinduced,circulation
    real(wp), allocatable, dimension(:,:)::liftdirection,rotationdirection

    real(wp), allocatable, dimension(:) :: Urelmag, Urelmag_mean
    real(wp), allocatable, dimension(:,:)::Urel_mean
    real(wp), allocatable, dimension(:,:)::Uinduced_mean, circulation_mean
    real(wp), allocatable, dimension(:,:)::liftdirection_mean

    real(wp), allocatable, dimension(:) :: AOA_mean
    
    !> color: used to identify which blade does a AL/AS element belongs to
    integer, allocatable, dimension(:) :: color


    real(wp), allocatable, dimension(:) :: angle_attack, angle_twist, chord_blade 

    !> CD2 and CL2 are 2D coefficients. CD might be equal to CD2, or a 3D corrected value.
    real(wp), allocatable, dimension(:) :: CD, CL, CD2, CL2
    !> maximum number of blades: 3
    real(wp), dimension(3) :: pitch

    !> U1_ref is instantaneous disk-average reference velocity upstream.
    !> U_ref is dynamic time average of U1_ref
    real(wp) :: U_ref, U1_ref, indf_axis, Tipspeedratio, CT, indf_tangent

    !> TSR_modified is used in Du's 3D correction
    real(wp) :: TSR_modified

    !> s2l: actuator line element index for each actuator surface element
    integer, allocatable, dimension(:) :: s2l
    integer, allocatable, dimension(:) :: s2l1, s2l2
    real(wp), allocatable, dimension(:) :: s2lc

    ! Fr_mean, Ft_mean, Fa_mean

    real(wp) :: friction_factor, dh
    real(wp), allocatable, dimension(:) :: dh_IP
    real(wp), allocatable, dimension(:,:) :: centIP, U_IPlagr
  end type
   
  !< **** variable read from FOIL%02d, CD%02d, CL%02d for actuator line model
  type ACLInfo
    !< num_AC: number of input points for angle of pitch and chord length
    integer :: num_AC, num_CD, num_CL
    !< pitch angle and chord length from input file
    real(wp), allocatable, dimension(:) :: r_AC, angle_twistInp, chord_bladeInp
    real(wp), allocatable, dimension(:) :: ang_CD, CDInp, ang_CL, CLInp
    !> this type airfoil begins at r=r_beg and ends at r=r_end
    real(wp) :: r_beg, r_end  
    !>
    real(wp) :: angle_zero_lift, CD_zero_angle
  end type

  !< Parameters for prescribed turbine motion
  type FSIPrePara
    integer, dimension(6) :: dof
    real(wp), dimension(6) :: refpos, amp, omega, phase0, phase
    real(wp), dimension(6) :: disp, vel, dd, disp_o, disp0 
  end type

  type(FSInfo), allocatable, dimension(:), public, save :: fsi_wt, fsi_acl2ref, fsi_nac
  type(IBMNodes), allocatable, dimension(:), public, save :: wtm, ibm_ACD, ibm_acl2ref,ibm_nac
  type(ACLInfo), allocatable, dimension(:), public, save :: acl
  type(FSIPrePara), allocatable, dimension(:), public, save :: pp_wt

  !> read input:
  public :: read_turbine_model_param
  public :: read_acl_param, read_acs_param, read_nacelle_param
  public :: ACL_read_ucd, disk_read_ucd, airfoil_ACL, surface_read_xpatch
  public :: read_turbine_control, read_nacelle_control

  !> top level subroutines:
  public :: rotor_model_acl, rotor_model_acs, nacelle_model_run
  
  !> runtime:
  public :: Pre_process, Calc_U_lagr, Uref_ACL, Calc_turbineangvel
  public :: rotor_Rot, Calc_F_lagr_ACL, Calc_forces_ACL, Calc_F_eul
  public :: Prescribe_para_update, Prescribe_FSI_update

  !> for AS model:
  public :: calc_s2l, ForceProjection_l2s

  !> for nacelle model:
  public :: Coordinates_IP, Calc_F_lagr_nacelle
  public :: Calc_forces_nacelle, nacelle_move, Calc_nacelle_angvel

  !> miscellaneous:
  public :: ArbitraryRotate
  public :: dfunc_s1h, dfunc_s2h, dfunc_s3h, dfunc_s4h
  public :: calc_vol_flux

  !> controller:
  public :: controller_pid
  real(wp) :: kc_p, kc_i, kc_d, sp_bforce
  real(wp) :: err_0, err_1, err_2
  integer :: controller_n

  !> output:
  public :: output_turbine, Export_Wake
  public :: Export_LineLocation, Export_SurfaceLocation
  public :: Export_ForceOnBlade
  public :: Read_LineLocation, Read_SurfaceLocation, Read_NacelleLocation
  
contains
  subroutine dotx(a,b,c)
    implicit none
    real(wp), dimension(3), intent(in) :: a,b
    real(wp), intent(out) :: c

    c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end subroutine dotx
  
  subroutine crossx(a,b,c)
  implicit none
    real(wp), dimension(3), intent(in) :: a,b
    real(wp), dimension(3), intent(out) :: c
!    real(wp) :: a,b
    
    if(size(a,1).eq.1 .or. size(b,1).eq.1) then
      print *, 'plyudebug, crossx: Incorrect size of array'
    endif
 !   if(size(a0,1).eq.1) a = reshape(a0,(/3,1/))
  !  if(size(b0,1).eq.1) b = reshape(b0,(/3,1/))
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end subroutine crossx

  !> interface for turbine model read
  subroutine read_turbine_model_param(iturbine_, time_)
    implicit none
    integer, intent(in) :: iturbine_
    real(wp), intent(in) :: time_

    !> enable the Smoothing for Discontinuity, see smooth.f90 and
    !  spectral.f90 for details. The turbines will brought a sharp jump to
    !  the fluid field, which will lead to Gibbs osscilation during the spectral
    !  difference process. We treat this jump as a discontinuity and apply a
    !  smoothing technique to it.
    ids = 0 
    if (ids .ne. 0 .and. myid .eq. 0) then
      print *, 'ids =', ids, '; Smoothing for pdfx, pdfy_x is enabled.'
    endif

    !> plyunote: as I can see, when turbine is inside the domain, it is necessary
    !! to set it as 2. This setting(idsd) is more effective than other two (ids,
    !! idsp) to alleviate the pressure oscillation. When idsd equals to 0, the
    !! original pressure oscillation spreads into the whole domain. After
    !! setting idsd to 2, the pressure oscillation is controlled in the near
    !! wake. How about setting idsd to 1? My test shows it would blow out.
    idsd = 0
    if (idsd .ne. 0 .and. myid .eq. 0) then
      print *, 'idsd =', idsd,'; Smoothing for dealiasxy is enabled.'
    endif

    idsp = 0 
    if (idsp .ne. 0 .and. myid .eq. 0) then
      print *, 'idsp =', idsp,'; Smoothing for poisson equation is enabled.'
    endif

    if (ids.ne.0 .or. idsd.ne.0 .or. idsp.ne.0) then
      call discontinuity_smooth_init
    endif
    
    if(iturbine_ .eq. 1) then
      call read_ad_dy_param
      call veldisc_initial(time_)
      call reread_admnr
    else if(iturbine_ .eq. 3) then
      call read_acl_param(time_)
    else if(iturbine_ .eq. 5) then
      call read_acs_param(time_)
    else if(iturbine_ .eq. 7) then
      call read_admr_param(time_)
    endif

  end subroutine read_turbine_model_param
   
 !> print info about node of wtm
  subroutine print_ibmnode1 (ibi, nv, tag)
    implicit none
    integer :: ibi, nv, tag
    
    if (ibi.eq.7 .and.  nv .le. 2) then
      print *, 'ibmnode1, ',tag,myid, wtm(ibi)%x_bp(nv), wtm(ibi)%y_bp(nv), &
        wtm(ibi)%z_bp(nv), ', origin:', fsi_wt(ibi)%x_c, fsi_wt(ibi)%y_c, &
        fsi_wt(ibi)%z_c, ', axis:', fsi_wt(ibi)%nx_tb, fsi_wt(ibi)%ny_tb, &
        fsi_wt(ibi)%nz_tb
    endif
  end subroutine print_ibmnode1
   
  subroutine output_turbine
    implicit none
    if (rotor_model .eq. 3) then
      call Export_LineLocation(wtm, NumberOfTurbines)
      call Export_ForceOnBlade(wtm, fsi_wt, NumberOfTurbines, 'ForceAL_')
      if (nacelle_model .ne. 0) then
        call Export_NacelleLocation(ibm_nac, NumberOfNacelle)
      endif
    else if (rotor_model .eq. 5 .or. rotor_model .eq. 7) then
      call Export_LineLocation(ibm_acl2ref, NumberOfTurbines)
      call Export_SurfaceLocation(wtm, NumberOfTurbines)
      call Export_ForceOnBlade(wtm, fsi_wt, NumberOfTurbines, 'ForceAS_')
      call Export_ForceOnBlade(ibm_acl2ref, fsi_acl2ref, &
        NumberOfTurbines,'ForceAL_')
      if (nacelle_model .ne. 0) then
        call Export_NacelleLocation(ibm_nac, NumberOfNacelle)
      endif
    endif
  end subroutine output_turbine
 
  !> plyunote: I have two plans
  !! 1. simple plan, output all grid points along streamwise direction at
  !certain locations. The line is in the same horizontal plane with rotating
  !centers.
  !! 2. complex plan, extend a ray from rotating center along its axis
  !direction. As every trubine might have different yaw angle, the monitoring
  !points of different turbine might be imparallel
  !! as I don't have time, I prefer to implement the simple plan now.
  subroutine Export_Wake (fsi, NumberOfObjects, uin, vin, win, level)
    use lib_array, only: interp1d
    implicit none
    type(FSInfo), dimension(:) :: fsi
    integer :: NumberOfObjects, level
    real(wp), dimension(:,:,:) :: uin(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: vin(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: win(xsz(1), xsz(2), 1-level:xsz(3)+level)
    
    character(64) :: filen
    real(wp) :: loctmp(3), r1, r2, r3, r0
    real(wp), dimension(:,:), allocatable :: wline, wline_local
    integer :: indextmp(3), ibi, i, i1, i2, i3, i4, i5, i6, nxwake


    do ibi = 1, NumberOfObjects
      if(ti .eq. ti_first+1 .and. myid.eq.0) then
        write(filen, '(a9,i0.3,a4)') 'wakeline_', ibi, '.dat'
        open(18009, file=trim(filen), action='write')
        close(18009)
      endif
      
      if (mod(ti, output_wake_step) .eq. 0) then
        r0 = fsi(ibi)%r_rotor
        i1 = floor((fsi(ibi)%x_c0 + wake_rel_xmin * r0)/(xl/nx)) 
        i2 = floor((fsi(ibi)%x_c0 + wake_rel_xmax * r0)/(xl/nx))
        nxwake = i2-i1+1

        !> strategy: initialize the array at every processor, then assign value
        !if it lies in that processor, then sync and write.
        allocate(wline(nxwake,6), wline_local(nxwake,6))
        wline = -1.0d50; wline_local = -1.0d50

        do i = 1, nxwake
          i3 = i1 + i - 1
          if (i3>nx) i3 = i3 - nx
          loctmp(1) = i3 * (xl/nx) 
          loctmp(2) = fsi(ibi)%y_c0 + wake_rel_y * r0
          
          i4 = i3 - xst(1) + 1
          i5 = floor(loctmp(2) / (yl/ny)) - xst(2) + 1
        
          if (i5>=1 .and. i5<=xsz(2)) then 
            if (i4<xst(1) .or. i4>xsz(1)) then 
              print *, 'i1,i2,i3,i4,i5,x_c0, loctmp(1:2)=', i1, i2, i3, &
                i4, i5, fsi(ibi)%x_c0, loctmp(1:2)
            endif
            
            loctmp(3) = (fsi(ibi)%z_c0 + hh(i4,i5)) / (hbar+hh(i4,i5))
            if (zz(0)<=loctmp(3) .and. loctmp(3)<zz(xsz(3))) then
              
              !!> locate the index of loctmp(3) in zz
              !i6 = locate(zz, loctmp(3)) ! zz(1-level:xsz(3)+level)
              !print*, loctmp(3), zz

              wline_local(i,1) = (i3-1)*(xl/nx)
              wline_local(i,2) = (xst(2)+i5-1-1)*(yl/ny)
              wline_local(i,3) = fsi(ibi)%z_c0
              wline_local(i,4) = interp1d(zz, uin(i4,i5,:), loctmp(3))
              wline_local(i,5) = interp1d(zz, vin(i4,i5,:), loctmp(3))
            endif
            if (zw(0)<=loctmp(3) .and. loctmp(3)<zw(xsz(3))) then
              !print*, loctmp(3), zw
              wline_local(i,6) = interp1d(zw, win(i4,i5,:), loctmp(3))
            endif

          endif

        enddo
        
        call mpi_allreduce(wline_local, wline, nxwake*6, &
          mpi_double_precision, mpi_max, mpi_comm_world, ierr_wt)

        if (myid .eq. 0) then
          write(filen, '(a9,i0.3,a4)') 'wakeline_', ibi, '.dat'
          open(18009, file=trim(filen), action='write', position='append')
          do i = 1, nxwake
            write(18009, *) wline(i, 1:6)
          enddo
          close(18009)
        endif

        deallocate(wline, wline_local)
      endif
    enddo

  end subroutine Export_Wake
  
  !> Pre_process: pair ibm elements to nearest fluid grid, and get minimal band
  !! for search 
  subroutine Pre_process(ibm, NumberOfObjects, pretype)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    integer :: NumberOfObjects
    integer :: pretype ! 1 for odinary ibm nodes, 2 for IP of ibm nodes
    ! IP: Interpolation Point
    
    integer :: xs, xe, ys, ye, zs, ze
    integer :: i, j, k, l, ibi, n_elmt
    integer :: lxs, lxe, lys, lye, lzs, lze

    integer, allocatable, dimension(:) :: sum_iclose, sum_jclose, sum_kclose
    integer, allocatable, dimension(:) :: iclose, jclose, kclose
    integer, allocatable, dimension(:,:,:) :: ijkbound

    real(wp) :: dmin, r1, r2, r3, d1, tmpx, tmpy, tmpz
    integer :: ic, jc, kc, ic1, ic2, jc1, jc2, kc1, kc2

    real(wp), dimension(3) :: cent

    xs = xst(1); xe = xst(1) + xsz(1) - 1
    ys = xst(2); ye = xst(2) + xsz(2) - 1
    zs = xst(3); ze = xst(3) + xsz(3) - 1
    !lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze
    lxs = xs-Itpwidth; lxe = xe+Itpwidth
    lys = ys-Itpwidth; lye = ye+Itpwidth
    lzs = zs-Itpwidth; lze = ze+Itpwidth

    if (size(ibm,1) .ne. NumberOfObjects) then
      print *, "Error in Pre_process, array size doesn't match"
      return
    endif

    do ibi = 1,NumberOfObjects
      n_elmt = ibm(ibi)%n_elmt
      allocate(iclose(n_elmt),jclose(n_elmt),kclose(n_elmt))
      allocate(sum_iclose(n_elmt), sum_jclose(n_elmt))
      allocate(sum_kclose(n_elmt))
      
      ! index 1 is number of elements, index 2 is coordinate index(x,y,z)
      ! index 3 is identifier of min or max (1 for min, 2 for max)
      allocate(ijkbound(n_elmt,3,2))

      iclose = 0; jclose = 0; kclose = 0
      sum_iclose = 0; sum_jclose = 0; sum_kclose = 0
      
      do l = 1,n_elmt
        dmin = 1.0e6_wp
        if (pretype .eq. 1) then
          cent(1) = ibm(ibi)%cent_x(l); cent(2) = ibm(ibi)%cent_y(l)
          cent(3) = ibm(ibi)%cent_z(l)
        else if (pretype .eq. 2) then
          cent(1) = ibm(ibi)%centIP(l,1); cent(2) = ibm(ibi)%centIP(l,2)
          cent(3) = ibm(ibi)%centIP(l,3)
        else
          print *, 'plyudebug, warning: pretype of Pre_process is incorrect'
        endif 

        !< plyunote: ic, jc is determined based on uniform grid assumption
        ic = floor(cent(1)/(xl/nx))
        jc = floor(cent(2)/(yl/ny))
        if(ic<1 .or. jc<1 .or. ic>nx .or. jc>ny) then
          !if (myid.eq.7 .and. ibi.eq.9 .and. l.eq.1) then
            print *,'Turbine outside domain, (myid,ibi,i_elmt, coordxyz=', &
              myid, ibi, l, cent(1:3)
          !endif
        endif
        if(ic<1) ic = 1
        if(jc<1) jc = 1
        if(ic>nx) ic = nx
        if(jc>ny) ic = ny
         
        do k = zs, ze
          r1 = cent(1) - cartx(ic)
          r2 = cent(2) - carty(jc)
          !r3 = cent(3) - cartz(ic,jc,k)
          r3 = cent(3) - (zzo(k)*(hbar+hho(ic,jc))-hho(ic,jc))
          d1 = sqrt(r1**2+r2**2+r3**2)
          if (d1<dmin) then
            dmin = d1
            kc = k
          endif 
        enddo
        
        r1 = cent(1) - cartx(ic)
        r2 = cent(2) - carty(jc)
        !r3 = cent(3) - cartz(ic,jc,kc)
        r3 = cent(3) - (zzo(kc)*(hbar+hho(ic,jc))-hho(ic,jc))
        dmin = sqrt(r1**2+r2**2+r3**2)

        iclose(l) = ic; jclose(l) = jc; kclose(l) = kc
        ibm(ibi)%iclose(l) = ic; ibm(ibi)%jclose(l) = jc
        ibm(ibi)%kclose(l) = kc
      enddo

      do l = 1, n_elmt
        ic = iclose(l); jc = jclose(l); kc = kclose(l)
              
        ic1 = ic - Itpwidth; ic2 = ic + Itpwidth
        if( (ic1>lxe) .or. (ic2<lxs) ) then
          !ibm(ibi)%i_min(l) = lxs
          !ibm(ibi)%i_max(l) = lxs
          ijkbound(l,1,1:2) = lxs
        else
          !ibm(ibi)%i_min(l) = max(ic1, lxs)
          !ibm(ibi)%i_max(l) = min(ic2, lxe)
          ijkbound(l,1,1) = max(ic1, lxs)
          ijkbound(l,1,2) = min(ic2, lxe)
        endif

        jc1 = jc - Itpwidth; jc2 = jc + Itpwidth 
        if( (jc1>lye) .or. (jc2<lys) ) then
          !ibm(ibi)%j_min(l) = lys
          !ibm(ibi)%j_max(l) = lys
          ijkbound(l,2,1:2) = lys
        else
          !ibm(ibi)%j_min(l) = max(jc1, lys)
          !ibm(ibi)%j_max(l) = min(jc2, lye)
          ijkbound(l,2,1) = max(jc1, lys)
          ijkbound(l,2,2) = min(jc2, lye)
        endif

        kc1 = kc - Itpwidth; kc2 = kc + Itpwidth
        if( (kc1>lze) .or. (kc2<lzs) ) then
          !ibm(ibi)%k_min(l) = lzs
          !ibm(ibi)%k_max(l) = lzs
          ijkbound(l,3,1:2) = lzs
        else
          !ibm(ibi)%k_min(l) = max(kc1, lzs)
          !ibm(ibi)%k_max(l) = min(kc2, lze)
          ijkbound(l,3,1) = max(kc1, lzs)
          ijkbound(l,3,2) = min(kc2, lze)
        endif

      !print*,'plyudebug,ic,jc,kc=',ic,jc,kc
      !print*,'plyudebug,i_min, imax=',ibm(ibi)%i_min(l),ibm(ibi)%i_max(l)
      !print*,'plyudebug,j_min, jmax=',ibm(ibi)%j_min(l),ibm(ibi)%j_max(l)
      !print*,'plyudebug,k_min, kmax=',ibm(ibi)%k_min(l),ibm(ibi)%k_max(l)
      enddo

      if (pretype .eq. 1) then
        ibm(ibi)%i_min(1:n_elmt) = ijkbound(1:n_elmt,1,1)
        ibm(ibi)%i_max(1:n_elmt) = ijkbound(1:n_elmt,1,2)
        ibm(ibi)%j_min(1:n_elmt) = ijkbound(1:n_elmt,2,1)
        ibm(ibi)%j_max(1:n_elmt) = ijkbound(1:n_elmt,2,2)
        ibm(ibi)%k_min(1:n_elmt) = ijkbound(1:n_elmt,3,1)
        ibm(ibi)%k_max(1:n_elmt) = ijkbound(1:n_elmt,3,2)
      elseif (pretype .eq. 2) then
        ibm(ibi)%iIP_min(1:n_elmt) = ijkbound(1:n_elmt,1,1)
        ibm(ibi)%iIP_max(1:n_elmt) = ijkbound(1:n_elmt,1,2)
        ibm(ibi)%jIP_min(1:n_elmt) = ijkbound(1:n_elmt,2,1)
        ibm(ibi)%jIP_max(1:n_elmt) = ijkbound(1:n_elmt,2,2)
        ibm(ibi)%kIP_min(1:n_elmt) = ijkbound(1:n_elmt,3,1)
        ibm(ibi)%kIP_max(1:n_elmt) = ijkbound(1:n_elmt,3,2)
      else
        print *, 'plyudebug, warning: pretype of Pre_process is incorrect'
      endif

      if(allocated(iclose)) deallocate(iclose)
      if(allocated(jclose)) deallocate(jclose)
      if(allocated(kclose)) deallocate(kclose)
      if(allocated(sum_iclose)) deallocate(sum_iclose)
      if(allocated(sum_jclose)) deallocate(sum_jclose)
      if(allocated(sum_kclose)) deallocate(sum_kclose)
      if(allocated(ijkbound)) deallocate(ijkbound)
    enddo

  end subroutine Pre_process
  
  subroutine ACL_read_ucd (ibm,ibi,fsi)
  implicit none
    type(IBMNodes) :: ibm
    type(FSinfo) :: fsi
    integer :: ibi

    integer :: n_v, n_elmt, nv_blade, nelmt_blade
    integer :: nb, i, j, ii, ncmp
    real(wp), allocatable, dimension(:) :: x_bp, y_bp, z_bp

    integer :: n1e, n2e, n3e
    real(wp) :: dr
    real(wp), allocatable, dimension(:) :: nv1, nv2, nv3, dA
    real(wp), allocatable, dimension(:) :: nf_x, nf_y, nf_z
    real(wp), allocatable, dimension(:) :: nt_x, nt_y, nt_z
    real(wp), allocatable, dimension(:) :: ns_x, ns_y, ns_z

    real(wp) :: dA_sum, dA_sum_local
    integer :: dA_n, dA_n_local
    
    dA_n =0; dA_n_local = 0
    dA_sum = 0.0; dA_sum_local = 0.0
    !print*,'plyudebug cpu=',myid,', ACL_read_ucd, 1'
    
    if (myid .eq. 0) then
      open(15003, file='acldata000', action='read', status='old')
      read(15003,*) n_v
    endif
    
    call mpi_bcast(n_v, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    
    n_elmt = n_v - 1
    nv_blade = n_v
    nelmt_blade = n_elmt
    ibm%n_v = n_v * num_blade
    ibm%n_elmt = n_elmt * num_blade
    n_v = ibm%n_v
    n_elmt = ibm%n_elmt

    allocate(x_bp(n_v), y_bp(n_v), z_bp(n_v))
    allocate(ibm%x_bp(n_v), ibm%y_bp(n_v), ibm%z_bp(n_v))
    allocate(ibm%x_bp_i(n_v), ibm%y_bp_i(n_v), ibm%z_bp_i(n_v))
    allocate(ibm%x_bp_o(n_v), ibm%y_bp_o(n_v), ibm%z_bp_o(n_v))
    allocate(ibm%x_bp0(n_v), ibm%y_bp0(n_v), ibm%z_bp0(n_v))  

    !print*,'plyudebug cpu=',myid,', ACL_read_ucd, 2'

    do nb=1,num_blade
      if (myid .eq. 0) then
        if(nb .ne. 1) read(15003,*) ii
      endif
      do j=1,nv_blade
        i = (nb-1)*nv_blade+j
        if(myid .eq. 0) then
          read(15003,*) x_bp(i), y_bp(i), z_bp(i)
        endif
    
        call mpi_bcast(x_bp(i), 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(y_bp(i), 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(z_bp(i), 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)


        ibm%x_bp_i(i) = x_bp(i)
        ibm%y_bp_i(i) = y_bp(i)
        ibm%z_bp_i(i) = z_bp(i)
        
        x_bp(i) = x_bp(i) + fsi%x_c
        y_bp(i) = y_bp(i) + fsi%y_c
        z_bp(i) = z_bp(i) + fsi%z_c

        ibm%x_bp0(i) = x_bp(i)
        ibm%y_bp0(i) = y_bp(i)
        ibm%z_bp0(i) = z_bp(i)
        
        ibm%x_bp(i) = x_bp(i)
        ibm%y_bp(i) = y_bp(i)
        ibm%z_bp(i) = z_bp(i)
        
        ibm%x_bp_o(i) = x_bp(i)
        ibm%y_bp_o(i) = y_bp(i)
        ibm%z_bp_o(i) = z_bp(i)
        
      enddo
    enddo
    
    !print*,'plyudebug cpu=',myid,', ACL_read_ucd, 3'
    
    ncmp = 3
    !> note that u, uold, and urm1 is not initialized
    allocate(ibm%u(n_v,ncmp),ibm%uold(n_v,ncmp),ibm%urm1(n_v,ncmp))
    allocate(ibm%nv1(n_elmt), ibm%nv2(n_elmt), ibm%nv3(n_elmt))
    allocate(ibm%nf_x(n_elmt), ibm%nf_y(n_elmt), ibm%nf_z(n_elmt))
    allocate(ibm%nt_x(n_elmt), ibm%nt_y(n_elmt), ibm%nt_z(n_elmt))
    allocate(ibm%ns_x(n_elmt), ibm%ns_y(n_elmt), ibm%ns_z(n_elmt))
    allocate(ibm%dA(n_elmt))
    
    allocate(nv1(n_elmt), nv2(n_elmt), nv3(n_elmt))
    allocate(nf_x(n_elmt), nf_y(n_elmt), nf_z(n_elmt))
    allocate(nt_x(n_elmt), nt_y(n_elmt), nt_z(n_elmt))
    allocate(ns_x(n_elmt), ns_y(n_elmt), ns_z(n_elmt))
    allocate(dA(n_elmt))

    allocate(ibm%cent_x(n_elmt), ibm%cent_y(n_elmt), ibm%cent_z(n_elmt))
    allocate(ibm%AOA_mean(n_elmt))
    allocate(ibm%Urelmag_mean(n_elmt), ibm%Urelmag(n_elmt))
    allocate(ibm%Urel_mean(n_elmt, ncmp), ibm%Uinduced_mean(n_elmt, ncmp))
    allocate(ibm%circulation_mean(n_elmt, ncmp))
    allocate(ibm%liftdirection_mean(n_elmt, ncmp))
    allocate(ibm%F_lagr_x(n_elmt), ibm%F_lagr_y(n_elmt), ibm%F_lagr_z(n_elmt))
    allocate(ibm%U_lagr_x(n_elmt), ibm%U_lagr_y(n_elmt), ibm%U_lagr_z(n_elmt))
    allocate(ibm%i_min(n_elmt), ibm%j_min(n_elmt), ibm%k_min(n_elmt))
    allocate(ibm%i_max(n_elmt), ibm%j_max(n_elmt), ibm%k_max(n_elmt))
    allocate(ibm%iclose(n_elmt), ibm%jclose(n_elmt), ibm%kclose(n_elmt))
    allocate(ibm%angle_attack(n_elmt), ibm%angle_twist(n_elmt), ibm%chord_blade(n_elmt))
    allocate(ibm%CD(n_elmt), ibm%CL(n_elmt), ibm%CD2(n_elmt), ibm%CL2(n_elmt))
    allocate(ibm%color(n_elmt), ibm%s2l(n_elmt))
    allocate(ibm%Urel(n_elmt,ncmp), ibm%Uinduced(n_elmt,ncmp), ibm%circulation(n_elmt,ncmp))
    allocate(ibm%liftdirection(n_elmt,ncmp), ibm%rotationdirection(n_elmt,ncmp))
    
    !print*,'plyudebug cpu=',myid,', ACL_read_ucd, 4'

    ! colored as blade index
    do nb = 1,num_blade
      do j = 1,nelmt_blade
        i = (nb-1)*nelmt_blade + j
        ii = (nb-1)*nv_blade + j
        
        ibm%color(i) = nb
        ibm%nv1(i) = ii; ibm%nv2(i) = ii+1; ibm%nv3(i) = ii+1 !> plyunote:nv2=nv3 here

        n1e = ibm%nv1(i); n2e = ibm%nv2(i);
        nf_x(i) = ibm%x_bp(n2e) - ibm%x_bp(n1e);
        nf_y(i) = ibm%y_bp(n2e) - ibm%y_bp(n1e);
        nf_z(i) = ibm%z_bp(n2e) - ibm%z_bp(n1e);

        dr = sqrt(nf_x(i)**2+nf_y(i)**2+nf_z(i)**2)
        
        ibm%nf_x(i) = nf_x(i)/dr; ibm%nf_y(i) = nf_y(i)/dr; ibm%nf_z(i) = nf_z(i)/dr
        ibm%dA(i) = dr
        ibm%cent_x(i) = (ibm%x_bp(n1e)+ibm%x_bp(n2e))/2.
        ibm%cent_y(i) = (ibm%y_bp(n1e)+ibm%y_bp(n2e))/2.
        ibm%cent_z(i) = (ibm%z_bp(n1e)+ibm%z_bp(n2e))/2.
        
        ibm%nt_x(i) =0.; ibm%nt_y(i) =0.; ibm%nt_z(i) = 0.
        ibm%ns_x(i) =0.; ibm%ns_y(i) =0.; ibm%ns_z(i) = 0.
              
        dA_sum = dA_sum + ibm%dA(i)
        dA_n = dA_n +1
      enddo
    enddo

    dA_mean = dA_sum/dA_n
    !> for actuator line
    dh_fsi = dA_mean
    !if(myid .eq.0) print*,'dh_fsi=',dh_fsi

    !print*,'plyudebug cpu=',myid,', ACL_read_ucd, 5'
    
    if(myid .eq. 0) then
      close(15003)
    endif

    !print *,'plyudebug, size(ibm%Urelmag_mean)=',size(ibm%Urelmag_mean,1)
    !ibm%Urelmag_mean(:) = 0.0_wp
    !print *,'plyudebug, array 1'
    ibm%Urelmag_mean(:) = 0.0_wp; 
    ibm%Urel_mean(:,:) = 0.0_wp;
    !print *, '2'
    ibm%Uinduced_mean(:,:) = 0.0_wp
    !print *, '3'
    ibm%circulation_mean(:,:) = 0.0_wp; 
    !print *,'4'
    ibm%liftdirection_mean(:,:) = 0.0_wp
    
    !print*,'plyudebug cpu=',myid,', ACL_read_ucd, 6'

  end subroutine ACL_read_ucd

  subroutine disk_read_ucd(ibm, ibi, fsi, fname, ioffset)
  implicit none
    type(IBMNodes) :: ibm
    integer :: ibi, ioffset
    type(FSInfo) :: fsi
    character(len=*) :: fname
    
    integer :: n_v, n_elmt, n1e, n2e, n3e, i, ii, ii1, ii2, ii3, ncmp
    real(wp) :: dr, dx12, dy12, dz12, dx13, dy13, dz13, rr

    character(64) :: ss

    integer, allocatable, dimension(:) :: nv1, nv2, nv3
    real(wp), allocatable, dimension(:) :: x_bp, y_bp, z_bp
    real(wp), allocatable, dimension(:) :: nf_x, nf_y, nf_z
    real(wp), allocatable, dimension(:) :: nt_x, nt_y, nt_z
    real(wp), allocatable, dimension(:) :: ns_x, ns_y, ns_z

    real(wp) :: dA_sum
    integer :: dA_n
    
    dA_sum = 0.0; dA_n = 0

    if (myid .eq. 0) then    
      open(15004, file=fname, action='read', status='old')
      read(15004, *)
      read(15004, *)
      read(15004, *)

      read(15004, *) n_v, n_elmt
    endif

    call  mpi_bcast(n_v, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call  mpi_bcast(n_elmt, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)

    ibm%n_v = n_v
    ibm%n_elmt = n_elmt

    allocate(x_bp(n_v), y_bp(n_v), z_bp(n_v))
    allocate(ibm%x_bp(n_v), ibm%y_bp(n_v), ibm%z_bp(n_v))
    allocate(ibm%x_bp_i(n_v), ibm%y_bp_i(n_v), ibm%z_bp_i(n_v))
    allocate(ibm%x_bp_o(n_v), ibm%y_bp_o(n_v), ibm%z_bp_o(n_v))
    allocate(ibm%x_bp0(n_v), ibm%y_bp0(n_v), ibm%z_bp0(n_v))

    do i = 1,n_v 
      if (myid .eq. 0) then
        read(15004,*) ii, x_bp(i), y_bp(i), z_bp(i)
      endif
      call  mpi_bcast(x_bp(i), 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call  mpi_bcast(y_bp(i), 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call  mpi_bcast(z_bp(i), 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)


      ibm%x_bp_i(i) = x_bp(i)
      ibm%y_bp_i(i) = y_bp(i)
      ibm%z_bp_i(i) = z_bp(i)
      
      x_bp(i) = x_bp(i) + fsi%x_c
      y_bp(i) = y_bp(i) + fsi%y_c
      z_bp(i) = z_bp(i) + fsi%z_c

      ibm%x_bp0(i) = x_bp(i)
      ibm%y_bp0(i) = y_bp(i)
      ibm%z_bp0(i) = z_bp(i)

      if(ioffset .ne. 0) then
        rr = loc_refvel * 2.0 * fsi%r_rotor
        x_bp(i) = x_bp(i) - rr * fsi%nx_tb
        y_bp(i) = y_bp(i) - rr * fsi%ny_tb
        z_bp(i) = z_bp(i) - rr * fsi%nz_tb
      endif
  
      ibm%x_bp(i) = x_bp(i)
      ibm%y_bp(i) = y_bp(i)
      ibm%z_bp(i) = z_bp(i)
      
      ibm%x_bp_o(i) = x_bp(i)
      ibm%y_bp_o(i) = y_bp(i)
      ibm%z_bp_o(i) = z_bp(i)
    enddo

    ncmp = 3
    allocate(ibm%u(n_v,ncmp),ibm%uold(n_v,ncmp),ibm%urm1(n_v,ncmp))
    allocate(ibm%nv1(n_elmt), ibm%nv2(n_elmt), ibm%nv3(n_elmt))
    allocate(ibm%nf_x(n_elmt), ibm%nf_y(n_elmt), ibm%nf_z(n_elmt))
    allocate(ibm%nt_x(n_elmt), ibm%nt_y(n_elmt), ibm%nt_z(n_elmt))
    allocate(ibm%ns_x(n_elmt), ibm%ns_y(n_elmt), ibm%ns_z(n_elmt))
    allocate(ibm%dA(n_elmt))
    
    allocate(nv1(n_elmt), nv2(n_elmt), nv3(n_elmt))
    allocate(nf_x(n_elmt), nf_y(n_elmt), nf_z(n_elmt))
    allocate(nt_x(n_elmt), nt_y(n_elmt), nt_z(n_elmt))
    allocate(ns_x(n_elmt), ns_y(n_elmt), ns_z(n_elmt))
    !allocate(dA(n_elmt))

    allocate(ibm%cent_x(n_elmt), ibm%cent_y(n_elmt), ibm%cent_z(n_elmt))
    allocate(ibm%AOA_mean(n_elmt))
    allocate(ibm%Urelmag_mean(n_elmt), ibm%Urelmag(n_elmt))
    allocate(ibm%Urel_mean(n_elmt, ncmp), ibm%Uinduced_mean(n_elmt, ncmp))
    allocate(ibm%circulation_mean(n_elmt, ncmp))
    allocate(ibm%liftdirection_mean(n_elmt, ncmp))
    allocate(ibm%F_lagr_x(n_elmt), ibm%F_lagr_y(n_elmt), ibm%F_lagr_z(n_elmt))
    allocate(ibm%U_lagr_x(n_elmt), ibm%U_lagr_y(n_elmt), ibm%U_lagr_z(n_elmt))
    allocate(ibm%i_min(n_elmt), ibm%j_min(n_elmt), ibm%k_min(n_elmt))
    allocate(ibm%i_max(n_elmt), ibm%j_max(n_elmt), ibm%k_max(n_elmt))
    allocate(ibm%iclose(n_elmt), ibm%jclose(n_elmt), ibm%kclose(n_elmt))
    allocate(ibm%angle_attack(n_elmt), ibm%angle_twist(n_elmt), ibm%chord_blade(n_elmt))
    allocate(ibm%CD(n_elmt), ibm%CL(n_elmt))
    allocate(ibm%color(n_elmt), ibm%s2l(n_elmt))
    allocate(ibm%s2l1(n_elmt), ibm%s2l2(n_elmt), ibm%s2lc(n_elmt))
    allocate(ibm%Urel(n_elmt,ncmp), ibm%Uinduced(n_elmt,ncmp), ibm%circulation(n_elmt,ncmp))
    allocate(ibm%liftdirection(n_elmt,ncmp), ibm%rotationdirection(n_elmt,ncmp))

    ! for nacelle wall model
    allocate(ibm%dh_IP(n_elmt))
    allocate(ibm%centIP(n_elmt,ncmp))
    allocate(ibm%U_IPlagr(n_elmt,ncmp))
    allocate(ibm%iIP_min(n_elmt), ibm%jIP_min(n_elmt), ibm%kIP_min(n_elmt))
    allocate(ibm%iIP_max(n_elmt), ibm%jIP_max(n_elmt), ibm%kIP_max(n_elmt))

    if(myid .eq. 0) then
      do i = 1,n_elmt
        read(15004, *) ii1, ii2, ss, nv1(i), nv2(i), nv3(i)
      enddo
    endif
    call  mpi_bcast(nv1, n_elmt, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call  mpi_bcast(nv2, n_elmt, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call  mpi_bcast(nv3, n_elmt, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)


    ibm%nv1 = nv1; ibm%nv2 = nv2; ibm%nv3 = nv3

    if(myid .eq. 0) then
      close(15004)
    endif

    do i = 1,n_elmt
      n1e = ibm%nv1(i); n2e = ibm%nv2(i); n3e = ibm%nv3(i)
      
      dx12 = x_bp(n2e) - x_bp(n1e)
      dy12 = y_bp(n2e) - y_bp(n1e)
      dz12 = z_bp(n2e) - z_bp(n1e)
      
      dx13 = x_bp(n3e) - x_bp(n1e)
      dy13 = y_bp(n3e) - y_bp(n1e)
      dz13 = z_bp(n3e) - z_bp(n1e)
      
      nf_x(i) = dy12 * dz13 - dz12 * dy13
      nf_y(i) = -dx12 * dz13 + dz12 * dx13
      nf_z(i) = dx12 * dy13 - dy12 * dx13

      dr = sqrt(nf_x(i)**2+nf_y(i)**2+nf_z(i)**2)

      nf_x(i) = nf_x(i)/dr; nf_y(i) = nf_y(i)/dr; nf_z(i) = nf_z(i)/dr
      
      if( ( (1.-nf_x(i))<=1.d-6 .and. (-1.+nf_x(i))<1.d-6 )  & 
        .or. ( (nf_x(i)+1.)<=1.d-6 .and. (-1.-nf_x(i))<1.d-6 ) ) then
        ns_x(i) = 0. ; ns_y(i) = 1. ; ns_z(i) = 0.
      else
        rr = sqrt(nf_y(i)**2 + nf_z(i)**2 )
        !> assume turbine mainly points to x direction
        !> ns is in y-z plane and perpendicular to projection of nf
        !> nt is cross multiplification of ns and nf
        !> the value need to be checked
        ns_x(i) = 0.
        ns_y(i) = nf_z(i) / rr
        ns_z(i) = -nf_y(i) / rr
        nt_x(i) = rr
        nt_y(i) = -nf_y(i)*nf_x(i) / rr
        nt_z(i) = -nf_z(i)*nf_x(i) / rr
      endif

      ibm%dA(i) = dr/2.

      dA_sum = dA_sum + ibm%dA(i)
      dA_n = dA_n + 1
      ! calc the center of the element
      ibm%cent_x(i) = (x_bp(n1e)+x_bp(n2e)+x_bp(n3e))/3.
      ibm%cent_y(i) = (y_bp(n1e)+y_bp(n2e)+y_bp(n3e))/3.
      ibm%cent_z(i) = (z_bp(n1e)+z_bp(n2e)+z_bp(n3e))/3.

      ibm%nf_x(i) = nf_x(i); ibm%nf_y(i) = nf_y(i); ibm%nf_z(i) = nf_z(i)
      ibm%nt_x(i) = nt_x(i); ibm%nt_y(i) = nt_y(i); ibm%nt_z(i) = nt_z(i)
      ibm%ns_x(i) = ns_x(i); ibm%ns_y(i) = ns_y(i); ibm%ns_z(i) = ns_z(i)

      !if(myid.eq.0) then
      !  print *, i, n1e, n2e, n3e, ibm%cent_x(i), ibm%cent_y(i), ibm%cent_z(i)
      !endif
    enddo

    dA_mean = dA_sum / dA_n
    !> for triangle elemnts
    dh_uref = sqrt(dA_mean)
    if(myid.eq.0)print*,'dh_uref=',dh_uref

    ibm%Urelmag_mean = 0.; ibm%Urel_mean = 0.; ibm%Uinduced_mean = 0.
    ibm%circulation_mean = 0.; ibm%liftdirection_mean = 0.

  end subroutine disk_read_ucd

  subroutine reread_uref (ibm, fsi, NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi
    integer :: NumberOfObjects
    
    integer :: ibi
    logical :: l_exist
     
    if (myid .eq. 0) then
      inquire(file='restart_uref.dat', exist=l_exist)
      if (l_exist) then
        print *, '!!!!!!! Reading restart_uref.dat !!!!!!'
        open(15010, file='restart_uref.dat', action='read', status='old')
        do ibi = 1, NumberOfObjects
          ibm(ibi)%U1_ref = 0.8_wp; ibm(ibi)%U_ref = 0.8_wp
          read(15010,*) ibm(ibi)%U1_ref, ibm(ibi)%U_ref
        enddo
        close(15010)
      else
        print *, '!!!!!!! Creating restart_uref.dat !!!!!!!'
        open(15010, file='restart_uref.dat', action='write', status='new')
        do ibi = 1, NumberOfObjects
          !> plyunote: 0.8 is estimated with the assumption that u_top=1.0
          ibm(ibi)%U1_ref = 0.8_wp; ibm(ibi)%U_ref = 0.8_wp
          write(15010,*) ibm(ibi)%U1_ref, ibm(ibi)%U_ref, ' # U1_ref, U_ref',&
            ' for turbine ', ibi
        enddo
        close(15010)
      endif 
    endif
    
    do ibi = 1, NumberOfTurbines 
      call mpi_bcast(ibm(ibi)%U1_ref, 1, mpi_double_precision, 0, &
        mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%U_ref, 1, mpi_double_precision, 0, &
        mpi_comm_2d_cart, ierr_wt)
    enddo

  end subroutine reread_uref

  recursive subroutine quicksort2(a, ncol, first, last)
    implicit none
    real(wp) :: a(ncol, 2), x, temp(2)
    integer :: ncol, first, last
    integer :: i, j

    x = a( (first+last) / 2, 1 )
    i = first
    j = last
    do
       do while (a(i, 1) < x)
          i=i+1
       end do
       do while (x < a(j,1))
          j=j-1
       end do
       if (i >= j) exit
       temp = a(i,:);  a(i,:) = a(j,:);  a(j,:) = temp
       i=i+1
       j=j-1
    end do
    if (first < i-1) call quicksort2(a, ncol, first, i-1)
    if (j+1 < last)  call quicksort2(a, ncol, j+1, last)
  end subroutine quicksort2

  recursive subroutine binary_search(idx, x, a, n, first, last)
    implicit none
    integer :: idx, first, last, n, mid
    real(wp) :: x, a(n)

    if (last >= first) then
      mid = first + (last - first)/2
      if (a(mid) == x) then
        idx = mid
        return
      elseif (mid == first) then
        idx = mid
        return
      elseif (mid >= (last - 1)) then
        idx = last - 1
        return
      elseif (a(mid) < x .and. a(mid+1) >= x) then
        idx = mid
        return
      elseif (a(mid) < x) then
        call binary_search(idx, x, a, n, mid, last)
        return
      else        
        call binary_search(idx, x, a, n, first, mid-1)
        return
      endif
    endif
  end subroutine binary_search

  !> Range adjustment and sort.
  !! Adjust the angle to the range of [-pi, pi]. Later ALM will use asin 
  !! to calculate attack angle, whose range is [-pi, pi].
  subroutine reorganize_cdcl_table(angle, coeff, n)
    implicit none
    real(wp), dimension(:), intent(inout) :: angle, coeff
    integer, intent(in) :: n

    integer :: i
    real(wp), dimension(:,:) :: table(n,2)

    !> subtract 2*pi
    do i = 1, n
      if (angle(i)>180.0) angle(i) = angle(i) - 360.0
    enddo

    !> sort
    table(:,1) = angle(:)
    table(:,2) = coeff(:)
    call quicksort2(table, n, 1, n)
    angle(:) = table(:,1)
    coeff(:) = table(:,2)
    
  end subroutine reorganize_cdcl_table

  !> Find the attack angle that has a lift coefficient of 0.0
  subroutine locate_y0(x_y0, ang_arr, coeff_arr, n)
    implicit none
    integer :: n
    real(wp) :: x_y0, ang_arr(n), coeff_arr(n)
    
    integer :: idx_x0, idx_y0, width, i, inrange

    idx_x0 = n/2
    call binary_search(idx_x0, 0.0_wp, ang_arr, n, 1, n)
    width = min(idx_x0-1, n - idx_x0)

    inrange = 0; idx_y0 = idx_x0
    do i = 0, width
      if (coeff_arr(idx_x0+i)<=0.0 .and. coeff_arr(idx_x0+i+1)>0.0 ) then
        idx_y0 = idx_x0+i; inrange = 1; exit
      elseif (coeff_arr(idx_x0-i)<=0.0 .and. coeff_arr(idx_x0-i+1)>0.0 ) then
        idx_y0 = idx_x0-i; inrange = 1; exit
      endif
    enddo
    if (inrange == 0) print *, "Warning: cannot find zero lift coefficient in table."

    x_y0 = ang_arr(idx_y0)
    call interp_2p(0.0_wp, x_y0, coeff_arr(idx_y0), ang_arr(idx_y0), &
      coeff_arr(idx_y0+1), ang_arr(idx_y0+1))

  end subroutine locate_y0

  !< Read information of ACL, like chord, twist angle, lift and drag 
  !> coefficients, for actuator line simulation
  subroutine airfoil_ACL(ibm, fsi)
    implicit none
    !type(ACL) :: acl
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi

    integer :: n_CD, n_CL, ifoil, i, num_AC, num_CD, num_CL, ibi, j
    integer :: kk, iset, inrange
    real(wp) :: r, fac1, fac2
    character(128) :: string, filen

    if (.not.(allocated(acl(1)%chord_bladeInp))) then
      do ifoil = 1,num_foiltype
        !> read FOILxx
        if (myid .eq. 0) then
          filen=''
          !> foil number starts from 1 in program
          !> while it starts from FOIL00 in input filelists
          write(filen, '(a4,i0.2)') 'FOIL', ifoil-1
          open(15005, file=trim(filen), action='read', status='old')
          read(15005,*)
          read(15005,*)
          
          read(15005,*) num_AC
        endif

        call mpi_bcast(num_AC, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)

        acl(ifoil)%num_AC = num_AC

        allocate(acl(ifoil)%r_AC(num_AC), acl(ifoil)%angle_twistInp(num_AC))
        allocate(acl(ifoil)%chord_bladeInp(num_AC))

        if (myid .eq. 0) then
        do i = 1,num_AC
          read(15005,*) acl(ifoil)%r_AC(i), acl(ifoil)%angle_twistInp(i), &
                  acl(ifoil)%chord_bladeInp(i)
        enddo
        endif

        call mpi_bcast(acl(ifoil)%r_AC, num_AC, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(acl(ifoil)%angle_twistInp, num_AC, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(acl(ifoil)%chord_bladeInp, num_AC, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)


        acl(ifoil)%r_beg = acl(ifoil)%r_AC(1)
        acl(ifoil)%r_end = acl(ifoil)%r_AC(num_AC)

        if (myid .eq. 0) then
          close(15005)
        endif

        !> read CDxx
        if (myid .eq. 0) then
          filen=''
          write(filen, '(a2,i0.2)') 'CD', ifoil-1
          open(15006, file=trim(filen), action='read', status='old')
          read(15006,*)
          read(15006,*)

          read(15006,*) num_CD
        endif

        call mpi_bcast(num_CD, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)

        acl(ifoil)%num_CD = num_CD

        allocate(acl(ifoil)%ang_CD(num_CD), acl(ifoil)%CDInp(num_CD))

        if(myid .eq. 0) then
          do i = 1,num_CD
            read(15006,*) acl(ifoil)%ang_CD(i), acl(ifoil)%CDInp(i) 
          enddo
               
          close(15006)

          !> plyunote: I left the draft of reorganization here. But I actually did 
          !! the adjustment of CD and CL database outside this program. Basically I
          !! assume the database use a range of attack angle of [-pi, pi]. I wish the
          !! range can cover possible stall angle in pitch motion.
          !call reorganize_cdcl_table(acl(ifoil)%ang_CD, acl(ifoil)%CDInp, num_CD)
        endif

        call mpi_bcast(acl(ifoil)%ang_CD, num_CD, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(acl(ifoil)%CDInp, num_CD, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)       

        !> for 3D correction: find C_D at alpha=0
        acl(ifoil)%CD_zero_angle = 0.0
        call interp_in_table(0.0_wp, acl(ifoil)%CD_zero_angle, &
          acl(ifoil)%ang_CD(:), acl(ifoil)%CDInp(:), acl(ifoil)%num_CD, &
          inrange)
        if (inrange == 0 .and. myid == 0) print *, "Warning: AOA=0 not found in CD table", &
          ' for foil type ', ifoil
        if (inrange == 1 .and. myid == 0) print *, "C_D = ",acl(ifoil)%CD_zero_angle, &
           " at AOA=0 for foil type ", ifoil

        !> read CLxx
        if(myid .eq. 0) then
          filen=''
          write(filen, '(a2,i0.2)') 'CL', ifoil-1
          open(15007, file=trim(filen), action='read', status='old')
          read(15007,*)
          read(15007,*)

          read(15007,*) num_CL
        endif

        call mpi_bcast(num_CL, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)

        acl(ifoil)%num_CL = num_CL

        allocate(acl(ifoil)%ang_CL(num_CL), acl(ifoil)%CLInp(num_CL))

        if(myid .eq. 0) then
          do i = 1,num_CL
            read(15007,*) acl(ifoil)%ang_CL(i), acl(ifoil)%CLInp(i) 
          enddo

          close(15007)
        endif

        call mpi_bcast(acl(ifoil)%ang_CL, num_CL, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(acl(ifoil)%CLInp, num_CL, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)

        !> for 3D correction: find zero-lift angle
        call locate_y0(acl(ifoil)%angle_zero_lift, acl(ifoil)%ang_CL, acl(ifoil)%CLInp, num_CL)
        if (myid .eq. 0) print *, "Zero-lift angle is ", acl(ifoil)%angle_zero_lift
      enddo
    endif

    do ibi = 1,NumberOfTurbines
      do i = 1,ibm(ibi)%n_elmt
        iset = 0
        r = sqrt( (ibm(ibi)%cent_x(i)-fsi(ibi)%x_c)**2 &
          + (ibm(ibi)%cent_y(i)-fsi(ibi)%y_c)**2 &
          + (ibm(ibi)%cent_z(i)-fsi(ibi)%z_c)**2 )

        do ifoil = 1,num_foiltype
          if(r>=acl(ifoil)%r_beg .and. r<=acl(ifoil)%r_end) then
            do j = 1,acl(ifoil)%num_AC-1
              if (r>=acl(ifoil)%r_AC(j) .and. r<=acl(ifoil)%r_AC(j+1)) then
                iset = 1
                fac1 = (acl(ifoil)%r_AC(j+1)-r)/( acl(ifoil)%r_AC(j+1) &
                  -acl(ifoil)%r_AC(j) )
                fac2 = (-acl(ifoil)%r_AC(j)+r)/( acl(ifoil)%r_AC(j+1) &
                  -acl(ifoil)%r_AC(j) )
                ibm(ibi)%angle_twist(i) = fac1*acl(ifoil)%angle_twistInp(j) &
                  +fac2*acl(ifoil)%angle_twistInp(j+1)
                ibm(ibi)%chord_blade(i) = fac1*acl(ifoil)%chord_bladeInp(j) &
                  +fac2*acl(ifoil)%chord_bladeInp(j+1)
              endif ! end of r
            enddo ! end of j =1, num_AC-1
          endif ! end of r:(r_beg,r_end)
        enddo ! end of ifoil=1,num_foiltype

        if (iset .eq. 0 .and. myid .eq. 0) then
          print *,'plyudebug, airfoil_ACL: element out of range of FOIL input,',&
            i,r,acl(1)%r_beg,acl(1)%r_end
        endif
      enddo ! end of i =1, n_elmt
    enddo ! end of ibi = 1, NumberOfTurbines


  end subroutine airfoil_ACL

  subroutine update_zgrid
    use utils
    use lib_array, only: is_nan
    implicit none
    
    integer :: i, j, k, kk
    logical :: inan

    do k = 1, nz
      zz1(k) = -1.0D16
      zzo(k) = 0.
      zw1(k) = -1.0D16
      zwo(k) = 0.
    enddo

    do k = 1, xsz(3)
      kk = xst(3) + k -1
      zz1(kk) = zz(k)
      zw1(kk) = zw(k)
      !print*,myid,xst(3),xsz(3),k,kk,zz(k)
    enddo

    call mpi_allreduce(zz1,zzo,nz,mpi_double_precision, &
      mpi_max, mpi_comm_world, ierr_wt)
    call mpi_allreduce(zw1,zwo,nz,mpi_double_precision, &
      mpi_max, mpi_comm_world, ierr_wt)

    !call mpi_barrier(mpi_comm_world, ierr_wt); print *, 11, myid
    !do k = 1, nz
      !inan = is_nan(zzo(k))
    !  if (inan) print *, "NaN in zzo,", k
    !  inan = is_nan(zwo(k))
    !  if (inan) print *, "NaN in zwo,", k
    !enddo
    !call mpi_barrier(mpi_comm_world, ierr_wt); print *, 12, myid

    !print *, size(eta,1), size(eta,2), size(eo,1), size(eo,2), xsz(1), xsz(2)

    !call alltoone(eta, eo)
    !call alltoone(hh, hho)
    !call alltoone(hx, hxo)
    !call alltoone(hy, hyo)
    call gather_2d_xy(eta, eo)
    !call mpi_barrier(mpi_comm_world, ierr_wt); print *, 13, myid
    call gather_2d_xy(hh, hho)
    !call mpi_barrier(mpi_comm_world, ierr_wt); print *, 14, myid
    call gather_2d_xy(hx, hxo)
    !call mpi_barrier(mpi_comm_world, ierr_wt); print *, 15, myid
    call gather_2d_xy(hy, hyo)
    !call mpi_barrier(mpi_comm_world, ierr_wt); print *, 16, myid
    !if (myid.eq.0) print *, eo(:,1)

    !call mpi_barrier(mpi_comm_world, ierr_wt)
    !if (myid .eq. 0) then
    !  do k = 1, nz
    !    print *, k, zzo(k), zwo(k)
    !  enddo
    !endif
    !call mpi_barrier(mpi_comm_world, ierr_wt)
    !if (myid .eq. 0) then
    !  do i = 1, nx
    !    do j = 1, ny
    !      print *, i, j, hho(i,j)
    !    enddo
    !  enddo
    !endif
    !call mpi_barrier(mpi_comm_world, ierr_wt)
    !print *, 17, myid
    !call mpi_barrier(mpi_comm_world, ierr_wt)
    !print *, myid, size(cartz), size(cartzw)
    !call mpi_barrier(mpi_comm_world, ierr_wt)


    !do i = 1, nx
    !  do j = 1, ny
    !    inan = is_nan(hho(i,j))
    !    if (inan) print *, "NaN in hh0,", i, j
    !    do k = 1, nz
    !      cartz(i,j,k) = zzo(k) * (hbar + hho(i,j)) - hho(i,j) 
    !      cartzw(i,j,k) = zwo(k) * (hbar + hho(i,j)) - hho(i,j) 
    !      !if(myid.eq.0) print*,i,j,k,zzo(k),eo(i,j),hho(i,j),cartz(i,j,k) 
    !    enddo
    !  enddo
    !enddo

    !call mpi_barrier(mpi_comm_world, ierr_wt); print *, 18, myid

  end subroutine update_zgrid

  subroutine collect_grid
    use spectral_hos
    use mpi
    implicit none

    integer :: i, j, k, kk
    integer :: mpisize
    integer, dimension(:,:), allocatable :: xstall, xszall
    integer, dimension(:), allocatable :: temps, tempr

    integer :: rankrow, rankcol, rank2d, coordrow, coordcol
    integer, dimension(2) :: coord2d

    if (.not. allocated(eo)) then
      allocate(eo(nx_global, ny_global), hho(nx_global, ny_global))
    endif
    allocate(hxo(nx_global, ny_global), hyo(nx_global, ny_global))
   
    allocate(zzo(nz), zz1(nz))
    allocate(zwo(nz), zw1(nz))
    
    allocate(cartx(nx), carty(ny))
    
    !> plyunote: in early version, I allocated memory for cartz,
    !!           it will lead to memory error when nx*ny*nz is too large,
    !!           so I choose to calculate temporary cartz before it is used.
    !allocate(cartz(nx,ny,nz), cartzw(nx, ny, nz))

    do i = 1, nx
      cartx(i) = (real(i)-0.5)*xl/real(nx)
    enddo

    do j = 1, ny
      carty(j) = (real(j)-0.5)*yl/real(ny)
    enddo

      !print *, 'plyudebug-',myid,': size(zz1,1)=',size(zz1,1)
      !print *, 'plyudebug-',myid,': size(zzo,1)=',size(zzo,1)

    call update_zgrid
    !print *, myid, 1; call mpi_barrier(mpi_comm_2d_cart, ierr_wt)
    
    !> plyunote: about two sub-communicater
    !> mpi_comm_2d_col limit the communications in y-direction division.
    !> mpi_comm_2d_row limit the communications in z-direction division.

    call mpi_comm_rank(mpi_comm_2d_cart, rank2d, ierr_wt) ; !print *, myid, 2
    call mpi_comm_rank(mpi_comm_2d_col, rankcol, ierr_wt); !print *, myid, 3
    call mpi_comm_rank(mpi_comm_2d_row, rankrow, ierr_wt); !print *, myid, 4

    call mpi_cart_coords(mpi_comm_2d_cart, rank2d, 2, coord2d, ierr_wt)
    !call mpi_cart_coords(mpi_comm_2d_col, rankcol, 1, coordcol, ierr_wt)
    !call mpi_cart_coords(mpi_comm_2d_row, rankrow, 1, coordrow, ierr_wt)

    if (myid.eq.0) then
      print *, 'nproc1=', nproc1
      print *, 'xsz(1:3)=', xsz
      print *, 'ysz(1:3)=', ysz
    endif
    
    call mpi_barrier(mpi_comm_world, ierr_wt) 
    if(myid.eq.0) then
      print *, 'plyunote: Domain decomposition info:'
      print *, 'rank2d, coord2d(1:2), rankcol, rankrow, xst(1:3)'
    endif
    call mpi_barrier(mpi_comm_world, ierr_wt)
    print *, rank2d, coord2d(1), coord2d(2), rankcol, rankrow, xst(1:3)
    call mpi_barrier(mpi_comm_world, ierr_wt)

    !> plyunote: mpi world info, typically decomposition in y and z
!    call mpi_comm_size(mpi_comm_world, mpisize, ierr_wt)
!    allocate(xstall(mpisize,3), xszall(mpisize,3))
!    allocate(temps(mpisize*3), tempr(mpisize*3))
!    
!    xstall(:,:) = 0; xszall(:,:) = 0
!    temps(:) = 0; tempr(:) = 0 
!    
!    temps((myid*3+1):(myid*3+3)) = xst(1:3)
!    call mpi_allreduce(temps,tempr,mpisize*3,mpi_integer,mpi_sum,&
!      mpi_comm_world,ierr_wt)
!    do i = 1, mpisize
!      xstall(i,:) = tempr((i*3-2):(i*3))
!    enddo
!    
!    temps((myid*3+1):(myid*3+3)) = xsz(1:3)
!    call mpi_allreduce(temps,tempr,mpisize*3,mpi_integer,mpi_sum,&
!      mpi_comm_world,ierr_wt)
!    do i = 1, mpisize
!      xszall(i,:) = tempr((i*3-2):(i*3))
!    enddo
!
!    if (myid.eq.1) then
!      print *, "plyudebug, mpi_comm_size:", mpisize
!      print *, "plyudebug, mpi_world, start index:"
!      print *, xstall(:,1)
!      print *, xstall(:,2)
!      print *, xstall(:,3)
!      print *, "plyudebug, mpi_world, size:"
!      print *, xszall(:,1)
!      print *, xszall(:,2)
!      print *, xszall(:,3)
!    endif

  end subroutine collect_grid

  subroutine dfunc_gaussian(r, epsilon, pdf)
    implicit none
    real(wp), intent(in) :: r, epsilon
    real(wp), intent(out) :: pdf

    !> epsilon should be > 0, I didn't check the input here.
    pdf = 1.0_wp / epsilon / sqrt_pi * exp(-(r/epsilon)**2)
  end subroutine dfunc_gaussian

  subroutine dfunc_s1h(r, dirac)
    implicit none
    real(wp) :: dirac
    real(wp) :: r

    dirac = 0.0_wp

    if(dabs(r)<=0.5_wp) then
      dirac = 0.75 - r**2
    else if (dabs(r)<=1.5_wp) then
      dirac = 1.125_wp - 1.5_wp * dabs(r) + 0.5_wp * r**2
    else
      dirac = 0.0_wp
    endif 
  end subroutine dfunc_s1h
  
  subroutine dfunc_s2h(r, dirac)
    implicit none
    real(wp) :: dirac
    real(wp) :: r

    dirac = 0.0_wp

    if(dabs(r)<=1.5_wp) then
      dirac = 0.5_wp/twopi*(twopi/2.0_wp+2.0_wp &
        * sin(twopi/8.0_wp*(2.0_wp*r+1.0_wp)) &
        - 2.0_wp*sin(twopi/8.0*(2.0_wp*r-1.0_wp)))
    else if (dabs(r)<=2.5_wp) then
      dirac = -0.25_wp/twopi*(-2.5_wp*twopi+twopi*dabs(r) &
        + 4.0_wp*sin(twopi/8.0_wp*(2.0_wp*dabs(r)-1.0_wp)))
    else
      dirac = 0.0_wp
    endif 
  end subroutine dfunc_s2h

  subroutine dfunc_s3h(r, dirac)
    implicit none
    real(wp) :: dirac
    real(wp) :: r

    dirac = 0.0_wp

    if(dabs(r)<=1.0_wp) then
      dirac = 17.0_wp/48.0_wp+sqrt(3.0_wp)*twopi/216.0_wp &
        + dabs(r)/4.0_wp - r*r/4.0_wp + (1.0_wp-2.0_wp*dabs(r)) &
        / 16.0_wp * sqrt(-12.0_wp*r*r+12.0_wp*dabs(r)+1.0_wp) &
        - sqrt(3.0_wp)/12.0_wp*asin(sqrt(3.0_wp)/2.0_wp &
        * (2.0_wp*dabs(r)-1.0_wp))
    else if (dabs(r)<=2.0_wp) then
      dirac = 55.0_wp/48.0_wp-sqrt(3.0_wp)*twopi/216.0_wp &
        - 13.0_wp*dabs(r)/12.0_wp + r*r/4.0_wp + (2.0_wp*dabs(r)-3.0_wp) &
        / 48.0_wp * sqrt(-12.0_wp*r*r+36.0_wp*dabs(r)-23.0_wp) &
        + sqrt(3.0_wp)/36.0_wp*asin(sqrt(3.0_wp)/2.0_wp &
        * (2.0_wp*dabs(r)-3.0_wp))
    else
      dirac = 0.0_wp
    endif 
  end subroutine dfunc_s3h

  subroutine dfunc_s4h(r, dirac)
    implicit none
    real(wp) :: dirac
    real(wp) :: r

    dirac = 0.0_wp

    if (dabs(r)<=0.5_wp) then
      dirac = 3.0_wp/8.0_wp+twopi/64.0_wp-r**2/4.0_wp
    else if(dabs(r)<=1.5_wp) then
      dirac = 0.25_wp+(1.0_wp-dabs(r))/8.0_wp*sqrt( &
        - 2.0_wp+8.0_wp*dabs(r)-4.0_wp*r**2) &
        - 0.125_wp*asin(sqrt(2.0_wp)*(dabs(r)-1.0_wp))
    else if(dabs(r)<=2.5_wp) then
      dirac = 17.0_wp/16.0_wp-twopi/128.0_wp-0.75_wp*dabs(r) &
        + 0.125_wp*r**2+(dabs(r)-2.0_wp)/16.0_wp &
        * sqrt(-14.0_wp+16.0_wp*dabs(r)-4.0_wp*dabs(r)**2) &
        + 1.0_wp/16.0_wp*asin(sqrt(2.0_wp)*(dabs(r)-2.0_wp))
    else
      dirac = 0.0_wp
    endif

    !if(isnan(dirac)) print*,'plyudebug,r=',r,',dirac=',dirac
  end subroutine dfunc_s4h

  subroutine surface_read_xpatch(ibm, ibi, fsi)
    implicit none
    type(IBMNodes) :: ibm
    type(FSInfo) :: fsi
    integer :: ibi

    integer :: n_v, n_elmt, i, ii, n1e, n2e, n3e, ncmp, color0
    real(wp) :: dr, rr
    real(wp), dimension(3) :: ds12, ds13, tmparray
  
    n_v = 0
    if (myid .eq. 0) then
      open(15008, file='acsdata000', action='read', status='old')
      read(15008,*)
      read(15008,*)
      read(15008,*)
      read(15008,*)

      read(15008,*) n_v
    endif
    call mpi_bcast(n_v, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    ibm%n_v = n_v

    allocate(ibm%x_bp(n_v), ibm%y_bp(n_v), ibm%z_bp(n_v))
    allocate(ibm%x_bp_i(n_v), ibm%y_bp_i(n_v), ibm%z_bp_i(n_v))
    allocate(ibm%x_bp_o(n_v), ibm%y_bp_o(n_v), ibm%z_bp_o(n_v))
    allocate(ibm%x_bp0(n_v), ibm%y_bp0(n_v), ibm%z_bp0(n_v))
    
    !> plyunote: u is defined at nodes, not center of elements
    !! so its size is n_v, not n_elmt
    allocate(ibm%u(n_v,3), ibm%uold(n_v,3), ibm%urm1(n_v,3))
    !allocate(ibm%tmprt(n_v))

    if (myid .eq. 0) then
      do i = 1, n_v
        read(15008,*) ibm%x_bp(i), ibm%y_bp(i), ibm%z_bp(i)
      enddo
    endif

    call mpi_bcast(ibm%x_bp, n_v, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(ibm%y_bp, n_v, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(ibm%z_bp, n_v, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)

    do i = 1, n_v
      !> plyunote: ignore reflength here

      ibm%x_bp_i(i)=ibm%x_bp(i); ibm%y_bp_i(i)=ibm%y_bp(i);
      ibm%z_bp_i(i)=ibm%z_bp(i)

      ibm%x_bp(i) = ibm%x_bp(i) + fsi%x_c
      ibm%y_bp(i) = ibm%y_bp(i) + fsi%y_c
      ibm%z_bp(i) = ibm%z_bp(i) + fsi%z_c

      ibm%x_bp0(i)=ibm%x_bp(i); ibm%y_bp0(i)=ibm%y_bp(i);
      ibm%z_bp0(i)=ibm%z_bp(i)

      ibm%x_bp_o(i)=ibm%x_bp(i); ibm%y_bp_o(i)=ibm%y_bp(i);
      ibm%z_bp_o(i)=ibm%z_bp(i)

      ibm%u(i,1:3)=0.0_wp; ibm%uold(i,1:3)=0.0_wp; ibm%urm1(i,1:3)=0.0_wp
    enddo

    if (myid .eq. 0) then
      read(15008,*)
      read(15008,*)
      read(15008,*) n_elmt, ii
    endif
    call mpi_bcast(n_elmt, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)

    ibm%n_elmt = n_elmt

    ncmp = 3
    allocate(ibm%nv1(n_elmt), ibm%nv2(n_elmt), ibm%nv3(n_elmt))
    allocate(ibm%nf_x(n_elmt), ibm%nf_y(n_elmt), ibm%nf_z(n_elmt))
    allocate(ibm%nt_x(n_elmt), ibm%nt_y(n_elmt), ibm%nt_z(n_elmt))
    allocate(ibm%ns_x(n_elmt), ibm%ns_y(n_elmt), ibm%ns_z(n_elmt))
    allocate(ibm%dA(n_elmt))
    
    !allocate(nv1(n_elmt), nv2(n_elmt), nv3(n_elmt))
    !allocate(nf_x(n_elmt), nf_y(n_elmt), nf_z(n_elmt))
    !allocate(nt_x(n_elmt), nt_y(n_elmt), nt_z(n_elmt))
    !allocate(ns_x(n_elmt), ns_y(n_elmt), ns_z(n_elmt))
    !allocate(dA(n_elmt))

    allocate(ibm%cent_x(n_elmt), ibm%cent_y(n_elmt), ibm%cent_z(n_elmt))
    allocate(ibm%AOA_mean(n_elmt))
    allocate(ibm%Urelmag_mean(n_elmt), ibm%Urelmag(n_elmt))
    allocate(ibm%Urel_mean(n_elmt, ncmp), ibm%Uinduced_mean(n_elmt, ncmp))
    allocate(ibm%circulation_mean(n_elmt, ncmp))
    allocate(ibm%liftdirection_mean(n_elmt, ncmp))
    allocate(ibm%F_lagr_x(n_elmt), ibm%F_lagr_y(n_elmt), ibm%F_lagr_z(n_elmt))
    allocate(ibm%U_lagr_x(n_elmt), ibm%U_lagr_y(n_elmt), ibm%U_lagr_z(n_elmt))
    allocate(ibm%i_min(n_elmt), ibm%j_min(n_elmt), ibm%k_min(n_elmt))
    allocate(ibm%i_max(n_elmt), ibm%j_max(n_elmt), ibm%k_max(n_elmt))
    allocate(ibm%iclose(n_elmt), ibm%jclose(n_elmt), ibm%kclose(n_elmt))
    allocate(ibm%angle_attack(n_elmt), ibm%angle_twist(n_elmt), ibm%chord_blade(n_elmt))
    allocate(ibm%CD(n_elmt), ibm%CL(n_elmt))
    allocate(ibm%color(n_elmt), ibm%s2l(n_elmt))
    allocate(ibm%s2l1(n_elmt), ibm%s2l2(n_elmt), ibm%s2lc(n_elmt))
    allocate(ibm%Urel(n_elmt,ncmp), ibm%Uinduced(n_elmt,ncmp), ibm%circulation(n_elmt,ncmp))
    allocate(ibm%liftdirection(n_elmt,ncmp), ibm%rotationdirection(n_elmt,ncmp))

    if (myid .eq. 0) then
      do i = 1, n_elmt
        read(15008,*) ibm%nv1(i), ibm%nv2(i), ibm%nv3(i), ii, ibm%color(i)
      enddo
      color0 = ibm%color(1) - 1
      do i = 1, n_elmt
        ibm%color(i) = ibm%color(i) - color0  
      enddo
      close(15008)
    endif

    call mpi_bcast(ibm%nv1, n_elmt, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(ibm%nv2, n_elmt, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(ibm%nv3, n_elmt, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(ibm%color, n_elmt, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)

    do i = 1, n_elmt
      n1e = ibm%nv1(i); n2e = ibm%nv2(i); n3e = ibm%nv3(i)
      ds12(1) = ibm%x_bp(n2e) - ibm%x_bp(n1e)
      ds12(2) = ibm%y_bp(n2e) - ibm%y_bp(n1e)
      ds12(3) = ibm%z_bp(n2e) - ibm%z_bp(n1e)
      
      ds13(1) = ibm%x_bp(n3e) - ibm%x_bp(n1e)
      ds13(2) = ibm%y_bp(n3e) - ibm%y_bp(n1e)
      ds13(3) = ibm%z_bp(n3e) - ibm%z_bp(n1e)

      call crossx(ds12, ds13, tmparray)
      dr = norm2(tmparray)
      ibm%nf_x(i) = tmparray(1)/dr; ibm%nf_y(i) = tmparray(2)/dr 
      ibm%nf_z(i) = tmparray(3)/dr

      !> plyunote: in original version, streamwise direction is z.
      !! now it is x. so the function is a little different
      !! the projection is z->x, x->y, y->z
      if (  (((1.0-ibm%nf_x(i))<=1.0e-6).and.((-1.0+ibm%nf_x(i))<1.0e-6)) .or.&
        (((1.0+ibm%nf_x(i))<=1.0e-6).and.((-1.0-ibm%nf_x(i))<1.0e-6)) ) then
        ibm%ns_y(i) = 1.0_wp; ibm%ns_z(i) = 0.0_wp; ibm%ns_x(i) = 0.0_wp
        ibm%nt_y(i) = 0.0_wp; ibm%nt_z(i) = 1.0_wp; ibm%nt_x(i) = 0.0_wp
      else 
        rr = sqrt(ibm%nf_y(i)**2+ibm%nf_z(i)**2)
        ibm%ns_y(i) = ibm%nf_y(i) / rr
        ibm%ns_z(i) = -ibm%nf_y(i) / rr
        ibm%ns_x(i) = 0.0_wp
        
        ibm%nt_x(i) = rr 
        ibm%nt_y(i) = -ibm%nf_x(i) * ibm%nf_y(i) / rr
        ibm%nt_z(i) = -ibm%nf_x(i) * ibm%nf_z(i) / rr
      endif

      !> plyunote: for triangular elements, dA is correct?
      !! answer: yes. dA=0.5*|ds12|*|ds13|*sin(A213)=0.5*|ds12 x ds13|
      ibm%dA(i) = dr/2.0_wp

      ibm%cent_x(i) = (ibm%x_bp(n1e)+ibm%x_bp(n2e)+ibm%x_bp(n3e))/3.0_wp
      ibm%cent_y(i) = (ibm%y_bp(n1e)+ibm%y_bp(n2e)+ibm%y_bp(n3e))/3.0_wp
      ibm%cent_z(i) = (ibm%z_bp(n1e)+ibm%z_bp(n2e)+ibm%z_bp(n3e))/3.0_wp
    enddo

    do i = 1, n_elmt
      ibm%Urelmag_mean(i) = 0.0_wp;
      !ibm%Fr_mean(i) = 0.0_wp; ibm%Ft_mean(i) = 0.0_wp; ibm%Fa_mean(i) = 0.0_wp
      !ibm%Ur_mean(i) = 0.0_wp; ibm%Ut_mean(i) = 0.0_wp; ibm%Ua_mean(i) = 0.0_wp
      ibm%AOA_mean(i) = 0.0_wp
      ibm%Urel_mean(i,1:3) = 0.0_wp; ibm%Uinduced_mean(i,1:3) = 0.0_wp
      ibm%circulation_mean(i,1:3) = 0.0_wp; ibm%liftdirection_mean(i,1:3)=0.0_wp
    enddo
  end subroutine surface_read_xpatch
 
  subroutine calc_s2l(ibm_surface, ibm_line, fsi_surface, fsi_line, &
    NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm_surface, ibm_line
    type(FSInfo), dimension(:) :: fsi_surface, fsi_line
    integer :: NumberOfObjects

    integer :: ibi, elmt_s, elmt_l
    real(wp) :: rs_pjt, rl_amp, dd_min, dd, proj_close, rl_close
    real(wp), dimension(3) :: rs, rl, nl
    integer :: colors, colorl

    integer :: i0, i1, i2, n1
    real(wp) :: r1, r2, c

    character(128) :: filen

    do ibi = 1, NumberOfObjects
      n1 = ibm_line(ibi)%n_elmt / num_blade
      
      !> output the pairing of line and surface elements
      filen = ''
      if (myid .eq. 0) then
        write(filen, '(a8,i0.3,a4)') 'pairing_', ibi,'.dat'
        open(18004, file=filen, action='write')
        write(18004,*) 'i_surface i_line color distance radius_s projection &
          &i1 i2 c_in_elmt'
      endif
      !print *,'plyudebug, calc_s2l: myid=',myid,',surface_nelmt=',&
      !  ibm_surface(ibi)%n_elmt,',line_nelmt=',ibm_line(ibi)%n_elmt
      do elmt_s = 1, ibm_surface(ibi)%n_elmt
        dd_min = 100000.0_wp; proj_close = 100000.0_wp; rl_close = 10000.0_wp
        rs(1) = ibm_surface(ibi)%cent_x(elmt_s) - fsi_surface(ibi)%x_c
        rs(2) = ibm_surface(ibi)%cent_y(elmt_s) - fsi_surface(ibi)%y_c
        rs(3) = ibm_surface(ibi)%cent_z(elmt_s) - fsi_surface(ibi)%z_c

        colors = ibm_surface(ibi)%color(elmt_s)
        
        do elmt_l = 1, ibm_line(ibi)%n_elmt
          rl(1) = ibm_line(ibi)%cent_x(elmt_l) - fsi_line(ibi)%x_c
          rl(2) = ibm_line(ibi)%cent_y(elmt_l) - fsi_line(ibi)%y_c
          rl(3) = ibm_line(ibi)%cent_z(elmt_l) - fsi_line(ibi)%z_c
          rl_amp = norm2(rl)
          nl(1:3) = rl(1:3) / rl_amp
          rs_pjt = dabs(rs(1)*nl(1)+rs(2)*nl(2)+rs(3)*nl(3))
          colorl = ibm_line(ibi)%color(elmt_l)

          dd = dabs(rs_pjt - rl_amp)
          if( (dd-dd_min)<1.0e-9_wp .and. colors == colorl) then
            dd_min = dd;
            proj_close = rs_pjt
            rl_close = rl_amp
            ibm_surface(ibi)%s2l(elmt_s) = elmt_l
          endif
        enddo ! end of elmt_l

        i0 = ibm_surface(ibi)%s2l(elmt_s)

        !> find two nearest line elements
        if (i0 .eq. 1) then
          i1 = 1; i2 = 2
        else if (i0 .eq. ibm_line(ibi)%n_elmt) then
          i2 = ibm_line(ibi)%n_elmt; i1 = i2 - 1
        else if (proj_close >= rl_close .and. mod((i0+1),n1).ne.0) then
          i1 = i0; i2 = i0 + 1
        else if (proj_close >= rl_close .and. mod((i0+1),n1).eq.0) then
          i1 = i0 - 1; i2 = i0
        else if (proj_close < rl_close .and. mod(i0,n1).ne.1) then
          i1 = i0 - 1; i2 = i0          
        else if (proj_close < rl_close .and. mod(i0,n1).eq.1) then
          i1 = i0; i2 = i0 + 1
        else
          print *, 'plyudebug, calc_s2l: any other case for i1, i2 '
        endif 

        rl(1) = ibm_line(ibi)%cent_x(i1) - fsi_line(ibi)%x_c
        rl(2) = ibm_line(ibi)%cent_y(i1) - fsi_line(ibi)%y_c
        rl(3) = ibm_line(ibi)%cent_z(i1) - fsi_line(ibi)%z_c
        r1 = norm2(rl)

        rl(1) = ibm_line(ibi)%cent_x(i2) - fsi_line(ibi)%x_c
        rl(2) = ibm_line(ibi)%cent_y(i2) - fsi_line(ibi)%y_c
        rl(3) = ibm_line(ibi)%cent_z(i2) - fsi_line(ibi)%z_c
        r2 = norm2(rl)

        c = (proj_close - r1) / (r2 - r1)
        ibm_surface(ibi)%s2l1(elmt_s) = i1
        ibm_surface(ibi)%s2l2(elmt_s) = i2
        ibm_surface(ibi)%s2lc(elmt_s) = c

        if (myid .eq. 0) then
          !> here, you may find large value of dd_min, maybe the color
          !! sequence in acldata000 and acsdata000 are different
          if ((proj_close/(norm2(rs(1:3))+1.0e-9_wp))<0.2) then
            print *, 'plyudebug, calc_s2l: bad pairing, maybe blade color&
              & in acldata000 and acsdata000 are different'
          endif
          write(18004,*) elmt_s, ibm_surface(ibi)%s2l(elmt_s),&
            colors, dd_min, norm2(rs(1:3)), proj_close, i1, i2, c
        endif
      enddo ! end of elmt_s

      if (myid .eq. 0) close(18004)
    enddo ! end of ibi
  end subroutine calc_s2l
  
  subroutine calc_d2l(ibm_surface, ibm_line, fsi_surface, fsi_line, &
    NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm_surface, ibm_line
    type(FSInfo), dimension(:) :: fsi_surface, fsi_line
    integer :: NumberOfObjects

    integer :: ibi, elmt_s, elmt_l
    real(wp) :: rs_pjt, rl_amp, dd_min, dd, proj_close, rl_close
    real(wp), dimension(3) :: rs, rl, nl
    integer :: colors, colorl

    integer :: i0, i1, i2, n1
    real(wp) :: r1, r2, c

    character(128) :: filen

    do ibi = 1, NumberOfObjects
      n1 = ibm_line(ibi)%n_elmt / num_blade
      
      !> output the pairing of line and surface elements
      filen = ''
      if (myid .eq. 0) then
        write(filen, '(a8,i0.3,a4)') 'pairing_', ibi,'.dat'
        open(18004, file=filen, action='write')
        write(18004,*) 'i_surface i_line color distance radius_s projection &
          &i1 i2 c_in_elmt'
      endif
      !print *,1 
      !print *,'plyudebug, calc_d2l: myid=',myid,',surface_nelmt=',&
      !  ibm_surface(ibi)%n_elmt,',line_nelmt=',ibm_line(ibi)%n_elmt
      do elmt_s = 1, ibm_surface(ibi)%n_elmt
        !print *,'elmt_s=',elmt_s
        dd_min = 100000.0_wp; proj_close = 100000.0_wp; rl_close = 10000.0_wp
        rs(1) = ibm_surface(ibi)%cent_x(elmt_s) - fsi_surface(ibi)%x_c
        rs(2) = ibm_surface(ibi)%cent_y(elmt_s) - fsi_surface(ibi)%y_c
        rs(3) = ibm_surface(ibi)%cent_z(elmt_s) - fsi_surface(ibi)%z_c

        rs_pjt = norm2(rs)
        !colors = ibm_surface(ibi)%color(elmt_s)
        
        !do elmt_l = 1, ibm_line(ibi)%n_elmt
        do elmt_l = 1, n1
          rl(1) = ibm_line(ibi)%cent_x(elmt_l) - fsi_line(ibi)%x_c
          rl(2) = ibm_line(ibi)%cent_y(elmt_l) - fsi_line(ibi)%y_c
          rl(3) = ibm_line(ibi)%cent_z(elmt_l) - fsi_line(ibi)%z_c
          rl_amp = norm2(rl)
          !print *,2,'elmt_l=',elmt_l
          !nl(1:3) = rl(1:3) / rl_amp
          !rs_pjt = dabs(rs(1)*nl(1)+rs(2)*nl(2)+rs(3)*nl(3))
          !colorl = ibm_line(ibi)%color(elmt_l)

          dd = dabs(rs_pjt - rl_amp)
          !if( (dd-dd_min)<1.0e-9_wp .and. colors == colorl) then
          if( (dd-dd_min)<1.0e-9_wp ) then
            dd_min = dd;
            proj_close = rs_pjt
            rl_close = rl_amp
            ibm_surface(ibi)%s2l(elmt_s) = elmt_l
          endif
        enddo ! end of elmt_l
        
        i0 = ibm_surface(ibi)%s2l(elmt_s)
        !print*,3,'s2l(elmt_s=',elmt_s,')=',i0
        !> find two nearest line elements
        if (i0 .eq. 1) then
          i1 = 1; i2 = 2
        else if (i0 .eq. ibm_line(ibi)%n_elmt) then
          i2 = ibm_line(ibi)%n_elmt; i1 = i2 - 1
        else if (proj_close >= rl_close .and. mod(i0,n1).ne.0) then
          i1 = i0; i2 = i0 + 1
        else if (proj_close >= rl_close .and. mod(i0,n1).eq.0) then
          i1 = i0 - 1; i2 = i0
        else if (proj_close < rl_close .and. mod(i0,n1).ne.1) then
          i1 = i0 - 1; i2 = i0          
        else if (proj_close < rl_close .and. mod(i0,n1).eq.1) then
          i1 = i0; i2 = i0 + 1
        else
          print *, 'plyudebug, calc_s2l: any other case for i1, i2 '
        endif 
        !print*,4,'i1,i2=',i1,i2
        rl(1) = ibm_line(ibi)%cent_x(i1) - fsi_line(ibi)%x_c
        rl(2) = ibm_line(ibi)%cent_y(i1) - fsi_line(ibi)%y_c
        rl(3) = ibm_line(ibi)%cent_z(i1) - fsi_line(ibi)%z_c
        r1 = norm2(rl)

        rl(1) = ibm_line(ibi)%cent_x(i2) - fsi_line(ibi)%x_c
        rl(2) = ibm_line(ibi)%cent_y(i2) - fsi_line(ibi)%y_c
        rl(3) = ibm_line(ibi)%cent_z(i2) - fsi_line(ibi)%z_c
        r2 = norm2(rl)

        c = (proj_close - r1) / (r2 - r1)
        ibm_surface(ibi)%s2l1(elmt_s) = i1
        ibm_surface(ibi)%s2l2(elmt_s) = i2
        ibm_surface(ibi)%s2lc(elmt_s) = c
        !print*,5
        if (myid .eq. 0) then
          !> here, you may find large value of dd_min, maybe the color
          !! sequence in acldata000 and acsdata000 are different
          if ((proj_close/(norm2(rs(1:3))+1.0e-9_wp))<0.2) then
            print *, 'plyudebug, calc_s2l: bad pairing, maybe blade color&
              & in acldata000 and acsdata000 are different'
          endif
          write(18004,*) elmt_s, ibm_surface(ibi)%s2l(elmt_s),&
            0, dd_min, norm2(rs(1:3)), proj_close, i1, i2, c
        endif
      enddo ! end of elmt_s

      if (myid .eq. 0) close(18004)
    enddo ! end of ibi
  end subroutine calc_d2l
  
  subroutine Coordinates_IP(ibm, NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    integer :: NumberOfObjects

    integer :: ibi, l
    real(wp) :: dh

    do ibi = 1, NumberOfObjects
      do l = 1, ibm(ibi)%n_elmt
        !> plyunote: make sure dh is read from Nacelle.inp
        dh = ibm(ibi)%dh
        ibm(ibi)%dh_IP(l) = dh
        ibm(ibi)%centIP(l,1) = ibm(ibi)%cent_x(l) + dh * ibm(ibi)%nf_x(l)
        ibm(ibi)%centIP(l,2) = ibm(ibi)%cent_y(l) + dh * ibm(ibi)%nf_y(l)
        ibm(ibi)%centIP(l,3) = ibm(ibi)%cent_z(l) + dh * ibm(ibi)%nf_z(l)
      enddo
    enddo
  end subroutine Coordinates_IP

  !subroutine ColorIB_tmp(ibm, NumberOfObjects)
    ! yes, it is empty now.
    !PetscInt        i, j, k, l, ibi;
    !for (ibi=0; ibi<NumberOfObjects; ibi++) {
    !  for (l=0; l<ibm[ibi].n_elmt; l++) {
    !    if(dabs(ibm[ibi].cent_z[l]-(-0.03051))<1.e-3 &&dabs(dabs(ibm[ibi].nf_z[l])-1.0)<1.e-4)
    !      ibm[ibi].color[l]=1000;
    !    if(dabs(ibm[ibi].cent_z[l]-(1.5-0.03051))<1.e-3&&dabs(dabs(ibm[ibi].nf_z[l])-1.0)<1.e-4)
    !      ibm[ibi].color[l]=1000;
    !    if(dabs(ibm[ibi].cent_z[l]-0.0352)<1.e-3&&dabs(dabs(ibm[ibi].nf_z[l])-1.0)<1.e-4)
    !      ibm[ibi].color[l]=1000;
    !    if(dabs(ibm[ibi].nf_y[l])<1.e-3)
    !      ibm[ibi].color[l]=1000;
    !    if(dabs(ibm[ibi].cent_z[l]-(1.4695))<1.e-3&&dabs(dabs(ibm[ibi].nf_z[l])-1.0)<1.e-4)
    !      ibm[ibi].color[l]=1000;
    !  }
    !}
  !end subroutine ColorIB_tmp

  subroutine read_turbine_control(time_)
    implicit none
    real(wp) :: time_

    real(wp) :: rr
    integer:: ibi, idof
    logical :: l_exist

    first_step_time = time_
    ti = nint(time_/dt)
    ti_first = nint(first_step_time/dt)

    sqrt_pi = sqrt(pi)

    if(myid.eq.0) print *, 'plyudebug: read_acl/acs_param 1: collect_grid'
    
    call collect_grid

    if(myid.eq.0) print *, 'plyudebug: read_acl/acs_param 2: read control.dat'
    
    rotor_model = 3; NumberOfTurbines = 1; num_foiltype = 1
    num_blade = 3; loc_refvel = 2.0; FixTurbineAngvel = 1
    FixTipspeedRatio = 0; halfwidth_dfunc = 1.0
    nacelle_model = 1; rotate_nacelle = 0
    NumberOfNacelle = 1; NumNacellePerLoc = 1
    fsitype = 0

    if(myid .eq. 0) then
      !> read control parameters for wind turbine model
      open(15001, file='control.dat', action='read', status='old')
      read(15001,*) rotor_model
      read(15001,*) NumberOfTurbines, num_foiltype, num_blade
      read(15001,*) loc_refvel
      read(15001,*) FixTurbineAngvel, FixTipspeedRatio
      read(15001,*) URefDynamicAverage
      read(15001,*) halfwidth_dfunc
      !> read control parameters for nacelle model
      !> cf_nacelle_fromfile, rotate_nacelle, nacelle_model, reflength_nacelle
      !> r_nacelle, L_nacelle, dh_nacelle
      read(15001,*) nacelle_model
      read(15001,*) rotate_nacelle
      read(15001,*) NumberOfNacelle, NumNacellePerLoc
      read(15001,*) deltafunc_U, deltafunc_F
      read(15001,*) Shen1_AL, Shen1_AL_tipcorrection, Shen1_AL_tipcorrectionratio_Fa
      read(15001,*) correction3D
      read(15001,*) InletRelaxation, IR_bw, IR_dataindex, IR_scale, IR_ti_shift
      read(15001,*) c_bforce
      read(15001,*) adjust_CT_ACL, adjust_CP_ACL, remove_lateral_ACL_force
      read(15001,*) adjust_CT_nacelle
      read(15001,*) fsitype
      !read(15001,*) IR_folder 
      read(15001,*) output_wake_switch, output_wake_step, wake_rel_xmin, &
        wake_rel_xmax, wake_rel_y
      close(15001)
    endif      
    call mpi_bcast(rotor_model, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(NumberOfTurbines, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(num_foiltype, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(num_blade, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(loc_refvel, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(FixTurbineAngvel, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(FixTipspeedRatio, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(URefDynamicAverage, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(halfwidth_dfunc, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)

    call mpi_bcast(nacelle_model, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(rotate_nacelle, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(NumberOfNacelle, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(NumNacellePerLoc, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)

    call mpi_bcast(deltafunc_U, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(deltafunc_F, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    
    call mpi_bcast(Shen1_AL, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(Shen1_AL_tipcorrection, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(Shen1_AL_tipcorrectionratio_Fa, 1, mpi_double_precision, 0, &
      mpi_comm_2d_cart, ierr_wt)
    
    call mpi_bcast(correction3D, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    
    call mpi_bcast(InletRelaxation, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(IR_bw, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(IR_dataindex, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(IR_scale, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(IR_ti_shift, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    
    call mpi_bcast(c_bforce, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(adjust_CT_ACL, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(adjust_CP_ACL, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(remove_lateral_ACL_force, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(adjust_CT_nacelle, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    
    call mpi_bcast(fsitype, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    
    call mpi_bcast(output_wake_switch, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(output_wake_step, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(wake_rel_xmin, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(wake_rel_xmax, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(wake_rel_y, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    
    !print *, 'CPU',myid,':rotor_model=',rotor_model,&
    !  'NumberOfTurbines=',NumberOfTurbines

    !> some other control parameters
    !halfwidth_dfunc = 3.0
    !if(myid.eq.0)print*,'dx=',xl/nx,',halfwidth_delta=',halfwidth_dfunc*xl/nx
    Itpwidth = floor(halfwidth_dfunc) + 2
    !if(myid.eq.0) print *, 'plyudebug, Itpwidth=', Itpwidth

    allocate(fsi_wt(NumberOfTurbines))
    do ibi = 1, NumberOfTurbines
      fsi_wt(ibi)%angvel(1:3) = 0.0_wp;
      fsi_wt(ibi)%angvel_axis = 0.0_wp;
      fsi_wt(ibi)%ang0_yaw = 0.0_wp;
      fsi_wt(ibi)%ang0_rot = 0.0_wp;
    enddo
    
    !> plyutodo: it seems to handle with single turbine when reading acldata00
    !>   need to check whether all NumberOfTurbines*num_blade blades is included 
    !>   in this acldata00, or in several files acldataxx
    allocate(wtm(NumberOfTurbines)) !< n_v = n_v_single_blade * num_blade
    
    if(myid.eq.0) print *, 'plyudebug: read_acl/acs_param 3: read Turbine.inp'

    if(myid .eq. 0) then
      open(15002, file='Turbine.inp', action='read', status='old')
      read(15002,*)
      do ibi=1,NumberOfTurbines
        read(15002,*) fsi_wt(ibi)%nx_tb, fsi_wt(ibi)%ny_tb, fsi_wt(ibi)%nz_tb,&
          fsi_wt(ibi)%x_c, fsi_wt(ibi)%y_c, fsi_wt(ibi)%z_c, &
          fsi_wt(ibi)%Tipspeedratio, fsi_wt(ibi)%angvel_fixed, &
          fsi_wt(ibi)%r_rotor, wtm(ibi)%pitch(1), fsi_wt(ibi)%CT0, &
          fsi_wt(ibi)%CP0
      enddo
      close(15002)
    endif 

    if (myid .eq. 0) then
      !> plyunote: read or create Turbine2.inp, which has initial rotation info
      inquire(file='Turbine2.inp', exist=l_exist)
      if (l_exist) then
        print *, '!!!!!! Reading Turbine2.inp !!!!!!'
        open(15011, file='Turbine2.inp', action='read', status='old')
        do ibi = 1, NumberOfTurbines
          read(15011, *) fsi_wt(ibi)%ang0_rot, fsi_wt(ibi)%ang0_yaw, &
            fsi_wt(ibi)%angvel0 
        enddo
        close(15011)
      else
        print *, '!!!!!! Creating Turbine2.inp !!!!!!'
        open(15011, file='Turbine2.inp', action='write', status='new')
        do ibi = 1, NumberOfTurbines
          fsi_wt(ibi)%ang0_rot = 0.0_wp
          fsi_wt(ibi)%ang0_yaw = 0.0_wp
          fsi_wt(ibi)%angvel0 = fsi_wt(ibi)%angvel_fixed
          write(15011, *) fsi_wt(ibi)%ang0_rot, fsi_wt(ibi)%ang0_yaw, &
            fsi_wt(ibi)%angvel0, ' # ang0_rot, ang0_yaw, angvel0 for ',&
            'turbine ', ibi 
        enddo
        close(15011)
      endif
    endif
    
    do ibi=1,NumberOfTurbines
      
      call mpi_bcast(fsi_wt(ibi)%nx_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%ny_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%nz_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%x_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%y_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%z_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%Tipspeedratio, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)    
      call mpi_bcast(fsi_wt(ibi)%angvel_fixed, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%r_rotor, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(wtm(ibi)%pitch(1), 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%CP0, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%CT0, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      
      call mpi_bcast(fsi_wt(ibi)%ang0_rot, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%ang0_yaw, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_wt(ibi)%angvel0, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    
      rr = sqrt(fsi_wt(ibi)%nx_tb**2+fsi_wt(ibi)%ny_tb**2+fsi_wt(ibi)%nz_tb**2)
      fsi_wt(ibi)%nx_tb = fsi_wt(ibi)%nx_tb / rr;
      fsi_wt(ibi)%ny_tb = fsi_wt(ibi)%ny_tb / rr;
      fsi_wt(ibi)%nz_tb = fsi_wt(ibi)%nz_tb / rr;
      fsi_wt(ibi)%nx_tb0 = fsi_wt(ibi)%nx_tb;
      fsi_wt(ibi)%ny_tb0 = fsi_wt(ibi)%ny_tb;
      fsi_wt(ibi)%nz_tb0 = fsi_wt(ibi)%nz_tb;

      wtm(ibi)%Tipspeedratio = fsi_wt(ibi)%Tipspeedratio

      fsi_wt(ibi)%x_c0 = fsi_wt(ibi)%x_c
      fsi_wt(ibi)%y_c0 = fsi_wt(ibi)%y_c
      fsi_wt(ibi)%z_c0 = fsi_wt(ibi)%z_c

      fsi_wt(ibi)%ang_axis = fsi_wt(ibi)%ang0_rot
      fsi_wt(ibi)%angvel_axis = fsi_wt(ibi)%angvel0
    enddo
    
    if (rotor_model .eq. 3) then
      allocate(acl(num_foiltype))
      allocate(ibm_ACD(NumberOfTurbines))
    endif
    
    if (rotor_model .eq. 5 .or. rotor_model .eq. 7) then
      allocate(ibm_acl2ref(NumberOfTurbines))
      allocate(fsi_acl2ref(NumberOfTurbines))
      do ibi = 1, NumberOfTurbines
        fsi_acl2ref(ibi)%angvel(1:3) = 0.0_wp;
        fsi_acl2ref(ibi)%angvel_axis = 0.0_wp;

        fsi_acl2ref(ibi)%x_c = fsi_wt(ibi)%x_c
        fsi_acl2ref(ibi)%y_c = fsi_wt(ibi)%y_c
        fsi_acl2ref(ibi)%z_c = fsi_wt(ibi)%z_c
        
        fsi_acl2ref(ibi)%x_c0 = fsi_wt(ibi)%x_c0
        fsi_acl2ref(ibi)%y_c0 = fsi_wt(ibi)%y_c0
        fsi_acl2ref(ibi)%z_c0 = fsi_wt(ibi)%z_c0

        fsi_acl2ref(ibi)%nx_tb = fsi_wt(ibi)%nx_tb
        fsi_acl2ref(ibi)%ny_tb = fsi_wt(ibi)%ny_tb
        fsi_acl2ref(ibi)%nz_tb = fsi_wt(ibi)%nz_tb

        ibm_acl2ref(ibi)%Tipspeedratio = wtm(ibi)%Tipspeedratio
        fsi_acl2ref(ibi)%r_rotor = fsi_wt(ibi)%r_rotor
        fsi_acl2ref(ibi)%angvel_fixed = fsi_wt(ibi)%angvel_fixed
        ibm_acl2ref(ibi)%pitch(1) = wtm(ibi)%pitch(1)

      enddo
      allocate(acl(num_foiltype))
      allocate(ibm_ACD(NumberOfTurbines))
    endif

    if(nacelle_model .ne. 0) then
      !> plyunote: just an assumption
      NumLoc = NumberOfNacelle / NumNacellePerLoc
      if (NumLoc .ne. NumberOfTurbines) then
        print *, 'plyudebug: NumberOfNacelle/NumNacellePerLoc is not &
          &equal to NumberOfTurbines'
      endif

      allocate(ibm_nac(NumberOfNacelle), fsi_nac(NumberOfNacelle))
      do ibi = 1, NumberOfNacelle
        fsi_nac(ibi)%angvel(1:3) = 0.0_wp
        fsi_nac(ibi)%angvel_axis = 0.0_wp 
      enddo

    endif

    !> Prescribed motion, these parameters will not be writen to FSInfo
    !! of each turbine part. They act as global variables and predict 
    !! real-time loaction and velocity for each FSI element.
    if (fsitype .eq. 1) then
      allocate(pp_wt(NumberOfTurbines))
      if (myid .eq. 0) then
        open(15015, file='FSI_Prescribe.inp', action='read', status='old')
        do ibi = 1, NumberOfTurbines
          do idof = 1, 6
            read(15015,*) pp_wt(ibi)%dof(idof), pp_wt(ibi)%refpos(idof), &
              pp_wt(ibi)%amp(idof), pp_wt(ibi)%omega(idof), &
              pp_wt(ibi)%phase0(idof)
          enddo
        enddo
        close(15015)
      endif
      
      do ibi = 1, NumberOfTurbines
        call mpi_bcast(pp_wt(ibi)%dof, 6, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(pp_wt(ibi)%refpos, 6, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(pp_wt(ibi)%amp, 6, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(pp_wt(ibi)%omega, 6, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        call mpi_bcast(pp_wt(ibi)%phase0, 6, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
        
        do idof = 1, 6        
          if (pp_wt(ibi)%dof(idof) .ne. 0 .and. abs(pp_wt(ibi)%amp(idof)).le.1e-9) then
              print *, 'Warning: prescribed amplitude for static DOF should be &
                &zero. idof, amp=', idof, pp_wt(ibi)%amp(idof)
          endif
          
          !if (pp_wt(ibi)%dof(idof) .ne. 0) then
            pp_wt(ibi)%phase(idof) = pp_wt(ibi)%phase0(idof)+pp_wt(ibi)%omega(idof)*time_
            pp_wt(ibi)%disp0(idof) = pp_wt(ibi)%amp(idof) * cos(pp_wt(ibi)%phase0(idof))
            pp_wt(ibi)%disp(idof) = pp_wt(ibi)%amp(idof) * cos(pp_wt(ibi)%phase(idof))
            pp_wt(ibi)%vel(idof) = - pp_wt(ibi)%omega(idof) * pp_wt(ibi)%amp(idof) &
              * sin(pp_wt(ibi)%phase(idof))
            
            pp_wt(ibi)%dd(idof) = 0.0_wp

            if(ibi.eq.7 .and. myid.eq.0) then
              print *, 'read fsipara, myid, idof, phase0, phase, amp, disp=', myid, idof, &
                pp_wt(ibi)%phase0(idof), pp_wt(ibi)%phase(idof), pp_wt(ibi)%amp(idof), &
                pp_wt(ibi)%disp(idof)
            endif
          !endif
        enddo
      enddo

     !!> for prescribed motion, not necessary to reread fsiposition of last step 
     ! if (myid .eq. 0) then
     !   open(15013, file='FSI2_rotor.inp', action='read', status='old')
     !   do ibi = 1, NumberOfTurbines
     !     read(15013,*) fsi_wt(ibi)%x_c, fsi_wt(ibi)%y_c, fsi_wt(ibi)%z_c, &
     !       fsi_wt(ibi)%nx_tb, fsi_wt(ibi)%ny_tb, fsi_wt(ibi)%nz_tb
     !   enddo
     !   close(15013)
     ! endif
     ! 
     ! do ibi = 1, NumberOfTurbines
     !   call mpi_bcast(fsi_wt(ibi)%x_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
     !   call mpi_bcast(fsi_wt(ibi)%y_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
     !   call mpi_bcast(fsi_wt(ibi)%z_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
     !   call mpi_bcast(fsi_wt(ibi)%nx_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
     !   call mpi_bcast(fsi_wt(ibi)%ny_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
     !   call mpi_bcast(fsi_wt(ibi)%nz_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
     ! enddo
    endif

  end subroutine read_turbine_control
  
  subroutine read_nacelle_control(time_)
    implicit none
    real(wp) :: time_

    real(wp) :: rr
    integer:: ibi

    if(myid.eq.0) print *, 'plyudebug: read_nacelle_control: read Nacelle.inp'

    if(myid .eq. 0) then
      open(15009, file='Nacelle.inp', action='read', status='old')
      read(15009,*)
      do ibi=1,NumberOfNacelle
        read(15009,*) fsi_nac(ibi)%nx_tb, fsi_nac(ibi)%ny_tb, fsi_nac(ibi)%nz_tb,&
          fsi_nac(ibi)%x_c, fsi_nac(ibi)%y_c, fsi_nac(ibi)%z_c, &
          fsi_nac(ibi)%angvel_fixed, fsi_nac(ibi)%rotate_alongaxis, &
          ibm_nac(ibi)%friction_factor, ibm_nac(ibi)%dh, &
          fsi_nac(ibi)%xnacelle_upstreamend,&
          fsi_nac(ibi)%ynacelle_upstreamend, fsi_nac(ibi)%znacelle_upstreamend,&
          fsi_nac(ibi)%CT0, fsi_nac(ibi)%area_CT, fsi_nac(ibi)%ratio_CT0
      enddo
      close(15009)
    endif 
    
    do ibi=1,NumberOfNacelle      
      call mpi_bcast(fsi_nac(ibi)%nx_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%ny_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%nz_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%x_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%y_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%z_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%angvel_fixed, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)    
      call mpi_bcast(fsi_nac(ibi)%rotate_alongaxis, 1, mpi_integer, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm_nac(ibi)%friction_factor, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm_nac(ibi)%dh, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%xnacelle_upstreamend, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%ynacelle_upstreamend, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%znacelle_upstreamend, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%CT0, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%area_CT, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(fsi_nac(ibi)%ratio_CT0, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    
      rr = sqrt(fsi_nac(ibi)%nx_tb**2+fsi_nac(ibi)%ny_tb**2+fsi_nac(ibi)%nz_tb**2)
      fsi_nac(ibi)%nx_tb = fsi_nac(ibi)%nx_tb / rr;
      fsi_nac(ibi)%ny_tb = fsi_nac(ibi)%ny_tb / rr;
      fsi_nac(ibi)%nz_tb = fsi_nac(ibi)%nz_tb / rr;
      fsi_nac(ibi)%nx_tb0 = fsi_nac(ibi)%nx_tb;
      fsi_nac(ibi)%ny_tb0 = fsi_nac(ibi)%ny_tb;
      fsi_nac(ibi)%nz_tb0 = fsi_nac(ibi)%nz_tb;

      fsi_nac(ibi)%x_c0 = fsi_nac(ibi)%x_c
      fsi_nac(ibi)%y_c0 = fsi_nac(ibi)%y_c
      fsi_nac(ibi)%z_c0 = fsi_nac(ibi)%z_c

      fsi_nac(ibi)%angvel_axis = 0.0_wp
    enddo

   ! if (myid .eq. 0) then
   !   if(fsitype .eq. 1) then
   !     open(15014, file='FSI2_nacelle.inp', action='read', status='old')
   !     do ibi = 1, NumberOfNacelle
   !       read(15014,*) fsi_nac(ibi)%x_c, fsi_nac(ibi)%y_c, fsi_nac(ibi)%z_c, &
   !         fsi_nac(ibi)%nx_tb, fsi_nac(ibi)%ny_tb, fsi_nac(ibi)%nz_tb
   !     enddo
   !     close(15014)
   !   endif
   ! endif
   ! 
   ! do ibi = 1, NumberOfNacelle
   !   call mpi_bcast(fsi_nac(ibi)%x_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
   !   call mpi_bcast(fsi_nac(ibi)%y_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
   !   call mpi_bcast(fsi_nac(ibi)%z_c, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
   !   call mpi_bcast(fsi_nac(ibi)%nx_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
   !   call mpi_bcast(fsi_nac(ibi)%ny_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
   !   call mpi_bcast(fsi_nac(ibi)%nz_tb, 1, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
   ! enddo

  end subroutine read_nacelle_control

  subroutine read_acl_param(time_)
    implicit none
    real(wp) :: time_

    real(wp) :: rr
    integer:: ibi

    call read_turbine_control(time_)

    do ibi = 1, NumberOfTurbines
      fsi_wt(ibi)%rotor_model_type = 3
    enddo
    
    !> plyutodo: sending out these parameters via MPI
    
    if(myid.eq.0) print *, 'plyudebug: read_acl_param 4: read ACL_ucd'

    !> read turbineActuator Line grid
    !> it is unclear that why there is only one acldata 
    do ibi = 1,NumberOfTurbines
      call ACL_read_ucd(wtm(ibi), ibi, fsi_wt(ibi))
    enddo

    !call print_ibmnode1(7, 1, 100)

    !> Pre_process
    if(myid.eq.0) print *, 'plyudebug: read_acl_param 5: preprocess(wtm)'
    call Pre_process(wtm, NumberOfTurbines, 1)

    !> read Uref grid
    if(myid.eq.0) print *, 'plyudebug: read_acl_param 6: read Uref_ucd)'
    do ibi = 1,NumberOfTurbines
      !> plyunote: last parameter is ioffset, which should be added to
      !control.dat
      call disk_read_ucd(ibm_ACD(ibi), ibi, fsi_wt(ibi), 'Urefdata000', 1)
    enddo
    call reread_uref(wtm, fsi_wt, NumberOfTurbines)

    dh_mean = sqrt(dh_fsi*dh_uref)

    if(myid.eq.0) print *, 'plyudebug: read_acl_param 7: preprocess(Uref)'
    call Pre_process(ibm_ACD, NumberOfTurbines, 1)

    if(myid.eq.0) print *, 'plyudebug: read_acl_param 8: read ACL_coeff'
    call airfoil_ACL(wtm, fsi_wt) 

    if(myid .eq. 0) then
      !call Test_ArbitraryRotate
    endif
    
    if(nacelle_model .ne. 0) then
      call read_nacelle_param(time_)
    endif

  end subroutine read_acl_param

  subroutine read_acs_param(time_)
  implicit none
    real(wp) :: time_

    real(wp) :: rr
    integer:: ibi

    call read_turbine_control(time_)
    
    do ibi = 1, NumberOfTurbines
      fsi_wt(ibi)%rotor_model_type = 5
      fsi_acl2ref(ibi)%rotor_model_type = 3
    enddo
        
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 4: read_surface(wtm)'
    do ibi = 1,NumberOfTurbines
      call surface_read_xpatch(wtm(ibi), ibi, fsi_wt(ibi))
    enddo

    !> Pre_process
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 5: preprocess(wtm)'
    call Pre_process(wtm, NumberOfTurbines, 1)

    !> read Uref grid
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 6: read Uref_ucd)'
    do ibi = 1,NumberOfTurbines
      !> plyunote: last parameter is ioffset, which should be added to
      !control.dat
      call disk_read_ucd(ibm_ACD(ibi), ibi, fsi_wt(ibi), 'Urefdata000',1)
    enddo
    
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 7: preprocess(Uref)'
    call Pre_process(ibm_ACD, NumberOfTurbines, 1)
    
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 8: read &
      &ACL_ucd(ibm_acl2ref)'

    !> read turbineActuator (reference) Line grid
    do ibi = 1,NumberOfTurbines
      call ACL_read_ucd(ibm_acl2ref(ibi), ibi, fsi_acl2ref(ibi))
    enddo

    !> Pre_process
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 9: preprocess(ibm_acl2ref)'
    call Pre_process(ibm_acl2ref, NumberOfTurbines, 1)

    if(myid.eq.0) print *, 'plyudebug: read_acs_param 10: pairing(s2l)'
    call calc_s2l(wtm, ibm_acl2ref, fsi_wt, fsi_acl2ref, NumberOfTurbines)

    if(myid.eq.0) print *, 'plyudebug: read_acs_param 11: read ACL_coeff'
    call airfoil_ACL(wtm, fsi_wt)
    call airfoil_ACL(ibm_acl2ref, fsi_acl2ref) 

    if(myid .eq. 0) then
      !call Test_ArbitraryRotate
    endif

    if(nacelle_model .ne. 0) then
      call read_nacelle_param(time_)
    endif

  end subroutine read_acs_param
  
  subroutine read_admr_param(time_)
    implicit none
    real(wp) :: time_

    real(wp) :: rr
    integer:: ibi

    call read_turbine_control(time_)
    
    do ibi = 1, NumberOfTurbines
      fsi_wt(ibi)%rotor_model_type = 7 
      fsi_acl2ref(ibi)%rotor_model_type = 3
    enddo
        
    if(myid.eq.0) print *, 'plyudebug: read_admr_param 4: read_surface(wtm)'
    !> plyunote: for rotor_model 5, Actuator Surface Model, geometry input file
    !  is xpatch because it has color property
    !> plyunote: for rotor_model 7, Actuator Disk Model with rotation, although
    !  many processes are same with ASM, the input file can be generated
    !  separately, or reuse reference disk input file.
    do ibi = 1,NumberOfTurbines
      call disk_read_ucd(wtm(ibi), ibi, fsi_wt(ibi), 'admrdata000', 0)
      !call surface_read_xpatch(wtm(ibi), ibi, fsi_wt(ibi))
    enddo
    call reread_uref(wtm, fsi_wt, NumberOfTurbines)
    
    !> Pre_process
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 5: preprocess(wtm)'
    call Pre_process(wtm, NumberOfTurbines, 1)

    !> read Uref grid
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 6: read Uref_ucd)'
    do ibi = 1,NumberOfTurbines
      !> plyunote: last parameter is ioffset, which should be added to
      !control.dat
      call disk_read_ucd(ibm_ACD(ibi), ibi, fsi_wt(ibi), 'Urefdata000',1)
    enddo
    
    
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 7: preprocess(Uref)'
    call Pre_process(ibm_ACD, NumberOfTurbines, 1)
    
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 8: read &
      &ACL_ucd(ibm_acl2ref)'

    !> read turbineActuator (reference) Line grid
    do ibi = 1,NumberOfTurbines
      call ACL_read_ucd(ibm_acl2ref(ibi), ibi, fsi_acl2ref(ibi))
    enddo
    call reread_uref(ibm_acl2ref, fsi_acl2ref, NumberOfTurbines)

    !> Pre_process
    if(myid.eq.0) print *, 'plyudebug: read_acs_param 9: preprocess(ibm_acl2ref)'
    call Pre_process(ibm_acl2ref, NumberOfTurbines, 1)

    if(myid.eq.0) print *, 'plyudebug: read_admr_param 10: pairing(s2l)'
    call calc_d2l(wtm, ibm_acl2ref, fsi_wt, fsi_acl2ref, NumberOfTurbines)

    if(myid.eq.0) print *, 'plyudebug: read_acs_param 11: read ACL_coeff'
    call airfoil_ACL(wtm, fsi_wt)
    call airfoil_ACL(ibm_acl2ref, fsi_acl2ref) 

    if(myid .eq. 0) then
      !call Test_ArbitraryRotate
    endif

    if(nacelle_model .ne. 0) then
      call read_nacelle_param(time_)
    endif

  end subroutine read_admr_param
  
  subroutine read_nacelle_param(time_)
    implicit none
    real(wp) :: time_

    real(wp) :: rr
    integer:: iloc, iper, inac
    character(128) :: fname

    if(myid.eq.0) print *, "plyudebug: read_nacelle_param, read_nacelle_control"
    call read_nacelle_control(time_)
    
    if(myid.eq.0) print *, "plyudebug: read_nacelle_param, disk_read_ucd"
    do iloc = 1, NumLoc
      do iper = 1, NumNacellePerLoc
        inac = (iloc - 1) * NumNacellePerLoc + iper
        write(fname, '(a7,i0.3)') 'nacelle', iper-1
        call disk_read_ucd(ibm_nac(inac), inac, fsi_nac(inac), fname, 0)
      enddo
    enddo  
    
    if(myid.eq.0) print *, "plyudebug: read_nacelle_param, Pre_process 1"
    call Pre_process(ibm_nac, NumberOfNacelle, 1)
    
    if(myid.eq.0) print *, "plyudebug: read_nacelle_param, Coordinates_IP"
    call Coordinates_IP(ibm_nac, NumberOfNacelle)
    
    !> Pre_process_IP is integrated into Pre_process
    if(myid.eq.0) print *, "plyudebug: read_nacelle_param, Pre_process 2"
    call Pre_process(ibm_nac, NumberOfNacelle, 2)

    !call ColorIB_tmp(ibm_nac, NumberOfNacelle)
    
    if(myid.eq.0) print*, "plyudebug: read_nacelle_param: completed" 
        
  end subroutine read_nacelle_param
  
  !> Interpolate the velocity at the Lagrangian points
  subroutine Calc_U_lagr(u, v, w, ibm, NumberOfObjects, level)
  implicit none
    real(wp), dimension(:,:,:) :: u(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: v(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: w(xsz(1), xsz(2), 1-level:xsz(3)+level)
    type(IBMNodes), dimension(:) :: ibm
    !> here fsi is need if handling with perioidc boundary condition
    !type(FSInfo) , dimension(:):: fsi
    integer :: NumberOfObjects, level

    integer :: ibi, l, i, j, k, ii, jj, kk
    integer :: xs, ys, zs, xe, ye, ze 
    real(wp) :: xc, yc, zc, zcw,  cent_x, cent_y, cent_z
    real(wp) :: dhx_, dhy_, dhz_
    real(wp) :: rxi, ryi, rzi, r1i, r2i, r3i
    real(wp) :: rxk, ryk, rzk, r1k, r2k, r3k
    real(wp) :: dfunci, dfunck, df1, df2, df3, dh1, dh2, dh3, dV
    !> dxdq(1,1:3) = DxDcsi, DxDphi, DxDzeta
    !> dxdq(2,1:3) = DyDcsi, DyDphi, DyDzeta
    !> dxdq(3,1:3) = DzDcsi, DzDphi, DzDzeta
    
    !> dqdx(1,1:3) = DcsiDx, DcsiDy, DcsiDz
    !> dqdx(2,1:3) = DphiDx, DphiDy, DphiDz
    !> dqdx(3,1:3) = DzetaDx, DzetaDy, DzetaDz
    
    real(wp), dimension(3,3) :: dxdq, dqdx

    real(wp), allocatable, dimension(:) :: u_local, u_sum
    real(wp) :: u_max, u_temp

    !dh0 = sqrt(xl/nx * dh_mean)
    !dh0 = 1.5*max(xl/nx, dh_mean)
    !dh0 = sqrt(xl/nx*yl/ny*zl/nz)
    !if(myid.eq.0) print*,'dh0=',dh0

    xs = xst(1); xe = xst(1) + xsz(1) - 1
    ys = xst(2); ye = xst(2) + xsz(2) - 1
    zs = xst(3); ze = xst(3) + xsz(3) - 1

    !if(myid.eq.0) print *, eo(:,1)

    do ibi = 1,NumberOfObjects
      ibm(ibi)%U_lagr_x=0.0; ibm(ibi)%U_lagr_y=0.0; ibm(ibi)%U_lagr_z=0.0
    enddo

    do ibi = 1,NumberOfObjects
      do l = 1, ibm(ibi)%n_elmt
        cent_x = ibm(ibi)%cent_x(l); cent_y = ibm(ibi)%cent_y(l)
        cent_z = ibm(ibi)%cent_z(l)
        
        do k = ibm(ibi)%k_min(l), ibm(ibi)%k_max(l)
        do j = ibm(ibi)%j_min(l), ibm(ibi)%j_max(l)
        do i = ibm(ibi)%i_min(l), ibm(ibi)%i_max(l)
          if (i>=xs .and. i<xe .and. j>=ys .and. j<=ye &
            .and. k>=zs .and. k<=ze) then 
            !> here xc,yc,zc should be location of velocity variables
            !> plyunote: check whether staggeered grid is used here
            xc = cartx(i)
            yc = carty(j)
            !cartz(i,j,k) = zzo(k) * (hbar + hho(i,j)) - hho(i,j) 
            !zc = cartz(i,j,k)
            !zcw = cartzw(i,j,k)
            zc = zzo(k) * (hbar + hho(i,j)) - hho(i,j)
            zcw = zwo(k) * (hbar + hho(i,j)) - hho(i,j)

            !> plyunote: this section is implemented in a different form from
            !> fotis's version
            !dxdq(1,1) = 1.; dxdq(1,2) = 0.; dxdq(1,3) = 0.
            !dxdq(2,1) = 0.; dxdq(2,2) = 1.; dxdq(2,3) = 0.
            !dxdq(3,1) = hxo(i,j)*(-1.0+zzo(k)); dxdq(3,2) = hyo(i,j)*(-1.0+zzo(k))
            !dxdq(3,3) = hbar + eo(i,j)

            !dhx_ = sqrt(dxdq(1,1)**2+dxdq(2,1)**2+dxdq(3,1)**2)
            !dhy_ = sqrt(dxdq(1,2)**2+dxdq(2,2)**2+dxdq(3,2)**2)
            !dhz_ = sqrt(dxdq(1,3)**2+dxdq(2,3)**2+dxdq(3,3)**2)
            

            rxi = xc - cent_x; ryi = yc - cent_y; rzi = zc - cent_z
            rxk = xc - cent_x; ryk = yc - cent_y; rzk = zcw - cent_z
            dh1 = xl/nx ! max(xl/nx, dh_mean) 
            !r1 = (rx*dxdq(1,1) + ry*dxdq(2,1) + rz*dxdq(3,1))/dhx_/dhx_/dh1
            dh2 = yl/ny !max(yl/ny, dh_mean)
            !r2 = (rx*dxdq(1,2) + ry*dxdq(2,2) + rz*dxdq(3,2))/dhy_/dhy_/dh2
            dh3 = dz(k-xst(3)+1) * (hbar+eo(i,j)) !max(zl/(nz-1),dh_mean) 
            !r3 = (rx*dxdq(1,3) + ry*dxdq(2,3) + rz*dxdq(3,3))/dhz_/dhz_/dh3

            r1i = rxi / dh1; r2i = ryi / dh2; r3i = rzi / dh3
            r1k = rxk / dh1; r2k = ryk / dh2; r3k = rzk / dh3

            !> delta function: use dfunc_s4h in Yang(JCP,2009)
            !! doi:10.1016/j.jcp.2009.07.023
            
            if (deltafunc_U .eq. 1) then
              call dfunc_s1h(r1i, df1); call dfunc_s1h(r2i,df2); call dfunc_s1h(r3i,df3)
              dfunci = df1 * df2 * df3
              call dfunc_s1h(r1k, df1); call dfunc_s1h(r2k,df2); call dfunc_s1h(r3k,df3)
              dfunck = df1 * df2 * df3
            else if (deltafunc_U .eq. 2) then
              call dfunc_s2h(r1i, df1); call dfunc_s2h(r2i,df2); call dfunc_s2h(r3i,df3)
              dfunci = df1 * df2 * df3
              call dfunc_s2h(r1k, df1); call dfunc_s2h(r2k,df2); call dfunc_s2h(r3k,df3)
              dfunck = df1 * df2 * df3
            else if (deltafunc_U .eq. 3) then
              call dfunc_s3h(r1i, df1); call dfunc_s3h(r2i,df2); call dfunc_s3h(r3i,df3)
              dfunci = df1 * df2 * df3
              call dfunc_s3h(r1k, df1); call dfunc_s3h(r2k,df2); call dfunc_s3h(r3k,df3)
              dfunck = df1 * df2 * df3
            else if (deltafunc_U .eq. 4) then
              call dfunc_s4h(r1i, df1); call dfunc_s4h(r2i,df2); call dfunc_s4h(r3i,df3)
              dfunci = df1 * df2 * df3
              call dfunc_s4h(r1k, df1); call dfunc_s4h(r2k,df2); call dfunc_s4h(r3k,df3)
              dfunck = df1 * df2 * df3
            else
              print *, 'deltafunc_U: should be 2, 3, or 4'
            endif

            dV = 1.0_wp
            !dV = dh1 * dh2 * dh3
                  
            ii = i - xst(1) + 1
            jj = j - xst(2) + 1
            kk = k - xst(3) + 1
            !> plyunote: please check U_lagr_x initialization
            ibm(ibi)%U_lagr_x(l) = ibm(ibi)%U_lagr_x(l) + u(ii,jj,kk)*dfunci*dV
            ibm(ibi)%U_lagr_y(l) = ibm(ibi)%U_lagr_y(l) + v(ii,jj,kk)*dfunci*dV  
            ibm(ibi)%U_lagr_z(l) = ibm(ibi)%U_lagr_z(l) + w(ii,jj,kk)*dfunck*dV
            
            !if (dfunc>1.0d-3 .and. ibi .eq. 2 .and. l .eq.506 ) then
            !print*, zc,cent_z,rz,dh3,r3
            !print*,'plyudebug,l_elmt=',l,'/',ibm(ibi)%n_elmt,',i,j,k=',i,j,k,&
            !  ',ii,jj,kk=',ii,jj,kk
            !print*,'plyudebug,xc=',xc,',yc=',yc,',zc=',zc
            !print*,'plyudebug,rx=',rx,',ry=',ry,',rz=',rz
            !print*,'plyudebug,dhx_=',dhx_,',dhy_=',dhy_,',dhz_=',dhz_
            !print*,'plyudebug,df1=',df1,',df2=',df2,',df3=',df3
            
            !print*,'plyudebug,r1,r2,r3:',r1,r2,r3,dfunc
            !print*,'plyudebug,u=',u(ii,jj,kk),',v=',v(ii,jj,kk),',w=',w(ii,jj,kk) 
            !print*,'plyudebug,dh1,dh2,dh3=',dh1,dh2,dh3,',dfunc=',dfunc
            !print*,'plyudebug,U_lagr=',ibm(ibi)%U_lagr_x(l),ibm(ibi)%U_lagr_y(l),ibm(ibi)%U_lagr_z(l)  
            !endif
          endif
        enddo
        enddo
        enddo
      enddo
    enddo

    !> mpi
    do ibi = 1, NumberOfObjects
      allocate(u_local(ibm(ibi)%n_elmt*3), u_sum(ibm(ibi)%n_elmt*3))
      do i = 1, ibm(ibi)%n_elmt
        u_local(3*i-2) = ibm(ibi)%U_lagr_x(i) 
        u_local(3*i-1) = ibm(ibi)%U_lagr_y(i) 
        u_local(3*i) = ibm(ibi)%U_lagr_z(i) 
      enddo

      !> plyunote: mpi_comm_world or mpi_comm_2d_cart?
      call mpi_allreduce(u_local, u_sum, ibm(ibi)%n_elmt*3, &
        mpi_double_precision, mpi_sum, mpi_comm_world, ierr_wt)
      
      u_max = 0.0_wp; u_temp = 0.0_wp
      do i = 1, ibm(ibi)%n_elmt
        ibm(ibi)%U_lagr_x(i) = u_sum(3*i-2)
        ibm(ibi)%U_lagr_y(i) = u_sum(3*i-1)
        ibm(ibi)%U_lagr_z(i) = u_sum(3*i)
        u_temp = sqrt(u_sum(3*i-2)**2+u_sum(3*i-1)**2+u_sum(3*i)**2)
        if (u_max<u_temp) then
          u_max = u_temp
        endif
      enddo
      !if(myid .eq. 0) print *, 'plyudebug, ti=', ti, ', u_max=', u_max

      if(allocated(u_local)) deallocate(u_local)
      if(allocated(u_sum)) deallocate(u_sum)
    enddo
  end subroutine Calc_U_lagr
  
  !> Interpolate the velocity at the Lagrangian interpolation points
  !> It has the same structure with Calc_U_lagr
  !> If possible, I would like to user pointer to implement this
  !> However, fortran pointer can not point to private member of derived type
  subroutine Calc_U_IPlagr(u, v, w, ibm, NumberOfObjects, level)
  implicit none
    real(wp), dimension(:,:,:) :: u(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: v(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: w(xsz(1), xsz(2), 1-level:xsz(3)+level)
    type(IBMNodes), dimension(:) :: ibm
    !> here fsi is need if handling with perioidc boundary condition
    !type(FSInfo) , dimension(:):: fsi
    integer :: NumberOfObjects, level

    integer :: ibi, l, i, j, k, ii, jj, kk
    real(wp) :: xc, yc, zc, cent_x, cent_y, cent_z
    real(wp) :: dhx_, dhy_, dhz_, rx, ry, rz, r1, r2, r3
    real(wp) :: dfunc, df1, df2, df3, dh1, dh2, dh3, dV
    !> dxdq(1,1:3) = DxDcsi, DxDphi, DxDzeta
    !> dxdq(2,1:3) = DyDcsi, DyDphi, DyDzeta
    !> dxdq(3,1:3) = DzDcsi, DzDphi, DzDzeta
    
    !> dqdx(1,1:3) = DcsiDx, DcsiDy, DcsiDz
    !> dqdx(2,1:3) = DphiDx, DphiDy, DphiDz
    !> dqdx(3,1:3) = DzetaDx, DzetaDy, DzetaDz
    
    real(wp), dimension(3,3) :: dxdq, dqdx

    real(wp), allocatable, dimension(:) :: u_local, u_sum
    real(wp) :: u_max, u_temp

    !dh0 = sqrt(xl/nx * dh_mean)
    !dh0 = 1.5*max(xl/nx, dh_mean)
    !dh0 = sqrt(xl/nx*yl/ny*zl/nz)
    !if(myid.eq.0) print*,'dh0=',dh0

    do ibi = 1,NumberOfObjects
      ibm(ibi)%U_IPlagr(:,:) = 0.0
    enddo

    do ibi = 1,NumberOfObjects
      do l = 1, ibm(ibi)%n_elmt
        cent_x = ibm(ibi)%centIP(l,1); cent_y = ibm(ibi)%centIP(l,2)
        cent_z = ibm(ibi)%centIP(l,3)
        
        do k = ibm(ibi)%kIP_min(l), ibm(ibi)%kIP_max(l)
        do j = ibm(ibi)%jIP_min(l), ibm(ibi)%jIP_max(l)
        do i = ibm(ibi)%iIP_min(l), ibm(ibi)%iIP_max(l)
                !> here xc,yc,zc should be location of velocity variables
          !> plyunote: check whether staggeered grid is used here
          xc = cartx(i)
          yc = carty(j)
          !zc = cartz(i,j,k)
          zc = zzo(k) * (hbar + hho(i,j)) - hho(i,j)

          !> plyunote: this section is implemented in a different form from
          !> fotis's version
          dxdq(1,1) = 1.; dxdq(1,2) = 0.; dxdq(1,3) = 0.
          dxdq(2,1) = 0.; dxdq(2,2) = 1.; dxdq(2,3) = 0.
          dxdq(3,1) = hxo(i,j)*(-1.0+zzo(k)); dxdq(3,2) = hyo(i,j)*(-1.0+zzo(k))
          dxdq(3,3) = hbar + eo(i,j)

          dhx_ = sqrt(dxdq(1,1)**2+dxdq(2,1)**2+dxdq(3,1)**2)
          dhy_ = sqrt(dxdq(1,2)**2+dxdq(2,2)**2+dxdq(3,2)**2)
          dhz_ = sqrt(dxdq(1,3)**2+dxdq(2,3)**2+dxdq(3,3)**2)
          

          rx = xc - cent_x; ry = yc - cent_y; rz = zc - cent_z
          dh1 = xl/nx ! max(xl/nx, dh_mean) 
          r1 = (rx*dxdq(1,1) + ry*dxdq(2,1) + rz*dxdq(3,1))/dhx_/dhx_/dh1
          dh2 = yl/ny !max(yl/ny, dh_mean)
          r2 = (rx*dxdq(1,2) + ry*dxdq(2,2) + rz*dxdq(3,2))/dhy_/dhy_/dh2
          dh3 = dz(k-xst(3)+1) !max(zl/(nz-1),dh_mean) 
          r3 = (rx*dxdq(1,3) + ry*dxdq(2,3) + rz*dxdq(3,3))/dhz_/dhz_/dh3

                !> delta function: use dfunc_s4h in Yang(JCP,2009)
          !! doi:10.1016/j.jcp.2009.07.023
          if (deltafunc_U .eq. 1) then
            call dfunc_s1h(r1, df1); call dfunc_s1h(r2,df2); call dfunc_s1h(r3,df3)
          else if (deltafunc_U .eq. 2) then
            call dfunc_s2h(r1, df1); call dfunc_s2h(r2,df2); call dfunc_s2h(r3,df3)
          else if (deltafunc_U .eq. 3) then
            call dfunc_s3h(r1, df1); call dfunc_s3h(r2,df2); call dfunc_s3h(r3,df3)
          else if (deltafunc_U .eq. 4) then
            call dfunc_s4h(r1, df1); call dfunc_s4h(r2,df2); call dfunc_s4h(r3,df3)
          else
            print *, 'deltafunc_U: should be 2, 3, or 4'
          endif
         
          dfunc = df1 * df2 * df3 
          dV = 1.0 !dV = dh1 * dh2 * dh3
                
          ii = i - xst(1) + 1
          jj = j - xst(2) + 1
          kk = k - xst(3) + 1
          !> plyunote: please check U_lagr_x initialization
          ibm(ibi)%U_IPlagr(l,1) = ibm(ibi)%U_IPlagr(l,1) + u(ii,jj,kk)*dfunc*dV
          ibm(ibi)%U_IPlagr(l,2) = ibm(ibi)%U_IPlagr(l,2) + v(ii,jj,kk)*dfunc*dV  
          ibm(ibi)%U_IPlagr(l,3) = ibm(ibi)%U_IPlagr(l,3) + w(ii,jj,kk)*dfunc*dV
          
          !if (dfunc>1.0d-3 .and. ibi .eq. 2 .and. l .eq.506 ) then
          !print*, zc,cent_z,rz,dh3,r3
          !print*,'plyudebug,l_elmt=',l,'/',ibm(ibi)%n_elmt,',i,j,k=',i,j,k,&
          !  ',ii,jj,kk=',ii,jj,kk
          !print*,'plyudebug,xc=',xc,',yc=',yc,',zc=',zc
          !print*,'plyudebug,rx=',rx,',ry=',ry,',rz=',rz
          !print*,'plyudebug,dhx_=',dhx_,',dhy_=',dhy_,',dhz_=',dhz_
          !print*,'plyudebug,df1=',df1,',df2=',df2,',df3=',df3
          
          !print*,'plyudebug,r1,r2,r3:',r1,r2,r3,dfunc
          !print*,'plyudebug,u=',u(ii,jj,kk),',v=',v(ii,jj,kk),',w=',w(ii,jj,kk) 
          !print*,'plyudebug,dh1,dh2,dh3=',dh1,dh2,dh3,',dfunc=',dfunc
          !print*,'plyudebug,U_lagr=',ibm(ibi)%U_lagr_x(l),ibm(ibi)%U_lagr_y(l),ibm(ibi)%U_lagr_z(l)  
          !endif
        enddo
        enddo
        enddo
      enddo
    enddo

    !> mpi
    do ibi = 1, NumberOfObjects
      allocate(u_local(ibm(ibi)%n_elmt*3), u_sum(ibm(ibi)%n_elmt*3))
      do i = 1, ibm(ibi)%n_elmt
        u_local(3*i-2) = ibm(ibi)%U_IPlagr(i,1) 
        u_local(3*i-1) = ibm(ibi)%U_IPlagr(i,2) 
        u_local(3*i) = ibm(ibi)%U_IPlagr(i,3) 
      enddo

      !> plyunote: mpi_comm_world or mpi_comm_2d_cart?
      call mpi_allreduce(u_local, u_sum, ibm(ibi)%n_elmt*3, &
        mpi_double_precision, mpi_sum, mpi_comm_world, ierr_wt)
      
      u_max = 0.0_wp; u_temp = 0.0_wp
      do i = 1, ibm(ibi)%n_elmt
        ibm(ibi)%U_IPlagr(i,1) = u_sum(3*i-2)
        ibm(ibi)%U_IPlagr(i,2) = u_sum(3*i-1)
        ibm(ibi)%U_IPlagr(i,3) = u_sum(3*i)
        u_temp = sqrt(u_sum(3*i-2)**2+u_sum(3*i-1)**2+u_sum(3*i)**2)
        if (u_max<u_temp) then
          u_max = u_temp
        endif
      enddo
      if(myid .eq. 0) print *, 'plyudebug, ti=', ti, ', u_max=', u_max

      if(allocated(u_local)) deallocate(u_local)
      if(allocated(u_sum)) deallocate(u_sum)
    enddo
  end subroutine Calc_U_IPlagr

  !> calculating the reference velocity for actuator line model
  subroutine Uref_ACL(u, v, w, ibm_ACL, NumberOfObjects, level)
  implicit none
    real(wp), dimension(:,:,:) :: u(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: v(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: w(xsz(1), xsz(2), 1-level:xsz(3)+level)
    type(IBMNodes), dimension(:) :: ibm_ACL
    integer :: NumberOfObjects, level

    integer :: l, ibi
    real(wp) :: A_Sum, U_Sum
    real(wp) :: nx, ny, nz, U_axis

    !< tc: exponential time coefficient for relaxation process
    !  tc=T*u_{*}/H=0.27 corresponds to T=10min using possible reference
    !  dimensional values u_{*}=0.45, and H=1000m
    real(wp) :: tc, fexp

    call Calc_U_lagr(u, v, w, ibm_ACD, NumberOfObjects, level) 

    do ibi = 1,NumberOfObjects
      nx = fsi_wt(ibi)%nx_tb; ny = fsi_wt(ibi)%ny_tb; nz = fsi_wt(ibi)%nz_tb
      U_sum = 0.0; A_sum = 0.0
      do l = 1,ibm_ACD(ibi)%n_elmt
        U_axis = ibm_ACD(ibi)%U_lagr_x(l)*nx + ibm_ACD(ibi)%U_lagr_y(l)*ny &
          + ibm_ACD(ibi)%U_lagr_z(l)*nz
        !> check value of dA, for disk, it is 0.5*dabs(ab x ac)
        U_Sum = U_Sum + U_axis * ibm_ACD(ibi)%dA(l)
        A_Sum = A_Sum + ibm_ACD(ibi)%dA(l)
        !if (myid .eq. 0 .and. ibi .eq. 1) then
        !  print *, 'plyudebug, Uref_ACL: ACD elmt ',l,', U_lagr=',&
        !    ibm_ACD(ibi)%U_lagr_x(l), ibm_ACD(ibi)%U_lagr_y(l),&
        !    ibm_ACD(ibi)%U_lagr_z(l), ', dA=', ibm_ACD(ibi)%dA(l)
        !endif
      enddo
      ibm_ACL(ibi)%U1_ref = U_Sum / A_Sum
    enddo

    if (URefDynamicAverage .eq. 0) then
      do ibi = 1, NumberOfObjects
        ibm_ACL(ibi)%U_ref = ibm_ACL(ibi)%U1_ref
      enddo
    elseif (URefDynamicAverage .eq. 1) then
      !< tc: relaxation time window for temporal averaging (calaf et al. 2010)
      if (ti < (0.15 * 0.27*hbar/usbot/dt)) then
        tc = 0.15 * 0.27
      elseif (ti < (0.45 * 0.27*hbar/usbot/dt)) then
        tc = 0.45 * 0.27
      else
        tc = 0.27
      endif
      fexp=dt/(tc*hbar/usbot)/(1.+dt/(tc*hbar/usbot))
      !> calculate ud: disc- and time-averaged velocity
      !  usbot = 1.0_wp / (2.5_wp * log(hbar/z0))
      do ibi = 1, NumberOfObjects
        ibm_ACL(ibi)%U_ref = ibm_ACL(ibi)%U_ref*(1.0-fexp) &
          + ibm_ACL(ibi)%U1_ref * fexp
      enddo
    else
      print *, "URefDynamicAverage mode unimplemented. Please use 0 or 1."
    endif

  end subroutine Uref_ACL

  subroutine Calc_turbineangvel(dt_, ibm, fsi)
  implicit none
    real(wp) :: dt_
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi

    integer :: ibi
    real(wp) :: R, TSR
    !< plyunote: please check angvel_fixed and Tipspeedratio initialization
    
    !< plyunote: is angvel0 necessary? seems not. Initial angular velocity
    !            can be get from following process.
    
    do ibi = 1, NumberOfTurbines
      R = fsi(ibi)%r_rotor
      if(FixTipSpeedRatio .eq. 1) then
        TSR = ibm(ibi)%Tipspeedratio
        fsi(ibi)%angvel_axis = TSR * ibm(ibi)%U_ref / R
      else if (FixTurbineAngVel .eq. 1) then
        fsi(ibi)%angvel_axis = fsi(ibi)%angvel_fixed
        ibm(ibi)%Tipspeedratio = fsi(ibi)%angvel_axis * R / ibm(ibi)%U_ref
      endif
      ibm(ibi)%TSR_modified = fsi(ibi)%angvel_axis * R / &
        sqrt((ibm(ibi)%U_ref)**2+(fsi(ibi)%angvel_axis*R)**2)
      if(myid.eq.0) then
        print*,'Turbine_',ibi,':angvel=',fsi(ibi)%angvel_axis,&
          ', TSR=',ibm(ibi)%Tipspeedratio,', U1ref=',ibm(ibi)%U1_ref, &
          ', Uref=', ibm(ibi)%U_ref
      endif
    enddo
  end subroutine Calc_turbineangvel
  
  subroutine Calc_nacelle_angvel(NumberOfObject, fsi, NumNacPerTurb, fsi_ref)
    implicit none
    integer :: NumberOfObject, NumNacPerTurb
    !type(IBMNodes), dimension(:) :: ibm, ibm_ref
    type(FSInfo), dimension(:) :: fsi, fsi_ref

    integer :: ibi, idx_turb
    
    do ibi = 1, NumberOfObject
      idx_turb = (ibi-1) / NumNacPerTurb + 1
      if (fsi(ibi)%rotate_alongaxis .eq. 0) then ! no rotation
        fsi(ibi)%angvel_axis = 0.0
      elseif (fsi(ibi)%rotate_alongaxis .eq. 1) then ! same angular velocity with rotor
        fsi(ibi)%angvel_axis = fsi_ref(idx_turb)%angvel_axis
      elseif (fsi(ibi)%rotate_alongaxis .eq. 2) then ! fixed angular velocity
        fsi(ibi)%angvel_axis = fsi(ibi)%angvel_fixed
      endif
    enddo
  end subroutine Calc_nacelle_angvel

  subroutine ArbitraryRotate(p, theta, r0, q)
  implicit none
    real(wp), dimension(3), intent(in) :: p, r0
    real(wp), intent(in) :: theta
    real(wp), dimension(3), intent(out) :: q

    real(wp), dimension(3) :: r
    real(wp) :: costheta, sintheta, rr
    integer :: nb, n_elmt_1

    rr = norm2(r0)+1.0e-16_wp
    r = r0 / rr

    costheta = cos(theta)
    sintheta = sin(theta)
    !print *,'plyudebug,test1.2,',theta,costheta,sintheta
    !print*,'plyudebug,test1.3,',r(1),r(2),r(3), rr
    q = 0.0_wp
    q(1) = (costheta+(1.0_wp-costheta)*r(1)**2)*p(1) + &
      ((1.0_wp-costheta)*r(1)*r(2)-r(3)*sintheta)*p(2) + &
      ((1.0_wp-costheta)*r(1)*r(3)+r(2)*sintheta)*p(3)

    q(2) = ((1.0_wp-costheta)*r(1)*r(2)+r(3)*sintheta)*p(1)+&
      (costheta+(1.0_wp-costheta)*r(2)**2)*p(2)+&
      ((1.0_wp-costheta)*r(2)*r(3)-r(1)*sintheta)*p(3)

    q(3) = ((1.0_wp-costheta)*r(1)*r(3)-r(2)*sintheta)*p(1)+&
      ((1.0_wp-costheta)*r(2)*r(3)+r(1)*sintheta)*p(2)+&
      (costheta+(1.0_wp-costheta)*r(3)*r(3))*p(3)
  end subroutine ArbitraryRotate
  
  subroutine Test_ArbitraryRotate
    real(wp), dimension(3) :: p,q,r
    real(wp) :: theta

    p(1) = 1.0; p(2) = 0.0; p(3) = 0.0
    r(1) = 0.0; r(2) = 1.0; r(3) = 0.0
    q = 0.0
    
    theta = twopi/4.0
    call ArbitraryRotate(p,theta,r,q)
    print*,'plyudebug,test1.1,',q(1),q(2),q(3)

    p = q
    theta = twopi/4.0
    call ArbitraryRotate(p,theta,r,q)
    print*,'plyudebug,test1.1,',q(1),q(2),q(3)
    
    p = q
    theta = twopi/4.0
    call ArbitraryRotate(p,theta,r,q)
    print*,'plyudebug,test1.1,',q(1),q(2),q(3)

    p(1) = 1.0; p(2) = 0.0; p(3) = 0.0
    theta = twopi*3.0/4.0
    call ArbitraryRotate(p,theta,r,q)
    print*,'plyudebug,test1.1,',q(1),q(2),q(3)

  end subroutine Test_ArbitraryRotate

  subroutine Export_LineLocation (ibm, NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    integer :: NumberOfObjects

    character(128) :: filen
    integer :: ibi, i
  
    if (myid .eq. 0) then
      do ibi = 1, NumberOfObjects
        filen=''
        !> foil number starts from 1 in program
        !> while it starts from FOIL00 in input filelists
        write(filen, '(a5,i0.6,a1,i0.3,a7)') 'line_', ti,'_',ibi,'_nf.dat'
        open(18001, file=trim(filen), action='write', status='replace')
        
        write(18001,*) 'TITLE = "Actuator line mesh"'
        write(18001,*) 'VARIABLES = "X", "Y", "Z", "ub_x", "ub_y", "ub_z", &
          &"color", "U_lagr_x", "U_lagr_y", "U_lagr_z", "F_lagr_x", &
          &"F_lagr_y", "F_lagr_z"'
        write(18001,*) 'ZONE T="Points", DATAPACKING=BLOCK, NODES=',&
          ibm(ibi)%n_v,', ELEMENTS=', ibm(ibi)%n_elmt, ', ZONETYPE=FELINESEG&
          &, VARLOCATION=([1-6]=NODAL,[7-13]=CELLCENTERED)'
        write(18001,*) 'STRANDID=0.1 SOLUTIONTIME=', ti*dt
        
        !do i = 1,ibm(ibi)%n_elmt
        !  write(18001,*) ibm(ibi)%cent_x(i), ibm(ibi)%cent_y(i), ibm(ibi)%cent_z(i)
        !enddo

        do i = 1, ibm(ibi)%n_v
          write(18001,*) ibm(ibi)%x_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18001,*) ibm(ibi)%y_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18001,*) ibm(ibi)%z_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18001,*) ibm(ibi)%u(i,1)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18001,*) ibm(ibi)%u(i,2)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18001,*) ibm(ibi)%u(i,3)
        enddo
        
        do i = 1, ibm(ibi)%n_elmt
          write(18001,*) ibm(ibi)%color(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18001,*) ibm(ibi)%U_lagr_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18001,*) ibm(ibi)%U_lagr_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18001,*) ibm(ibi)%U_lagr_z(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18001,*) ibm(ibi)%F_lagr_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18001,*) ibm(ibi)%F_lagr_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18001,*) ibm(ibi)%F_lagr_z(i)
        enddo

        do i = 1, ibm(ibi)%n_elmt
          write(18001,*) ibm(ibi)%nv1(i), ibm(ibi)%nv2(i)
        enddo

        close(18001)
      enddo
    endif
  end subroutine Export_LineLocation
  
  subroutine Read_LineLocation (ibm, NumberOfObjects, ti_)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    integer :: NumberOfObjects

    character(128) :: filen
    integer :: ti_, ibi, i
  
    if (myid .eq. 0) then
      do ibi = 1, NumberOfObjects
        filen=''
        !> foil number starts from 1 in program
        !> while it starts from FOIL00 in input filelists
        write(filen, '(a5,i0.6,a1,i0.3,a7)') 'line_', ti_,'_',ibi,'_nf.dat'
        open(18001, file=trim(filen), action='read', status='old')
        
        read(18001,*) !'TITLE = "Actuator line mesh"'
        read(18001,*) !'VARIABLES = "X", "Y", "Z", "ub_x", "ub_y", "ub_z", &
         ! &"color", "U_lagr_x", "U_lagr_y", "U_lagr_z", "F_lagr_x", &
         ! &"F_lagr_y", "F_lagr_z"'
        read(18001,*) !'ZONE T="Points", DATAPACKING=BLOCK, NODES=',&
         ! ibm(ibi)%n_v,', ELEMENTS=', ibm(ibi)%n_elmt, ', ZONETYPE=FELINESEG&
         ! &, VARLOCATION=([1-6]=NODAL,[7-13]=CELLCENTERED)'
        read(18001,*) !'STRANDID=0.1 SOLUTIONTIME=', ti*dt
        
        !do i = 1,ibm(ibi)%n_elmt
        !  write(18001,*) ibm(ibi)%cent_x(i), ibm(ibi)%cent_y(i), ibm(ibi)%cent_z(i)
        !enddo

        do i = 1, ibm(ibi)%n_v
          read(18001,*) ibm(ibi)%x_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18001,*) ibm(ibi)%y_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18001,*) ibm(ibi)%z_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18001,*) ibm(ibi)%u(i,1)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18001,*) ibm(ibi)%u(i,2)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18001,*) ibm(ibi)%u(i,3)
        enddo
        
        do i = 1, ibm(ibi)%n_elmt
          read(18001,*) ibm(ibi)%color(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18001,*) ibm(ibi)%U_lagr_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18001,*) ibm(ibi)%U_lagr_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18001,*) ibm(ibi)%U_lagr_z(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18001,*) ibm(ibi)%F_lagr_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18001,*) ibm(ibi)%F_lagr_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18001,*) ibm(ibi)%F_lagr_z(i)
        enddo

        do i = 1, ibm(ibi)%n_elmt
          read(18001,*) ibm(ibi)%nv1(i), ibm(ibi)%nv2(i)
        enddo

        close(18001)
      enddo
    endif

    do ibi = 1, NumberOfObjects
      call mpi_bcast(ibm(ibi)%x_bp, ibm(ibi)%n_v, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%y_bp, ibm(ibi)%n_v, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%z_bp, ibm(ibi)%n_v, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%F_lagr_x, ibm(ibi)%n_elmt, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%F_lagr_y, ibm(ibi)%n_elmt, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%F_lagr_z, ibm(ibi)%n_elmt, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)

      !> plyutodo: some variables are not sent to other cpus.
    enddo

  end subroutine Read_LineLocation
  
  subroutine Export_SurfaceLocation (ibm, NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    integer :: NumberOfObjects

    character(128) :: filen
    integer :: ibi, i
  
    if (myid .eq. 0) then
      do ibi = 1, NumberOfObjects
        filen=''
        !> foil number starts from 1 in program
        !> while it starts from FOIL00 in input filelists
        write(filen, '(a8,i0.6,a1,i0.3,a7)') 'surface_', ti,'_',ibi,'_nf.dat'
        open(18005, file=trim(filen), action='write', status='replace')
        
        write(18005,*) 'TITLE = "Actuator surface mesh"'
        write(18005,*) 'VARIABLES = "X", "Y", "Z", "ub_x", "ub_y", "ub_z", &
          &"color", "F_lagr_x", "F_lagr_y", "F_lagr_z"'
        write(18005,*) 'ZONE T="Triangles", NODES=', ibm(ibi)%n_v,', ELEMENTS=',&
          ibm(ibi)%n_elmt, ', ZONETYPE=FETRIANGLE, &
          &VARLOCATION=([1-6]=NODAL,[7-10]=CELLCENTERED)'
        write(18005,*) 'STRANDID=0.1 SOLUTIONTIME=', ti*dt
        
        !do i = 1,ibm(ibi)%n_elmt
        !  write(18001,*) ibm(ibi)%cent_x(i), ibm(ibi)%cent_y(i), ibm(ibi)%cent_z(i)
        !enddo

        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%x_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%y_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%z_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%u(i,1)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%u(i,2)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%u(i,3)
        enddo
        
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%color(i)
        enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  write(18005,*) ibm(ibi)%U_lagr_x(i)
        !enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  write(18005,*) ibm(ibi)%U_lagr_y(i)
        !enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  write(18005,*) ibm(ibi)%U_lagr_z(i)
        !enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%F_lagr_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%F_lagr_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%F_lagr_z(i)
        enddo

        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%nv1(i), ibm(ibi)%nv2(i), ibm(ibi)%nv3(i)
        enddo

        close(18005)
      enddo
    endif
    call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
  end subroutine Export_SurfaceLocation
  
  subroutine Read_SurfaceLocation (ibm, NumberOfObjects, ti_)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    integer :: NumberOfObjects, ti_

    character(128) :: filen
    integer :: ibi, i
  
    if (myid .eq. 0) then
      do ibi = 1, NumberOfObjects
        filen=''
        !> foil number starts from 1 in program
        !> while it starts from FOIL00 in input filelists
        write(filen, '(a8,i0.6,a1,i0.3,a7)') 'surface_', ti_,'_',ibi,'_nf.dat'
        open(18005, file=trim(filen), action='read', status='old')
        
        read(18005,*) !'TITLE = "Actuator surface mesh"'
        read(18005,*) !'VARIABLES = "X", "Y", "Z", "ub_x", "ub_y", "ub_z", &
        !  &"color", "F_lagr_x", "F_lagr_y", "F_lagr_z"'
        read(18005,*) !'ZONE T="Triangles", NODES=', ibm(ibi)%n_v,', ELEMENTS=',&
        !  ibm(ibi)%n_elmt, ', ZONETYPE=FETRIANGLE, &
        !  &VARLOCATION=([1-6]=NODAL,[7-10]=CELLCENTERED)'
        read(18005,*) !'STRANDID=0.1 SOLUTIONTIME=', ti*dt

        !do i = 1,ibm(ibi)%n_elmt
        !  write(18001,*) ibm(ibi)%cent_x(i), ibm(ibi)%cent_y(i), ibm(ibi)%cent_z(i)
        !enddo

        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%x_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%y_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%z_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%u(i,1)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%u(i,2)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%u(i,3)
        enddo
        
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%color(i)
        enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  read(18005,*) ibm(ibi)%U_lagr_x(i)
        !enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  read(18005,*) ibm(ibi)%U_lagr_y(i)
        !enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  read(18005,*) ibm(ibi)%U_lagr_z(i)
        !enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%F_lagr_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%F_lagr_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%F_lagr_z(i)
        enddo

        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%nv1(i), ibm(ibi)%nv2(i), ibm(ibi)%nv3(i)
        enddo

        close(18005)
      enddo
    endif
    call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
    
    do ibi = 1, NumberOfObjects
      call mpi_bcast(ibm(ibi)%x_bp, ibm(ibi)%n_v, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%y_bp, ibm(ibi)%n_v, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%z_bp, ibm(ibi)%n_v, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%F_lagr_x, ibm(ibi)%n_elmt, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%F_lagr_y, ibm(ibi)%n_elmt, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%F_lagr_z, ibm(ibi)%n_elmt, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
    enddo
  end subroutine Read_SurfaceLocation

  subroutine Export_NacelleLocation (ibm, NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    integer :: NumberOfObjects

    character(128) :: filen
    integer :: ibi, i
  
    if (myid .eq. 0) then
      do ibi = 1, NumberOfObjects
        filen=''
        !> foil number starts from 1 in program
        !> while it starts from FOIL00 in input filelists
        write(filen, '(a8,i0.6,a1,i0.3,a7)') 'surfnac_', ti,'_',ibi,'_nf.dat'
        open(18005, file=trim(filen), action='write', status='replace')
        
        write(18005,*) 'TITLE = "Nacelle surface mesh"'
        write(18005,*) 'VARIABLES = "X", "Y", "Z", "ub_x", "ub_y", "ub_z", &
          &"color", "F_lagr_x", "F_lagr_y", "F_lagr_z", "nf_x", "nf_y", "nf_z"'
        write(18005,*) 'ZONE T="Triangles", NODES=', ibm(ibi)%n_v,', ELEMENTS=',&
          ibm(ibi)%n_elmt, ', ZONETYPE=FETRIANGLE, &
          &VARLOCATION=([1-6]=NODAL,[7-13]=CELLCENTERED)'
        write(18005,*) 'STRANDID=0.1 SOLUTIONTIME=', ti*dt
        
        !do i = 1,ibm(ibi)%n_elmt
        !  write(18001,*) ibm(ibi)%cent_x(i), ibm(ibi)%cent_y(i), ibm(ibi)%cent_z(i)
        !enddo

        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%x_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%y_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%z_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%u(i,1)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%u(i,2)
        enddo
        do i = 1, ibm(ibi)%n_v
          write(18005,*) ibm(ibi)%u(i,3)
        enddo
        
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%color(i)
        enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  write(18005,*) ibm(ibi)%U_lagr_x(i)
        !enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  write(18005,*) ibm(ibi)%U_lagr_y(i)
        !enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  write(18005,*) ibm(ibi)%U_lagr_z(i)
        !enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%F_lagr_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%F_lagr_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%F_lagr_z(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%nf_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%nf_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%nf_z(i)
        enddo

        do i = 1, ibm(ibi)%n_elmt
          write(18005,*) ibm(ibi)%nv1(i), ibm(ibi)%nv2(i), ibm(ibi)%nv3(i)
        enddo

        close(18005)
      enddo
    endif
  end subroutine Export_NacelleLocation

  subroutine Read_NacelleLocation (ibm, NumberOfObjects, ti_)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    integer :: NumberOfObjects, ti_

    character(128) :: filen
    integer :: ibi, i
  
    if (myid .eq. 0) then
      do ibi = 1, NumberOfObjects
        filen=''
        !> foil number starts from 1 in program
        !> while it starts from FOIL00 in input filelists
        write(filen, '(a8,i0.6,a1,i0.3,a7)') 'surfnac_', ti_,'_',ibi,'_nf.dat'
        open(18005, file=trim(filen), action='read', status='old')
        
        read(18005,*) !'TITLE = "Nacelle surface mesh"'
        read(18005,*) !'VARIABLES = "X", "Y", "Z", "ub_x", "ub_y", "ub_z", &
        !  &"color", "F_lagr_x", "F_lagr_y", "F_lagr_z", "nf_x", "nf_y", "nf_z"'
        read(18005,*) !'ZONE T="Triangles", NODES=', ibm(ibi)%n_v,', ELEMENTS=',&
        !  ibm(ibi)%n_elmt, ', ZONETYPE=FETRIANGLE, &
        !  &VARLOCATION=([1-6]=NODAL,[7-13]=CELLCENTERED)'
        read(18005,*) !'STRANDID=0.1 SOLUTIONTIME=', ti*dt
        
        !do i = 1,ibm(ibi)%n_elmt
        !  write(18001,*) ibm(ibi)%cent_x(i), ibm(ibi)%cent_y(i), ibm(ibi)%cent_z(i)
        !enddo

        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%x_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%y_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%z_bp(i)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%u(i,1)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%u(i,2)
        enddo
        do i = 1, ibm(ibi)%n_v
          read(18005,*) ibm(ibi)%u(i,3)
        enddo
        
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%color(i)
        enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  write(18005,*) ibm(ibi)%U_lagr_x(i)
        !enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  write(18005,*) ibm(ibi)%U_lagr_y(i)
        !enddo
        !do i = 1, ibm(ibi)%n_elmt
        !  write(18005,*) ibm(ibi)%U_lagr_z(i)
        !enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%F_lagr_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%F_lagr_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%F_lagr_z(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%nf_x(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%nf_y(i)
        enddo
        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%nf_z(i)
        enddo

        do i = 1, ibm(ibi)%n_elmt
          read(18005,*) ibm(ibi)%nv1(i), ibm(ibi)%nv2(i), ibm(ibi)%nv3(i)
        enddo

        close(18005)
      enddo
    endif
    
    do ibi = 1, NumberOfObjects
      call mpi_bcast(ibm(ibi)%x_bp, ibm(ibi)%n_v, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%y_bp, ibm(ibi)%n_v, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%z_bp, ibm(ibi)%n_v, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%F_lagr_x, ibm(ibi)%n_elmt, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%F_lagr_y, ibm(ibi)%n_elmt, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
      call mpi_bcast(ibm(ibi)%F_lagr_z, ibm(ibi)%n_elmt, mpi_double_precision, &
        0, mpi_comm_2d_cart, ierr_wt)
    enddo
  end subroutine Read_NacelleLocation
  
  !> exporting the forces on the turbines, actuator line model
  subroutine Calc_forces_ACL(ibm, fsi, NumberOfObjects, fname)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi
    integer :: NumberOfObjects
    real(wp) :: time_
    character(len=*) :: fname

    integer :: ibi, l
    real(wp) :: A_Sum, P_Sum, U_Sum, F_Sum_amp, M_Sum
    real(wp) :: r_amp, F_axis, U_axis
    real(wp), dimension(3) :: F_Sum_vec, M_vec, r_vec, n_vec, temp_vec
    character(128) :: filen
    real(wp) :: C_Thrust, Power_aero, C_Power, ratio_CT, ratio_CP

    !> variables for remove_lateral_ACL_force
    real(wp) :: F_sum_y, F_mean_y

    do ibi = 1, NumberOfObjects
      A_Sum=0.0_wp; P_Sum=0.0_wp; U_Sum=0.0_wp; F_Sum_amp=0.0_wp; M_Sum=0.0_wp
      F_Sum_vec=0.0_wp; M_vec=0.0_wp; F_sum_y=0.0_wp
      n_vec(1:3) = (/ fsi(ibi)%nx_tb, fsi(ibi)%ny_tb, fsi(ibi)%nz_tb /)
      do l = 1, ibm(ibi)%n_elmt
        A_Sum = A_Sum + ibm(ibi)%dA(l)
        F_axis = ibm(ibi)%F_lagr_x(l)*n_vec(1)+ibm(ibi)%F_lagr_y(l) &
          *n_vec(2)+ibm(ibi)%F_lagr_z(l)*n_vec(3)

        !> plyunote: ibmsurface%U_lagr_x is not assigned.
        !> therefore U_axis, P_Sum, U_Sum are meaningless for ACS
        U_axis = ibm(ibi)%U_lagr_x(l)*n_vec(1)+ibm(ibi)%U_lagr_y(l) &
          *n_vec(2)+ibm(ibi)%U_lagr_z(l)*n_vec(3)
        !> plyunote: F = F_lagr_x * dA
        F_Sum_amp = F_Sum_amp + F_axis * ibm(ibi)%dA(l)
        P_Sum = P_Sum + F_axis * ibm(ibi)%dA(l) * U_axis
        U_Sum = U_Sum + U_axis * ibm(ibi)%dA(l)
        F_Sum_vec(1) = F_Sum_vec(1) + ibm(ibi)%F_lagr_x(l) * ibm(ibi)%dA(l)
        F_Sum_vec(2) = F_Sum_vec(2) + ibm(ibi)%F_lagr_y(l) * ibm(ibi)%dA(l)
        F_Sum_vec(3) = F_Sum_vec(3) + ibm(ibi)%F_lagr_z(l) * ibm(ibi)%dA(l)
        F_sum_y = F_sum_y + ibm(ibi)%F_lagr_y(l) * ibm(ibi)%dA(l)
        r_vec(1) = ibm(ibi)%cent_x(l) - fsi(ibi)%x_c
        r_vec(2) = ibm(ibi)%cent_y(l) - fsi(ibi)%y_c
        r_vec(3) = ibm(ibi)%cent_z(l) - fsi(ibi)%z_c
        temp_vec(1:3) =(/ibm(ibi)%F_lagr_x(l),ibm(ibi)%F_lagr_y(l),ibm(ibi)%F_lagr_z(l)/)
        temp_vec(1:3) = temp_vec(1:3) * ibm(ibi)%dA(l)
        call crossx(r_vec, temp_vec, M_vec)
        M_Sum = M_Sum + M_vec(1)*n_vec(1)+M_vec(2)*n_vec(2)+M_vec(3)*n_vec(3)

        !> debug:
        !if (norm2(r_vec)>0.29 .and. norm2(r_vec)<0.30 .and. ibi.eq.1  &
        !  .and. myid.eq.0) then
        !  print *,'plyudebug info, force_ACL: rotor_model=',fsi(ibi)%rotor_model_type
        !  print *, '  cent coord:', ibm(ibi)%cent_x(l), ibm(ibi)%cent_y(l), ibm(ibi)%cent_z(l)
        !  print *, '  r_vec(1:3), dA:', r_vec(1:3), ibm(ibi)%dA(l)
        !  print *, '  F_lagr_xyz:', ibm(ibi)%F_lagr_x(l), ibm(ibi)%F_lagr_y(l),&
        !    ibm(ibi)%F_lagr_z(l)
        !  print *, '  F_axis, U_axis, M_vec(1:3):', F_axis, U_axis, M_vec(1:3)
        !endif
      enddo

      U_Sum = U_Sum / A_Sum

      fsi(ibi)%Torque_aero = -M_Sum
      fsi(ibi)%Force_axis = -F_Sum_amp
      C_Thrust = dabs(fsi(ibi)%Force_axis)/(0.5*ibm(ibi)%U_ref**2*twopi/2.0&
        *fsi(ibi)%r_rotor**2)
      Power_aero = dabs(fsi(ibi)%Torque_aero * fsi(ibi)%angvel_axis)
      C_Power = Power_aero / (0.25*twopi*(fsi(ibi)%r_rotor)**2*ibm(ibi)%U_ref**3)
      
      if(myid .eq. 0) then
        !if (ibi .eq. 1) then
          print *, 'Thrust=', fsi(ibi)%Force_axis, ', Torque=', &
            fsi(ibi)%Torque_aero, ', C_Thrust=', C_Thrust, &
            ', Power=', Power_aero,', C_Power=',C_Power
        !endif
        
        filen=''
        !> foil number starts from 1 in program
        !> while it starts from FOIL00 in input filelists
        write(filen, '(a11,i0.3,a4)') fname, ibi,'.dat'
        
        if(ti .eq. (ti_first+1)) then
          open(18003, file=trim(filen), action='write')
          write(18003,*) 'Variables="time","angle","angvel_axis","Force_axis",&
            &"Torque_fluid","Uref","Ud","TSR","Force_x","Force_y","Force_z"'
        else
          open(18003, file=trim(filen), action='write', position='append')
        endif

        write(18003,"(11(e14.5))") ti*dt, fsi(ibi)%ang_axis, fsi(ibi)%angvel_axis,&
          fsi(ibi)%Force_axis, fsi(ibi)%Torque_aero, ibm(ibi)%U_ref, &
          U_Sum, ibm(ibi)%Tipspeedratio, F_Sum_vec(1),F_Sum_vec(2),F_Sum_vec(3)
              
        close(18003)
      endif

      if(adjust_CT_ACL .eq. 1) then
        ratio_CT = fsi(ibi)%CT0 / (C_Thrust+1.0e-10)
        if(myid.eq.0) then
          print *, "Calc_forces_ACL: scale CT from ", C_Thrust,& 
            " to ", fsi(ibi)%CT0, " , ratio_CT_ACL=", ratio_CT
        endif
        
        do l = 1, ibm(ibi)%n_elmt
          ibm(ibi)%F_lagr_x(l) = ibm(ibi)%F_lagr_x(l) * ratio_CT
          ibm(ibi)%F_lagr_y(l) = ibm(ibi)%F_lagr_y(l) * ratio_CT
          ibm(ibi)%F_lagr_z(l) = ibm(ibi)%F_lagr_z(l) * ratio_CT
        enddo
      endif

      !> the total force of the three blades might have a spanwise component,
      !! which may lead to a lateral mean flow in wind farm case after long-time
      !! development. 
      !> therefore, for wind farm case, I suggest to remove the mean lateral
      !! force.
      !> the removing is applied to blade force only. The nacelle force is not
      !! affected.
      !> But, the implementation is limited to the situation that rotation axis
      !! parallel to x direction.
      if (remove_lateral_ACL_force .eq. 1) then
        F_mean_y = F_sum_y / ibm(ibi)%n_elmt
        do l = 1, ibm(ibi)%n_elmt
          ibm(ibi)%F_lagr_y(l) = ibm(ibi)%F_lagr_y(l) - F_mean_y/ibm(ibi)%dA(l)
        enddo
        if (myid.eq.0) then
          print *, 'plyudebug: remove_lateral_ACL_force, turbine i=', &
            ibi, ', mean=',F_mean_y
        endif
      endif

    enddo
  end subroutine Calc_forces_ACL
  
  !> exporting the forces on the nacelle model
  subroutine Calc_forces_nacelle(ibm, fsi, NumberOfObjects, fname)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi
    integer :: NumberOfObjects
    real(wp) :: time_
    character(len=*) :: fname

    integer :: ibi, l
    real(wp) :: A_Sum, P_Sum, U_Sum, F_Sum_amp, M_Sum
    real(wp) :: r_amp, F_axis, U_axis
    real(wp), dimension(3) :: F_Sum_vec, M_vec, r_vec, n_vec, temp_vec
    character(128) :: filen
    real(wp) :: C_Thrust, Power_aero, C_Power, ratio_CT, ratio_CP

    do ibi = 1, NumberOfObjects
      A_Sum=0.0; P_Sum=0.0; U_Sum=0.0; F_Sum_amp=0.0; M_Sum=0.0
      F_Sum_vec=0.0; M_vec=0.0
      n_vec(1:3) = (/ fsi(ibi)%nx_tb, fsi(ibi)%ny_tb, fsi(ibi)%nz_tb /)
      do l = 1, ibm(ibi)%n_elmt
        A_Sum = A_Sum + ibm(ibi)%dA(l)
        F_axis = ibm(ibi)%F_lagr_x(l)*n_vec(1)+ibm(ibi)%F_lagr_y(l) &
          *n_vec(2)+ibm(ibi)%F_lagr_z(l)*n_vec(3)

        !> plyunote: ibmsurface%U_lagr_x is not assigned.
        !> therefore U_axis, P_Sum, U_Sum are meaningless for ACS
        U_axis = ibm(ibi)%U_lagr_x(l)*n_vec(1)+ibm(ibi)%U_lagr_y(l) &
          *n_vec(2)+ibm(ibi)%U_lagr_z(l)*n_vec(3)
        !> plyunote: F = F_lagr_x * dA
        F_Sum_amp = F_Sum_amp + F_axis * ibm(ibi)%dA(l)
        !P_Sum = P_Sum + F_axis * ibm(ibi)%dA(l) * U_axis
        U_Sum = U_Sum + U_axis * ibm(ibi)%dA(l)
        F_Sum_vec(1) = F_Sum_vec(1) + ibm(ibi)%F_lagr_x(l) * ibm(ibi)%dA(l)
        F_Sum_vec(2) = F_Sum_vec(2) + ibm(ibi)%F_lagr_y(l) * ibm(ibi)%dA(l)
        F_Sum_vec(3) = F_Sum_vec(3) + ibm(ibi)%F_lagr_z(l) * ibm(ibi)%dA(l)
        !r_vec(1) = ibm(ibi)%cent_x(l) - fsi(ibi)%x_c
        !r_vec(2) = ibm(ibi)%cent_y(l) - fsi(ibi)%y_c
        !r_vec(3) = ibm(ibi)%cent_z(l) - fsi(ibi)%z_c
        !temp_vec(1:3) =(/ibm(ibi)%F_lagr_x(l),ibm(ibi)%F_lagr_y(l),ibm(ibi)%F_lagr_z(l)/)
        !temp_vec(1:3) = temp_vec(1:3) * ibm(ibi)%dA(l)
        !call crossx(r_vec, temp_vec, M_vec)
        !M_Sum = M_Sum + M_vec(1)*n_vec(1)+M_vec(2)*n_vec(2)+M_vec(3)*n_vec(3)

       ! !> debug:
       ! if (norm2(r_vec)>0.0 .and. norm2(r_vec)<20.0 .and. ibi.eq.1  &
       !   .and. myid.eq.0) then
       !   print *,'plyudebug info, force_ACL: rotor_model=',fsi(ibi)%rotor_model_type
       !   print *, '  cent coord:', ibm(ibi)%cent_x(l), ibm(ibi)%cent_y(l), ibm(ibi)%cent_z(l)
       !   print *, '  r_vec(1:3), dA:', r_vec(1:3), ibm(ibi)%dA(l)
       !   print *, '  F_lagr_xyz:', ibm(ibi)%F_lagr_x(l), ibm(ibi)%F_lagr_y(l),&
       !     ibm(ibi)%F_lagr_z(l)
       !   print *, '  F_axis, U_axis, M_vec(1:3):', F_axis, U_axis, M_vec(1:3)
       ! endif
      enddo

      U_Sum = U_Sum / A_Sum

      !fsi(ibi)%Torque_aero = -M_Sum
      fsi(ibi)%Force_axis = -F_Sum_amp
      C_Thrust = dabs(fsi(ibi)%Force_axis)/(0.5*ibm(ibi)%U_ref**2 * fsi(ibi)%area_CT)
      !Power_aero = dabs(fsi(ibi)%Torque_aero * fsi(ibi)%angvel_axis)
      !C_Power = Power_aero / (0.25*twopi*(fsi(ibi)%r_rotor)**2*ibm(ibi)%U_ref**3)

      ratio_CT = -1.0
      if(adjust_CT_nacelle .ne. 0) then
        if(adjust_CT_nacelle .eq. 1) then
          ratio_CT = fsi(ibi)%CT0 / (C_Thrust+1.0e-10)
        endif
        if(adjust_CT_nacelle .eq. 2) then
          ratio_CT = fsi(ibi)%ratio_CT0
        endif
        !if(myid.eq.0) then
        !  print *, "Calc_forces_nacelle: scale CD from ", C_Thrust,& 
        !    " to ", fsi(ibi)%CT0, " , ratio_CD_nacelle=", ratio_CT
        !endif

        do l = 1, ibm(ibi)%n_elmt
          ibm(ibi)%F_lagr_x(l) = ibm(ibi)%F_lagr_x(l) * ratio_CT
          ibm(ibi)%F_lagr_y(l) = ibm(ibi)%F_lagr_y(l) * ratio_CT
          ibm(ibi)%F_lagr_z(l) = ibm(ibi)%F_lagr_z(l) * ratio_CT
        enddo
      endif
      
      if(myid .eq. 0) then
        !if (ibi .eq. 1) then
        !  print *, 'nacelle F_Drag=', fsi(ibi)%Force_axis, ', C_Drag=', C_Thrust
        !endif
        
        filen=''
        !> foil number starts from 1 in program
        !> while it starts from FOIL00 in input filelists
        write(filen, '(a,i0.3,a4)') fname, ibi,'.dat'
        
        if(ti .eq. (ti_first+1)) then
          open(18003, file=trim(filen), action='write')
          write(18003,*) 'Variables="time","Force_axis","Uref","Ud","Area","Force_x","Force_y","Force_z","RatioCT"'
        else
          open(18003, file=trim(filen), action='write', position='append')
        endif

        write(18003,"(9(e14.5))") ti*dt, fsi(ibi)%Force_axis, ibm(ibi)%U_ref, &
          U_Sum, A_Sum, F_Sum_vec(1),F_Sum_vec(2),F_Sum_vec(3),ratio_CT
              
        close(18003)
      endif
    enddo
  end subroutine Calc_forces_nacelle

  subroutine rotor_Rot(time_, fsi, ibm, dt_, ibi, dimension_)
    implicit none
    real(wp), intent(in) :: time_
    type(FSInfo) :: fsi
    type(IBMNodes) :: ibm
    !type(FSIPrePara) :: para
    real(wp), intent(in) :: dt_
    integer, intent(in) :: ibi, dimension_

    integer :: n_v, n_elmt, i, istart, j, n1e, n2e, n3e
    integer :: nb, n_elmt_1
    real(wp), dimension(3) :: p, q, nr, na, nt, tmparray, ds12, ds13
    real(wp) :: ang_axis, theta, rx, ry, rz, rr, Ut
    real(wp) :: dr, dx12, dy12, dz12

    character(128) :: filen
    integer :: idof
    real(wp) :: rot_ang
    real(wp), dimension(3) :: rot_in, rot_out, rot_axis

    type(FSIPrePara) :: para

    if (fsitype .eq. 1) then
      para = pp_wt(ibi)
    endif

    !< plyunote: parameter rstart_turbinerotation is ignored
    n_v = ibm%n_v
    !print *,'plyudebug, rotor_Rot:n_v=',n_v,',size(x_bp0)=',&
    !  size(ibm%x_bp0)
    
    !< several cases for ti and ti_first:
    !! 1. run for first time, ti=0, ti_first=0
    !! 2. restart from save of ti=0, ti_first=0, ti=0 when read in, ti=1 when run
    !! 3. restart from save of ti=N, same
    if (ti .eq. (ti_first+1)) then
      !call print_ibmnode1(ibi, 1, 1)
      na(1) = fsi%nx_tb
      na(2) = fsi%ny_tb
      na(3) = fsi%nz_tb
      theta = fsi%ang0_rot
      
     ! do i = 1, n_v
     !   p(1) = ibm%x_bp0(i) - fsi%x_c0
     !   p(2) = ibm%y_bp0(i) - fsi%y_c0
     !   p(3) = ibm%z_bp0(i) - fsi%z_c0
     !   
     !  ! if(fsitype .eq. 1) then
     !  !   rot_in(1:3) = p(1:3)
     !  !   rot_out(1:3) = rot_in(1:3)
     !  !   do idof = 6, 4, -1
     !  !     if(para%dof(idof).ne.0) then
     !  !       !> here rotation angle dd is got from para. Maybe better from fsi.
     !  !       rot_ang = para%disp(idof)
     !  !       rot_axis = 0.0; rot_axis(idof-3) = 1.0
     !  !       call ArbitraryRotate(rot_in, rot_ang, rot_axis, rot_out)
     !  !       rot_in = rot_out
     !  !     endif
     !  !   enddo
     !  !   p(1:3) = rot_out(1:3)
     !  ! endif

     !   !call print_ibmnode1(ibi, i, 2)

     !   !> plyunote: fsi%ang_axis should be added to turbine.inp or elsewhere
     !   !> It is a temporary solution to use fsi%angvel_fixed and ti_first to
     !   !! guess initial rotation angle. For cases of variable angular speed,
     !   !! fsi%ang_axis should be reread from a save
     !   
     !   !fsi%ang_axis = ti_first * dt * fsi%angvel_fixed
     !   !theta = fsi%ang_axis

     !   call ArbitraryRotate(p, theta, na, q)

     !   !> this position considers only rotor rotation. No FSI motion.
     !   ibm%x_bp(i) = q(1) + fsi%x_c0
     !   ibm%y_bp(i) = q(2) + fsi%y_c0
     !   ibm%z_bp(i) = q(3) + fsi%z_c0

     !   !call print_ibmnode1(ibi, i, 3)
     ! enddo
      
      fsi%ang_axis = theta
      if (myid.eq.0) print *, 'plyudebug, initial_ang_axis=', fsi%ang_axis
    endif

    !> plyuinfo: fsi%angvel_axis should be added to turbine.inp or elsewhere
    theta = fsi%angvel_axis * dt_
    fsi%ang_axis = fsi%ang_axis + theta
    if (myid.eq.0) then
      print *, 'Angleinfo: ang_axis, ang0_yaw, angvel_axis', &
        fsi%ang_axis, 0.0_wp, fsi%angvel_axis
    endif

    do i = 1, n_v
      !> plyunote: even though following printing is unnecessary for compuation,
      !! it had to stay there in case of a strange error.
      !! You can try to comment these lines.
      !! I guess we had to call or use some variables before further usages.
      if (ibi .eq. 1 .and. i .eq. 10 .and. myid.eq.0) then
        print *, 'plyudebug, rotor_Rot, point before rot:', &
          ibm%x_bp(i), ibm%y_bp(i), ibm%z_bp(i)
      endif

      !p(1) = ibm%x_bp(i) - fsi%x_co
      p(1) = ibm%x_bp0(i) - fsi%x_c0
      p(2) = ibm%y_bp0(i) - fsi%y_c0
      p(3) = ibm%z_bp0(i) - fsi%z_c0
      na(1) = fsi%nx_tb
      na(2) = fsi%ny_tb
      na(3) = fsi%nz_tb

      if(fsitype .eq. 1) then
        rot_in(1:3) = p(1:3)
        rot_out(1:3) = rot_in(1:3)
        do idof = 6, 4, -1
          if(para%dof(idof).ne.0) then
            !> here rotation angle dd is got from para. Maybe better from fsi.
            rot_ang = para%disp(idof)
            rot_axis = 0.0; rot_axis(idof-3) = 1.0
            call ArbitraryRotate(rot_in, rot_ang, rot_axis, rot_out)
            rot_in = rot_out
          endif
        enddo
        p(1:3) = rot_out(1:3)
      endif

      !call print_ibmnode1(ibi, i, 4)
      
      !call ArbitraryRotate(p, theta, na, q)
      call ArbitraryRotate(p, fsi%ang_axis, na, q)

      !call print_ibmnode1(ibi, i, 5)
      
      !> here translation from fsi is included in fsi%x_c
      ibm%x_bp(i) = q(1) + fsi%x_c
      ibm%y_bp(i) = q(2) + fsi%y_c
      ibm%z_bp(i) = q(3) + fsi%z_c
      
      !> plyunote: even though following printing is unnecessary for compuation,
      !! it had to stay there in case of a strange error.
      !! You can try to comment these lines.
      !! I guess we had to call or use some variables before further usages.
      if (ibi .eq. 1 .and. i .eq. 10 .and. myid.eq.0) then
        print *, 'plyudebug, rotor_Rot, point after rot:', &
          ibm%x_bp(i), ibm%y_bp(i), ibm%z_bp(i)
      endif

      rx=ibm%x_bp(i)-fsi%x_c; ry=ibm%y_bp(i)-fsi%y_c; rz=ibm%z_bp(i)-fsi%z_c
      rr = sqrt(rx**2 + ry**2 + rz**2)+1.0e-19_wp
      nr(1) = rx/rr; nr(2) = ry/rr; nr(3) = rz/rr

      nt(1) = na(2)*nr(3) - na(3)*nr(2)
      nt(2) = na(3)*nr(1) - na(1)*nr(3)
      nt(3) = na(1)*nr(2) - na(2)*nr(1)
      !print*,'plyudebug,rotot_Rot,na,nr,nt'
      !print*,na(1),na(2),na(3)
      !print*,nr(1),nr(2),nr(3)
      !print*,nt(1),nt(2),nt(3)
      
      Ut = fsi%angvel_axis * rr
      
      !> here u is pure rotation velocity. Translational velocity not included.
      ibm%u(i, 1:3) = Ut * nt(1:3)

      !if (fsitype .eq. 0) then
      !  ibm%u(i,1:3) = Ut * nt(1:3)
      !elseif (fsitype .eq. 1) then
      !  ibm%u(i,1:3) = Ut * nt(1:3) + fsi%vel(1:3)
      !else
      !  print *, "This fsitype-rotor_Rot has not been implemented."
      !endif
      
      !print*,'plyudebug,angvel_axis=',fsi%angvel_axis,',rr=',rr
    enddo

    !> plyunote: case dimension==2 is ignored here

    !> for ACL, the nf_x, nf_y, nf_z denote the direction of the actuator line
    !> dA denotes the length of each element
    n_elmt_1 = ibm%n_elmt / num_blade
    if (dimension_ .eq. 1) then 
      do nb = 1, num_blade
      do j = 1, n_elmt_1
        i = (nb-1) * n_elmt_1 + j

        n1e = ibm%nv1(i); n2e = ibm%nv2(i)
        dx12 = ibm%x_bp(n2e) - ibm%x_bp(n1e)
        dy12 = ibm%y_bp(n2e) - ibm%y_bp(n1e)
        dz12 = ibm%z_bp(n2e) - ibm%z_bp(n1e)

        dr = sqrt(dx12**2+dy12**2+dz12**2)
        dx12 = dx12/dr; dy12 = dy12/dr; dz12 = dz12/dr
        ibm%nf_x(i) = dx12; ibm%nf_y(i) = dy12; ibm%nf_z(i) = dz12
        !dr = sqrt(ibm%nf_x(i)**2 + ibm%nf_y(i)**2 + ibm%nf_z(i)**2)
        
        ibm%ns_x(i) = 0.; ibm%ns_y(i) = 0.; ibm%ns_z(i) = 0.
        !> plyunote: why nt_x is set to zero here?
        !! answer: in process of reading point info, ACL set ns and nt to 0,
        !! AD and AS set ns and nt to nonzero. 
        !! So here ns and nt should be handle based on element type.
        ibm%nt_x(i) = 0.; ibm%nt_y(i) = 0.; ibm%nt_z(i) = 0.

        ibm%dA(i) = dr

        !> cent_x of ACL, denotes center point of this two-vertex line element
        ibm%cent_x(i) = (ibm%x_bp(n1e)+ibm%x_bp(n2e))/2.
        ibm%cent_y(i) = (ibm%y_bp(n1e)+ibm%y_bp(n2e))/2.
        ibm%cent_z(i) = (ibm%z_bp(n1e)+ibm%z_bp(n2e))/2.
        
        !if (ibi .eq. 1 .and. i .eq. 1) then
        !  print *, 'plyudebug, rotor_Rot, line elmt after rot:', &
        !    ibm%cent_x(i), ibm%cent_y(i), ibm%cent_z(i)
        !  print *, 'PID and its node index:', myid, n1e, n2e
        !  print *, 'its 1st node:', ibm%x_bp(n1e),ibm%y_bp(n1e),ibm%z_bp(n1e)
        !  print *, 'its 2nd node:', ibm%x_bp(n2e),ibm%y_bp(n2e),ibm%z_bp(n2e)
        !endif
      enddo
      enddo
    endif

    if (dimension_ .eq. 2) then
    do i = 1, ibm%n_elmt
      n1e = ibm%nv1(i); n2e = ibm%nv2(i); n3e = ibm%nv3(i)
      ds12(1) = ibm%x_bp(n2e) - ibm%x_bp(n1e)
      ds12(2) = ibm%y_bp(n2e) - ibm%y_bp(n1e)
      ds12(3) = ibm%z_bp(n2e) - ibm%z_bp(n1e)
      
      ds13(1) = ibm%x_bp(n3e) - ibm%x_bp(n1e)
      ds13(2) = ibm%y_bp(n3e) - ibm%y_bp(n1e)
      ds13(3) = ibm%z_bp(n3e) - ibm%z_bp(n1e)

      call crossx(ds12, ds13, tmparray)
      dr = norm2(tmparray)
      ibm%nf_x(i) = tmparray(1)/dr; ibm%nf_y(i) = tmparray(2)/dr 
      ibm%nf_z(i) = tmparray(3)/dr

      !> plyunote: in original version, streamwise direction is z.
      !! now it is x. so the function is a little different
      !! the projection is z->x, x->y, y->z
      if (  (((1.0-ibm%nf_x(i))<=1.0e-6).and.((-1.0+ibm%nf_x(i))<1.0e-6)) .or.&
        (((1.0+ibm%nf_x(i))<=1.0e-6).and.((-1.0-ibm%nf_x(i))<1.0e-6)) ) then
        ibm%ns_y(i) = 1.0_wp; ibm%ns_z(i) = 0.0_wp; ibm%ns_x(i) = 0.0_wp
        ibm%nt_y(i) = 0.0_wp; ibm%nt_z(i) = 1.0_wp; ibm%nt_x(i) = 0.0_wp
      else 
        rr = sqrt(ibm%nf_y(i)**2+ibm%nf_z(i)**2)
        ibm%ns_y(i) = ibm%nf_y(i) / rr
        ibm%ns_z(i) = -ibm%nf_y(i) / rr
        ibm%ns_x(i) = 0.0_wp
        
        ibm%nt_x(i) = rr 
        ibm%nt_y(i) = -ibm%nf_x(i) * ibm%nf_y(i) / rr
        ibm%nt_z(i) = -ibm%nf_x(i) * ibm%nf_z(i) / rr
      endif

      !> plyunote: for triangular elements, dA is correct?
      !! answer: yes. dA=0.5*|ds12|*|ds13|*sin(A213)=0.5*|ds12 x ds13|
      ibm%dA(i) = dr/2.0_wp

      ibm%cent_x(i) = (ibm%x_bp(n1e)+ibm%x_bp(n2e)+ibm%x_bp(n3e))/3.0_wp
      ibm%cent_y(i) = (ibm%y_bp(n1e)+ibm%y_bp(n2e)+ibm%y_bp(n3e))/3.0_wp
      ibm%cent_z(i) = (ibm%z_bp(n1e)+ibm%z_bp(n2e)+ibm%z_bp(n3e))/3.0_wp
      
      !if (ibi .eq. 1 .and. i .eq. 10 .and. myid.eq.0) then
      !  print *, 'plyudebug, rotor_Rot, surface elmt after rot:', &
      !    ibm%cent_x(i), ibm%cent_y(i), ibm%cent_z(i)
      !endif
    enddo
    endif

   ! if (myid .eq. 0) then
   !   filen=''
   !   !> foil number starts from 1 in program
   !   !> while it starts from FOIL00 in input filelists
   !   write(filen, '(a4,i0.6,a1,i0.3,a7)') 'line', ti,'_',ibi,'_nf.dat'
   !   open(18001, file=trim(filen), action='write', status='replace')
   !   do i = 1,ibm%n_elmt
   !  write(18001,*) ibm%cent_x(i), ibm%cent_y(i), ibm%cent_z(i)
   !   enddo
   !   close(18001)
   ! endif

  end subroutine rotor_Rot

  subroutine nacelle_move(ibm, fsi, para)
    implicit none
    type(FSInfo) :: fsi
    type(IBMNodes) :: ibm
    type(FSIPrepara) :: para

    integer :: i, ibi, idof
    real(wp), dimension(3) :: rot_in, rot_out, rot_axis
    real(wp) :: rot_ang
    
    integer :: n1e, n2e, n3e
    real(wp), dimension(3) :: tmparray, ds12, ds13
    real(wp) :: dr, rr

    !> Restart. See rotor_Rot() for its meaning.
   ! if (ti .eq. (ti_first+1)) then
   !   do i = 1, ibm%n_v
   !     rot_in(1) = ibm%x_bp0(i) - fsi%x_c0
   !     rot_in(2) = ibm%y_bp0(i) - fsi%y_c0
   !     rot_in(3) = ibm%z_bp0(i) - fsi%z_c0
   !     rot_out(1:3) = rot_in(1:3)
   !     do idof = 6, 4, -1
   !       if (para%dof(idof).ne.0) then
   !         rot_ang = para%disp(idof)
   !         rot_axis = 0.0; rot_axis(idof-3) = 1.0
   !         call ArbitraryRotate(rot_in, rot_ang, rot_axis, rot_out)
   !         rot_in = rot_out
   !       endif
   !     enddo
   !     ibm%x_bp(i) = fsi%x_c + rot_out(1)
   !     ibm%y_bp(i) = fsi%y_c + rot_out(2)
   !     ibm%z_bp(i) = fsi%z_c + rot_out(3)        
   !   enddo
   ! endif
    
    do i = 1, ibm%n_v
      !rot_in(1) = ibm%x_bp(i) - fsi%x_co
      rot_in(1) = ibm%x_bp0(i) - fsi%x_c0
      rot_in(2) = ibm%y_bp0(i) - fsi%y_c0
      rot_in(3) = ibm%z_bp0(i) - fsi%z_c0
      rot_out(1:3) = rot_in(1:3)
      do idof = 6, 4, -1
        if (para%dof(idof).ne.0) then
          rot_ang = para%disp(idof)
          rot_axis = 0.0; rot_axis(idof-3) = 1.0
          call ArbitraryRotate(rot_in, rot_ang, rot_axis, rot_out)
          rot_in = rot_out
        endif
      enddo
      ibm%x_bp(i) = fsi%x_c + rot_out(1)
      ibm%y_bp(i) = fsi%y_c + rot_out(2)
      ibm%z_bp(i) = fsi%z_c + rot_out(3)
    enddo

    do i = 1, ibm%n_elmt
      n1e = ibm%nv1(i); n2e = ibm%nv2(i); n3e = ibm%nv3(i)
      ds12(1) = ibm%x_bp(n2e) - ibm%x_bp(n1e)
      ds12(2) = ibm%y_bp(n2e) - ibm%y_bp(n1e)
      ds12(3) = ibm%z_bp(n2e) - ibm%z_bp(n1e)

      ds13(1) = ibm%x_bp(n3e) - ibm%x_bp(n1e)
      ds13(2) = ibm%y_bp(n3e) - ibm%y_bp(n1e)
      ds13(3) = ibm%z_bp(n3e) - ibm%z_bp(n1e)

      call crossx(ds12, ds13, tmparray)
      dr = norm2(tmparray)
      ibm%nf_x(i) = tmparray(1)/dr; ibm%nf_y(i) = tmparray(2)/dr
      ibm%nf_z(i) = tmparray(3)/dr

      !> plyunote: in original version, streamwise direction is z.
      !! now it is x. so the function is a little different
      !! the projection is z->x, x->y, y->z
      if (  (((1.0-ibm%nf_x(i))<=1.0e-6).and.((-1.0+ibm%nf_x(i))<1.0e-6)) .or.&
        (((1.0+ibm%nf_x(i))<=1.0e-6).and.((-1.0-ibm%nf_x(i))<1.0e-6)) ) then
        ibm%ns_y(i) = 1.0_wp; ibm%ns_z(i) = 0.0_wp; ibm%ns_x(i) = 0.0_wp
        ibm%nt_y(i) = 0.0_wp; ibm%nt_z(i) = 1.0_wp; ibm%nt_x(i) = 0.0_wp
      else
        rr = sqrt(ibm%nf_y(i)**2+ibm%nf_z(i)**2)
        ibm%ns_y(i) = ibm%nf_y(i) / rr
        ibm%ns_z(i) = -ibm%nf_y(i) / rr
        ibm%ns_x(i) = 0.0_wp

        ibm%nt_x(i) = rr
        ibm%nt_y(i) = -ibm%nf_x(i) * ibm%nf_y(i) / rr
        ibm%nt_z(i) = -ibm%nf_x(i) * ibm%nf_z(i) / rr
      endif

      !> plyunote: for triangular elements, dA is correct?
      !! answer: yes. dA=0.5*|ds12|*|ds13|*sin(A213)=0.5*|ds12 x ds13|
      ibm%dA(i) = dr/2.0_wp

      ibm%cent_x(i) = (ibm%x_bp(n1e)+ibm%x_bp(n2e)+ibm%x_bp(n3e))/3.0_wp
      ibm%cent_y(i) = (ibm%y_bp(n1e)+ibm%y_bp(n2e)+ibm%y_bp(n3e))/3.0_wp
      ibm%cent_z(i) = (ibm%z_bp(n1e)+ibm%z_bp(n2e)+ibm%z_bp(n3e))/3.0_wp

      !if (ibi .eq. 1 .and. i .eq. 10 .and. myid.eq.0) then
      !  print *, 'plyudebug, rotor_Rot, surface elmt after rot:', &
      !    ibm%cent_x(i), ibm%cent_y(i), ibm%cent_z(i)
      !endif
    enddo
  end subroutine nacelle_move

  subroutine Prescribe_para_update(time_, para, NumberOfObjects)
    implicit none
    real(wp) :: time_
    type(FSIPrePara), dimension(:) :: para
    integer :: NumberOfObjects

    integer :: ibi, idof

    do ibi = 1, NumberOfObjects
      do idof = 1, 6        
        !if (pp_wt(ibi)%dof(idof) .ne. 0) then
          pp_wt(ibi)%disp_o(idof) = pp_wt(ibi)%disp(idof)
          pp_wt(ibi)%phase(idof) = pp_wt(ibi)%phase0(idof) + &
            pp_wt(ibi)%omega(idof) * time_
          pp_wt(ibi)%disp(idof) = pp_wt(ibi)%amp(idof) * cos(pp_wt(ibi)%phase(idof))
          pp_wt(ibi)%dd(idof) = pp_wt(ibi)%disp(idof) - pp_wt(ibi)%disp_o(idof)
          pp_wt(ibi)%vel(idof) = - pp_wt(ibi)%omega(idof) * pp_wt(ibi)%amp(idof) &
            * sin(pp_wt(ibi)%phase(idof))
        !endif
      enddo
    enddo
  end subroutine Prescribe_para_update

  subroutine Prescribe_FSI_update(fsi, ibm, para, ibi)
    implicit none
    type(FSInfo) :: fsi
    type(IBMNodes) :: ibm
    type(FSIPrePara) :: para 
    integer :: ibi

    integer :: idof
    real(wp) :: rot_ang
    real(wp), dimension(3) :: rot_in, rot_out, rot_axis

    !> For large angle rotation, the sequence of rotation along three axises
    !! is important. Here we write the formula of Yaw-Pitch-Roll method
    !! based on a set of  Body (Z-Y-X) rotations
    !> Usually, we only turn on the rotation along a single axis. So the
    !! sequences will have no effect on final result.
    
    fsi%x_co = fsi%x_c
    fsi%y_co = fsi%y_c
    fsi%z_co = fsi%z_c

    !rot_in(1) = fsi%x_c0 - (para%refpos(1) + para%disp0(1))
    rot_in(1) = fsi%x_c0 - (para%refpos(1))
    rot_in(2) = fsi%y_c0 - (para%refpos(2))
    rot_in(3) = fsi%z_c0 - (para%refpos(3))

    !if (ibi .eq. 7 .and. myid.eq.0) then 
      !print *, 'FSIUpdate1 rot_in:myid, x_c=', myid, fsi%x_c, fsi%y_c, fsi%z_c
      !print *, 'calc rot_in:refpos=', para%refpos(1:6)
      !print *, 'calc rot_in:disp_o=', para%disp_o(1:6)
    !endif
    
    !> here we use fsi%disp from last timestep, so make sure it was updated
    !rot_in(1:3) = fsi(ibi)%disp(1:3)
    rot_out(1:3) = rot_in(1:3)

    ! disp_i = xc_i - refpos; disp(t=0)=0; (difficult in restart)
    !rot_in(1:3) = fsi(ibi)%disp_i(1:3)

    !rot_ang = para(ibi)%disp(6)
    !rot_ang = para(ibi)%disp(6) - para(ibi)%disp_o(6)
    !rot_axis = (/0.0, 0.0, 1.0/)
    !call ArbitraryRotate(rot_in, rot_ang, rot_axis, rot_out)
    
    !rot_in = rot_out
    !rot_ang = para(ibi)%disp(5) - para(ibi)%disp_o(5)
    !rot_axis = (/0.0, 1.0, 0.0/)
    !call ArbitraryRotate(rot_in, rot_ang, rot_axis, rot_out)
     
    !rot_in = rot_out
    !rot_ang = para(ibi)%disp(4) - para(ibi)%disp_o(4)
    !rot_axis = (/1.0, 0.0, 0.0/)
    !call ArbitraryRotate(rot_in, rot_ang, rot_axis, rot_out)


    do idof = 6, 4, -1
      if (para%dof(idof) .ne. 0) then
        !if (ibi.eq.7 .and. myid.eq.0) print *, 'idof, rot_in:',idof, rot_in(1:3)

        rot_ang = para%disp(idof)
        rot_axis = 0.0; rot_axis(idof-3) = 1.0
        call ArbitraryRotate(rot_in, rot_ang, rot_axis, rot_out)
        rot_in = rot_out
      endif
    enddo

    !if (ibi.eq.7 .and. myid.eq.0) print *, 'idof, rot_out:',rot_out(1:3)

    fsi%disp(1:3) = rot_out(1:3)
    !fsi%x_c = para%refpos(1) + para%disp(1) + rot_out(1) 
    fsi%x_c = para%refpos(1) + rot_out(1) 
    fsi%y_c = para%refpos(2) + rot_out(2) 
    fsi%z_c = para%refpos(3) + rot_out(3)
    fsi%disp(4:6) = para%disp(4:6)
    
    !if (ibi.eq.7 .and. myid.eq.0) print *, 'idof, final fsi x_c', fsi%x_c, &
    !  fsi%y_c, fsi%z_c
    
    !> update velocity, here rot_in is omega, rot_axis is r
    rot_in(1:3) = 1.0 * para%vel(4:6)
    rot_axis(1:3) = rot_out(1:3)
    call crossx(rot_in, rot_axis, rot_out)
    fsi%vel(1:3) = para%vel(1:3) + rot_out(1:3)
    fsi%vel(4:6) = para%vel(4:6)
    
    !> update turbine rotor axis direction
    rot_in(1:3) = (/fsi%nx_tb0, fsi%ny_tb0, fsi%nz_tb0/)
    rot_out(1:3) = rot_in(1:3)
     
    do idof = 6, 4, -1
      if (para%dof(idof) .ne. 0) then
        !rot_ang = para%disp(idof) - para%disp_o(idof)
        rot_ang = para%disp(idof)
        rot_axis = 0.0; rot_axis(idof-3) = 1.0
        call ArbitraryRotate(rot_in, rot_ang, rot_axis, rot_out)
        rot_in = rot_out
      endif
    enddo

    fsi%nx_tb = rot_out(1); fsi%ny_tb = rot_out(2); fsi%nz_tb = rot_out(3)
    
  end subroutine Prescribe_FSI_update

  subroutine interp_2p(x0, y0, x1, y1, x2, y2)
    implicit none
    real(wp), intent(in) :: x0, x1, y1, x2, y2
    real(wp), intent(out) :: y0
    
    real(wp) :: fac, epsi
    
    epsi = 1.0e-6_wp
    if ((x2 - x1) < epsi) then
      y0 = y1
    else
      fac = (x0 - x1)/(x2 - x1)
      y0 = y1*(1.0 - fac) + y2*fac
    endif
  end subroutine interp_2p
  
  subroutine interp_in_table(x0, y0, xarr, yarr, narr, inrange)
    implicit none
    real(wp), intent(in) :: x0
    real(wp), intent(out) :: y0
    integer, intent(in) :: narr
    real(wp), dimension(:), intent(in) :: xarr(narr), yarr(narr)
    integer, intent(out) :: inrange
    
    integer :: i
    
    inrange = 0 
    do i = 1, (narr-1)
      if (x0.ge.xarr(i) .and. x0.le.xarr(i+1)) then
        call interp_2p(x0, y0, xarr(i), yarr(i), xarr(i+1), yarr(i+1))
        inrange = 1
        exit
      endif
      !if (myid.eq.0) print *, i, inrange, x0, xarr(i), xarr(i+1), yarr(i), yarr(i+1)
    enddo
  end subroutine interp_in_table

  !> calculate force at lagrangian points using actuator line model
  !! called in solver.c before solving the momentum equation
  subroutine Calc_F_lagr_ACL(ibm, fsi, NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi
    integer :: NumberOfObjects

    real(wp), dimension(3,3) :: dxdq, dqdx

    integer :: n_elmt, ibi, l, j, nv1, nv2, nv3, ifoil, ll
    real(wp) :: A_sum, U_ref, rr, rx, ry, rz, tmp, u_tangent
    real(wp), dimension(3) :: U_b, n_relV, n_L, n_blade, n_rot
    real(wp) :: fac1, fac2, r, r1, r2

    real(wp) :: dmax_AOA, isgn, theta, refangle_AL
    integer :: maxiteration_rotormodel, icount, i
    real(wp), allocatable, dimension(:,:) :: coor_refblade
    real(wp), dimension(3) :: p, q, nr, na, nt

    real(wp) :: ux, uy, uz, fac, Ublade, Ur
    integer :: inrange, num_CD, num_CL

    real(wp) :: factor, Fa, Ft, angle_twist, angle_AOA, angle_inflow
    real(wp) :: nx, F_correc
    real(wp), dimension(3) :: tmparray

    real(wp) :: UU, Taumag, r_amp
    real(wp), allocatable, dimension(:) :: AOA_old 
    real(wp), dimension(3) :: Ureldirection, circulationdirection
    real(wp), dimension(3) :: Urel_noinduce, r_vec, Tau_vec

    !> 3D correctoin
    real(wp) :: CD_zero_angle, angle_zero_lift, &
      c_v_r, R_v_r, corr_3d_fl, corr_3d_fd

    real(wp) :: shen_R, shen_pi, shen_phi1, shen_Lambda, shen_g1
    real(wp) :: shen_g1_Fa, shen_F1, shen_F1_Fa
    do ibi = 1,NumberOfObjects
      n_elmt = ibm(ibi)%n_elmt
      do l = 2, n_elmt
        rx = ibm(ibi)%cent_x(l-1) - fsi(ibi)%x_c
        ry = ibm(ibi)%cent_y(l-1) - fsi(ibi)%y_c
        rz = ibm(ibi)%cent_z(l-1) - fsi(ibi)%z_c
        r1 = sqrt(rx**2 + ry**2 + rz**2)+1.d-20
        
        rx = ibm(ibi)%cent_x(l) - fsi(ibi)%x_c
        ry = ibm(ibi)%cent_y(l) - fsi(ibi)%y_c
        rz = ibm(ibi)%cent_z(l) - fsi(ibi)%z_c
        r2 = sqrt(rx**2 + ry**2 + rz**2)+1.d-20

        !> plyunote: here near_hub_inflow_correction is ignored   
      enddo

      dmax_AOA = 10000.0
      icount = 0

      allocate(coor_refblade(n_elmt,3))

      do i = 1, n_elmt
        p(1) = ibm(ibi)%cent_x(i) - fsi(ibi)%x_c
        p(2) = ibm(ibi)%cent_y(i) - fsi(ibi)%y_c
        p(3) = ibm(ibi)%cent_z(i) - fsi(ibi)%z_c
        na(1) = fsi(ibi)%nx_tb
        na(2) = fsi(ibi)%ny_tb
        na(3) = fsi(ibi)%nz_tb
              
        isgn = (1.0e-16_wp+fsi(ibi)%angvel_axis) &
          /(dabs(fsi(ibi)%angvel_axis)+1.0e-16_wp)
        if(isgn > 0.0_wp) then
          isgn = 1.0_wp
        else
          isgn = -1.0_wp
        endif
  
        !> plyuinfo: refangle_AL should be added to turbine.inp or elsewhere
        refangle_AL = 20.0_wp
        theta = isgn * refangle_AL * twopi/2.0_wp/180.0_wp

        if (rotor_model .eq. 5) theta = -theta

        call ArbitraryRotate(p, theta, na, q)

        coor_refblade(i,1) = q(1) + fsi(ibi)%x_c
        coor_refblade(i,2) = q(2) + fsi(ibi)%y_c
        coor_refblade(i,3) = q(3) + fsi(ibi)%z_c
      enddo

      allocate(AOA_old(n_elmt))
      ibm(ibi)%Uinduced(:,:) = 0.0_wp
      AOA_old(:) = 1000.0_wp
      
      maxiteration_rotormodel = 1
      if (rotor_model .eq. 5 .or. rotor_model .eq. 6 .or. rotor_model.eq.7) then
        maxiteration_rotormodel = 100
      endif

      do icount = 1, maxiteration_rotormodel
        dmax_AOA = 0.0_wp;
        do l = 1, n_elmt
          dmax_AOA = dmax_AOA + dabs(AOA_old(l)-ibm(ibi)%angle_attack(l))/n_elmt
        enddo
        
        !if(maxiteration_rotormodel > 1 .and. myid.eq.0) then
        !  print *, 'Rotor model, iter:',icount,', Convergence of AOA', dmax_AOA
        !endif
        
        if(dabs(dmax_AOA)<0.01_wp) exit

        do l = 1, n_elmt
          AOA_old(l) = ibm(ibi)%angle_attack(l)
        enddo

        do l = 1, n_elmt
          nv1 = ibm(ibi)%nv1(l); nv2 = ibm(ibi)%nv2(l); nv3 = ibm(ibi)%nv3(l)
          rx = ibm(ibi)%cent_x(l) - fsi(ibi)%x_c
          ry = ibm(ibi)%cent_y(l) - fsi(ibi)%y_c
          rz = ibm(ibi)%cent_z(l) - fsi(ibi)%z_c

          r = sqrt(rx**2+ry**2+rz**2)+1.0e-16_wp
          n_blade(1) = rx/r; n_blade(2) = ry/r; n_blade(3) = rz/r

          if (rotor_model.eq.3 .or. rotor_model.eq.5 .or. rotor_model.eq.6 .or. &
            rotor_model .eq. 7) then
            !> rotor_model: 3,6 for AL; 5 for AS
            ux = 0.5_wp*(ibm(ibi)%u(nv1,1)+ibm(ibi)%u(nv2,1))
            uy = 0.5_wp*(ibm(ibi)%u(nv1,2)+ibm(ibi)%u(nv2,2))
            uz = 0.5_wp*(ibm(ibi)%u(nv1,3)+ibm(ibi)%u(nv2,3))
          else if (rotor_model .eq. 2) then
            ux = 1.0_wp/3.0_wp*(ibm(ibi)%u(nv1,1)+ibm(ibi)%u(nv2,1)+ibm(ibi)%u(nv3,1))
            uy = 1.0_wp/3.0_wp*(ibm(ibi)%u(nv1,2)+ibm(ibi)%u(nv2,2)+ibm(ibi)%u(nv3,2))
            uz = 1.0_wp/3.0_wp*(ibm(ibi)%u(nv1,3)+ibm(ibi)%u(nv2,3)+ibm(ibi)%u(nv3,3))
          endif

          Ublade = sqrt(ux**2+uy**2+uz**2)
          !> n_rot = Omega x r
          n_rot(1) = ux/(Ublade+1.d-20)
          n_rot(2) = uy/(Ublade+1.d-20)
          n_rot(3) = uz/(Ublade+1.d-20)

          !print *,'plyudebug,Ublade=',Ublade,',n_rot=',n_rot(1),n_rot(2),n_rot(3)

          do ll = 1,3
            ibm(ibi)%rotationdirection(l,ll) = n_rot(ll)
          enddo

          ibm(ibi)%Urel(l,1) = ibm(ibi)%U_lagr_x(l) - ux
          ibm(ibi)%Urel(l,2) = ibm(ibi)%U_lagr_y(l) - uy
          ibm(ibi)%Urel(l,3) = ibm(ibi)%U_lagr_z(l) - uz
          if (fsitype .eq. 1) then
            ibm(ibi)%Urel(l,1:3) = ibm(ibi)%Urel(l,1:3) - fsi(ibi)%vel(1:3)
          endif
          
          if (max(ibm(ibi)%Urel(l,1),ibm(ibi)%Urel(l,2),ibm(ibi)%Urel(l,3))>1.0e5_wp) then
            print *,'----------plyudebug, Calc_F_lagr, myid=',myid,', l=',l,'--------'
            print *,'plyudebug,U_lagr_xyz=',ibm(ibi)%U_lagr_x(l),ibm(ibi)%U_lagr_y(l)&
              ,ibm(ibi)%U_lagr_z(l)
            print *,'plyudebug,uxyz=', ux, uy, uz
            print *,'plyudebug,Urel=',ibm(ibi)%Urel(l,1),ibm(ibi)%Urel(l,2),&
              ibm(ibi)%Urel(l,3)
          endif

          !> plyunote: nearhub is ignored here

          do ll = 1,3
            ibm(ibi)%Urel(l,ll) = ibm(ibi)%Urel(l,ll) - ibm(ibi)%Uinduced(l,ll)
          enddo

          Ur = ibm(ibi)%Urel(l,1)*n_blade(1) +  ibm(ibi)%Urel(l,2)*n_blade(2) &
            + ibm(ibi)%Urel(l,3)*n_blade(3)

          !> Get the cross-section Urel
          do ll = 1,3
            ibm(ibi)%Urel(l,ll) = ibm(ibi)%Urel(l,ll) - Ur*n_blade(ll)
          enddo

          !> here U_ref is for local element
          !> there is another U_ref defined for entire turbine
          U_ref = sqrt(ibm(ibi)%Urel(l,1)**2+ibm(ibi)%Urel(l,2)**2+ibm(ibi)%Urel(l,3)**2)
          if(U_ref > 1.0e5_wp) then
            print*,'U_ref out of range for element ',l,' of turbine ',ibi
          endif

          do ll = 1,3
            n_relV(ll) = ibm(ibi)%Urel(l,ll)/(U_ref+1.0e-20_wp)
          enddo

          !> plyunote: acos cannot tell positive or negative of AOA
          !tmp = -n_relV(1)*n_rot(1) - n_relV(2)*n_rot(2) - n_relV(3)*n_rot(3)
          !tmp = tmp / (1.0_wp+1.0e-9_wp)
          !ibm(ibi)%angle_attack(l) = acos(tmp) * 360.0_wp / twopi &
          !  - ibm(ibi)%angle_twist(l) - ibm(ibi)%pitch(1)

          !print *,'plyudebug,tmp=',tmp,',acos=',acos(tmp),',deg=',acos(tmp)*360.0/twopi
          !print*,'plyudebug,twist=',ibm(ibi)%angle_twist(l),',pitch1=',ibm(ibi)%pitch(1)
          
          !> asin can tell.
          !> plyunote: here n_blade is used, only valid for straight blade,
          !  without cone.
          tmparray(1:3) = 0.0_wp
          call crossx(n_relV, -n_rot, tmparray)
          tmp = tmparray(1)*n_blade(1)+tmparray(2)*n_blade(2)+tmparray(3)*n_blade(3)
          ibm(ibi)%angle_attack(l) = asin(tmp) * 360.0_wp / twopi &
            - ibm(ibi)%angle_twist(l) - ibm(ibi)%pitch(1)

          !> estimate reference incidence angle AOA_infty by multiplying 
          !> local incidence angle with factor, the factor is 0.75 for NACA0012.
          !> See: PhD-Schito-2011-Large eddy simulation of wind turbines: interaction
          !  with turbulent flow
          !  chapter 2.1.2 Effective velocity model
          if (dabs(ibm(ibi)%angle_attack(l))<10.0_wp) then
            !ibm(ibi)%angle_attack(l) = ibm(ibi)%angle_attack(l) / 0.75_wp
          endif

          !> plyunote: we make a convention here: angles provided in input file
          !  should be between -45 to 90 
          
          !if(ibm(ibi)%angle_attack(l) < 0.0_wp) then
          !  ibm(ibi)%angle_attack(l) = ibm(ibi)%angle_attack(l) + 360.0_wp
          !else if(ibm(ibi)%angle_attack(l) > 360.0_wp) then
          !  ibm(ibi)%angle_attack(l) = ibm(ibi)%angle_attack(l) - 360.0_wp
          !endif

          tmparray(1:3) = 0.0_wp
          !> tmparray is used because liftdirection is (n_elmt,3), crossx is (3,1)
          call crossx(n_relV, n_blade, tmparray)
          do ll = 1,3
            ibm(ibi)%liftdirection(l,ll) = tmparray(ll)
          enddo
          
          tmp = fsi(ibi)%nx_tb * ibm(ibi)%liftdirection(l,1) + fsi(ibi)%ny_tb &
            *ibm(ibi)%liftdirection(l,2) +fsi(ibi)%nz_tb *ibm(ibi)%liftdirection(l,3)

          if (tmp .lt. 0.0) then
            do ll = 1,3
            ibm(ibi)%liftdirection(l,ll) = -ibm(ibi)%liftdirection(l,ll)
            enddo
          endif

          CD_zero_angle = 0.0; angle_zero_lift = 0.0
          do ifoil = 1, num_foiltype
            num_CD = acl(ifoil)%num_CD
            num_CL = acl(ifoil)%num_CL
            if ( r .ge. acl(ifoil)%r_beg .and. r .le. acl(ifoil)%r_end) then
              call interp_in_table(ibm(ibi)%angle_attack(l), ibm(ibi)%CD(l), &
                acl(ifoil)%ang_CD(:), acl(ifoil)%CDInp(:), acl(ifoil)%num_CD, &
                inrange)    

              if (inrange.eq.0) then
                !> The AOA is not in the range of database, then guess it based on
                !! emperical formula.
                if (ibm(ibi)%angle_attack(l) < acl(ifoil)%ang_CD(1) .and. &
                  ibm(ibi)%angle_attack(l) > -45.0) then
                  ibm(ibi)%CD(l) = (-1.2-acl(ifoil)%CDInp(1)) &
                    * (ibm(ibi)%angle_attack(l) - acl(ifoil)%ang_CD(1)) &
                    / (-45.0-acl(ifoil)%ang_CD(1)) + acl(ifoil)%CDInp(1)
                elseif (ibm(ibi)%angle_attack(l) <= -45.0) then
                  ibm(ibi)%CD(l) = -1.2
                elseif (ibm(ibi)%angle_attack(l) > acl(ifoil)%ang_CD(num_CD) .and. &
                  ibm(ibi)%angle_attack(l) < 45.0) then
                  ibm(ibi)%CD(l) = (1.2-acl(ifoil)%CDInp(num_CD)) &
                    * (ibm(ibi)%angle_attack(l) - acl(ifoil)%ang_CD(num_CD)) &
                    / (45.0-acl(ifoil)%ang_CD(num_CD)) + acl(ifoil)%CDInp(num_CD)
                elseif (ibm(ibi)%angle_attack(l) >= 45.0) then
                  ibm(ibi)%CD(l) = 1.2
                endif

                if(myid.eq.0) print *,'attack angle is out of range of CD table',&
                  l,ibm(ibi)%angle_attack(l), ibm(ibi)%CD(l)
              endif

              call interp_in_table(ibm(ibi)%angle_attack(l), ibm(ibi)%CL(l), &
                acl(ifoil)%ang_CL(:), acl(ifoil)%CLInp(:), acl(ifoil)%num_CL, &
                inrange)              

              if (inrange.eq.0) then
                if (ibm(ibi)%angle_attack(l) < acl(ifoil)%ang_CL(1) .and. &
                  ibm(ibi)%angle_attack(l) > -45.0) then
                  ibm(ibi)%CL(l) = (-1.05-acl(ifoil)%CLInp(1)) &
                    * (ibm(ibi)%angle_attack(l) - acl(ifoil)%ang_CL(1)) &
                    / (-45.0-acl(ifoil)%ang_CL(1)) + acl(ifoil)%CLInp(1)
                elseif (ibm(ibi)%angle_attack(l) <= -45.0 .and. & 
                  ibm(ibi)%angle_attack(l) >= -90.0) then
                  ibm(ibi)%CL(l) = 1.05*sin(2.0*ibm(ibi)%angle_attack(l)*twopi/360.0)
                elseif (ibm(ibi)%angle_attack(l) < -90.0) then
                  ibm(ibi)%CL(l) = 0.0

                elseif (ibm(ibi)%angle_attack(l) > acl(ifoil)%ang_CL(num_CL) .and. &
                  ibm(ibi)%angle_attack(l) < 45.0) then
                  ibm(ibi)%CL(l) = (1.05-acl(ifoil)%CLInp(num_CL)) &
                    * (ibm(ibi)%angle_attack(l) - acl(ifoil)%ang_CL(num_CL)) &
                    / (45.0-acl(ifoil)%ang_CL(num_CL)) + acl(ifoil)%CLInp(num_CL)
                elseif (ibm(ibi)%angle_attack(l) >= 45.0 .and. & 
                  ibm(ibi)%angle_attack(l) <= 90.0) then
                  ibm(ibi)%CL(l) = 1.05*sin(2.0*ibm(ibi)%angle_attack(l)*twopi/360.0)
                elseif (ibm(ibi)%angle_attack(l) > 90.0) then
                  ibm(ibi)%CL(l) = 0.0
                endif

               !> ignore this case now
                if(myid.eq.0) print *,'attack angle is out of range of CL table',&
                  l, ibm(ibi)%angle_attack(l), ibm(ibi)%CL(l)
              endif

              CD_zero_angle = acl(ifoil)%CD_zero_angle
              angle_zero_lift = acl(ifoil)%angle_zero_lift              
            endif

          enddo

          !> save 2D coefficients
          ibm(ibi)%CD2(l) = ibm(ibi)%CD(l)
          ibm(ibi)%CL2(l) = ibm(ibi)%CL(l)

          correction3D_CH = 0
          correction3D_DS = 0
          if (correction3D .eq. 1) then
            correction3D_DS = 1
          endif
          
          !> add 3D and rotational effects, Chariaropoulos PK
          correction3D_CH = 0
          if (correction3D_CH .ne. 0) then

          endif

          !> add 3D effects, see Du and Selig, https://arc.aiaa.org/doi/pdf/10.2514/6.1998-21
          ! correction3D_DS = 1
          if (correction3D_DS .eq. 1) then
            c_v_r = ibm(ibi)%chord_blade(l) / r
            R_v_r = fsi(ibi)%r_rotor / r
            corr_3d_fl = 1.0_wp/twopi*(1.6_wp*c_v_r/0.1267_wp &
              *(1.0_wp-c_v_r**(1.0_wp/ibm(ibi)%TSR_modified*R_v_r)) &
              /(1.0_wp+c_v_r**(1.0_wp/ibm(ibi)%TSR_modified*R_v_r))-1.0_wp)
            corr_3d_fd = 1.0_wp/twopi*(1.6_wp*c_v_r/0.1267_wp &
              *(1.0_wp-c_v_r**(0.5_wp/ibm(ibi)%TSR_modified*R_v_r)) &
              /(1.0_wp+c_v_r**(0.5_wp/ibm(ibi)%TSR_modified*R_v_r))-1.0_wp)
            ibm(ibi)%CL(l) = ibm(ibi)%CL2(l) + corr_3d_fl * & 
              (twopi*(ibm(ibi)%angle_attack(l)-angle_zero_lift)/360.0_wp*twopi - ibm(ibi)%CL2(l))
            ibm(ibi)%CD(l) = ibm(ibi)%CD2(l) - corr_3d_fd * (ibm(ibi)%CD2(l) - CD_zero_angle)
          endif

          !Shen1_AL = 0
          shen_F1 = 1.0_wp; shen_F1_Fa = 1.0_wp
          !> tip loss correction
          if (Shen1_AL .ne. 0) then
            shen_R = fsi(ibi)%r_rotor
            shen_pi = twopi / 2.0_wp
            shen_phi1 = (ibm(ibi)%angle_attack(l)+ibm(ibi)%angle_twist(l) &
              +ibm(ibi)%pitch(1))*twopi/360.0_wp
            shen_Lambda = dabs(ibm(ibi)%Tipspeedratio)
            shen_g1 = 0.0_wp

            !> here a relationship between lambda and g should be determined
            !! from comparison of simulation and experiment
            !> firstly, set g as 1.0
            !> run several omega case, compare the result, calculate g

            !> for lambda=8.6, Shen estimate g=0.6488
            !> VFS estimate g=0.176

            if (shen_Lambda < 7.0636) then
              shen_g1 = -1.7636_wp*shen_Lambda+11.6112_wp
            else if (shen_Lambda >=7.0636 .and. shen_Lambda < 9.8146) then
              shen_g1 = -0.5805 * shen_Lambda + 3.2542
            else
              shen_g1 = -0.5077 * shen_Lambda + 2.5397
            endif

            shen_g1 = shen_g1 * Shen1_AL_tipcorrection
            shen_g1_Fa = shen_g1 * Shen1_AL_tipcorrectionratio_Fa

            tmp = exp(shen_g1)
            shen_F1 = 2.0_wp *acos(exp(-tmp*num_blade*(max(shen_R-r,0.0_wp)) &
              /2.0_wp/r/dabs(sin(shen_phi1))))/shen_pi
            
            tmp = exp(shen_g1_Fa)
            shen_F1_Fa = 2.0_wp *acos(exp(-tmp*num_blade*(max(shen_R-r,0.0_wp))&
              /2.0_wp/r/dabs(sin(shen_phi1))))/shen_pi

          endif

          factor = 1.0_wp
          if(rotor_model.eq.3 .or. rotor_model.eq.5 .or. rotor_model.eq.6 &
            .or. rotor_model.eq.7) then
            !> plyunote: where is chord_blade assigned?
            !> answer: read from FOILxx in airfoil_ACL, then interpolate
            factor = ibm(ibi)%chord_blade(l)
          endif

          ibm(ibi)%Urelmag(l) = U_ref

          !> plyunote: F_lagr_x is the force from blade to fluid
          !> force from drag coefficient C_D
          ibm(ibi)%F_lagr_x(l) = -0.5 * U_ref**2 * dabs(ibm(ibi)%CD(l)) &
            * n_relV(1) * factor
          ibm(ibi)%F_lagr_y(l) = -0.5 * U_ref**2 * dabs(ibm(ibi)%CD(l)) &
            * n_relV(2) * factor
          ibm(ibi)%F_lagr_z(l) = -0.5 * U_ref**2 * dabs(ibm(ibi)%CD(l)) &
            * n_relV(3) * factor

          !> force from lift coefficient C_L
          ibm(ibi)%F_lagr_x(l) = ibm(ibi)%F_lagr_x(l) - 0.5 * U_ref**2 * &
            dabs(ibm(ibi)%CL(l)) * ibm(ibi)%liftdirection(l,1) * factor 
          ibm(ibi)%F_lagr_y(l) = ibm(ibi)%F_lagr_y(l) - 0.5 * U_ref**2 * &
            dabs(ibm(ibi)%CL(l)) * ibm(ibi)%liftdirection(l,2) * factor 
          ibm(ibi)%F_lagr_z(l) = ibm(ibi)%F_lagr_z(l) - 0.5 * U_ref**2 * &
            dabs(ibm(ibi)%CL(l)) * ibm(ibi)%liftdirection(l,3) * factor 

          !> Fa is thrust force, Ft is torque force,
          !! adjustment coefficient is ignored
          Fa = (ibm(ibi)%F_lagr_x(l)*fsi(ibi)%nx_tb + ibm(ibi)%F_lagr_y(l) &
            *fsi(ibi)%ny_tb + ibm(ibi)%F_lagr_z(l)*fsi(ibi)%nz_tb) &
            * shen_F1_Fa
          Ft = (ibm(ibi)%F_lagr_x(l)*n_rot(1)+ibm(ibi)%F_lagr_y(l)*n_rot(2) &
            +ibm(ibi)%F_lagr_z(l)*n_rot(3)) * shen_F1
          
          ibm(ibi)%F_lagr_x(l) = Fa * fsi(ibi)%nx_tb + Ft * n_rot(1)
          ibm(ibi)%F_lagr_y(l) = Fa * fsi(ibi)%ny_tb + Ft * n_rot(2)
          ibm(ibi)%F_lagr_z(l) = Fa * fsi(ibi)%nz_tb + Ft * n_rot(3)

          !> plyunote: correction3D_CL is ignored here

          !> plyunote: nacelle is ignored here

        enddo ! end of n_elmt

        ! circulation on the blade
        do l = 1, n_elmt
          r = sqrt((ibm(ibi)%cent_x(l)-fsi(ibi)%x_c)**2+(ibm(ibi)%cent_y(l) &
            -fsi(ibi)%y_c)**2+(ibm(ibi)%cent_z(l)-fsi(ibi)%z_c)**2)

          Urel_noinduce(1:3) = ibm(ibi)%Urel(l,1:3) + ibm(ibi)%Uinduced(l,1:3)
          UU = norm2(Urel_noinduce)
          Ureldirection(1:3) = Urel_noinduce / (UU+1.0e-16_wp)

          do ll = 1, 3
            tmparray(ll) = ibm(ibi)%liftdirection(l,ll)
          enddo
          call crossx(tmparray, Ureldirection, circulationdirection)

          Taumag = -(ibm(ibi)%F_lagr_x(l)*ibm(ibi)%liftdirection(l,1) &
            + ibm(ibi)%F_lagr_y(l)*ibm(ibi)%liftdirection(l,2) &
            + ibm(ibi)%F_lagr_z(l)*ibm(ibi)%liftdirection(l,3) )/UU

          !do ll = 1, 3
          ibm(ibi)%circulation(l,1:3) = Taumag * circulationdirection(1:3)
          !enddo
        enddo

        if (rotor_model .eq. 5 .or. rotor_model.eq.7) then
          do l = 1, n_elmt
            ibm(ibi)%Uinduced(l,1:3) = 0.0

            do ll = 1, n_elmt
              r_vec(1) = -coor_refblade(ll,1) + ibm(ibi)%cent_x(l)
              r_vec(2) = -coor_refblade(ll,2) + ibm(ibi)%cent_y(l)
              r_vec(3) = -coor_refblade(ll,3) + ibm(ibi)%cent_z(l)
              r_amp = norm2(r_vec)

              Tau_vec(1:3) = ibm(ibi)%dA(ll) * ibm(ibi)%circulation(ll,1:3)

              fac = 1.0_wp/(2.0_wp*twopi*r_amp*r_amp*r_amp+1.0e-16_wp)

              !> plyunote: r_nacelle 
              call crossx(Tau_vec, r_vec, tmparray)
              ibm(ibi)%Uinduced(l,1:3) = ibm(ibi)%Uinduced(l,1:3) &
                + tmparray(1:3)*fac
            enddo ! end of ll = 1, n_elmt
          enddo ! end of l=1,n_elmt
        endif

        if (rotor_model .eq. 6 ) then
          do l = 1, n_elmt
            ibm(ibi)%Uinduced(l,1:3) = 0.0

            do ll = 1, n_elmt
              r_vec(1) = coor_refblade(l,1) - ibm(ibi)%cent_x(ll)
              r_vec(2) = coor_refblade(l,2) - ibm(ibi)%cent_y(ll)
              r_vec(3) = coor_refblade(l,3) - ibm(ibi)%cent_z(ll)
              r_amp = norm2(r_vec)

              Tau_vec(1:3) = ibm(ibi)%dA(ll) * ibm(ibi)%circulation(ll,1:3)

              fac = 1.0_wp/(2.0_wp*twopi*r_amp*r_amp*r_amp+1.0e-16_wp)

              !> plyunote: r_nacelle 
              call crossx(Tau_vec, r_vec, tmparray)
              ibm(ibi)%Uinduced(l,1:3) = ibm(ibi)%Uinduced(l,1:3) &
                + tmparray(1:3)*fac
            enddo ! end of ll = 1, n_elmt
          enddo ! end of l=1,n_elmt
        endif

      enddo ! end of icount
      
      deallocate(coor_refblade)
      deallocate(AOA_old)

    enddo ! end of numberofturbines
  end subroutine Calc_F_lagr_ACL

  subroutine Calc_F_lagr_nacelle (ibm, fsi, NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi
    integer :: NumberOfObjects

    integer :: l, ibi, nv1, nv2, nv3, n_elmt, ll
    ! sumf_p: sumforce_pressure, sumf_f: sumforce_friction
    real(wp) :: sumf_p, sumf_f, sumf1_p, sumf1_f, sumf2_f
    real(wp), dimension(3) :: Ut, Un, Ubt, Ubn, Ub, Ut_IP, Un_IP
    real(wp), dimension(3) :: sumf, sumf1, nf, fn, ft, Urel, U_lagr
    real(wp) :: sum_dh, dh, Ubn_amp, Un_amp


    do ibi = 1, NumberOfObjects
      sumf_p = 0.0_wp; sumf_f = 0.0_wp; sumf1_p = 0.0_wp
      sumf1_f = 0.0_wp; sumf2_f = 0.0_wp
      sumf(1:3) = 0.0_wp; sumf1(1:3) = 0.0_wp
      sum_dh = 0.0_wp
      n_elmt = ibm(ibi)%n_elmt

      do l = 1, n_elmt
        nv1=ibm(ibi)%nv1(l);nv2=ibm(ibi)%nv2(l);nv3=ibm(ibi)%nv3(l)
        U_lagr(1) = ibm(ibi)%U_lagr_x(l)
        U_lagr(2) = ibm(ibi)%U_lagr_y(l)
        U_lagr(3) = ibm(ibi)%U_lagr_z(l)
        
        dh = ibm(ibi)%dh 
        !> plyunote: dh is strange here. from equations below, 
        !! dh should change with l (index of surface elements)
        !! answer: dh is only a factor. Whatever dh is, the force will
        !! try to change fluid velocity to target boundary velocity.
        !! dh provide a flexible way to adjust the strength of this process

        sum_dh = sum_dh + dh

        nf(1) = ibm(ibi)%nf_x(l); nf(2) = ibm(ibi)%nf_y(l)
        nf(3) = ibm(ibi)%nf_z(l)

        !> Ub: boundary condition, or so-called target velocity
        do ll = 1, 3
          Ub(ll) = (ibm(ibi)%u(nv1,ll) + ibm(ibi)%u(nv2,ll) &
            + ibm(ibi)%u(nv3,ll))/3.0_wp
        enddo
        if (fsitype .eq. 1) then
          if (myid.eq.0 .and. ibi.eq.1 .and. l.eq.1) then
            print *, 'Ub=', Ub(1:3)
            print *, 'pp_wt%vel=', pp_wt((ibi-1)/NumNacellePerLoc+1)%vel(1:3)
          endif
          Ub(1:3) = Ub(1:3) + pp_wt((ibi-1)/NumNacellePerLoc+1)%vel(1:3)
        endif
        Ubn_amp = Ub(1)*nf(1)+Ub(2)*nf(2)+Ub(3)*nf(3)

        Ubn(1:3) = Ubn_amp * nf(1:3)
        Ubt(1:3) = Ub(1:3) - Ubn(1:3)

        if(nacelle_model .eq. 7) then
          ft(1:3) = 0.0_wp
        else if (nacelle_model .eq. 1) then
          Un_amp = U_lagr(1)*nf(1) + U_lagr(2)*nf(2) + U_lagr(3)*nf(3)
          Un(1:3) = Un_amp * nf(1:3)
          Ut(1:3) = U_lagr(1:3) - Un(1:3)

          !> normal force
          Urel(1:3) = Ubn(1:3) - Un(1:3)
          fn(1:3) = dh*Urel(1:3)/dt

          !> tangential force
          Urel(1:3) = Ubt(1:3) - Ut(1:3)
          ft(1:3) = dh*Urel(1:3)/dt
        endif
        !> plyunote: here wall model is ignored

        ibm(ibi)%F_lagr_x(l) = fn(1) + ft(1)
        ibm(ibi)%F_lagr_y(l) = fn(2) + ft(2)
        ibm(ibi)%F_lagr_z(l) = fn(3) + ft(3)

        if(myid.eq.0 .and. ibi.eq.1 .and. l.eq.1) print *, 'nacelle F_lagr,', ibm(ibi)%F_lagr_x(l)
        
        !> plyunote: output of force infomation is ignored here        

      enddo
    enddo
  end subroutine Calc_F_lagr_nacelle

  subroutine Export_ForceOnBlade(ibm, fsi, NumberOfObjects, fname)
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi
    integer :: NumberOfObjects
    character(len=*) :: fname

    character(128) :: filen
    integer :: ibi, l, i

    integer :: nv1, nv2, nv3
    real(wp) :: r_amp, Ublade, Ft, Fa, AOA
    real(wp), dimension(3) :: n_blade, r_vec, u_vec, n_rot

    filen = ''

    if (myid .eq. 0) then
      do ibi = 1, NumberOfObjects
        filen = ''
        write(filen,'(a8,i0.6,a1,i0.3,a4)') fname, ti,'_',ibi,'.dat'
        open(18002,file=filen,action='write',status='replace')
        
        write(18002,*) 'Variables = r Ft Fa AOA C CL CD CL2 CD2 dA &
          &Urel Urel.x Urel.y Urel.z &
          &liftdirection.x liftdirection.y liftdirection.z &
          &u_ib v_ib w_ib U_lagr_x U_lagr_y U_lagr_z Fx Fy Fz'
        
        do l = 1, ibm(ibi)%n_elmt
          nv1=ibm(ibi)%nv1(l); nv2=ibm(ibi)%nv2(l); nv3=ibm(ibi)%nv3(l)
          
          r_vec(1) = ibm(ibi)%cent_x(l) - fsi(ibi)%x_c
          r_vec(2) = ibm(ibi)%cent_y(l) - fsi(ibi)%y_c
          r_vec(3) = ibm(ibi)%cent_z(l) - fsi(ibi)%z_c
          r_amp = norm2(r_vec)

          n_blade = r_vec / (r_amp+1.d-20)

          do i = 1,3
            u_vec(i) = 0.5*(ibm(ibi)%u(nv1,i)+ibm(ibi)%u(nv2,i))
          enddo
          Ublade = norm2(u_vec)
          n_rot = u_vec / (Ublade + 1.d-20)

          Ft = -(ibm(ibi)%F_lagr_x(l)*n_rot(1)+ibm(ibi)%F_lagr_y(l)&
            *n_rot(2)+ibm(ibi)%F_lagr_z(l)*n_rot(3))
          Fa = -(ibm(ibi)%F_lagr_x(l)*fsi(ibi)%nx_tb+&
            ibm(ibi)%F_lagr_y(l)*fsi(ibi)%ny_tb+&
            ibm(ibi)%F_lagr_z(l)*fsi(ibi)%nz_tb)
          AOA = ibm(ibi)%angle_attack(l)

          !print *,'plyudebug, Export_ForceOnBlade:size(ibm%u)',&
          !  size(ibm(ibi)%u,1),',n_elmt=',ibm(ibi)%n_elmt

          write(18002,*) r_amp, Ft, Fa, AOA, ibm(ibi)%chord_blade(l),&
            ibm(ibi)%CL(l), ibm(ibi)%CD(l), ibm(ibi)%CL2(l), ibm(ibi)%CD2(l),&
            ibm(ibi)%dA(l), ibm(ibi)%Urelmag(l), ibm(ibi)%Urel(l,1), &
            ibm(ibi)%Urel(l,2), ibm(ibi)%Urel(l,3), &
            ibm(ibi)%liftdirection(l,1), ibm(ibi)%liftdirection(l,2),&
            ibm(ibi)%liftdirection(l,3),&
            u_vec(1), u_vec(2), u_vec(3), &
            ibm(ibi)%U_lagr_x(l), ibm(ibi)%U_lagr_y(l), ibm(ibi)%U_lagr_z(l),&
            ibm(ibi)%F_lagr_x(l),&
            ibm(ibi)%F_lagr_y(l), ibm(ibi)%F_lagr_z(l)
        enddo

        close(18002)
      enddo
    endif

    
  end subroutine Export_ForceOnBlade
  
  !> distribute the force on the turbines to background grid
  subroutine Calc_F_eul(wtforce_, wtforce_y_, wtforce_z_, fturbx, fturby, &
     level, ibm, fsi, NumberOfObjects, dh, mtype, flag_de)
    use spectral
    use utils
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi
    integer :: NumberOfObjects, mtype
    real(wp) :: dh
    
    !real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: wtforce_(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_y_(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_z_(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: flag_de(xsz(1), xsz(2), xsz(3))
    integer :: level
    real(wp) :: fturbx, fturby

    integer :: ioutputforce
    real(wp), allocatable, dimension(:,:) :: wtfx, wtfy, wtfz, ftemp2
    !real(wp), allocatable, dimension(:,:) :: ux, uy, uz 
    integer, allocatable, dimension(:,:) :: pid, pidtemp
    real(wp), allocatable, dimension(:,:,:) :: ftemp
    character(len=128) :: filen

    real(wp), dimension(3,3) :: dxdq, dqdx
    real(wp) :: dhx_, dhy_, dhz_, h1, h2, h3, dh1, dh2, dh3, dV, dv0, J_det

    integer :: i, j, k, l, ibi, ii, jj, kk
    integer :: xs, ys, zs, xe, ye, ze
    real(wp) :: xi, yi, zi, xk, yk, zk 
    real(wp) :: xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr
    real(wp) :: dfunci, dfunck, df1, df2, df3, dh0, dA
    real(wp) :: df1k, df2k, df3k, epsi_r
    real(wp) :: r1i, r2i, r3i, r1k, r2k, r3k, rxi, ryi, rzi, rxk, ryk, rzk
    real(wp) :: f_temp, f_max, dh_max, dh_min
    real(wp), dimension(3) :: sum_f_eul, sum_f_lag, tmparray
    
    dh0 = sqrt(xl/nx*dh_mean)
    dh1 = xl/nx
    dh2 = yl/ny
    dh3 = zl/(nz-1)
         
    xs = xst(1); xe = xst(1) + xsz(1) - 1
    ys = xst(2); ye = xst(2) + xsz(2) - 1
    zs = xst(3); ze = xst(3) + xsz(3) - 1

    !print *, 'plyudebug, myid=',myid,', dh123=',dh1,dh2,dh3

    if (mtype .eq. 1) then
      ! only do initialization in 1st call from rotor model
      ! ignore this process in 2nd call from nacelle model
      !wtforce(:,:,(1-level):(xsz(3)+level)) = 0.0_wp
      !wtforce_y(:,:,(1-level):(xsz(3)+level)) = 0.0_wp
      !wtforce_z(:,:,(1-level):(xsz(3)+level)) = 0.0_wp
      wtforce_(:,:,:) = 0.0_wp
      wtforce_y_(:,:,:) = 0.0_wp
      wtforce_z_(:,:,:) = 0.0_wp
      flag_de(:,:,:) = 0.0_wp 
    endif

!print *, "nac_eul, -1", mtype, myid; call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
    
    f_temp = 0.0_wp; f_max = 0.0_wp
    do ibi = 1, NumberOfObjects
    do l = 1, ibm(ibi)%n_elmt
      dA = ibm(ibi)%dA(l)
      do k = ibm(ibi)%k_min(l), ibm(ibi)%k_max(l)
      do j = ibm(ibi)%j_min(l), ibm(ibi)%j_max(l)
      do i = ibm(ibi)%i_min(l), ibm(ibi)%i_max(l)
        if (i>=xs .and. i<xe .and. j>=ys .and. j<=ye &
          .and. k>=zs .and. k<=ze) then 
          !> plyunote: xi,yi,zi should be overlapped with u, v nodes
          !> plyunote: xk,yk,zk should be overlapped with w nodes
          xi = cartx(i)
          yi = carty(j)
          !zi = cartz(i,j,k)
          zi = zzo(k) * (hbar + hho(i,j)) - hho(i,j)
          
          xk = cartx(i)
          yk = carty(j)
          !zk = cartzw(i,j,k)
          zk = zwo(k) * (hbar + hho(i,j)) - hho(i,j)
          
          xc = cartx(i)
          yc = carty(j)
          !zc = cartz(i,j,k)
          zc = zzo(k) * (hbar + hho(i,j)) - hho(i,j)

          rxi = xi - ibm(ibi)%cent_x(l);  ryi = yi - ibm(ibi)%cent_y(l)
          rzi = zi - ibm(ibi)%cent_z(l)
          rxk = xk - ibm(ibi)%cent_x(l);  ryk = yk - ibm(ibi)%cent_y(l)
          rzk = zk - ibm(ibi)%cent_z(l)

          !> plyunote: this section is implemented in a different form from
          !> fotis's version
          !dxdq(1,1) = 1._wp; dxdq(1,2) = 0._wp; dxdq(1,3) = 0._wp
          !dxdq(2,1) = 0._wp; dxdq(2,2) = 1._wp; dxdq(2,3) = 0._wp
          !dxdq(3,1) = hxo(i,j)*(-1.0+zzo(k)); dxdq(3,2) = hyo(i,j)*(-1.0+zzo(k))
          !dxdq(3,3) = hbar + eo(i,j)
          !J_det = hbar + eo(i,j)

          !dqdx(1,1) = 1.0_wp; dqdx(1,2) = 0.0_wp; dqdx(1,3) = 0.0_wp;
          !dqdx(2,1) = 0.0_wp; dqdx(2,2) = 1.0_wp; dqdx(2,3) = 0.0_wp;
          !dqdx(3,1) = hxo(i,j)/(hbar+eo(i,j))*(1.0_wp-zzo(k))
          !dqdx(3,2) = hyo(i,j)/(hbar+eo(i,j))*(1.0_wp-zzo(k))
          !dqdx(3,3) = 1.0_wp/(hbar+eo(i,j))

          !dhx_ = sqrt(dxdq(1,1)**2+dxdq(2,1)**2+dxdq(3,1)**2)
          !dhy_ = sqrt(dxdq(1,2)**2+dxdq(2,2)**2+dxdq(3,2)**2)
          !dhz_ = sqrt(dxdq(1,3)**2+dxdq(2,3)**2+dxdq(3,3)**2)
          !h1 = sqrt(dqdx(1,1)**2+dqdx(2,1)**2+dqdx(3,1)**2)
          !h2 = sqrt(dqdx(1,2)**2+dqdx(2,2)**2+dqdx(3,2)**2)
          !h3 = sqrt(dqdx(1,3)**2+dqdx(2,3)**2+dqdx(3,3)**2)
          
          !dh_max = max(dhx_, dhy_, dhz_)
          !dh_min = min(dhx_, dhy_, dhz_)
          !if(max(dabs(dh_max),dabs(dh_min))<1.0d-6) then
          !  print *, 'plyudebug, myid=',myid,',i,j,k=',i,j,k,',dhxyz=',dhx_,dhy_,dhz_
          !endif
          
          dh3 = dz(k-xst(3)+1) * (hbar+eo(i,j))
          !dh3 = dz(k-xst(3)+1) !> for testing 
          dV = dh1 * dh2 * dh3 
          !dV = dh1 * dh2 * dh3 / her(i-xst(1)+1,j-xst(2)+1)
          
          !rx = xc - cent_x; ry = yc - cent_y; rz = zc - cent_z
          !r1i = (rxi*dxdq(1,1) + ryi*dxdq(2,1) + rzi*dxdq(3,1))/dhx_/dhx_/dh1
          !r2i = (rxi*dxdq(1,2) + ryi*dxdq(2,2) + rzi*dxdq(3,2))/dhy_/dhy_/dh2
          !r3i = (rxi*dxdq(1,3) + ryi*dxdq(2,3) + rzi*dxdq(3,3))/dhz_/dhz_/dh3
          !r1k = (rxk*dxdq(1,1) + ryk*dxdq(2,1) + rzk*dxdq(3,1))/dhx_/dhx_/dh1
          !r2k = (rxk*dxdq(1,2) + ryk*dxdq(2,2) + rzk*dxdq(3,2))/dhy_/dhy_/dh2
          !r3k = (rxk*dxdq(1,3) + ryk*dxdq(2,3) + rzk*dxdq(3,3))/dhz_/dhz_/dh3
          
          !r1i = rxi*dqdx(1,1)/h1+ryi*dqdx(1,2)/h2+rzi*dqdx(1,3)/h3/dh1
          !r2i = rxi*dqdx(2,1)/h1+ryi*dqdx(2,2)/h2+rzi*dqdx(2,3)/h3/dh2
          !r3i = rxi*dqdx(3,1)/h1+ryi*dqdx(3,2)/h2+rzi*dqdx(3,3)/h3/dh3
          !r1k = rxk*dqdx(1,1)/h1+ryk*dqdx(1,2)/h2+rzk*dqdx(1,3)/h3/dh1
          !r2k = rxk*dqdx(2,1)/h1+ryk*dqdx(2,2)/h2+rzk*dqdx(2,3)/h3/dh2
          !r3k = rxk*dqdx(3,1)/h1+ryk*dqdx(3,2)/h2+rzk*dqdx(3,3)/h3/dh3

          r1i = rxi/dh1; r2i = ryi/dh2;  r3i = rzi/dh3
          r1k = rxk/dh1; r2k = ryk/dh2;  r3k = rzk/dh3
          
          !> delta function: use dfunc_s4h in Yang(JCP,2009)
          !! doi:10.1016/j.jcp.2009.07.023
          !> Or Gaussian function: 2D and 3D version see Martinez-Tossas(Wind Energy, 2017)
          !! doi: 10.1002/we.2081

          !> plyunote: as a temporary solution, the epsilon for Gaussian distrubtion is set as 
          !!           0.17*c, which is good for a flat plate, see Martinez-Tossas's paper.
          !!           For other type of airfoil section, the coefficient might need to be set
          !!           by a parameter in user control dictionary.
          

          if (deltafunc_F .eq. 0) then
            !> for Gaussian distribution, sigma=0.75 have a similar effect with dfunc_s1h
            !> we add a max func here to avoid that 0.17*r is too small and make
            !  the aerodynamic force concentrate at very few grids.
            epsi_r = max(0.17 * ibm(ibi)%chord_blade(l), 0.75*dh1)
            call dfunc_gaussian(rxi/dh1, epsi_r/dh1, df1)
            call dfunc_gaussian(rxk/dh1, epsi_r/dh1, df1k)
            epsi_r = max(0.17 * ibm(ibi)%chord_blade(l), 0.75*dh2)
            call dfunc_gaussian(ryi/dh2, epsi_r/dh2, df2)
            call dfunc_gaussian(ryk/dh2, epsi_r/dh2, df2k)
            epsi_r = max(0.17 * ibm(ibi)%chord_blade(l), 0.75*dh3)
            call dfunc_gaussian(rzi/dh3, epsi_r/dh3, df3)
            call dfunc_gaussian(rzk/dh3, epsi_r/dh3, df3k)
            dfunci = df1 * df2 * df3
            dfunck = df1k * df2k * df3k
          else if (deltafunc_F .eq. 1) then
            call dfunc_s1h(r1i, df1); 
            call dfunc_s1h(r2i,df2); call dfunc_s1h(r3i,df3)
            dfunci = df1 * df2 * df3
            call dfunc_s1h(r1k, df1); 
            call dfunc_s1h(r2k,df2); call dfunc_s1h(r3k,df3)
            dfunck = df1 * df2 * df3
          else if (deltafunc_F .eq. 2) then
            call dfunc_s2h(r1i, df1); 
            call dfunc_s2h(r2i,df2); call dfunc_s2h(r3i,df3)
            dfunci = df1 * df2 * df3
            call dfunc_s2h(r1k, df1); 
            call dfunc_s2h(r2k,df2); call dfunc_s2h(r3k,df3)
            dfunck = df1 * df2 * df3
          else if (deltafunc_F .eq. 3) then
            call dfunc_s3h(r1i, df1); 
            call dfunc_s3h(r2i,df2); call dfunc_s3h(r3i,df3)
            dfunci = df1 * df2 * df3
            call dfunc_s3h(r1k, df1); 
            call dfunc_s3h(r2k,df2); call dfunc_s3h(r3k,df3)
            dfunck = df1 * df2 * df3
          else if (deltafunc_F .eq. 4) then
            call dfunc_s4h(r1i, df1); 
            call dfunc_s4h(r2i,df2); call dfunc_s4h(r3i,df3)
            dfunci = df1 * df2 * df3
            call dfunc_s4h(r1k, df1); 
            call dfunc_s4h(r2k,df2); call dfunc_s4h(r3k,df3)
            dfunck = df1 * df2 * df3
          else
            print *, 'deltafunc_F: should be 2, 3, or 4'
          endif         
          
          !dfunc = df1 * df2 * df3
          
          ii = i - xst(1) + 1
          jj = j - xst(2) + 1
          kk = k - xst(3) + 1
          !> plyunote: please check U_lagr_x initialization
          !wtforce(ii,jj,kk) = wtforce(ii,jj,kk) + ibm(ibi)%F_lagr_x(l)*dA*dfunci/dV
          !wtforce_y(ii,jj,kk) = wtforce_y(ii,jj,kk) + ibm(ibi)%F_lagr_y(l)*dA*dfunci/dV
          !wtforce_z(ii,jj,kk) = wtforce_z(ii,jj,kk) + ibm(ibi)%F_lagr_z(l)*dA*dfunck/dV
          
          wtforce_(ii,jj,kk) = wtforce_(ii,jj,kk) + ibm(ibi)%F_lagr_x(l)*dA*dfunci/dV
          wtforce_y_(ii,jj,kk) = wtforce_y_(ii,jj,kk) + ibm(ibi)%F_lagr_y(l)*dA*dfunci/dV
          wtforce_z_(ii,jj,kk) = wtforce_z_(ii,jj,kk) + ibm(ibi)%F_lagr_z(l)*dA*dfunck/dV

          !> flag_de is used to identify discontinuity in fluid field.
          !> you can limit the flag to blade force only
          !if (mtype.eq.1 .and. (abs(dfunci)>1e-3 .or. abs(dfunck)>1e-3)) flag_de(ii,jj,kk) = 1.0_wp 
          !> or apply on both the blade force and nacelle force
          if ( (abs(dfunci)>1e-3 .or. abs(dfunck)>1e-3)) flag_de(ii,jj,kk) = 1.0_wp 

          !if (abs(dfunck)>1e-3 .and. (kk.eq.1 .or. kk.eq.2 &
          !  .or. kk.eq.xsz(3) .or. kk.eq.(xsz(3)-1)) .and. i.eq.52) then
          !if (abs(dfunci)>1e-2 .and. mtype.eq.2) then
          !  print *, myid, l, i, j, k, ii, jj, kk, ibm(ibi)%k_min(l),&
          !    ibm(ibi)%k_max(l), ibm(ibi)%F_lagr_x(l), dA, dV, &
          !    wtforce_(ii,jj,kk)
          !endif
          
          !> plyunote: I think it needn't any mpi communication here.

          f_temp = sqrt(wtforce_(ii,jj,kk)**2+wtforce_y_(ii,jj,kk)**2&
            +wtforce_z_(ii,jj,kk)**2)
          if(f_max < f_temp) f_max = f_temp
        endif
      enddo
      enddo
      enddo
    enddo
    enddo

!print *, "nac_eul, 0", mtype, myid; call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)

    if (nacelle_model .eq. 0 .or. (nacelle_model .ne.0 .and. mtype .ne. 1)) then
      do i = 1, xsz(1)
        do j = 1, xsz(2)
          do k = 1, xsz(3)
            if (isnan(wtforce_(i,j,k))) print *, "wtf_x nan at ", myid, xst(1:3), i, j, k
            if (isnan(wtforce_y_(i,j,k))) print *, "wtf_y nan at ", myid, xst(1:3), i, j, k
            if (isnan(wtforce_z_(i,j,k))) print *, "wtf_z nan at ", myid, xst(1:3), i, j, k
          enddo
        enddo
      enddo
    endif
   
   !> plyunote: As Di Yang did in the early code, the dealiasxy is not done in
   !!the turbine force calculation, instead it is done after velocity prediction
   !!step.
   !print *, "nac_eul, 0.1, myid" 
   ! call dealiasxy(wtforce_(:,:,1:xsz(3)))
   !print *, "nac_eul, 0.2, myid" 
   ! call dealiasxy(wtforce_y_(:,:,1:xsz(3)))
   !print *, "nac_eul, 0.3, myid" 
   ! call dealiasxy(wtforce_z_(:,:,1:xsz(3)))

!print *, "nac_eul, 1", mtype, myid; call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
    ! export wtforce to tecplot file
    ioutputforce = 0
    if(ioutputforce .ne. 0 .and. (nacelle_model .eq. 0 .or. &
      (nacelle_model .ne.0 .and. mtype .ne. 1)) .and. mod(ti, 20).eq.1) then
      !print *, 'temp0, allocate'
      allocate(wtfx(ny, nz), wtfy(ny, nz), wtfz(ny, nz))
      !allocate(ux(ny, nz), uy(ny, nz), uz(ny, nz))
      allocate(pid(ny, nz))
       
      !print *, 'temp1, gather'
      
      i=52
      ii = i - xst(1) + 1
      allocate(ftemp2(ny, nz), pidtemp(ny, nz))
      ftemp2 = 0.0_wp
      pidtemp = 0

      !print *, 'temp1.1, wtfx'
      ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
        wtforce_(ii,:,1:xsz(3))
      call MPI_Allreduce(ftemp2, wtfx, ny*nz, mpi_double_precision, &
        mpi_sum, mpi_comm_2d_cart,ierr_wt)
      
      !print *, 'temp1.2, wtfy'
      ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
        wtforce_y_(ii,:,1:xsz(3))
      call MPI_Allreduce(ftemp2, wtfy, ny*nz, mpi_double_precision, &
        mpi_sum, mpi_comm_2d_cart,ierr_wt)
      
      !print *, 'temp1.3, wtfz'
      ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
        wtforce_z_(ii,:,1:xsz(3))
      call MPI_Allreduce(ftemp2, wtfz, ny*nz, mpi_double_precision, &
        mpi_sum, mpi_comm_2d_cart,ierr_wt)
      
      !ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
      !  u(ii,:,1:xsz(3))
      !call MPI_Allreduce(ftemp2, ux, ny*nz, mpi_double_precision, &
      !  mpi_sum, mpi_comm_2d_cart,ierr_wt)
      !
      !ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
      !  v(ii,:,1:xsz(3))
      !call MPI_Allreduce(ftemp2, uy, ny*nz, mpi_double_precision, &
      !  mpi_sum, mpi_comm_2d_cart,ierr_wt)
      !
      !ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
      !  w(ii,:,1:xsz(3))
      !call MPI_Allreduce(ftemp2, uz, ny*nz, mpi_double_precision, &
      !  mpi_sum, mpi_comm_2d_cart,ierr_wt)
      
      !print *, 'temp1.4, pid'
      pidtemp(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
        myid
      call MPI_Allreduce(pidtemp, pid, ny*nz, mpi_integer, &
        mpi_sum, mpi_comm_2d_cart,ierr_wt)



      
      !allocate(ftemp(nx_global, ny_global, nz_global))
      !call gather_3d_xyz(wtforce(:,:,1:xsz(3)), ftemp, 0)
      !wtfx = ftemp(1, 1:ny, 1:nz)
      !call gather_3d_xyz(wtforce_y(:,:,1:xsz(3)), ftemp, 0)
      !wtfy = ftemp(1, 1:ny, 1:nz)
      !call gather_3d_xyz(wtforce_z(:,:,1:xsz(3)), ftemp, 0)
      !wtfz = ftemp(1, 1:ny, 1:nz)
      
      if (myid.eq.0) then
        write(filen, '(a,i0.10,a)') 'wtforce_', ti, '.dat'
        !print *, 'temp2, prepare filename ', filen
        open(18008, file=filen, status='replace')
        write(18008, *) "VARIABLES = X, Y, Z, FX, FY, FZ, PID"
        write(18008, *) 'ZONE T="', ti*dt, '" I=1 J=', ny, ' K=',nz,&
          ' SOLUTIONTIME=', ti*dt
        !i = 52
        do k = 1, nz
          do j = 1, ny
            !do i = 52, 52
            zc = zzo(k) * (hbar + hho(i,j)) - hho(i,j)
              write(18008, *) cartx(i), carty(j), zc, &
                wtfx(j,k), wtfy(j,k), wtfz(j,k), pid(j,k)
            !enddo
          enddo
        enddo
        close(18008)
      endif

      deallocate(wtfx, wtfy, wtfz, pid)
      deallocate(ftemp2, pidtemp)
      !print *, 'temp3, end'
    endif


    
    !print*,'plyudebug, ti=',ti,', myid=', myid, ', Calc_F_eul, f_max=', f_max

    !> check the sum of wtforce
    sum_f_eul(1:3) = 0.0_wp; sum_f_lag(1:3) = 0.0_wp; tmparray(1:3) = 0.0_wp 
    do i = 1, xsz(1)
    do j = 1, xsz(2)
    do k = 1, xsz(3)
      sum_f_eul(1) = sum_f_eul(1) + wtforce_(i,j,k)*(xl/nx*yl/ny*dz(k)/her(i,j))
      sum_f_eul(2) = sum_f_eul(2) + wtforce_y_(i,j,k)*(xl/nx*yl/ny*dz(k)/her(i,j))
      sum_f_eul(3) = sum_f_eul(3) + wtforce_z_(i,j,k)*(xl/nx*yl/ny*dz(k)/her(i,j))

      !if (abs(wtforce(i,j,k))>1e-6 .and. (k.eq.1 .or. k.eq.2 .or. k.eq.xsz(3) &
      !  .or. k.eq.(xsz(3)-1)) ) then
      !  print *, i, j, k, xst(1)+i-1, xst(2)+j-1, xst(3)+k-1, &
      !    wtforce(i,j,k)*(xl/nx*yl/ny*dz(k)/her(i,j))
      !endif

    enddo
    enddo
    enddo
call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
!print *, "nac_eul, 2", mtype, myid
    call mpi_allreduce(sum_f_eul,tmparray,3,mpi_double_precision,mpi_sum,&
      mpi_comm_2d_cart,ierr_wt)

    sum_f_eul(1:3) = tmparray(1:3)

    do ibi = 1, NumberOfObjects
      do l = 1, ibm(ibi)%n_elmt
        dA = ibm(ibi)%dA(l)
        sum_f_lag(1) = sum_f_lag(1) + ibm(ibi)%F_lagr_x(l)*dA
        sum_f_lag(2) = sum_f_lag(2) + ibm(ibi)%F_lagr_y(l)*dA
        sum_f_lag(3) = sum_f_lag(3) + ibm(ibi)%F_lagr_z(l)*dA
      enddo
    enddo
call MPI_Barrier(mpi_comm_2d_cart, ierr_wt)
!print *, "nac_eul, 3", mtype, myid
    !> plyunote: need zl*hbar?
    if (mtype .eq. 1) then
      if (InletRelaxation .eq. 0) then
        !fturbx = -sum_f_eul(1) / (xl*yl*zl*hbar)
        fturbx = 0.0_wp
        do ibi = 1, NumberOfObjects
          !fturbx = (fturbx + 4.0_wp/3.0_wp * (ibm(ibi)%U_ref)**2 * twopi/2.0_wp *
          !  (fsi(ibi)%r_rotor) **2) / (xl*yl*zl*hbar)
          
          !> Q: What's this used for?
          !> A: To compensate the drag of turbines in the flow field. But later
          !     I found that this force may be several times larger than
          !     pressure gradient, which might lead to a larger top velocity and
          !     need a long running time for a steady flow. 

          !>    Therefore, c_bforce is recommended to be zero (no compensation), 
          !     or thrust coefficient. likely to be 3/4
          
          fturbx = fturbx + (c_bforce * (fsi(ibi)%angvel_fixed * fsi(ibi)%r_rotor &
            / fsi(ibi)%Tipspeedratio)**2 * twopi/2.0_wp * (fsi(ibi)%r_rotor) **2) &
            / (xl*yl*zl*hbar)
        enddo
        if(myid.eq.0 .and. ti .eq. (ti_first+1)) then
          print *, 'bforce=',bforce,', fturb_x=', fturbx
        endif
        !fturby = -sum_f_eul(2) / (xl*yl*zl*hbar)
      else
        fturbx = 0.0_wp
        fturby = 0.0_wp
      endif
    endif
    
!print *, "nac_eul, 4", mtype, myid
    if (myid .eq. 0) then
      print *, 'Sum of eulerian force:', sum_f_eul(1:3)
      print *, 'Sum of lagrangian force:', sum_f_lag(1:3)
    endif
    
  end subroutine Calc_F_eul

  !> distribute the force on the turbines to background grid
  subroutine Calc_F_eul_old1(wtforce, wtforce_y, wtforce_z, fturbx, fturby, &
     level, ibm, fsi, NumberOfObjects, dh, mtype)
    use spectral
    use utils
    implicit none
    type(IBMNodes), dimension(:) :: ibm
    type(FSInfo), dimension(:) :: fsi
    integer :: NumberOfObjects, mtype
    real(wp) :: dh
    
    !real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), xsz(3))
    integer :: level
    real(wp) :: fturbx, fturby

    integer :: ioutputforce
    real(wp), allocatable, dimension(:,:) :: wtfx, wtfy, wtfz, ftemp2
    !real(wp), allocatable, dimension(:,:) :: ux, uy, uz 
    integer, allocatable, dimension(:,:) :: pid, pidtemp
    real(wp), allocatable, dimension(:,:,:) :: ftemp
    character(len=128) :: filen

    real(wp), dimension(3,3) :: dxdq, dqdx
    real(wp) :: dhx_, dhy_, dhz_, dh1, dh2, dh3, dV

    integer :: i, j, k, l, ibi, ii, jj, kk
    integer :: xs, ys, zs, xe, ye, ze
    real(wp) :: xi, yi, zi, xk, yk, zk 
    real(wp) :: xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr
    real(wp) :: dfunci, dfunck, df1, df2, df3, dh0, dA
    real(wp) :: r1i, r2i, r3i, r1k, r2k, r3k, rxi, ryi, rzi, rxk, ryk, rzk
    real(wp) :: f_temp, f_max, dh_max, dh_min
    real(wp), dimension(3) :: sum_f_eul, sum_f_lag, tmparray
    
    dh0 = sqrt(xl/nx*dh_mean)
    dh1 = xl/nx
    dh2 = yl/ny
    dh3 = zl/(nz-1)
         
    xs = xst(1); xe = xst(1) + xsz(1) - 1
    ys = xst(2); ye = xst(2) + xsz(2) - 1
    zs = xst(3); ze = xst(3) + xsz(3) - 1

    !print *, 'plyudebug, myid=',myid,', dh123=',dh1,dh2,dh3

    if (mtype .eq. 1) then
      ! only do initialization in 1st call from rotor model
      ! ignore this process in 2nd call from nacelle model
      !wtforce(:,:,(1-level):(xsz(3)+level)) = 0.0_wp
      !wtforce_y(:,:,(1-level):(xsz(3)+level)) = 0.0_wp
      !wtforce_z(:,:,(1-level):(xsz(3)+level)) = 0.0_wp
      wtforce(:,:,:) = 0.0_wp
      wtforce_y(:,:,:) = 0.0_wp
      wtforce_z(:,:,:) = 0.0_wp
    endif
    
    f_temp = 0.0_wp; f_max = 0.0_wp
    do ibi = 1, NumberOfObjects
    do l = 1, ibm(ibi)%n_elmt
      dA = ibm(ibi)%dA(l)
      do k = ibm(ibi)%k_min(l), ibm(ibi)%k_max(l)
      do j = ibm(ibi)%j_min(l), ibm(ibi)%j_max(l)
      do i = ibm(ibi)%i_min(l), ibm(ibi)%i_max(l)
        if (i>=xs .and. i<xe .and. j>=ys .and. j<=ye &
          .and. k>=zs .and. k<=ze) then 
          !> plyunote: xi,yi,zi should be overlapped with u, v nodes
          !> plyunote: xk,yk,zk should be overlapped with w nodes
          xi = cartx(i)
          yi = carty(j)
          !zi = cartz(i,j,k)
          zi = zzo(k) * (hbar + hho(i,j)) - hho(i,j)
          
          xk = cartx(i)
          yk = carty(j)
          !zk = cartzw(i,j,k)
          zk = zwo(k) * (hbar + hho(i,j)) - hho(i,j)
          
          xc = cartx(i)
          yc = carty(j)
          !zc = cartz(i,j,k)
          zc = zzo(k) * (hbar + hho(i,j)) - hho(i,j)

          rxi = xi - ibm(ibi)%cent_x(l);  ryi = yi - ibm(ibi)%cent_y(l)
          rzi = zi - ibm(ibi)%cent_z(l)
          rxk = xk - ibm(ibi)%cent_x(l);  ryk = yk - ibm(ibi)%cent_y(l)
          rzk = zk - ibm(ibi)%cent_z(l)

          !> plyunote: this section is implemented in a different form from
          !> fotis's version
          dxdq(1,1) = 1.; dxdq(1,2) = 0.; dxdq(1,3) = 0.
          dxdq(2,1) = 0.; dxdq(2,2) = 1.; dxdq(2,3) = 0.
          dxdq(3,1) = hxo(i,j)*(-1.0+zzo(k)); dxdq(3,2) = hyo(i,j)*(-1.0+zzo(k))
          dxdq(3,3) = hbar + eo(i,j)

          dhx_ = sqrt(dxdq(1,1)**2+dxdq(2,1)**2+dxdq(3,1)**2)
          dhy_ = sqrt(dxdq(1,2)**2+dxdq(2,2)**2+dxdq(3,2)**2)
          dhz_ = sqrt(dxdq(1,3)**2+dxdq(2,3)**2+dxdq(3,3)**2)
          
          dh_max = max(dhx_, dhy_, dhz_)
          dh_min = min(dhx_, dhy_, dhz_)
          if(max(dabs(dh_max),dabs(dh_min))<1.0d-6) then
            print *, 'plyudebug, myid=',myid,',i,j,k=',i,j,k,',dhxyz=',dhx_,dhy_,dhz_
          endif
          
          dh3 = dz(k-xst(3)+1) 
          dV = dh1 * dh2 * dh3 / her(i-xst(1)+1,j-xst(2)+1)
          !rx = xc - cent_x; ry = yc - cent_y; rz = zc - cent_z
          r1i = (rxi*dxdq(1,1) + ryi*dxdq(2,1) + rzi*dxdq(3,1))/dhx_/dhx_/dh1
          r2i = (rxi*dxdq(1,2) + ryi*dxdq(2,2) + rzi*dxdq(3,2))/dhy_/dhy_/dh2
          r3i = (rxi*dxdq(1,3) + ryi*dxdq(2,3) + rzi*dxdq(3,3))/dhz_/dhz_/dh3
          
          r1k = (rxk*dxdq(1,1) + ryk*dxdq(2,1) + rzk*dxdq(3,1))/dhx_/dhx_/dh1
          r2k = (rxk*dxdq(1,2) + ryk*dxdq(2,2) + rzk*dxdq(3,2))/dhy_/dhy_/dh2
          r3k = (rxk*dxdq(1,3) + ryk*dxdq(2,3) + rzk*dxdq(3,3))/dhz_/dhz_/dh3

          !> delta function: use dfunc_s4h in Yang(JCP,2009)
          !! doi:10.1016/j.jcp.2009.07.023
          if (deltafunc_F .eq. 1) then
            call dfunc_s1h(r1i, df1); call dfunc_s1h(r2i,df2); call dfunc_s1h(r3i,df3)
            dfunci = df1 * df2 * df3
            call dfunc_s1h(r1k, df1); call dfunc_s1h(r2k,df2); call dfunc_s1h(r3k,df3)
            dfunck = df1 * df2 * df3
          else if (deltafunc_F .eq. 2) then
            call dfunc_s2h(r1i, df1); call dfunc_s2h(r2i,df2); call dfunc_s2h(r3i,df3)
            dfunci = df1 * df2 * df3
            call dfunc_s2h(r1k, df1); call dfunc_s2h(r2k,df2); call dfunc_s2h(r3k,df3)
            dfunck = df1 * df2 * df3
          else if (deltafunc_F .eq. 3) then
            call dfunc_s3h(r1i, df1); call dfunc_s3h(r2i,df2); call dfunc_s3h(r3i,df3)
            dfunci = df1 * df2 * df3
            call dfunc_s3h(r1k, df1); call dfunc_s3h(r2k,df2); call dfunc_s3h(r3k,df3)
            dfunck = df1 * df2 * df3
          else if (deltafunc_F .eq. 4) then
            call dfunc_s4h(r1i, df1); call dfunc_s4h(r2i,df2); call dfunc_s4h(r3i,df3)
            dfunci = df1 * df2 * df3
            call dfunc_s4h(r1k, df1); call dfunc_s4h(r2k,df2); call dfunc_s4h(r3k,df3)
            dfunck = df1 * df2 * df3
          else
            print *, 'deltafunc_F: should be 2, 3, or 4'
          endif         
          
          !dfunc = df1 * df2 * df3
          
          ii = i - xst(1) + 1
          jj = j - xst(2) + 1
          kk = k - xst(3) + 1
          !> plyunote: please check U_lagr_x initialization
          wtforce(ii,jj,kk) = wtforce(ii,jj,kk) + ibm(ibi)%F_lagr_x(l)*dA*dfunci/dV
          wtforce_y(ii,jj,kk) = wtforce_y(ii,jj,kk) + ibm(ibi)%F_lagr_y(l)*dA*dfunci/dV
          wtforce_z(ii,jj,kk) = wtforce_z(ii,jj,kk) + ibm(ibi)%F_lagr_z(l)*dA*dfunck/dV

          !if (abs(dfunck)>1e-3 .and. (kk.eq.1 .or. kk.eq.2 &
          !  .or. kk.eq.xsz(3) .or. kk.eq.(xsz(3)-1)) .and. i.eq.52) then
          !  print *, myid, l, i, j, k, ii, jj, kk, ibm(ibi)%k_min(l),&
          !    ibm(ibi)%k_max(l), &
          !    wtforce(ii,jj,kk)*(xl/nx*yl/ny*dz(kk)/her(ii,jj))
          !endif
          
          !> plyunote: I think it needn't any mpi communication here.

          f_temp = sqrt(wtforce(ii,jj,kk)**2+wtforce_y(ii,jj,kk)**2&
            +wtforce_z(ii,jj,kk)**2)
          if(f_max < f_temp) f_max = f_temp
        endif
      enddo
      enddo
      enddo
    enddo
    enddo

    if (nacelle_model .eq. 0 .or. (nacelle_model .ne.0 .and. mtype .ne. 1)) then
      call dealiasxy(wtforce(:,:,1:xsz(3)))
      call dealiasxy(wtforce_y(:,:,1:xsz(3)))
      call dealiasxy(wtforce_z(:,:,1:xsz(3)))
    endif

    ! export wtforce to tecplot file
    ioutputforce = 0
    if(ioutputforce .ne. 0 .and. (nacelle_model .eq. 0 .or. &
      (nacelle_model .ne.0 .and. mtype .ne. 1)) .and. mod(ti, 20).eq.1) then
      !print *, 'temp0, allocate'
      allocate(wtfx(ny, nz), wtfy(ny, nz), wtfz(ny, nz))
      !allocate(ux(ny, nz), uy(ny, nz), uz(ny, nz))
      allocate(pid(ny, nz))
       
      !print *, 'temp1, gather'
      
      i=52
      ii = i - xst(1) + 1
      allocate(ftemp2(ny, nz), pidtemp(ny, nz))
      ftemp2 = 0.0_wp
      pidtemp = 0

      !print *, 'temp1.1, wtfx'
      ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
        wtforce(ii,:,1:xsz(3))
      call MPI_Allreduce(ftemp2, wtfx, ny*nz, mpi_double_precision, &
        mpi_sum, mpi_comm_2d_cart,ierr_wt)
      
      !print *, 'temp1.2, wtfy'
      ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
        wtforce_y(ii,:,1:xsz(3))
      call MPI_Allreduce(ftemp2, wtfy, ny*nz, mpi_double_precision, &
        mpi_sum, mpi_comm_2d_cart,ierr_wt)
      
      !print *, 'temp1.3, wtfz'
      ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
        wtforce_z(ii,:,1:xsz(3))
      call MPI_Allreduce(ftemp2, wtfz, ny*nz, mpi_double_precision, &
        mpi_sum, mpi_comm_2d_cart,ierr_wt)
      
      !ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
      !  u(ii,:,1:xsz(3))
      !call MPI_Allreduce(ftemp2, ux, ny*nz, mpi_double_precision, &
      !  mpi_sum, mpi_comm_2d_cart,ierr_wt)
      !
      !ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
      !  v(ii,:,1:xsz(3))
      !call MPI_Allreduce(ftemp2, uy, ny*nz, mpi_double_precision, &
      !  mpi_sum, mpi_comm_2d_cart,ierr_wt)
      !
      !ftemp2(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
      !  w(ii,:,1:xsz(3))
      !call MPI_Allreduce(ftemp2, uz, ny*nz, mpi_double_precision, &
      !  mpi_sum, mpi_comm_2d_cart,ierr_wt)
      
      !print *, 'temp1.4, pid'
      pidtemp(xst(2):(xst(2)+xsz(2)-1), xst(3):(xst(3)+xsz(3)-1)) = &
        myid
      call MPI_Allreduce(pidtemp, pid, ny*nz, mpi_integer, &
        mpi_sum, mpi_comm_2d_cart,ierr_wt)



      
      !allocate(ftemp(nx_global, ny_global, nz_global))
      !call gather_3d_xyz(wtforce(:,:,1:xsz(3)), ftemp, 0)
      !wtfx = ftemp(1, 1:ny, 1:nz)
      !call gather_3d_xyz(wtforce_y(:,:,1:xsz(3)), ftemp, 0)
      !wtfy = ftemp(1, 1:ny, 1:nz)
      !call gather_3d_xyz(wtforce_z(:,:,1:xsz(3)), ftemp, 0)
      !wtfz = ftemp(1, 1:ny, 1:nz)
      
      if (myid.eq.0) then
        write(filen, '(a,i0.10,a)') 'wtforce_', ti, '.dat'
        !print *, 'temp2, prepare filename ', filen
        open(18008, file=filen, status='replace')
        write(18008, *) "VARIABLES = X, Y, Z, FX, FY, FZ, PID"
        write(18008, *) 'ZONE T="', ti*dt, '" I=1 J=', ny, ' K=',nz,&
          ' SOLUTIONTIME=', ti*dt
        !i = 52
        do k = 1, nz
          do j = 1, ny
            !do i = 52, 52
            zc = zzo(k) * (hbar + hho(i,j)) - hho(i,j)
            write(18008, *) cartx(i), carty(j), zc, &
              wtfx(j,k), wtfy(j,k), wtfz(j,k), pid(j,k)
            !enddo
          enddo
        enddo
        close(18008)
      endif

      deallocate(wtfx, wtfy, wtfz, pid)
      deallocate(ftemp2, pidtemp)
      !print *, 'temp3, end'
    endif


    
    !print*,'plyudebug, ti=',ti,', myid=', myid, ', Calc_F_eul, f_max=', f_max

    !> check the sum of wtforce
    sum_f_eul(1:3) = 0.0_wp; sum_f_lag(1:3) = 0.0_wp; tmparray(1:3) = 0.0_wp 
    do i = 1, xsz(1)
    do j = 1, xsz(2)
    do k = 1, xsz(3)
      sum_f_eul(1) = sum_f_eul(1) + wtforce(i,j,k)*(xl/nx*yl/ny*dz(k)/her(i,j))
      sum_f_eul(2) = sum_f_eul(2) + wtforce_y(i,j,k)*(xl/nx*yl/ny*dz(k)/her(i,j))
      sum_f_eul(3) = sum_f_eul(3) + wtforce_z(i,j,k)*(xl/nx*yl/ny*dz(k)/her(i,j))

      !if (abs(wtforce(i,j,k))>1e-6 .and. (k.eq.1 .or. k.eq.2 .or. k.eq.xsz(3) &
      !  .or. k.eq.(xsz(3)-1)) ) then
      !  print *, i, j, k, xst(1)+i-1, xst(2)+j-1, xst(3)+k-1, &
      !    wtforce(i,j,k)*(xl/nx*yl/ny*dz(k)/her(i,j))
      !endif

    enddo
    enddo
    enddo

    call mpi_allreduce(sum_f_eul,tmparray,3,mpi_double_precision,mpi_sum,&
      mpi_comm_2d_cart,ierr_wt)

    sum_f_eul(1:3) = tmparray(1:3)

    do ibi = 1, NumberOfObjects
      do l = 1, ibm(ibi)%n_elmt
        dA = ibm(ibi)%dA(l)
        sum_f_lag(1) = sum_f_lag(1) + ibm(ibi)%F_lagr_x(l)*dA
        sum_f_lag(2) = sum_f_lag(2) + ibm(ibi)%F_lagr_y(l)*dA
        sum_f_lag(3) = sum_f_lag(3) + ibm(ibi)%F_lagr_z(l)*dA
      enddo
    enddo

    !> plyunote: need zl*hbar?
    if (mtype .eq. 1) then
      if (InletRelaxation .eq. 0) then
        !fturbx = -sum_f_eul(1) / (xl*yl*zl*hbar)
        fturbx = 0.0_wp
        do ibi = 1, NumberOfObjects
          !fturbx = (fturbx + 4.0_wp/3.0_wp * (ibm(ibi)%U_ref)**2 * twopi/2.0_wp *
          !  (fsi(ibi)%r_rotor) **2) / (xl*yl*zl*hbar)
          
          !> c_bforce is recommended to be thrust coefficient. likely to be 3/4
          
          fturbx = fturbx + (c_bforce * (fsi(ibi)%angvel_fixed * fsi(ibi)%r_rotor &
            / fsi(ibi)%Tipspeedratio)**2 * twopi/2.0_wp * (fsi(ibi)%r_rotor) **2) &
            / (xl*yl*zl*hbar)
        enddo
        if(myid.eq.0 .and. ti .eq. (ti_first+1)) then
          print *, 'bforce=',bforce,', fturb_x=', fturbx
        endif
        !fturby = -sum_f_eul(2) / (xl*yl*zl*hbar)
      else
        fturbx = 0.0_wp
        fturby = 0.0_wp
      endif
    endif
    

    if (myid .eq. 0) then
      print *, 'Sum of eulerian force:', sum_f_eul(1:3)
      print *, 'Sum of lagrangian force:', sum_f_lag(1:3)
    endif
    
  end subroutine Calc_F_eul_old1

  subroutine ForceProjection_l2s(ibm_surface, ibm_line, fsi_surface, fsi_line,&
    NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm_surface, ibm_line
    type(FSInfo), dimension(:) :: fsi_surface, fsi_line
    integer :: NumberOfObjects

    integer :: ibi, elmt_s, elmt_l, l

    integer :: n_elmt_line
    real(wp), allocatable, dimension(:,:) :: force_l, force_s
    real(wp), allocatable, dimension(:) :: area_l, area_s
    real(wp), dimension(3) :: sum_force_l, sum_force_s 
    real(wp) :: sum_area_l, sum_area_s

    integer :: n1e, n2e, n3e
    integer :: db

    integer :: i1, i2
    real(wp) :: c, f1, f2, c1, c2

    do ibi = 1, NumberOfObjects
      do elmt_s = 1, ibm_surface(ibi)%n_elmt
        elmt_l = ibm_surface(ibi)%s2l(elmt_s)
        i1 = ibm_surface(ibi)%s2l1(elmt_s)
        i2 = ibm_surface(ibi)%s2l2(elmt_s)
        c = ibm_surface(ibi)%s2lc(elmt_s)
        
        if (ibm_line(ibi)%chord_blade(i1)<1.0e-9_wp .or. &
          ibm_line(ibi)%chord_blade(i2)<1.0e-9_wp ) then
          print *,'plyudebug, ForceProjection_l2s, chord length too small',&
            elmt_l
        endif

        c1 = 1.0_wp / ibm_line(ibi)%chord_blade(i1) * (1-c)
        c2 = 1.0_wp / ibm_line(ibi)%chord_blade(i2) * c
        
        !> plyunote: linear interpolation is used here. VFS not.
        ibm_surface(ibi)%F_lagr_x(elmt_s) = ibm_line(ibi)%F_lagr_x(i1) &
          * c1 + ibm_line(ibi)%F_lagr_x(i2) * c2
        ibm_surface(ibi)%F_lagr_y(elmt_s) = ibm_line(ibi)%F_lagr_y(i1) &
          * c1 + ibm_line(ibi)%F_lagr_y(i2) * c2
        ibm_surface(ibi)%F_lagr_z(elmt_s) = ibm_line(ibi)%F_lagr_z(i1) &
          * c1 + ibm_line(ibi)%F_lagr_z(i2) * c2
      enddo
    enddo

    !> plyunote: check the force distribution
    if(myid .eq. 0) then
    do ibi = 1, 1 
      n_elmt_line = ibm_line(ibi)%n_elmt
      allocate(force_l(n_elmt_line,3), force_s(n_elmt_line,3))
      allocate(area_l(n_elmt_line), area_s(n_elmt_line))
      force_l(:,:)=0.0_wp; force_s(:,:)=0.0_wp
      area_l(:)=0.0_wp; area_s(:)=0.0_wp
      do elmt_s = 1, ibm_surface(ibi)%n_elmt
        elmt_l = ibm_surface(ibi)%s2l(elmt_s)

        db = 0 
        if (db .gt. 0) then
        if (elmt_l.eq.10) then
          n1e=ibm_surface(ibi)%nv1(elmt_s)
          n2e=ibm_surface(ibi)%nv2(elmt_s)
          n3e=ibm_surface(ibi)%nv3(elmt_s)

          print *,'plyudebug, s2l_s: p1, p2, p3, dA, F'
          print*, ibm_surface(ibi)%x_bp(n1e), &
            ibm_surface(ibi)%y_bp(n1e),ibm_surface(ibi)%z_bp(n1e)
          print*, ibm_surface(ibi)%x_bp(n2e), &
            ibm_surface(ibi)%y_bp(n2e),ibm_surface(ibi)%z_bp(n2e)
          print*, ibm_surface(ibi)%x_bp(n3e), &
            ibm_surface(ibi)%y_bp(n3e),ibm_surface(ibi)%z_bp(n3e)
          print*, ibm_surface(ibi)%dA(elmt_s)
          print*, ibm_surface(ibi)%F_lagr_x(elmt_s),&
            ibm_surface(ibi)%F_lagr_y(elmt_s),&
            ibm_surface(ibi)%F_lagr_z(elmt_s)
        endif
        endif

        force_s(elmt_l,1) = force_s(elmt_l,1) + ibm_surface(ibi)%F_lagr_x(elmt_s)&
          *ibm_surface(ibi)%dA(elmt_s)
        force_s(elmt_l,2) = force_s(elmt_l,2) + ibm_surface(ibi)%F_lagr_y(elmt_s)&
          *ibm_surface(ibi)%dA(elmt_s)
        force_s(elmt_l,3) = force_s(elmt_l,3) + ibm_surface(ibi)%F_lagr_z(elmt_s)&
          *ibm_surface(ibi)%dA(elmt_s)
        area_s(elmt_l) = area_s(elmt_l) + ibm_surface(ibi)%dA(elmt_s)
      enddo

      if (db .gt. 0) then
      n1e=ibm_line(ibi)%nv1(10); n2e=ibm_line(ibi)%nv2(10);
      print *, 'plyudebug, s2l_l: p1, p2, dA, c, F'
      print *, ibm_line(ibi)%x_bp(n1e),&
        ibm_line(ibi)%y_bp(n1e), ibm_line(ibi)%z_bp(n1e)
      print *, ibm_line(ibi)%x_bp(n2e),&
        ibm_line(ibi)%y_bp(n2e), ibm_line(ibi)%z_bp(n2e)
      print *, ibm_line(ibi)%dA(10), ibm_line(ibi)%chord_blade(10)
      print *, ibm_line(ibi)%F_lagr_x(10),&
        ibm_line(ibi)%F_lagr_y(10),ibm_line(ibi)%F_lagr_z(10)
      endif

      do elmt_l = 1, n_elmt_line
        force_l(elmt_l,1) = ibm_line(ibi)%F_lagr_x(elmt_l)*ibm_line(ibi)%dA(elmt_l)
        force_l(elmt_l,2) = ibm_line(ibi)%F_lagr_y(elmt_l)*ibm_line(ibi)%dA(elmt_l)
        force_l(elmt_l,3) = ibm_line(ibi)%F_lagr_z(elmt_l)*ibm_line(ibi)%dA(elmt_l)
        area_l(elmt_l) = ibm_line(ibi)%dA(elmt_l) &
          *ibm_line(ibi)%chord_blade(elmt_l)
      enddo
      sum_force_l(:)=0.0_wp; sum_force_s(:)=0.0_wp
      sum_area_l=0.0_wp; sum_area_s=0.0_wp
      do elmt_l = 1, n_elmt_line
        sum_force_l(1:3)=sum_force_l(1:3)+force_l(elmt_l,1:3)
        sum_force_s(1:3)=sum_force_s(1:3)+force_s(elmt_l,1:3)
        sum_area_l=sum_area_l+area_l(elmt_l)
        sum_area_s=sum_area_s+area_s(elmt_l)
      enddo     
      
      if (db .gt. 0) then
        print *,'plyudebug, s2linfo (ielemt area_line area_surface F_line,F_surface)'
        do elmt_l = 1, n_elmt_line
          print *, elmt_l, area_l(elmt_l), area_s(elmt_l), &
            norm2(force_l(elmt_l,1:3)), norm2(force_s(elmt_l,1:3))
        enddo
      endif
      
      print*,'Sum_force of line and surface:',sum_force_l(1:3),', ',sum_force_s(1:3)
      print*,'Sum_area of line and surface:', sum_area_l, sum_area_s

      deallocate(force_l, force_s)
      deallocate(area_l, area_s)
    enddo
    endif
  end subroutine ForceProjection_l2s
  
  subroutine ForceProjection_l2d(ibm_surface, ibm_line, fsi_surface, fsi_line,&
    NumberOfObjects)
    implicit none
    type(IBMNodes), dimension(:) :: ibm_surface, ibm_line
    type(FSInfo), dimension(:) :: fsi_surface, fsi_line
    integer :: NumberOfObjects

    integer :: ibi, elmt_s, elmt_l, l

    integer :: n_elmt_line, n1
    !real(wp), allocatable, dimension(:,:) :: force_l, force_s
    real(wp), allocatable, dimension(:) :: area_l, area_s
    real(wp) :: sum_ft_l, sum_ft_s, sum_fa_l, sum_fa_s 
    real(wp) :: sum_area_l, sum_area_s

    integer :: n1e, n2e, n3e
    integer :: db

    integer :: i1, i2, iblade
    real(wp) :: c, f1, f2, c1, c2, tmp

    real(wp), dimension(:), allocatable :: Ft, Fa
    real(wp), dimension(3) :: r_vec, n_blade, n_axis, n_rot, Flagr
    real(wp) :: r_amp, Ft1, Fa1

    do ibi = 1, NumberOfObjects
      !> For ADM-R, firstly decompose force of AL elements to thrust and
      !tangential direction
      n1 = ibm_line(ibi)%n_elmt / num_blade
      n_axis(1)=fsi_line(ibi)%nx_tb; n_axis(2)=fsi_line(ibi)%ny_tb
      n_axis(3)=fsi_line(ibi)%nz_tb

      allocate(Ft(n1), Fa(n1))
      Ft(:)=0.0_wp; Fa(:)=0.0_wp
      do elmt_l = 1, n1
        
        do iblade = 1, num_blade
          i1 = elmt_l + (iblade -1)*n1
          r_vec(1) = ibm_line(ibi)%cent_x(i1) - fsi_line(ibi)%x_c
          r_vec(2) = ibm_line(ibi)%cent_y(i1) - fsi_line(ibi)%y_c
          r_vec(3) = ibm_line(ibi)%cent_z(i1) - fsi_line(ibi)%z_c
          r_amp = norm2(r_vec) + 1.0e-16_wp
          n_blade = r_vec / r_amp

          if (ibm_line(ibi)%chord_blade(i1) < 1.0e-9_wp ) then
            print *, 'plyudebug, ForceProjection_l2d, chord length too small'
          endif

          Flagr(1) = ibm_line(ibi)%F_lagr_x(i1)
          Flagr(2) = ibm_line(ibi)%F_lagr_y(i1)
          Flagr(3) = ibm_line(ibi)%F_lagr_z(i1)
          Flagr = Flagr / (twopi*r_amp)

          call crossx(n_axis, n_blade, n_rot)

          Ft(elmt_l) = Ft(elmt_l) + dot_product(Flagr, n_axis) 
          Fa(elmt_l) = Fa(elmt_l) + dot_product(Flagr, n_rot)  
          
          !if(myid.eq.0) then
          !  print *, 'plyudebug, l2d, i1=',i1,',F_lagr=',Flagr
          !  print *, 'n_blade=',n_blade,',n_rot=',n_rot,',Ft,Fa=',Ft(elmt_l),Fa(elmt_l)       
          !endif
        enddo
      enddo

      do elmt_s = 1, ibm_surface(ibi)%n_elmt
        elmt_l = ibm_surface(ibi)%s2l(elmt_s)
        i1 = ibm_surface(ibi)%s2l1(elmt_s)
        i2 = ibm_surface(ibi)%s2l2(elmt_s)
        c = ibm_surface(ibi)%s2lc(elmt_s)
        
        if (ibm_line(ibi)%chord_blade(i1)<1.0e-9_wp .or. &
          ibm_line(ibi)%chord_blade(i2)<1.0e-9_wp ) then
          print *,'plyudebug, ForceProjection_l2d, chord length too small',&
            elmt_l
        endif

        !c1 = 1.0_wp / ibm_line(ibi)%chord_blade(i1) * (1-c)
        !c2 = 1.0_wp / ibm_line(ibi)%chord_blade(i2) * c

        Ft1 = Ft(i1)*(1.0_wp-c)+Ft(i2)*c; Fa1 = Fa(i1)*(1.0_wp-c)+Fa(i2)*c
        
        r_vec(1) = ibm_surface(ibi)%cent_x(elmt_s) - fsi_surface(ibi)%x_c
        r_vec(2) = ibm_surface(ibi)%cent_y(elmt_s) - fsi_surface(ibi)%y_c
        r_vec(3) = ibm_surface(ibi)%cent_z(elmt_s) - fsi_surface(ibi)%z_c
        r_amp = norm2(r_vec) + 1.0e-16_wp
        n_blade = r_vec / r_amp

        call crossx(n_axis, n_blade, n_rot)

        Flagr = Ft1 * n_axis + Fa1 * n_rot

        !if (myid.eq.0 .and. elmt_s .eq. 100) then
        !  print *,'i1,i2,c=',i1,i2,c,',Ft1,Fa1=',Ft1,Fa1
        !  print *,'n_blade=',n_blade,', n_rot=',n_rot,', Flagr=',Flagr
        !endif

        ibm_surface(ibi)%F_lagr_x(elmt_s) = Flagr(1) 
        ibm_surface(ibi)%F_lagr_y(elmt_s) = Flagr(2) 
        ibm_surface(ibi)%F_lagr_z(elmt_s) = Flagr(3) 
        !ibm_surface(ibi)%F_lagr_x(elmt_s) = ibm_line(ibi)%F_lagr_x(i1) &
        !  * c1 + ibm_line(ibi)%F_lagr_x(i2) * c2
        !ibm_surface(ibi)%F_lagr_y(elmt_s) = ibm_line(ibi)%F_lagr_y(i1) &
        !  * c1 + ibm_line(ibi)%F_lagr_y(i2) * c2
        !ibm_surface(ibi)%F_lagr_z(elmt_s) = ibm_line(ibi)%F_lagr_z(i1) &
        !  * c1 + ibm_line(ibi)%F_lagr_z(i2) * c2
      enddo
    enddo
  end subroutine ForceProjection_l2d
  
  !subroutine Add_Inlet_Forcing(u, v, w, wtforce, wtforce_y, wtforce_z, level)
  !  use param
  !  use decomp
  !  implicit none
  !  real(wp), dimension(:,:,:) :: u(xsz(1), xsz(2), 1-level:xsz(3)+level)
  !  real(wp), dimension(:,:,:) :: v(xsz(1), xsz(2), 1-level:xsz(3)+level)
  !  real(wp), dimension(:,:,:) :: w(xsz(1), xsz(2), 1-level:xsz(3)+level)
  !  !real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), 1-level:xsz(3)+level)
  !  !real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), 1-level:xsz(3)+level)
  !  !real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), 1-level:xsz(3)+level)
  !  real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), xsz(3))
  !  real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), xsz(3))
  !  real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), xsz(3))
  !  integer :: level

  !  !> u_tgt: target velocity
  !  real(wp) :: r, xi, f_max, utgt_fmax, u_fmax, urelax_fmax, ir_vol, ir_turb_vol
  !  real(wp) :: ir_rms, dvol, dturb, ir_c_rms, ir_c1, ir_usum
  !  real(wp), allocatable, dimension(:) :: rvec, su
  !  real(wp) :: u_relax(3), temp_out(3), temp_in(3)
  !  !real(wp), allocatable, dimension(:,:) :: f_relax, f_temp
  !  real(wp), allocatable, dimension(:,:,:,:) :: du_relax, du_temp
  !  real(wp) :: zlbot, z, zs, u1, u2
  !  integer :: i, j, k, j1, imin, imax, ismooth, i_fmax, j_fmax, k_fmax
  !  integer :: ir_num, DynamicIRCoeff
  !  character(len=64) :: filen

  !  f_max = 0.0_wp; i_fmax=1; j_fmax=1; k_fmax=1
  !  u_fmax=0.0_wp; utgt_fmax=0.0_wp; urelax_fmax=0.0_wp

  !  ir_num = 0; ir_vol = 0.0_wp; ir_turb_vol = 0.0_wp; ir_rms = 0.0_wp
  !  ir_usum = 0.0_wp
  ! 
  !  !nsmooth = 2  
  !  DynamicIRCoeff = 0
  !   
  !  allocate(su(1-level:xsz(3)+level))
  !  allocate(rvec(xsz(1)))
  !  allocate(du_relax(xsz(1),xsz(2),1-level:xsz(3)+level,3))
  !  allocate(du_temp(xsz(1),xsz(2),1-level:xsz(3)+level,3))
  !  
  !  du_relax(:,:,:,:) = 0.0_wp

  !  imin = floor(IR_start/(xl/nx))
  !  imax = floor(IR_end/(xl/nx))+1

  !  !print *, 'IR_1'    
  !  if (restop .le. 1.e-6) then
  !    zlbot = hbar * resbot/(resbot+restop)
  !    !zlsbot = 1.0_wp / resbot
  !    !zlstop = 0.0_wp
  !    
  !    !> usbot is imported from param module
  !    !usbot = 1.0_wp / (2.5_wp * log(hbar/z0))
  !    
  !    !ustop = 0.0_wp
  !    !print *, 'IR_2'

  !    !> Get su(xsz(3))
  !    !utop = usbot * (2.5_wp*log(hbar/z0))
  !    do k = 1-level, xsz(3)+level
  !      z = zz(k) * hbar
  !      ! for wave cases, z might less than 0, hard to use log law.
  !      if (z .le. zlbot .and. z>0.0) then
  !        zs = z/zlsbot
  !        su(k) = 0.0_wp
  !        u1 = usbot * zs
  !        u2 = usbot * (2.5_wp * log(z/z0))
  !        if (u1 .lt. u2) su(k) = u1
  !        if (u1 .ge. u2) su(k) = u2
  !        if (u2 .lt. 0.0_wp) su(k) = u1
  !      endif        
  !    enddo

  !    !> Get rvec(xsz(1))
  !    rvec(:) = 0.0_wp
  !    if (xst(1)<=imax .and. (xst(1)+xsz(1)-1)>=imin) then
  !      do i = 1, xsz(1)
  !        xi = (xst(1)+i-1)*(xl/nx)
  !        if (xi>=IR_start .and. xi <= IR_end) then
  !          rvec(i) = (IR_end - xi)/(IR_end - IR_start)
  !        endif
  !      enddo
  !    endif

  !    !> Get du_relax(xsz(1),xsz(2),1-level:xsz(3)+level,3)
  !    
  !    if (xst(1)<=imax .and. (xst(1)+xsz(1)-1)>=imin) then
  !      do k = 1-level, xsz(3)+level
  !        z = zz(k) * hbar
  !        if (isbot .and. (k==1 .or. k==0)) then

  !        elseif(istop .and. k==xsz(3)+level) then

  !        ! for wave cases, z might less than 0, hard to use log law.
  !        elseif(z>0.0_wp) then
  !          !print *, 'myid, k_local, k_global, dz, z=',myid,k,xst(3)+k-1,dz(k),z
  !          do j = 1, xsz(2)
  !            do i = 1, xsz(1)
  !              xi = (xst(1)+i-1)*(xl/nx)
  !              if (xi>=IR_start .and. xi <= IR_end) then
  !                u_relax(1) = (1.0_wp - rvec(i))*su(k)+rvec(i)*u(i,j,k)
  !                u_relax(2) = rvec(i)*v(i,j,k)
  !                u_relax(3) = rvec(i)*w(i,j,k)
  !                du_relax(i,j,k,1) = u_relax(1) - u(i,j,k)
  !                du_relax(i,j,k,2) = u_relax(2) - v(i,j,k)
  !                du_relax(i,j,k,3) = u_relax(3) - w(i,j,k)
  !                    
  !                !if (norm2(du_relax(i,j,k,:))>0.2) then
  !                !  print *, 'Error in du_relax: local i,j,k=',i,j,k,&
  !                !    ',global i,j,k=',xst(1)+i-1,xst(2)+j-1,xst(3)+k-1
  !                !  print *, 'rvec,su=',rvec(i),su(k),',u,v,w,u_relax,v_relax,w_relax=',&
  !                !    u(i,j,k), v(i,j,k), w(i,j,k), u_relax(1:3)                      
  !                !endif
  !                
  !                dturb = du_relax(i,j,k,1)**2+du_relax(i,j,k,2)**2+du_relax(i,j,k,3)**2
  !                dvol = xl/nx * yl/ny * dz(k)*hbar

  !                ir_num = ir_num+1
  !                ir_vol = ir_vol + dvol
  !                ir_usum = ir_usum + dvol * su(k)
  !                ir_turb_vol = ir_turb_vol + dturb * dvol 
  !                
  !                !if ((xst(1)+i-1).eq. floor((IR_start+IR_end)/2.0/(xl/nx))) then
  !                !  j1 = floor(ny/2.0) - xst(2) +1
  !                !  if(j1>=1 .and. j1<=xsz(2) .and. j1.eq.j) then
  !                !    print *, 'k,u,su,r,u_relax,du_relax=',&
  !                !      xst(3)+k-1, u(i,j1,k), su(k), rvec(i), u_relax(1), du_relax(i,j1,k,1)
  !                !  endif
  !                !endif
  !                             
  !              endif
  !            enddo
  !          enddo

  !        endif
  !      enddo
  !        
  !      do ismooth = 1, IR_NSmooth          
  !        if (isbot) then
  !          du_relax(:,:,1:4,:) = 0.0_wp
  !        endif

  !        if(istop) then
  !          du_relax(:,:,xsz(3)-3:xsz(3),:) = 0.0_wp
  !        endif

  !        du_temp(:,:,:,:) = 0.0_wp

  !        do j = 2, xsz(2)-1
  !          du_temp(:,j,:,:)=(du_relax(:,j-1,:,:)+2.0_wp*du_relax(:,j,:,:) &
  !            +du_relax(:,j+1,:,:))/4.0_wp
  !          du_temp(:,1,:,:)=(du_relax(:,1,:,:)+du_relax(:,2,:,:))/2.0_wp
  !          du_temp(:,xsz(2),:,:)=(du_relax(:,xsz(2)-1,:,:)+du_relax(:,xsz(2),:,:))/2.0_wp
  !        enddo
  !        du_relax(:,:,:,:) = du_temp(:,:,:,:)

  !        do k = 1, xsz(3)
  !          ! level = 1
  !          du_temp(:,:,k,:)=(du_relax(:,:,k-1,:)+2.0_wp*du_relax(:,:,k,:) &
  !            +du_relax(:,:,k+1,:))/4.0_wp
  !        enddo
  !        
  !        call update_ghost(du_temp(:,:,:,1),level)
  !        call update_ghost(du_temp(:,:,:,2),level)
  !        call update_ghost(du_temp(:,:,:,3),level)
  !        du_relax(:,:,:,:) = du_temp(:,:,:,:)
  !      enddo

  !      !print *, 'myid=',myid,', ir_num, ir_vol, ir_usum, ir_turb_vol=',&
  !      !  ir_num,ir_vol, ir_usum, ir_turb_vol

  !      
  !      if (DynamicIRCoeff .ne. 0) then
  !        temp_out(:) = (/ir_vol, ir_usum, ir_turb_vol/)
  !        temp_in(:) = 0.0_wp

  !        call mpi_allreduce(temp_out,temp_in,3,mpi_double_precision,&
  !          mpi_sum, mpi_comm_world, ierr_wt)

  !        ir_vol = temp_in(1); ir_usum = temp_in(2); ir_turb_vol = temp_in(3)

  !        ir_rms = sqrt(ir_turb_vol / (ir_vol+1.0e-16_wp)) &
  !          / (ir_usum/(ir_vol+1.0e-16_wp))
  !        ir_c_rms = 1.0_wp - arm / (ir_rms+1.0e-16_wp)
  !        ir_c1 = max(0.0_wp, ir_c_rms)
  !        ir_c1 = min(IR_cx, ir_c1)
  !        !print *, 'myid=',myid,'ir_rms, arm ,ir_c_rms, IR_c, ir_c1=',&
  !        !  ir_rms, arm, ir_c_rms, IR_c, ir_c1
  !      endif

  !      !ir_c1 = 0.05

  !      wtforce(:,:,1:xsz(3)) = wtforce(:,:,1:xsz(3)) &
  !        + du_relax(:,:,1:xsz(3),1) / dt * ir_cx
  !      wtforce_y(:,:,1:xsz(3)) = wtforce_y(:,:,1:xsz(3)) &
  !        + du_relax(:,:,1:xsz(3),2) / dt * ir_cy
  !      wtforce_z(:,:,1:xsz(3)) = wtforce_z(:,:,1:xsz(3)) &
  !        + du_relax(:,:,1:xsz(3),3) / dt * ir_cz

  !      !if (mod(ti,noutc).eq.0) then
  !      !  write(filen, '(a8,i0.3,a1,i0.6,a4)') 'wtforce_', myid,'_',ti,'.dat'
  !      !  open(18006, file=filen, action='write')
  !      !  
  !      !  write(18006,*) 'i_g, j_g, k_g, wtforce_x, wtforce_y, wtforce_z'
  !      !  do k = 1-level, xsz(3)+level
  !      !    do j = 1, xsz(2)
  !      !      do i = 1, xsz(1)
  !      !        write(18006,*) xst(1)+i-1, xst(2)+j-1, xst(3)+k-1,&
  !      !          wtforce(i,j,k), wtforce_y(i,j,k), wtforce_z(i,j,k)
  !      !      enddo
  !      !    enddo
  !      !  enddo
  !      !  close(18006)
  !      !endif
  !    endif
  !  else
  !    print *, "plyudebug: This case is not yet included in Add_Inlet_Forcing"
  !  endif
  !    
  !  deallocate(su, rvec)
  !  deallocate(du_relax, du_temp)

  !    !do j=1, xsz(2)
  !    !  if (abs(f_relax(j))>f_max) then
  !    !    f_max = abs(f_relax(j))
  !    !    i_fmax = xst(1)+i-1; j_fmax=xst(2)+j-1; k_fmax=xst(3)+j-1
  !    !    u_fmax = u(i,j,k); utgt_fmax=u_tgt(j); urelax_fmax=u_relax(j)
  !    !  endif
  !    !enddo
  !    !

  !    !print *, 'plyudebug, cpu ', myid, 'i, j, k=',i_fmax,j_fmax,k_fmax,&
  !
  !end subroutine Add_Inlet_Forcing

      
  subroutine rotor_model_acs(time_, u, v, w, wtforce, wtforce_y, wtforce_z,&
    fturbx, fturby, level, flag_de)
    implicit none
    real(wp) :: time_
    real(wp), dimension(:,:,:) :: u(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: v(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: w(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: flag_de(xsz(1), xsz(2), xsz(3))
    integer :: level
    real(wp) :: fturbx, fturby
    
    integer :: ibi
    real(wp) :: volflux(3)
 
    ti = nint(time_/dt)

    call update_zgrid 
    !call calc_vol_flux(u(:,:,1:xsz(3)), volflux)
    
    if(myid.eq.0) print *, 'plyudebug, acs_run 1: Uref_ACL'
    !> in Uref_ACL, Calc_U_lagr is called for ref_vel disk
    call Uref_ACL(u, v, w, ibm_acl2ref, NumberOfTurbines, level)

    if(myid.eq.0) print *, 'plyudebug, acs_run 2: turbineangvel'
    call Calc_turbineangvel(dt, ibm_acl2ref, fsi_acl2ref)

    if(myid.eq.0) print *, 'plyudebug, acs_run 3: copy Uref'
    do ibi = 1, NumberOfTurbines
      wtm(ibi)%U_ref = ibm_acl2ref(ibi)%U_ref
      fsi_wt(ibi)%angvel_axis = fsi_acl2ref(ibi)%angvel_axis
      fsi_wt(ibi)%ang_axis = fsi_acl2ref(ibi)%ang_axis
    enddo

    if(myid.eq.0) print *, 'plyudebug, acs_run 3.1: rotor_Rot(wtm)'
    do ibi = 1,NumberOfTurbines
      call rotor_Rot(time_, fsi_wt(ibi), wtm(ibi), dt, ibi, 2)
    enddo
    if(myid.eq.0) print *, 'plyudebug, acs_run 3.2: rotor_Rot(acl2ref)'
    do ibi = 1,NumberOfTurbines
      call rotor_Rot(time_, fsi_acl2ref(ibi), ibm_acl2ref(ibi), dt, ibi, 1)
    enddo
    !call Export_SurfaceLocation(wtm, NumberOfTurbines)
    if(myid.eq.0) print *, 'plyudebug, acs_run 4: Pre_process(wtm)'
    call Pre_process(wtm, NumberOfTurbines, 1)
    call Pre_process(ibm_acl2ref, NumberOfTurbines, 1)

    !> here Calc_U_lagr is called for ACL elements
    if(myid.eq.0) print *, 'plyudebug, acs_run 5: Calc_U_lagr'
    call Calc_U_lagr(u, v, w, ibm_acl2ref, NumberOfTurbines, level)

    if(myid.eq.0) print *, 'plyudebug, acs_run 6: Calc_F_lagr_ACL'
    call Calc_F_lagr_ACL(ibm_acl2ref, fsi_acl2ref, NumberOfTurbines)

    !> plyunote: force exportation moved to main 
    !call Export_ForceOnBlade(wtm, fsi_wt, NumberOfTurbines)
    
    call ForceProjection_l2s(wtm, ibm_acl2ref, fsi_wt, fsi_acl2ref,&
      NumberOfTurbines)

    if(myid.eq.0) print *, 'plyudebug, acs_run 7: Calc_forces_ACL'
    call Calc_forces_ACL(ibm_acl2ref, fsi_acl2ref, NumberOfTurbines,&
      'Turbine_AL_')
    call Calc_forces_ACL(wtm, fsi_wt, NumberOfTurbines,&
      'Turbine_AS_')

    if(myid.eq.0) print *, 'plyudebug, acs_run 8: Calc_F_eul'
    call Calc_F_eul(wtforce, wtforce_y, wtforce_z, fturbx, fturby, level, &
      wtm, fsi_wt, NumberOfTurbines, 1.0_wp, 1, flag_de)
    
    !if(InletRelaxation .eq. 2) call Add_Inlet_Forcing(u, v, w, wtforce, wtforce_y, wtforce_z, level)

    !> plyunote: wtforce initialization should be moved outside Calc_F_eul
    !> because both rotor_model_acs and nacelle_model_run will call Calc_F_eul

    if (nacelle_model .ne. 0) then
      call nacelle_model_run(time_, u, v, w, wtforce, wtforce_y, wtforce_z,&
        fturbx, fturby, level, flag_de)
    endif


  end subroutine rotor_model_acs
 
  subroutine rotor_model_admr(time_, u, v, w, wtforce, wtforce_y, wtforce_z,&
    fturbx, fturby, level, flag_de)
    implicit none
    real(wp) :: time_
    real(wp), dimension(:,:,:) :: u(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: v(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: w(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: flag_de(xsz(1), xsz(2), xsz(3))
    integer :: level
    real(wp) :: fturbx, fturby
    
    integer :: ibi
    real(wp) :: volflux(3)

    ti = nint(time_/dt)

    call update_zgrid 
    !call calc_vol_flux(u(:,:,1:xsz(3)), volflux)
    
    if(myid.eq.0) print *, 'plyudebug, admr_run 1: Uref_ACL'
    !> in Uref_ACL, Calc_U_lagr is called for ref_vel disk
    call Uref_ACL(u, v, w, ibm_acl2ref, NumberOfTurbines, level)

    if(myid.eq.0) print *, 'plyudebug, admr_run 2: turbineangvel'
    call Calc_turbineangvel(dt, ibm_acl2ref, fsi_acl2ref)

    if(myid.eq.0) print *, 'plyudebug, admr_run 3: copy Uref'
    do ibi = 1, NumberOfTurbines
      wtm(ibi)%U_ref = ibm_acl2ref(ibi)%U_ref
      fsi_wt(ibi)%angvel_axis = fsi_acl2ref(ibi)%angvel_axis
      fsi_wt(ibi)%ang_axis = fsi_acl2ref(ibi)%ang_axis
    enddo

    if(myid.eq.0) print *, 'plyudebug, admr_run 3.1: rotor_Rot(wtm)'
    do ibi = 1,NumberOfTurbines
      call rotor_Rot(time_, fsi_wt(ibi), wtm(ibi), dt, ibi, 2)
    enddo
    if(myid.eq.0) print *, 'plyudebug, admr_run 3.2: rotor_Rot(acl2ref)'
    do ibi = 1,NumberOfTurbines
      call rotor_Rot(time_, fsi_acl2ref(ibi), ibm_acl2ref(ibi), dt, ibi, 1)
    enddo

    if(myid.eq.0) print *, 'plyudebug, admr_run 4: Pre_process(wtm)'
    call Pre_process(wtm, NumberOfTurbines, 1)
    call Pre_process(ibm_acl2ref, NumberOfTurbines, 1)

    !> here Calc_U_lagr is called for ACL elements
    if(myid.eq.0) print *, 'plyudebug, admr_run 5: Calc_U_lagr'
    call Calc_U_lagr(u, v, w, ibm_acl2ref, NumberOfTurbines, level)

    if(myid.eq.0) print *, 'plyudebug, admr_run 6: Calc_F_lagr_ACL'
    call Calc_F_lagr_ACL(ibm_acl2ref, fsi_acl2ref, NumberOfTurbines)

    !> plyunote: force exportation moved to main 
    !call Export_ForceOnBlade(wtm, fsi_wt, NumberOfTurbines)
    
    call ForceProjection_l2d(wtm, ibm_acl2ref, fsi_wt, fsi_acl2ref,&
      NumberOfTurbines)

    if(myid.eq.0) print *, 'plyudebug, admr_run 7: Calc_forces_ACL'
    call Calc_forces_ACL(ibm_acl2ref, fsi_acl2ref, NumberOfTurbines,&
      'Turbine_AL_')
    call Calc_forces_ACL(wtm, fsi_wt, NumberOfTurbines,&
      'Turbine_AS_')

    if(myid.eq.0) print *, 'plyudebug, admr_run 8: Calc_F_eul'
    call Calc_F_eul(wtforce, wtforce_y, wtforce_z, fturbx, fturby, level, &
      wtm, fsi_wt, NumberOfTurbines, 1.0_wp, 1, flag_de)

    !if(InletRelaxation .eq. 2) call Add_Inlet_Forcing(u, v, w, wtforce, wtforce_y, wtforce_z, level)

    !> plyunote: wtforce initialization should be moved outside Calc_F_eul
    !> because both rotor_model_acs and nacelle_model_run will call Calc_F_eul

    if (nacelle_model .ne. 0) then
      call nacelle_model_run(time_, u, v, w, wtforce, wtforce_y, wtforce_z,&
        fturbx, fturby, level, flag_de)
    endif


  end subroutine rotor_model_admr
  
  subroutine rotor_model_acl(time_, u, v, w, wtforce, wtforce_y, wtforce_z, &
    fturbx, fturby, level, flag_de)
    implicit none
    real(wp) :: time_
    real(wp), dimension(:,:,:) :: u(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: v(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: w(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: flag_de(xsz(1), xsz(2), xsz(3))
    integer :: level
    real(wp) :: fturbx, fturby
    
    integer :: ibi
    real(wp) :: volflux(3)
 
    ti = nint(time_/dt)
    if (myid.eq.0) print *, 'ti=', ti
    call update_zgrid
    !call calc_vol_flux(u(:,:,1:xsz(3)), volflux)
    if(output_wake_switch .ne. 0) then
      call Export_Wake(fsi_wt, NumberOfTurbines, u, v, w, level)
    endif
    
    !call print_ibmnode1(7, 1, 10)
    
    if(fsitype .eq. 1) then
      call Prescribe_para_update(time_, pp_wt, NumberOfTurbines)
      do ibi = 1, NumberOfTurbines
        call Prescribe_FSI_update(fsi_wt(ibi), wtm(ibi), pp_wt(ibi), ibi)
      enddo
    endif 

    !call print_ibmnode1(7,1,11)
    
    if(myid.eq.0) print *, 'plyudebug, acl_run 1: Uref_ACL'
    !> in Uref_ACL, Calc_U_lagr is called for ref_vel disk
    !> this subroutine get value of Uref for reference disk
    !> and this value will only be used for angvel update.
    !> In Calc_F_lagr_ACL, Urefl(n_elmt,3) is calculated at the blade location
    call Uref_ACL(u, v, w, wtm, NumberOfTurbines, level)

    if(myid.eq.0) print *, 'plyudebug, acl_run 2: turbineangvel'
    call Calc_turbineangvel(dt, wtm, fsi_wt)

    if(myid.eq.0) print *, 'plyudebug, acl_run 3: rotor_Rot'
    do ibi = 1,NumberOfTurbines
      call rotor_Rot(time_, fsi_wt(ibi), wtm(ibi), dt, ibi, 1)
    enddo

    if(myid.eq.0) print *, 'plyudebug, acl_run 4: Pre_process(wtm)'
    call Pre_process(wtm, NumberOfTurbines, 1)

    !> here Calc_U_lagr is called for ACL elements
    if(myid.eq.0) print *, 'plyudebug, acl_run 5: Calc_U_lagr'
    call Calc_U_lagr(u, v, w, wtm, NumberOfTurbines, level)

    if(myid.eq.0) print *, 'plyudebug, acl_run 6: Calc_F_lagr_ACL'
    call Calc_F_lagr_ACL(wtm, fsi_wt, NumberOfTurbines)

    !> plyunote: force exportation moved to main 
    !call Export_ForceOnBlade(wtm, fsi_wt, NumberOfTurbines)

    !> plyunote: detail calculation is actally skipped now
    if(myid.eq.0) print *, 'plyudebug, acl_run 7: Calc_forces_ACL'
    call Calc_forces_ACL(wtm, fsi_wt, NumberOfTurbines, 'Turbine_AL_')

    if(myid.eq.0) print *, 'plyudebug, acl_run 8: Calc_F_eul'
    call Calc_F_eul(wtforce, wtforce_y, wtforce_z, fturbx, fturby, level, &
      wtm, fsi_wt, NumberOfTurbines, 1.0_wp, 1, flag_de)

    !if(InletRelaxation .eq. 2) call Add_Inlet_Forcing(u, v, w, wtforce, wtforce_y, wtforce_z, level)
    
    if (nacelle_model .ne. 0) then
      call nacelle_model_run(time_, u, v, w, wtforce, wtforce_y, wtforce_z,&
        fturbx, fturby, level, flag_de)
    endif

  end subroutine rotor_model_acl
 
  subroutine nacelle_model_run(time_, u, v, w, wtforce, wtforce_y, wtforce_z, &
    fturbx, fturby, level, flag_de)
    implicit none
    real(wp) :: time_
    real(wp), dimension(:,:,:) :: u(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: v(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: w(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), 1-level:xsz(3)+level)
    !real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp), dimension(:,:,:) :: wtforce(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_y(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: wtforce_z(xsz(1), xsz(2), xsz(3))
    real(wp), dimension(:,:,:) :: flag_de(xsz(1), xsz(2), xsz(3))
    integer :: level
    real(wp) :: fturbx, fturby

    ! nr is direction from nacelle body center to a surface element center
    ! na is direction of rotation axis of nacelle, maybe same with rotor's
    ! nt is circular velocity direction of a surface element
    real(wp), dimension(3) :: nr, na, nt, r_vec
    integer :: i, ibi, ipt, inac
    real(wp) :: r_amp, Ut

    !call Comput_actualshear_nacelle(ibm_nac, fsi_nac, NumberOfNacelle)
    !> set U_ref to 0 for t=0,1,tstart,tstart+1?

    !call Comput_modelcoef_nacelle(ibm_nac, fsi_nac, NumberOfNacelle)

    !call Comput_desireshear_nacelle(ibm_nac, fsi_nac, NumberOfNacelle)

    if (myid.eq.0) print *, "plyudebug: nacelle_model_run"

    if (rotate_nacelle .ne. 0) then
      call Calc_nacelle_angvel(NumberOfNacelle, fsi_nac, NumNacellePerLoc, fsi_wt)
    endif

    if (nacelle_model .eq. 7) then
     ! Comput_nut_nacelle(ibm_nac, fsi_nac, NumberOfNacelle)
    endif
   
    !> For prescribed motion, update nacelle element's location
    if (fsitype .eq. 1) then
      do ibi = 1, NumberOfNacelle
        call Prescribe_FSI_update(fsi_nac(ibi), ibm_nac(ibi), &
          pp_wt((ibi-1)/NumNacellePerLoc+1), (ibi-1)/NumNacellePerLoc+1)
        call nacelle_move(ibm_nac(ibi), fsi_nac(ibi), & 
          pp_wt((ibi-1)/NumNacellePerLoc+1))
      enddo
      
      !if(myid.eq.0) print *, "plyudebug: read_nacelle_param, Pre_process 1"
      call Pre_process(ibm_nac, NumberOfNacelle, 1)
      
      !if(myid.eq.0) print *, "plyudebug: read_nacelle_param, Coordinates_IP"
      call Coordinates_IP(ibm_nac, NumberOfNacelle)
      
      !> Pre_process_IP is integrated into Pre_process
      !if(myid.eq.0) print *, "plyudebug: read_nacelle_param, Pre_process 2"
      call Pre_process(ibm_nac, NumberOfNacelle, 2)
    endif
    
    call Calc_U_lagr(u, v, w, ibm_nac, NumberOfNacelle, level) 
    
    !> even though rotation of nacelle is permited here,
    !! the geometry is assumed to have no change during rotation
    do ibi = 1, NumberOfNacelle
      if (fsi_nac(ibi)%rotate_alongaxis .eq. 1) then
        do i = 1, ibm_nac(ibi)%n_v
          na(1) = fsi_nac(ibi)%nx_tb
          na(2) = fsi_nac(ibi)%ny_tb
          na(3) = fsi_nac(ibi)%nz_tb

          r_vec(1) = ibm_nac(ibi)%x_bp(i) - fsi_nac(ibi)%x_c
          r_vec(2) = ibm_nac(ibi)%y_bp(i) - fsi_nac(ibi)%y_c
          r_vec(3) = ibm_nac(ibi)%z_bp(i) - fsi_nac(ibi)%z_c
          r_amp = norm2(r_vec) + 1.0e-16_wp
          
          nr(1:3) = r_vec(1:3) / r_amp

          call crossx(na, nr, nt)

          fsi_nac(ibi)%angvel_axis = fsi_wt((ibi-1)/NumNacellePerLoc+1)%angvel_axis
          Ut = fsi_nac(ibi)%angvel_axis * r_amp

          ibm_nac(ibi)%u(i,1:3) = Ut * nt(1:3)

          if(abs(ibm_nac(ibi)%u(i,1))>1.0e10 .or. abs(ibm_nac(ibi)%u(i,2))>1.0e10) then
            print *, 'Error in nacelle u:', 'myid, ibi, i=', myid, ibi, i, &
              ', r_vec=', r_vec,', na=', na, ', nt=', nt, ', Ut=', Ut, ', u=', Ut*nt(1:3)
          endif
        enddo
      else
        do i = 1, ibm_nac(ibi)%n_v
          ibm_nac(ibi)%u(i,1:3) = 0.0_wp
        enddo   
      endif
    enddo

    do ibi = 1, NumLoc 
      do ipt = 1, NumNacellePerLoc
        inac = (ibi-1)*NumNacellePerLoc+ipt
        !if (rotor_model .ne. 0) then
          ibm_nac(inac)%U_ref = wtm(ibi)%U_ref
        !else
          ! plyunote: assume 1.0 is average inflow velocity
        !  ibm_nac(ibi)%U_ref = 1.0_wp
        !endif
      enddo
    enddo

    call Calc_F_lagr_nacelle(ibm_nac, fsi_nac, NumberOfNacelle)

    !> plyunote: Calc_Nut_eul is ignored here
    
    !if(myid.eq.0) print *, 'plyudebug, nacelle_run : Calc_forces_nac'
    call Calc_forces_nacelle(ibm_nac, fsi_nac, NumberOfNacelle,&
      'Nacelle_')
    !print *, "nacdbg, 1"
    !call Calc_forces_rotor(ibm_nac, fsi_nac, NumberOfNacelle)
    call Calc_F_eul(wtforce, wtforce_y, wtforce_z, fturbx, fturby, &
      level, ibm_nac, fsi_nac, NumberOfNacelle, 1.0_wp, 2, flag_de)
    !print *, "nacdbg, 2"
  end subroutine nacelle_model_run
  
  subroutine read_ad_dy_param
  ! Read parameters for Actuator Disk model, originally written by Di Yang
    implicit none
    
    !> plyunote: the read process need to be modified to parallel form
    open(15000, file="turbine_dy_param.inp", action="READ", status="OLD")
    read(15000,*) nxwt, nywt
    read(15000,*) rdisc, hdisc
    read(15000,*) c_bforce
    read(15000,*) i_thrust
    read(15000,*) c_thrust, c_thrust_induced
    close(15000) 
    
    !> plyunote: following lines are used for inlet Relaxation.
    !> read control parameters for wind turbine model
    open(15001, file='control.dat', action='read', status='old')
    read(15001,*) rotor_model
    read(15001,*) !NumberOfTurbines, num_foiltype, num_blade
    read(15001,*) !loc_refvel
    read(15001,*) !FixTurbineAngvel, FixTipspeedRatio
    read(15001,*) !halfwidth_dfunc
    !> read control parameters for nacelle model
    !> cf_nacelle_fromfile, rotate_nacelle, nacelle_model, reflength_nacelle
    !> r_nacelle, L_nacelle, dh_nacelle
    read(15001,*) !nacelle_model
    read(15001,*) !rotate_nacelle
    read(15001,*) !NumberOfNacelle, NumNacellePerLoc
    read(15001,*) !deltafunc_U, deltafunc_F
    read(15001,*) !Shen1_AL, Shen1_AL_tipcorrection, Shen1_AL_tipcorrectionratio_Fa
    
    read(15001,*) InletRelaxation, IR_bw, IR_dataindex, IR_scale, IR_ti_shift
    
    read(15001,*) !c_bforce
    read(15001,*) !adjust_CT_ACL, adjust_CP_ACL
    read(15001,*) !adjust_CT_nacelle
    read(15001,*) !fsitype
    !read(15001,*) IR_folder 
    close(15001)
     
  end subroutine read_ad_dy_param
  
  subroutine fraction_v2(ic,jp,kp,yc,zc, wtgamma_in, idisc_in, wndp_in)
    implicit none
    !>use mpi

    integer :: j, k, ic, jp, kp
    real(wp) :: yc, zc, wndp_in
    real(wp), dimension(:,:,:) :: wtgamma_in
    integer, dimension(:,:,:) :: idisc_in

    real(wp) :: rcell, rcc, aa, acell
    real(wp) :: ylt, zlt, yrt, zrt, ylb, zlb, yrb, zrb, rlt, rlb, rrt, rrb
    real(wp) :: yyl, zzl, ly, lz, ly2, lz2

    !< calculate cell boundary (value need be checked, plyu)
    ylt = (jp-1)*dy - dy/2.
    zlt = (zzo(kp)+zzo(kp+1))/2.*HBAR
    rlt = sqrt((ylt-yc)**2+(zlt-zc)**2)
    !< (strange: rlt is calculated in transformed coord
    !< while it is later compared with ground coord rdisc, plyu)

    ylb = (jp-1)*dy - dy/2.
    zlb = (zzo(kp)+zzo(kp-1))/2.*HBAR
    rlb = sqrt((ylb-yc)**2+(zlb-zc)**2)

    yrt = (jp-1)*dy + dy/2.
    zrt = (zzo(kp)+zzo(kp+1))/2.*HBAR
    rrt = sqrt((yrt-yc)**2+(zrt-zc)**2)

    yrb = (jp-1)*dy + dy/2.
    zrb = (zzo(kp)+zzo(kp-1))/2.*HBAR
    rrb = sqrt((yrb-yc)**2+(zrb-zc)**2)

    !< check which grid cell is complete in or out of disc
    !< grid cell is centered at grid point
    if(rlt .lt. rdisc .and. rlb .lt. rdisc .and. rrt .lt. rdisc &
      .and. rrb .lt. rdisc) then
      wtgamma_in(ic,jp,kp) = 1.
      idisc_in(ic,jp,kp) = 1
    endif

    rcell = sqrt((ylt-yrb)**2 + (zlt-zrb)**2)/2.
    rcc = sqrt(((jp-1)*dy-yc)**2+(zzo(kp)*hbar-zc)**2)
    !< this condition looks strange, plyu
    if(rcc .gt. (rcell+rdisc)) then
      wtgamma_in(ic,jp,kp) = 0.
    endif


    !> calculate fraction when 3 corners are inside disc

    !>     case 1
    !>     left top corner out of disc
    if(rlt.gt.rdisc.and.rlb.lt.rdisc.and.rrt.lt.rdisc &
      .and.rrb.lt.rdisc) then

      yyl=yc-ylt
      zzl=sqrt(rdisc**2-yyl**2)
      lz=zlt-(zc+zzl)

      zzl=zlt-zc
      yyl=sqrt(rdisc**2-zzl**2)
      ly=(yc-ylt)-yyl

      if(ly.lt.0.or.lz.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 1.'
        print*, 'ly=',ly
        print*, 'lz=',lz
        stop
      endif

      aa=0.5*ly*lz
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 1.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=(acell-aa)/acell

      idisc_in(ic,jp,kp)=1
    !> case 1 end here

    !>      case 2
    !>      right top corner out of disc
    else if(rrt.gt.rdisc.and.rlt.lt.rdisc.and.rlb.lt.rdisc &
      .and.rrb.lt.rdisc) then

      yyl=yrt-yc
      zzl=sqrt(rdisc**2-yyl**2)
      lz=zrt-(zc+zzl)

      zzl=zrt-zc
      yyl=sqrt(rdisc**2-zzl**2)
      ly=(yrt-yc)-yyl

      if(ly.lt.0.or.lz.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 2.'
        print*, 'ly=',ly
        print*, 'lz=',lz
        stop
      endif

      aa=0.5*ly*lz
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 2.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=(acell-aa)/acell

      idisc_in(ic,jp,kp)=1
      !>  case 2 end here

    !>      case 3
    !>      left bottom corner out of disc
    else if(rlb.gt.rdisc.and.rlt.lt.rdisc.and.rrt.lt.rdisc &
      .and.rrb.lt.rdisc) then

      yyl=yc-ylb
      zzl=sqrt(rdisc**2-yyl**2)
      lz=(zc-zlb)-zzl

      zzl=zc-zlb
      yyl=sqrt(rdisc**2-zzl**2)
      ly=(yc-ylb)-yyl

      if(ly.lt.0.or.lz.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 3.'
        print*, 'ly=',ly
        print*, 'lz=',lz
        stop
      endif

      aa=0.5*ly*lz
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 3.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=(acell-aa)/acell

      idisc_in(ic,jp,kp)=1
      !> case 3 end here

    !>     case 4
    !>      right bottom corner out of disc
    else if(rrb.gt.rdisc.and.rlt.lt.rdisc.and.rlb.lt.rdisc &
      .and.rrt.lt.rdisc) then

      yyl=yrb-yc
      zzl=sqrt(rdisc**2-yyl**2)
      lz=(zc-zrb)-zzl

      zzl=zc-zrb
      yyl=sqrt(rdisc**2-zzl**2)
      ly=(yrb-yc)-yyl

      if(ly.lt.0.or.lz.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 4.'
        print*, 'ly=',ly
        print*, 'lz=',lz
        stop
      endif

      aa=0.5*ly*lz
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 4.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=(acell-aa)/acell

      idisc_in(ic,jp,kp)=1
      !>  case 4 end here

    !>  end here

    !>  calculate fraction when 2 corners are inside disc

    !>     case 5
    !>     two top corners are out of disc


    else if(rlt.gt.rdisc.and.rrt.gt.rdisc.and.rlb.lt.rdisc &
      .and.rrb.lt.rdisc) then

      yyl=yc-ylt
      zzl=sqrt(rdisc**2-yyl**2)
      lz=zzl-(zlb-zc)

      yyl=yrt-yc
      zzl=sqrt(rdisc**2-yyl**2)
      lz2=zzl-(zrb-zc)

      if(lz.lt.0.or.lz2.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 5.'
        print*, 'lz=',lz
        print*, 'lz2=',lz2
        stop
      endif

      aa=0.5*(lz+lz2)*(yrb-ylb)
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 5.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=aa/acell

      idisc_in(ic,jp,kp)=1
      !> case 5 end here

    !>      case 6
    !>      two bottom corners are out of disc
    else if(rlb.gt.rdisc.and.rrb.gt.rdisc.and.rlt.lt.rdisc &
      .and.rrt.lt.rdisc) then

      yyl=yc-ylb
      zzl=sqrt(rdisc**2-yyl**2)
      lz=zzl-(zc-zlt)

      yyl=yrb-yc
      zzl=sqrt(rdisc**2-yyl**2)
      lz2=zzl-(zc-zrt)

      if(lz.lt.0.or.lz2.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 6.'
        print*, 'lz=',lz
        print*, 'lz2=',lz2
        stop
      endif

      aa=0.5*(lz+lz2)*(yrb-ylb)
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 6.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=aa/acell

      idisc_in(ic,jp,kp)=1
      !> case 6 end here

    !>      case 7
    !>      two left corners are out of disc
    else if(rlt.gt.rdisc.and.rlb.gt.rdisc.and.rrt.lt.rdisc &
      .and.rrb.lt.rdisc) then

      zzl=zlt-zc
      yyl=sqrt(rdisc**2-zzl**2)
      ly=yyl-(yc-yrt)

      zzl=zc-zlb
      yyl=sqrt(rdisc**2-zzl**2)
      ly2=yyl-(yc-yrb)

      if(ly.lt.0.or.ly2.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 7.'
        print*, 'ly=',ly
        print*, 'ly2=',ly2
        stop
      endif

      aa=0.5*(ly+ly2)*(zrt-zrb)
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 7.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=aa/acell

      idisc_in(ic,jp,kp)=1
      !> case 7 end here

    !>     case 8
    !>     two right corners are out of disc
    else if(rrt.gt.rdisc.and.rrb.gt.rdisc.and.rlt.lt.rdisc &
      .and.rlb.lt.rdisc) then

      zzl=zrt-zc
      yyl=sqrt(rdisc**2-zzl**2)
      ly=yyl-(ylt-yc)

      zzl=zc-zrb
      yyl=sqrt(rdisc**2-zzl**2)
      ly2=yyl-(ylb-yc)

      if(ly.lt.0.or.ly2.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 8.'
        print*, 'ly=',ly
        print*, 'ly2=',ly2
        stop
      endif

      aa=0.5*(ly+ly2)*(zlt-zlb)
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 8.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=aa/acell

      idisc_in(ic,jp,kp)=1
      !> case 8 end here

    !> end here

    !>      calculate fraction when 1 corners are inside disc

    !>     case 9
    !>     only right bottom corner is inside disc
    else if(rlt.gt.rdisc.and.rlb.gt.rdisc.and.rrt.gt.rdisc &
      .and.rrb.lt.rdisc) then

      yyl=yc-yrt
      zzl=sqrt(rdisc**2-yyl**2)
      lz=zzl-(zrb-zc)

      zzl=zrb-zc
      yyl=sqrt(rdisc**2-zzl**2)
      ly=yyl-(yc-yrb)

      if(ly.lt.0.or.lz.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 9.'
        print*, 'ly=',ly
        print*, 'lz=',lz
        stop
      endif

      aa=0.5*ly*lz
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 9.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=aa/acell

      idisc_in(ic,jp,kp)=1
      !> case 9 end here

    !>      case 10
    !>      only left bottom corner is inside disc
    else if(rlt.gt.rdisc.and.rrt.gt.rdisc.and.rrb.gt.rdisc &
      .and.rlb.lt.rdisc) then

      yyl=ylt-yc
      zzl=sqrt(rdisc**2-yyl**2)
      lz=zzl-(zlb-zc)

      zzl=zlb-zc
      yyl=sqrt(rdisc**2-zzl**2)
      ly=yyl-(ylb-yc)

      if(ly.lt.0.or.lz.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 10.'
        print*, 'ly=',ly
        print*, 'lz=',lz
        stop
      endif

      aa=0.5*ly*lz
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 10.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=aa/acell

      idisc_in(ic,jp,kp)=1
      !> case 10 end here

    !>     case 11
    !>     only right top corner is inside disc
    else if(rlt.gt.rdisc.and.rlb.gt.rdisc.and.rrb.gt.rdisc &
      .and.rrt.lt.rdisc) then

      yyl=yc-yrt
      zzl=sqrt(rdisc**2-yyl**2)
      lz=zzl-(zc-zrt)

      zzl=zc-zrt
      yyl=sqrt(rdisc**2-zzl**2)
      ly=yyl-(yc-yrt)

      if(ly.lt.0.or.lz.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 11.'
        print*, 'ly=',ly
        print*, 'lz=',lz
        stop
      endif

      aa=0.5*ly*lz
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 11.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=aa/acell

      idisc_in(ic,jp,kp)=1
      !> case 11 end here

    !>     case 12
    !>     only left bottom corner is inside disc
    else if(rrt.gt.rdisc.and.rlb.gt.rdisc.and.rrb.gt.rdisc &
      .and.rlt.lt.rdisc) then

      yyl=ylt-yc
      zzl=sqrt(rdisc**2-yyl**2)
      lz=zzl-(zc-zlt)

      zzl=zc-zlt
      yyl=sqrt(rdisc**2-zzl**2)
      ly=yyl-(ylt-yc)

      if(ly.lt.0.or.lz.lt.0) then
        print*, 'error in fraction!'
        print*, 'case 12.'
        print*, 'ly=',ly
        print*, 'lz=',lz
        stop
      endif

      aa=0.5*ly*lz
      acell=(zlt-zlb)*(yrt-ylt)
      if(aa.le.0.or.acell.le.0) then
        print*, 'error in fraction area!'
        print*, 'case 12.'
        print*, 'aa=',aa
        print*, 'acell=',acell
        stop
      endif
      wtgamma_in(ic,jp,kp)=aa/acell

      idisc_in(ic,jp,kp)=1

      endif
      !>  case 12 end here

    wndp_in=wndp_in+wtgamma_in(ic,jp,kp)

  end subroutine fraction_v2

  subroutine discloca_v3
    use spectral_hos
    implicit none
    
    integer :: i, j, k, ii, jj, kk, ic, jc, kc
    integer :: kl, kh, jr, kp, jp, nxr, nyr

    real(wp) :: xc, yc, zc, yp, zp, rp

    !real(wp), allocatable, dimension(:,:) :: eo, hho
    real(wp) :: ztb(nx, ny)

    ! eo, hho is allocated to size of hos, 
    ! while later computation use size of nx,ny,nz
    ! therefore nxhos should be equal to nx
    !allocate(eo(nxhos,nyhos), hho(nxhos,nyhos))
    !allocate(eo(xsz(1), xsz(2)), hho(xsz(1), xsz(2)))

    !> move the allocation of eo and hho to reading_part
    !allocate(eo(nx_global, ny_global), hho(nx_global, ny_global))

    do k=1,nz
      zz1(k) = -1.0D16
      zzo(k) = 0.
    enddo
    do k=1,xsz(3)
      kk = xst(3)+k-1
      zz1(kk)=zz(k)
    enddo
    call mpi_allreduce(zz1,zzo,nz,mpi_double_precision,&
      mpi_max,mpi_comm_world,ierr_wt)
 !   do k=1,nz
     ! print *, 'plyudebug: myid,k,zz1,zzo,xst,xsz',&
     !   myid,k,zz1(k),zzo(k),xst(3),xsz(3)
 !   enddo
    ! dx=xl/nx
    nxr = int(rdisc/dx) + 1
    nyr = int(rdisc/dy) + 1

    !< check spanwise domain size
    if (2.0 * rdisc * nywt .gt. yl) then
      print *, 'Insufficient domain size for wind turbine!'
      print *, 'Radius of wind turbine:', rdisc
      print *, 'Number of columns:', nywt
      print *, 'Minimum spanwise space required:', 2.*rdisc*nywt
      print *, 'Available spanwise space:', yl
      stop
    endif

    !< send surface information to one cpu
    call alltoone(eta, eo)
    call alltoone(hh, hho)

!    do j=1,ny
!      do i=1,nx
! print *,'plyudebug: myid,i,j,eta,hh,eo,hho',&
!   myid,i,j,eta(i,j),hh(i,j),eo(i,j),hho(i,j)
!      enddo
!    enddo

    !< initialize parameter
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          idisc(i,j,k) = 0
          ndisc(i,j,k) = 0
          wtgamma0(i,j,k) = 0.
          idisc_ref(i,j,k) = 0
          ndisc_ref(i,j,k) = 0
        enddo
      enddo
    enddo

    do i = 1,nxwt
      icut(i) = 0
      icut_ref(i) = 0
    enddo
    do j = 1,nywt
      jcut(j) = 0
      jcut_ref(j) = 0
    enddo
    kcut=0

    !< ztb might be modified for special case
    ztb = 0.

    !< calculate position of wind turbine
    !< example, if ny=32:
    !< if nywt=1, then jcut=16
    !< if nywt=2, then jcut=8,24
    !< if nywt=3, then jcut=6,16,26
    iwt = 0
    nxmod0 = nx -1
    nymod0 = ny -1
    do ii = 1,nxwt
      !icut(ii) = (2*ii-1)*nxmod0/nxwt/2+1
      icut(ii) = nint(0.375/(xl/nx))
      if (myid.eq.0) print *, 'Attention, xc is at ',icut(ii)*(xl/nx) 

      icut_ref(ii) = int((icut(ii)*xl/nx - 4.0*rdisc)/xl*nx)
      if (icut_ref(ii)<1) then
        icut_ref(ii) = icut_ref(ii) + nx
      endif
    enddo
    do jj = 1,nywt
      jcut(jj) = (2*jj-1)*nymod0/nywt/2+1
      jcut_ref(jj) = jcut(jj) 
    enddo

    do ii = 1,nxwt
      do jj = 1,nywt
        iwt = (ii-1)*nywt+jj 
        ic = icut(ii)
        jc = jcut(jj)
        yc = (jc-1)*dy
        wndp(iwt) = 0.
        !print *, 'plyudebug:myid,ii,jj,icut,jcut=',&
        !  myid,ii,jj,ic,jc
        do k = 1,nz
!> plyunote: in Di's code, z=zz(k)*(hbar+etas)-hhs
!! he also use hh=-eta*(1.-fexp)
          if( (zzo(k)*(hbar+hho(ic,jc)) - hho(ic,jc)) .gt. &
            (ztb(ic,jc)+hdisc) ) then
            kc = k
            exit
          endif
        enddo
        zc = hdisc
        !kcut = kcut + kc
        kcut = kc
        !print *,'plyudebug:myid,ii,jj,ic,jc,kc,eo,hho',&
        !  myid,ii,jj,ic,jc,kc,eo(ic,jc),hho(ic,jc)

        do k = 1,nz
          if( (zzo(k)*(hbar+hho(ic,jc)) - hho(ic,jc)) .gt. &
            (ztb(ic,jc)+hdisc-rdisc) ) then
            kl = k-1
            exit
          endif
        enddo
        do k = 1,nz
          if( (zzo(k)*(hbar+hho(ic,jc)) - hho(ic,jc)) .gt. &
            (ztb(ic,jc)+hdisc+rdisc) ) then
            kh = k
            exit
          endif
        enddo

        !print *, 'plyudebug_tm1: myid,ii,jj,ic,jc,kc,kl,kh=',&
  !  myid,ii,jj,ic,jc,kc,kl,kh
        do kp = kl,kh
          do jr = -nyr,nyr
            jp = jc+jr
            if(jp .lt. 1) jp = jp+nymod0
            if(jp .gt. nymod0) jp = jp-nymod0
            yp = (jp-1)*dy
            zp = zzo(kp)*(hbar+hho(ic,jp))-hho(ic,jp)
            rp = sqrt((yp-yc)**2+(zp-zc)**2)
            !< (yp, zp, rp are calculated but not used, plyu)

            ndisc(ic,jp,kp)=iwt
            call fraction_v2(ic,jp,kp,yc,zc,wtgamma0, idisc, wndp(iwt))
          enddo
        enddo
      enddo
    enddo
    
    do ii = 1,nxwt
      do jj = 1,nywt
        iwt = (ii-1)*nywt+jj 
        ic = icut_ref(ii)
        jc = jcut_ref(jj)
        yc = (jc-1)*dy
        wndp_ref(iwt) = 0.
        !print *, 'plyudebug:myid,ii,jj,icut,jcut=',&
        !  myid,ii,jj,ic,jc
        do k = 1,nz
          if( (zzo(k)*(hbar+hho(ic,jc)) - hho(ic,jc)) .gt. &
            (ztb(ic,jc)+hdisc) ) then
            kc = k
            exit
          endif
        enddo
        zc = hdisc
        !kcut = kcut + kc
        kcut = kc
        !print *,'plyudebug:myid,ii,jj,ic,jc,kc,eo,hho',&
        !  myid,ii,jj,ic,jc,kc,eo(ic,jc),hho(ic,jc)

        do k = 1,nz
          if( (zzo(k)*(hbar+hho(ic,jc)) - hho(ic,jc)) .gt. &
            (ztb(ic,jc)+hdisc-rdisc) ) then
            kl = k-1
            exit
          endif
        enddo
        do k = 1,nz
          if( (zzo(k)*(hbar+hho(ic,jc)) - hho(ic,jc)) .gt. &
            (ztb(ic,jc)+hdisc+rdisc) ) then
            kh = k
            exit
          endif
        enddo

        !print *, 'plyudebug_tm1: myid,ii,jj,ic,jc,kc,kl,kh=',&
  !  myid,ii,jj,ic,jc,kc,kl,kh
        do kp = kl,kh
          do jr = -nyr,nyr
            jp = jc+jr
            if(jp .lt. 1) jp = jp+nymod0
            if(jp .gt. nymod0) jp = jp-nymod0
            yp = (jp-1)*dy
            zp = zzo(kp)*(hbar+hho(ic,jp))-hho(ic,jp)
            rp = sqrt((yp-yc)**2+(zp-zc)**2)
            !< (yp, zp, rp are calculated but not used, plyu)

            ndisc_ref(ic,jp,kp)=iwt
            call fraction_v2(ic,jp,kp,yc,zc, wtgamma0_ref, idisc_ref,&
               wndp_ref(iwt))
          enddo
        enddo
      enddo
    enddo

    !KCUT=INT(KCUT/FLOAT(NXWT*NYWT))

   ! plyudebug 
   !do k = 1,nz
   !   do j = 1,ny
   !     do i = 1,nx
   !       if(idisc(i,j,k) .ne. 0) then
   !         !print *, 'plyudebug_tm2: i,j,k,ndisc,wtgamma0', &
   !         !  i,j,k,ndisc(i,j,k),wtgamma0(i,j,k)
   !       endif
   !     enddo
   !   enddo
   ! enddo

    !deallocate(eo, hho)

  end subroutine discloca_v3

  subroutine veldisc_initial (time_)
    !< disc-averaged and time-averaged velocity
    implicit none
    !use mpi
    !< use navier for u(nx,ny/ncpu_hos,*)
    !< use navier, only : u
    real(wp) :: time_


    integer i, j, k, jj

    real(wp) :: u0(nx,ny,nz), u1(nx,ny,nz)
    !real(wp) :: ud1(nxwt*nywt)

    !< tc: exponential time coefficient for relaxation process
    real(wp) :: tc, fexp
    
    first_step_time = time_
    ti = nint(time_/dt)
    ti_first = nint(first_step_time/dt)
    
    !< initialize variables
    tc = 14.14

    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr_wt)

    allocate(ndisc(nx, ny, nz), idisc(nx, ny, nz))
    allocate(wtgamma0(nx, ny, nz))
    allocate(icut(nxwt), jcut(nywt))
    allocate(ud(nxwt*nywt), ud1(nxwt*nywt))
    allocate(wndp(nxwt*nywt))
    
    allocate(ndisc_ref(nx, ny, nz), idisc_ref(nx, ny, nz))
    allocate(wtgamma0_ref(nx, ny, nz))
    allocate(icut_ref(nxwt), jcut_ref(nywt))
    allocate(ud_ref(nxwt*nywt), ud1_ref(nxwt*nywt))
    allocate(wndp_ref(nxwt*nywt))
    
    !> plyunote: zzo, zz1, eo, hho are also allocated in collect_grid of ACL model
    allocate(zzo(nz), zz1(nz))
    allocate(eo(nx_global, ny_global), hho(nx_global, ny_global))
    allocate(hxo(nx_global, ny_global), hyo(nx_global, ny_global))
    
    do k=1,nz
      do j=1,ny
        do i=1,nx
          u0(i,j,k)=0.
          u1(i,j,k)=0.
        enddo
      enddo
    enddo

    do i=1,nxwt
      do j=1,nywt
        ud1((i-1)*nywt+j)=0.  
        ud1_ref((i-1)*nywt+j)=0.  
      enddo
    enddo
    
    !> here a collection of U from all cpus are omitted by Di
    !> maybe it will happend later, plyu
    !> also, it should be noted that ndp is not assigned.
    
    !> calculate ud: disc- and time-averaged velocity
    do i=1,nxwt
      do j=1,nywt
        ud((i-1)*nywt+j)=ud1((i-1)*nywt+j)
        ud_ref((i-1)*nywt+j)=ud1_ref((i-1)*nywt+j)
      enddo
    enddo

  end subroutine veldisc_initial

  subroutine reread_admnr
    implicit none
    integer :: i, j
    if(myid .eq. 0) then
      open(15010, file='restart_ad.dat', action='read', status='old')
      do i = 1, nxwt
        do j = 1, nywt
          read(15010,*) ud1((i-1)*nywt+j), ud((i-1)*nywt+j), &
            ud1_ref((i-1)*nywt+j), ud_ref((i-1)*nywt+j)
        enddo
      enddo
      close(15010)
    endif 
    
    call mpi_bcast(ud1, nxwt*nywt, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(ud, nxwt*nywt, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(ud1_ref, nxwt*nywt, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)
    call mpi_bcast(ud_ref, nxwt*nywt, mpi_double_precision, 0, mpi_comm_2d_cart, ierr_wt)

  end subroutine reread_admnr

  subroutine veldisc_v2(u, level, time_)
    !> disc-averaged and time-averaged velocity
    implicit none
    
    integer :: level
    integer :: i, j, k, ii, jj, kk
    real(wp) :: time_
    
    real(wp) :: u(xsz(1),xsz(2),1-level:xsz(3)+level)
    real(wp) :: u0(nx,ny,nz), u1(nx,ny,nz)
    !real(wp) :: ud1(nxwt*nywt)

    !< tc: exponential time coefficient for relaxation process
    !  tc=T*u_{*}/H=0.27 corresponds to T=10min using possible reference
    !  dimensional values u_{*}=0.45, and H=1000m
    real(wp) :: tc, fexp
    
    ti = nint(time_ / dt)
    !print *, 1, ti

    !< tc: relaxation time window for temporal averaging (calaf et al. 2010)
    if (ti < (0.15 * 0.27*hbar/usbot/dt)) then
      tc = 0.15 * 0.27
    elseif (ti < (0.45 * 0.27*hbar/usbot/dt)) then
      tc = 0.45 * 0.27
    else
      tc = 0.27
    endif
    
    do k=1,nz
      do j=1,ny
        do i=1,nx
          u0(i,j,k)=0.
          u1(i,j,k)=0.
        enddo
      enddo
    enddo

    do i=1,nxwt
      do j=1,nywt
        ud1((i-1)*nywt+j)=0.  
        ud1_ref((i-1)*nywt+j)=0.  
      enddo
    enddo
    
    !> collect u from all cpus
    do k=1,xsz(3)
      kk = xst(3)+k-1
      do j=1,xsz(2)
        jj = xst(2)+j-1
        do i=1,xsz(1)
          ii = xst(1)+i-1
          u1(ii,jj,kk)=u(i,j,k)
        enddo
      enddo
    enddo

    call mpi_allreduce(u1,u0,nx*ny*nz, &
      mpi_double_precision,mpi_sum,mpi_comm_world,ierr_wt)

    !> calculate ud1: disc-averaged velocity
    !  In Di's version, wndp is sum of all turbine's wtgamma0, 
    !  so wtgamma0/wndp is a weight for average among all turbines.
    !  I don't really agree with this, especially in wind farm case,
    !  disk-average velocity for different locations in the farm should be
    !  different.
    !  I change wndp to an array corresponding every seperate turbine
    do k=1,nz
      do j=1,ny
        do i=1,nx
          if(idisc(i,j,k).eq.1) then
            ud1(ndisc(i,j,k))=ud1(ndisc(i,j,k)) &
              +u0(i,j,k)*wtgamma0(i,j,k)/wndp(ndisc(i,j,k))
          endif
          if(idisc_ref(i,j,k).eq.1) then
            ud1_ref(ndisc_ref(i,j,k))=ud1_ref(ndisc_ref(i,j,k)) &
              +u0(i,j,k)*wtgamma0_ref(i,j,k)/wndp_ref(ndisc_ref(i,j,k))
          endif
        enddo
      enddo
    enddo
    

    !> calculate ud: disc- and time-averaged velocity
    !  usbot = 1.0_wp / (2.5_wp * log(hbar/z0))
    do i=1,nxwt
       do j=1,nywt
          fexp=dt/(tc*hbar/usbot)/(1.+dt/(tc*hbar/usbot))
          !fexp=1.
          ud((i-1)*nywt+j)=ud((i-1)*nywt+j)*(1-fexp) &
              +ud1((i-1)*nywt+j)*fexp
          ud_ref((i-1)*nywt+j)=ud_ref((i-1)*nywt+j)*(1-fexp) &
              +ud1_ref((i-1)*nywt+j)*fexp
       enddo
    enddo

    if (myid.eq.0) then
      do i = 1, nxwt
        do j = 1, nywt
          print *, (i-1)*nywt+j, 'ud_instant=',ud1((i-1)*nywt+j), ', ud_mean=', &
            ud((i-1)*nywt+j), 'udref_instant=', ud1_ref((i-1)*nywt+j), &
            'udref_mean=', ud_ref((i-1)*nywt+j)
        enddo
      enddo
    endif
    
  end subroutine veldisc_v2
  
  subroutine wind_turbine_force(wtforce, level, fturbx)
    
    !> disc model for wind turbine force
    implicit none
  
    !real(wp) :: wtforce(xsz(1), xsz(2), 1-level:xsz(3)+level)
    real(wp) :: wtforce(xsz(1), xsz(2), xsz(3))
    real(wp) :: fturbx
    integer :: level
    
    !> ct: thrust coefficient (ct' in calaf_meneveau_meyers_2010)
    real(wp) :: ct
    
    integer :: i, j, k, ii, jj, kk

    ct = 1.33
    
    ! dx=twopi/pex/float(nxmod)
    !do k=1-level,xsz(3)+level
    do k=1,xsz(3)
      do j=1,xsz(2)
        do i=1,xsz(1)
          wtforce(i,j,k)=0.
        enddo
      enddo
    enddo

    do k=1,xsz(3)
      kk = xst(3)+k-1
      do j=1,xsz(2)
        ! old code
  !jj=j+myid*ny/ncpu_hos 

  ! new code in new domain decomposition
        jj = xst(2)+j-1
        do i=1,xsz(1)
          ii = xst(1)+i-1
          if(idisc(ii,jj,kk) .ne. 0) then
            !print *, 'plyudebug',i,j,k,ii,jj,kk,ndisc(ii,jj,kk),&
            !  wtgamma0(ii,jj,kk)
            if (i_thrust .eq. 1) then
              wtforce(i,j,k)=-0.5*c_thrust*ud_ref(ndisc(ii,jj,kk))**2 &
                *wtgamma0(ii,jj,kk)/dx
            elseif (i_thrust .eq. 2) then
              wtforce(i,j,k)=-0.5*c_thrust*ud1_ref(ndisc(ii,jj,kk))**2 &
                *wtgamma0(ii,jj,kk)/dx
            elseif (i_thrust .eq. 3) then
              wtforce(i,j,k)=-0.5*c_thrust_induced*ud(ndisc(ii,jj,kk))**2 &
                *wtgamma0(ii,jj,kk)/dx
            elseif (i_thrust .eq. 4) then
              wtforce(i,j,k)=-0.5*c_thrust_induced*ud1(ndisc(ii,jj,kk))**2 &
                *wtgamma0(ii,jj,kk)/dx
            else
              print *, "i_thrust value error, not in any categary"
            endif
          endif
        enddo
      enddo
    enddo

    !> c_bforce is recommended to be thrust coefficient. likely to be 3/4
    fturbx = 0.0_wp
    do i = 1, nxwt
      do j = 1, nywt
        fturbx = fturbx + (c_bforce * (ud_ref((i-1)*nywt+j) )&
              **2 * twopi/2.0_wp * rdisc **2) / (xl*yl*zl*hbar)
      enddo
    enddo

  end subroutine wind_turbine_force

  subroutine controller_pid (in_param, out_response, setpoint, ncontr)
    implicit none
    real(wp) :: in_param, out_response, setpoint
    integer :: ncontr

    real(wp) :: dp

!http://homepages.ed.ac.uk/jwp/control06/controlcourse/restricted/course/advanced/module4.html
!http://www.cnblogs.com/cv-pr/p/4785195.html
  
  !real(wp) :: kc_p, kc_i, kc_d, sp_bforce
  !real(wp) :: err_0, err_1, err_2
  !integer :: controller_n
    kc_p = 5.0; kc_i = 0.1; kc_d = 5.0

    ncontr = ncontr + 1

    if (ncontr .eq. 1) then
      err_2 = 0.0_wp; err_1 = 0.0_wp
      err_0 = setpoint - out_response
      dp = 0.0_wp 
    else if (ncontr .eq. 2) then
      err_2 = 0.0_wp
      err_1 = err_0
      err_0 = setpoint - out_response
      dp = kc_p * (err_0-err_1) + kc_i * err_0
    else
      err_2 = err_1
      err_1 = err_0
      err_0 = setpoint - out_response
      dp = kc_p * (err_0-err_1) + kc_i * err_0 + kc_d *(err_0-2.0*err_1+err_2)
    endif
     
    in_param = in_param + dp    

  end subroutine controller_pid

  subroutine calc_vol_flux(u_, volflux_)
    implicit none
    real(wp), dimension(:,:,:), intent(in) :: u_
    real(wp), intent(out) :: volflux_(3)

    real(wp) :: vol_, area_, volr_, arear_
    integer :: i, j, k
    real(wp) :: temps(2), tempr(2)

    area_ = 0.0
    vol_ = 0.0
    i = 1
    do j = 1, xsz(2)
      do k = 1, xsz(3)-1
        area_ = area_ + (zz(k+1)-zz(k))*(hbar+hh(i,j))
        vol_ = vol_ + (u_(i,j,k)+u_(i,j,k+1))/2.0* &
          (zz(k+1)-zz(k))*(hbar+hh(i,j))
      enddo
    enddo
    temps(1) = area_; temps(2) = vol_
    tempr = 0.0
    call mpi_allreduce(temps, tempr, 2, mpi_double_precision, mpi_sum, &
      mpi_comm_2d_cart, ierr_wt)
    arear_ = tempr(1); volr_ = tempr(2)
    
    volflux_(1) = volr_; volflux_(2) = arear_; volflux_(3) = volr_ / arear_

    if (myid.eq.0) print *, 'VolumeFluxInfo: vol, area, velocity = ', volflux_
  end subroutine calc_vol_flux

end module turbine_model
