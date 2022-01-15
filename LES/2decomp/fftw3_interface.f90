  integer(C_INT), parameter :: FFTW_MEASURE = 0
  integer(C_INT), parameter :: FFTW_DESTROY_INPUT = 1
  integer(C_INT), parameter :: FFTW_UNALIGNED = 2
  integer(C_INT), parameter :: FFTW_EXHAUSTIVE = 8
  integer(C_INT), parameter :: FFTW_PRESERVE_INPUT = 16
  integer(C_INT), parameter :: FFTW_PATIENT = 32
  integer(C_INT), parameter :: FFTW_ESTIMATE = 64
  integer(C_INT), parameter :: FFTW_ESTIMATE_PATIENT = 128

  interface
     type(C_PTR) function fftw_plan_many_dft_r2c(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,flags) &
                   bind(C, name='fftw_plan_many_dft_r2c')
       import
       integer(C_INT), value :: rank
       integer(C_INT), dimension(*), intent(in) :: n
       integer(C_INT), value :: howmany
       real(C_DOUBLE), dimension(*), intent(out) :: in
       integer(C_INT), dimension(*), intent(in) :: inembed
       integer(C_INT), value :: istride
       integer(C_INT), value :: idist
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
       integer(C_INT), dimension(*), intent(in) :: onembed
       integer(C_INT), value :: ostride
       integer(C_INT), value :: odist
       integer(C_INT), value :: flags
     end function fftw_plan_many_dft_r2c

     type(C_PTR) function fftw_plan_many_dft_c2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,flags) &
                   bind(C, name='fftw_plan_many_dft_c2r')
       import
       integer(C_INT), value :: rank
       integer(C_INT), dimension(*), intent(in) :: n
       integer(C_INT), value :: howmany
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
       integer(C_INT), dimension(*), intent(in) :: inembed
       integer(C_INT), value :: istride
       integer(C_INT), value :: idist
       real(C_DOUBLE), dimension(*), intent(out) :: out
       integer(C_INT), dimension(*), intent(in) :: onembed
       integer(C_INT), value :: ostride
       integer(C_INT), value :: odist
       integer(C_INT), value :: flags
     end function fftw_plan_many_dft_c2r

     subroutine fftw_destroy_plan(p) bind(C, name='fftw_destroy_plan')
       import
       type(C_PTR), value :: p
     end subroutine fftw_destroy_plan

     subroutine my_fftw_execute_dft_r2c(p,in,out) bind(C, name='fftw_execute_dft_r2c')
       import
       type(C_PTR), value :: p
       real(C_DOUBLE), dimension(*), intent(in) :: in
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
     end subroutine my_fftw_execute_dft_r2c

     subroutine my_fftw_execute_dft_c2r(p,in,out) bind(C, name='fftw_execute_dft_c2r')
       import
       type(C_PTR), value :: p
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
       real(C_DOUBLE), dimension(*), intent(out) :: out
     end subroutine my_fftw_execute_dft_c2r
  end interface

