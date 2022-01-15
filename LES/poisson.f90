module poisson
   use decomp, only : wp, real_type

   implicit none

   private

   public :: tdma_para
   public :: tridiag_lu, tridiag_lu_solve

contains
   subroutine tdma_para(d,l,u,x,nz,nsys,comm,rank,nproc)
     use MPI

     implicit none

     integer,  intent(IN) :: nsys,nz,comm,rank,nproc
     real(wp), dimension(nsys,nz) :: d,l,u,x

     integer :: k, error
     real(wp), dimension(nsys) :: mp
     real(wp), dimension(nsys,2) :: sendbuf,recvbuf !< send and receive buffer

     ! forward substitution
     if (rank /= 0) then 
        ! upper interface for cpus below the top cpu 
        call MPI_Recv(recvbuf,nsys*2,real_type,rank-1,0,comm,MPI_STATUS_IGNORE,error)
        d(:,1) = d(:,1) - l(:,1)*recvbuf(:,1)
        x(:,1) = x(:,1) - l(:,1)*recvbuf(:,2)
     end if
     ! internal points for cpus above the bottom cpus
     do k = 2, nz
        mp = l(:,k)/d(:,k-1)
        d(:,k) = d(:,k) - u(:,k-1)*mp
        x(:,k) = x(:,k) - x(:,k-1)*mp
     end do
     if (rank /= nproc-1) then
        sendbuf(:,1) = u(:,nz)/d(:,nz)
        sendbuf(:,2) = x(:,nz)/d(:,nz)
        call MPI_Send(sendbuf,nsys*2,real_type,rank+1,0,comm,error)
     end if

     ! backward substitution
     if (rank == nproc-1) then
        x(:,nz) = x(:,nz)/d(:,nz)
     else
        ! all cpus above the bottom should receive a message about solutions
        call MPI_Recv(recvbuf,nsys,real_type,rank+1,1,comm,MPI_STATUS_IGNORE,error)
        x(:,nz) = (x(:,nz) - u(:,nz)*recvbuf(:,1))/d(:,nz)
     end if
     do k = nz-1, 1, -1
        x(:,k) = (x(:,k) - u(:,k)*x(:,k+1))/d(:,k)
     end do
     if (rank /= 0) then
        call MPI_Send(x(:,1),nsys,real_type,rank-1,1,comm,error)
     end if

   end subroutine tdma_para
   
   !-----------------------------------------------------------
   !> @brief Perform LU decomposition of mutiple independent
   !! tridiagonal systems
   !
   !> The input matrix A is tridiagonal with decomposition A=LU.
   !! L is a lower triangle matrix with diagonal elements 1. U is
   !! a upper triangle matrix.
   !! See http://www.cfm.brown.edu/people/gk/chap6/node13.html
   !
   !> The matrix A can be distributed over several processors.
   !! For each tridiagonal system, each processors has nz equations.
   !! The input array d, l, u is of dimension [nsys, nz].
   !
   !> @param[in,out]  d   input main diagonal of A, output main
   !!                     diagonal of U
   !> @param[in,out]  u   input sup-diagonal of A, output super-diagonal
   !!                     of U
   !> @param[in,out]  l   input sub-diagonal of A, output sub-diagonal
   !!                     of L
   !> @param[in]      nz  number of equations distributed on local
   !!                     processor
   !> @param[in]    nsys  number of independent systems to solve
   !> @param[in]    comm  MPI communicator
   !-----------------------------------------------------------
   subroutine tridiag_lu(d,l,u,nz,nsys,comm)
     use MPI

     implicit none

     integer,  intent(IN) :: nsys,nz,comm
     real(wp), dimension(nsys,nz), intent(INOUT) :: d,l,u

     integer :: rank, nproc, k, error
     real(wp), dimension(nsys,2) :: sendbuf,recvbuf !< send and receive buffer

     call MPI_Comm_rank(comm, rank, error)
     call MPI_Comm_size(comm, nproc, error)

     ! LU factorization
     if (rank /= 0) then 
        call MPI_Recv(recvbuf,nsys*2,real_type,rank-1,0,comm,MPI_STATUS_IGNORE,error)
        l(:,1) = l(:,1)/recvbuf(:,1)
        d(:,1) = d(:,1) - l(:,1)*recvbuf(:,2)
     end if
     do k = 2, nz
        l(:,k) = l(:,k)/d(:,k-1)
        d(:,k) = d(:,k) - l(:,k)*u(:,k-1)
     end do
     if (rank /= nproc-1) then
        sendbuf(:,1) = d(:,nz)
        sendbuf(:,2) = u(:,nz)
        call MPI_Send(sendbuf,nsys*2,real_type,rank+1,0,comm,error)
     end if

   end subroutine tridiag_lu

   !-----------------------------------------------------------
   !> @brief Solve mutiple independent tridiagonal LU systems.
   !
   !> The input the LU decomposition of tridiagonal matrix A.
   !! L is a lower triangle matrix with diagonal elements 1. U is
   !! a upper triangle matrix.
   !! See http://www.cfm.brown.edu/people/gk/chap6/node13.html
   !
   !> The matrix L/U can be distributed over several processors.
   !! For each tridiagonal system, each processors has nz equations.
   !! The input array d, l, u is of dimension [nsys, nz].
   !
   !> @param[in]       d  main diagonal of U
   !> @param[in]       u  super-diagonal of U
   !> @param[in]       l  output sub-diagonal of L
   !> @param[in,out]   x  rhs and solution
   !> @param[in]      nz  number of equations distributed on local
   !!                     processor
   !> @param[in]    nsys  number of independent systems to solve
   !> @param[in]    comm  MPI communicator
   !-----------------------------------------------------------
   subroutine tridiag_lu_solve(d,l,u,x,nz,nsys,comm)
     use MPI

     implicit none

     integer,  intent(IN) :: nsys,nz,comm
     real(wp), dimension(nsys,nz), intent(IN) :: d,l,u
     real(wp), dimension(nsys,nz), intent(INOUT) :: x

     integer :: rank, nproc, k, error
     real(wp), dimension(nsys) :: recvbuf !< receive buffer

     call MPI_Comm_rank(comm, rank, error)
     call MPI_Comm_size(comm, nproc, error)

     ! forward substitution
     if (rank /= 0) then 
        call MPI_Recv(recvbuf,nsys,real_type,rank-1,0,comm,MPI_STATUS_IGNORE,error)
        x(:,1) = x(:,1) - l(:,1)*recvbuf
     end if
     do k = 2, nz
        x(:,k) = x(:,k) - l(:,k)*x(:,k-1)
     end do
     if (rank /= nproc-1) then
        call MPI_Send(x(:,nz),nsys,real_type,rank+1,0,comm,error)
     end if

     ! backward substitution
     if (rank == nproc-1) then
        x(:,nz) = x(:,nz)/d(:,nz)
     else
        call MPI_Recv(recvbuf,nsys,real_type,rank+1,1,comm,MPI_STATUS_IGNORE,error)
        x(:,nz) = (x(:,nz) - u(:,nz)*recvbuf)/d(:,nz)
     end if
     do k = nz-1, 1, -1
        x(:,k) = (x(:,k) - u(:,k)*x(:,k+1))/d(:,k)
     end do
     if (rank /= 0) then
        call MPI_Send(x(:,1),nsys,real_type,rank-1,1,comm,error)
     end if

   end subroutine tridiag_lu_solve

end module poisson
