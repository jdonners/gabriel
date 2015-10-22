program ex1d

  use gabriel
  use MPI

  interface operator(//)
    function string_and_integer(s,i)
      character(:),allocatable :: string_and_integer
      character(*),intent(in) :: s
      integer,intent(in) :: i
    end function
  end interface

  integer,parameter :: n=10

!  real,dimension(:,:,:),allocatable :: a,b
  real,dimension(:),allocatable :: a,b
  real,dimension(:,:),allocatable :: c
  integer ierr,rank,right,left,mpisize,i,j,k

  type(halo) :: h(10,2)
  type(decomposition) :: d
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mpisize,ierr)
  right=rank+1
  left=rank-1
  if (right.ge.mpisize)right=0
  if (left.lt.0)left=mpisize-1
  print*,'left,rank,right=',left,rank,right
  
  allocate(a(0:n+1))
  a=rank

  write(*,'(a,i3,a,12f13.3)')'BEFORE Rank',rank,' data=',a(:)

  print*,'Define subarray halos..'
  call h(1,1)%subarray(a,(/n,n/),(/n,n/),err=ierr)
  print*,'Errorcode: ',ierr

  call MPI_Finalize(ierr)

end

