program ex1d

  use gabriel
  use MPI_F08
  implicit none

  integer,parameter :: n=10

  logical :: correct = .true.
!  real,dimension(:,:,:),allocatable :: a,b
  real,dimension(:),allocatable :: a,b
  real,dimension(:,:),allocatable :: c
  integer ierr,rank,right,left,mpisize,i,j,k

  type(parcel) :: h(10,2)

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

  print*,'Define subarray parcels..'
  call h(1,1)%subarray(a,(/n,n/),(/n,n/),err=ierr)
  print*,'Errorcode: ',ierr
  if (ierr.ne.6) correct=.false.
  call MPI_Finalize(ierr)
  if (.not.correct) stop 1

end

