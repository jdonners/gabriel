program ex

  use gabriel
  use MPI_F08
  implicit none

  integer,parameter :: n=5

  real,dimension(:,:,:),allocatable :: a,a2,b
  integer ierr,rank,right,left,mpisize,i,j,k

  type(parcel) :: h(10,2)
  type(distribution) :: d
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mpisize,ierr)
  right=rank+1
  left=rank-1
  if (right.ge.mpisize)right=0
  if (left.lt.0)left=mpisize-1
  print*,'left,rank,right=',left,rank,right
  
  allocate(a(0:n+1,5,6))
  a=rank
  allocate(a2(0:n+1,5,6))
  a2=rank+100

  write(*,'(a,i3,a,12f13.3)')'BEFORE Rank',rank,' data=',a(:,1,1)

  call gabriel_init
  print*,'Define subarray parcels..'
  call h(1,1)%subarray(a,(/n,1,1/),(/n,5,6/))
  call h(1,2)%subarray(a,(/0,1,1/),(/0,5,6/))
  write(*,'(a)',advance='no')'Define joined parcels..'
  call h(1,1)%joined(2)
  call h(1,2)%joined(2)
  write(*,'(a)',advance='no')'1'
  call h(1,1)%add_joined(a)
  write(*,'(a)',advance='no')'2'
  call h(1,1)%add_joined(a2)
  write(*,'(a)',advance='no')'3'
  call h(1,1)%commit
  call h(1,2)%add_joined(a)
  call h(1,2)%add_joined(a2)
  call h(1,2)%commit
  write(*,'(a)',advance='yes')'4'

  print*,'Initialize distribution..'
  call d%init(1,1,MPI_COMM_WORLD)
  print*,'Add send..'
  call d%add_send(right,h(1,1))
  print*,'Add recv..'
  call d%add_recv(left,h(1,2))
  print*,'Create distribution..'
  call d%create

!  print*,'Update decomposition..'
!  print*,'Result should be valid'
!  if (h(1	,1)%is_valid_parcel(a)) then
!    print*,' Valid'
!  else
!    print*,'NOT valid'
!  endif

  print*,'Apply distribution..'
  call d%update

  write(*,'(a,i3,a,12f13.3)')'AFTER  Rank',rank,' a=',a(:,1,1)
  write(*,'(a,i3,a,12f13.3)')'AFTER  Rank',rank,' a2=',a2(:,1,1)
  deallocate(a,a2)

  call MPI_Finalize(ierr)

end

