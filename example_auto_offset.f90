program ex

  use gabriel
  use MPI

  integer,parameter :: n=10
  integer,parameter :: s=2

  real,dimension(:,:,:),allocatable :: a,b
  integer ierr,rank,right,left,mpisize,i,j,k
  integer hor,ver

  type(decomposition) :: d
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mpisize,ierr)
  right=rank+1
  left=rank-1
  if (right.ge.mpisize)right=0
  if (left.lt.0)left=mpisize-1
  print*,'left,rank,right=',left,rank,right
  
  hor=mod(rank,2)
  ver=rank/2
  allocate(a(0:s+1,0:s+1,6))
  a=rank

  write(*,'(a,4i3)')'Size,Rank,hor,ver=',mpisize,rank,hor,ver

  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',0,a(:,0,1),lbound(a,1),ubound(a,1)
  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',1,a(:,1,1),lbound(a,1),ubound(a,1)
  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',2,a(:,2,1),lbound(a,1),ubound(a,1)
  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',3,a(:,3,1),lbound(a,1),ubound(a,1)

  call d%autocreate(a,(/1,1,1/),(/s,s,6/),MPI_COMM_WORLD,offset=(/hor*s,ver*s,0/))

  print*,'Update decomposition..'
  call d%update(a,a)

  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',0,a(:,0,1)
  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',1,a(:,1,1)
  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',2,a(:,2,1)
  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',3,a(:,3,1)
  deallocate(a)

  call MPI_Finalize(ierr)

end

