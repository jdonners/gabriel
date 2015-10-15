program ex

  use gabriel
  use MPI

  integer,parameter :: n=10
  integer,parameter :: s=2

  real,dimension(:,:),allocatable :: a,b
  integer ierr,rank,right,left,mpisize,i,j,k
  integer sz,r

  type(decomposition) :: d
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mpisize,ierr)
  right=rank+1
  left=rank-1
  if (right.ge.mpisize)right=0
  if (left.lt.0)left=mpisize-1
  print*,'left,rank,right=',left,rank,right

  sz=mpisize*(mpisize+1)/2
  r=rank+1
  allocate(a((r-1)*r/2+1:(r-1)*r/2+r,sz))
  allocate(b(sz,(r-1)*r/2+1:(r-1)*r/2+r))
  print*,'shape(a),shape(b)=',shape(a),shape(b)
  print*,'lb(a),lb(b)=',lbound(a),lbound(b)
  print*,'ub(a),ub(b)=',ubound(a),ubound(b)
  a=r
  do i=lbound(a,1),ubound(a,1)
  do j=lbound(a,2),ubound(a,2)
    a(i,j)=a(i,j)+i*0.1+j*0.01
  enddo
  enddo

  write(*,'(a,4i3)')'Size,Rank=',mpisize,rank

!  write(*,'(a,i3,a,4f13.3,2i5)')'BEFORE Rank',rank,' data=',a(:,1),lbound(a,1),ubound(a,1)
 print*,'BEFORE Rank',rank,' data=',a(lbound(a,1),:),lbound(a,2),ubound(a,2)

  print*,'Create transform..'
  call d%transform(a,b,(/lbound(a,1),lbound(a,2)/),(/ubound(a,1),ubound(a,2)/),MPI_COMM_WORLD)

  print*,'Transform decomposition..'
  call d%update(a,b)

!  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',b(1,:)
  print*,' AFTER Rank',rank,' data=',b(:,lbound(b,2))
  deallocate(a,b)

  call MPI_Finalize(ierr)

end

