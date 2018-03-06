program ex

  use gabriel
  use MPI_F08
  implicit none

  integer,parameter :: n=10
  integer,parameter :: s=2

  logical :: correct=.true.
  
  real,dimension(:,:),allocatable :: a,b
  integer ierr,rank,right,left,mpisize,i,j,k
  integer sz,r
  character :: szchar,aszchar
  type(distribution) :: d
  type(box) :: b1, b2
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mpisize,ierr)
  if (mpisize.gt.3)call MPI_Abort(MPI_COMM_WORLD,1,ierr)

  right=rank+1
  left=rank-1
  if (right.ge.mpisize)right=0
  if (left.lt.0)left=mpisize-1
  print*,'left,rank,right=',left,rank,right

  sz=mpisize*(mpisize+1)/2
  szchar=char(ichar('0')+sz)
  
  print*,sz,'=',szchar
  r=rank+1
  allocate(a((r-1)*r/2+1:(r+1)*r/2,sz))
  allocate(b(sz,(r-1)*r/2+1:(r+1)*r/2))
  aszchar=char(ichar('0')+size(a,1))
  print*,'shape(a),shape(b)=',shape(a),shape(b)
  print*,'lb(a),lb(b)=',lbound(a),lbound(b)
  print*,'ub(a),ub(b)=',ubound(a),ubound(b)
  a=0.0
  b=0.0
  do i=lbound(a,1),ubound(a,1)
  do j=lbound(a,2),ubound(a,2)
    a(i,j)=a(i,j)+i*0.1+j*0.01
  enddo
  enddo

  write(*,'(a,4i3)')'Size,Rank=',mpisize,rank

  write(*,'(a,i3,a,'//aszchar//'f8.3,2i5)')'BEFORE Rank',rank,' data=',a(:,1),lbound(a,1),ubound(a,1)
 print*,'BEFORE Rank',rank,' data=',a(lbound(a,1),:),lbound(a,2),ubound(a,2)

  print*,'Create compositions..'
  call b1%init(a,lbound(a),ubound(a),MPI_COMM_WORLD)
  call b2%init(b,lbound(b),ubound(b),MPI_COMM_WORLD)
  print*,'Create transform..'
  call d%transform(b1,b2)

  print*,'Transform decomposition..'
  call d%update(a,b)

!  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',b(1,:)
  write(*,'(a,i3,a,'//szchar//'f8.3)')' AFTER Rank',rank,' data=',b(:,lbound(b,2))
  
  do i=lbound(b,1),ubound(b,1)
  do j=lbound(b,2),ubound(b,2)
    if (b(i,j).ne.i*0.1+j*0.01) correct=.false.
  enddo
  enddo
  
  deallocate(a,b)

  call MPI_Finalize(ierr)
  if (.not.correct) stop 1

end

