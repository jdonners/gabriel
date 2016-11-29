program ex

  use gabriel
  use MPI
  implicit none

  integer,parameter :: n=10
  integer,parameter :: s=2

  logical :: correct=.true.
  
  real,dimension(:),allocatable :: a,b
  integer ierr,rank,right,left,mpisize,i,j,k
  integer sz,r,oi,oo,e           ! offset in, offset out
  character :: szchar,aszchar
  type(distribution) :: d
  type(box) :: b1, b2
  
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mpisize,ierr)
  if (mpisize.gt.9)call MPI_Abort(MPI_COMM_WORLD,1,ierr)

  right=rank+1
  left=rank-1
  if (right.ge.mpisize)right=0
  if (left.lt.0)left=mpisize-1
  print*,'left,rank,right=',left,rank,right

  sz=mpisize*(mpisize+1)/2
  szchar=char(ichar('0')+sz)

  print*,sz,'=',szchar
  r=rank+1
  oi=0
  oo=0
  do i=1,rank
    oi=oi+i
    oo=oo+mpisize-i+1
  enddo
  allocate(a(1:r))
  allocate(b(1:mpisize-r+1))
  call b1%init(a,lbound(a),ubound(a),MPI_COMM_WORLD,offset=(/oi/),err=e)
  call b2%init(b,lbound(b),ubound(b),MPI_COMM_WORLD,offset=(/oo/),err=e)
  do i=1,r
    a(i)=i+oi
  enddo
  print*,'shape(a),shape(b)=',shape(a),shape(b)
  print*,'lb(a),lb(b)=',lbound(a),lbound(b)
  print*,'ub(a),ub(b)=',ubound(a),ubound(b)
  b=0.0

  write(*,'(a,4i3)')'Size,Rank=',mpisize,rank

  aszchar=char(ichar('0')+size(a,1))
  write(*,'(a,i3,a,'//aszchar//'f8.3)')'BEFORE Rank',rank,' data=',a(:)

  print*,'Create compositions..'
  print*,'Create transform..'
  call d%transform(b1,b2,err=e)

  print*,'Transform decomposition..'
  call d%update(a,b,err=e)

!  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',b(1,:)
  aszchar=char(ichar('0')+size(b,1))
  write(*,'(a,i3,a,'//aszchar//'f5.1)')' AFTER Rank',rank,' data=',b(:)
  
!  do i=lbound(b,1),ubound(b,1)
!  do j=lbound(b,2),ubound(b,2)
!    if (b(i,j).ne.i*0.1+j*0.01) correct=.false.
!  enddo
!  enddo
  
  deallocate(a,b)

  call MPI_Finalize(ierr)
  if (.not.correct) stop 1

end

