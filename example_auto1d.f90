program ex

  use gabriel
  use MPI_F08
  implicit none

  integer,parameter :: n=10
  integer,parameter :: s=3
  real,parameter :: delta=0.01

  logical :: correct = .true.

  real,dimension(:),allocatable :: a,b
  integer ierr,rank,right,left,mpisize,i,j,k
  integer hor

  type(box) :: c
  type(distribution) :: d
  

  call MPI_Init(ierr)
  print*,'Initialization...'
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mpisize,ierr)
  right=rank+1
  left=rank-1
  if (right.ge.mpisize)right=0
  if (left.lt.0)left=mpisize-1
  print*,'left,rank,right=',left,rank,right

  hor=rank
  allocate(a(hor*s-1:hor*s+s))
  a=rank
  do i=lbound(a,1),ubound(a,1)
    a(i)=a(i)+i*delta
  enddo

  write(*,'(a,4i4)')'Size,Rank,hor=',mpisize,rank,hor

  write(*,'(a,i3,a,5f13.3,2i5)')'BEFORE Rank',rank,' data=',a,lbound(a,1),ubound(a,1)
  

  print*,'Gabriel setup...'
  call gabriel_init
  call c%init(a,(/hor*s/),(/hor*s+s-1/),MPI_COMM_WORLD,periodic=(/.true./))
  call d%halo(c)
  call d%create
!  call d%autocreate(a,(/hor*s/),(/hor*s+s-1/),MPI_COMM_WORLD,periodic=(/.true./))
!  call d%autocreate(a,(/hor*s,ver*s+1,1/),(/hor*s+s-1,ver*s+s,6/),MPI_COMM_WORLD,periodic=(/.true.,.false.,.false./))

  print*,'Update decomposition..'
  call d%update(a)

  write(*,'(a,i4,a,5f13.3)')' AFTER Rank',rank,' data=',a

  if (rank.ne.0 .and. a(hor*s-1).ne.left+(hor*s-1)*delta) correct=.false.
  if (rank.ne.mpisize-1 .and. a(hor*s+s).ne.right+(hor*s+s)*delta) correct=.false.
  if (rank.eq.0 .and. a(hor*s-1).ne.left+(left*s+s-1)*delta) correct=.false.
  if (rank.eq.mpisize-1 .and. a(hor*s+s).ne.right+(right*s)*delta) correct=.false.
  deallocate(a)

  call MPI_Finalize(ierr)
  if (.not.correct) stop 1
end

