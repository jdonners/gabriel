program ex

  use gabriel
  use MPI_F08
  implicit none

  integer,parameter :: n=10
  integer,parameter :: s=2

  real,dimension(:,:),allocatable :: a,b
  integer ierr,rank,right,left,mpisize,i,j,k
  integer hor,ver
  logical :: correct=.true.

  type(distribution) :: d
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mpisize,ierr)
  right=rank+1
  left=rank-1
  if (right.ge.mpisize)right=0
  if (left.lt.0)left=mpisize-1
  print*,'left,rank,right=',left,rank,right

  hor=mod(rank,3)
  ver=rank/3
  allocate(a(hor*s-1:hor*s+s,ver*s:ver*s+s+1))
  a=rank
  do i=lbound(a,1),ubound(a,1)
  do j=lbound(a,2),ubound(a,2)
    a(i,j)=a(i,j)+i*0.1+j*0.01
  enddo
  enddo

  write(*,'(a,4i3)')'Size,Rank,hor,ver=',mpisize,rank,hor,ver

  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',ver*s,a(:,ver*s),lbound(a,1),ubound(a,1)
  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',ver*s+1,a(:,ver*s+1),lbound(a,1),ubound(a,1)
  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',ver*s+2,a(:,ver*s+2),lbound(a,1),ubound(a,1)
  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',ver*s+3,a(:,ver*s+3),lbound(a,1),ubound(a,1)

  call gabriel_init
  call d%autocreate(a,(/hor*s,ver*s+1/),(/hor*s+s-1,ver*s+s/),MPI_COMM_WORLD,periodic=(/.true.,.true./))
!  call d%autocreate(a,(/hor*s,ver*s+1,1/),(/hor*s+s-1,ver*s+s,6/),MPI_COMM_WORLD,periodic=(/.true.,.false.,.false./))

  print*,'Update decomposition..'
  call d%update(a,a)

  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s,a(:,ver*s)
  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s+1,a(:,ver*s+1)
  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s+2,a(:,ver*s+2)
  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s+3,a(:,ver*s+3)

  do i=lbound(a,1),ubound(a,1)
  do j=lbound(a,2),ubound(a,2)
    if (a(i,j).ne.value(i,j,s,mpisize)) then
      correct=.false.
      print*,'rank,i,j,a,value=',rank,i,j,a(i,j),value(i,j,s,mpisize)
    endif
  enddo
  enddo
  
  deallocate(a)

  call MPI_Finalize(ierr)

  if (.not.correct) stop 1

contains

  function value(i,j,s,mpisize)
    integer, intent(in) ::  i,j,s,mpisize
    real :: value
  
    integer rank
    integer hor,ver

    hor=i/s
!    if (i.lt.0) hor=min(mpisize-1,2)
    ver=(j-1)/s
    rank=ver*3+hor
    
    value=rank+i*0.1+j*0.01

  end function value
end
