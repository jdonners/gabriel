program ex

  use halos
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
  allocate(a(hor*s-1:hor*s+s,ver*s:ver*s+s+1,0:7))
  a=rank
  do i=lbound(a,1),ubound(a,1)
  do j=lbound(a,2),ubound(a,2)
  do k=lbound(a,3),ubound(a,3)
    a(i,j,k)=a(i,j,k)+i*0.1+j*0.01+k*0.001
  enddo
  enddo
  enddo

  write(*,'(a,4i3)')'Size,Rank,hor,ver=',mpisize,rank,hor,ver

  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',ver*s,a(:,ver*s,1),lbound(a,1),ubound(a,1)
  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',ver*s+1,a(:,ver*s+1,1),lbound(a,1),ubound(a,1)
  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',ver*s+2,a(:,ver*s+2,1),lbound(a,1),ubound(a,1)
  write(*,'(a,i3,a,i4,4f13.3,2i5)')'BEFORE Rank',rank,' data=',ver*s+3,a(:,ver*s+3,1),lbound(a,1),ubound(a,1)

  call d%autocreate(a,(/hor*s,ver*s+1,1/),(/hor*s+s-1,ver*s+s,6/),MPI_COMM_WORLD,periodic=(/.true.,.true.,.true./))
!  call d%autocreate(a,(/hor*s,ver*s+1,1/),(/hor*s+s-1,ver*s+s,6/),MPI_COMM_WORLD,periodic=(/.true.,.false.,.false./))

  print*,'Update decomposition..'
  call d%update(a,a)

  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s,a(:,ver*s,1)
  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s+1,a(:,ver*s+1,1)
  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s+2,a(:,ver*s+2,1)
  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s+3,a(:,ver*s+3,1)
  deallocate(a)

  call MPI_Finalize(ierr)

end

