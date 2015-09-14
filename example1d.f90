function string_and_integer(s,i)
  character(*) :: s
  integer i
  string_and_integer="a"
end function
function integer_and_string(i,s)
  character(*) :: s
  integer i
  integer_and_string="b"
end function

program ex

  use halos
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
  call h(1,1)%subarray(a,(/n/),(/n/))
  call h(1,2)%subarray(a,(/0/),(/0/))
  call h(2,1)%subarray(a,(/n-1/),(/n-1/))
  call h(2,2)%subarray(a,(/1/),(/1/))
  print*,'Define combined halos..'
  call h(3,1)%combined(h(1:2,1))
  call h(3,2)%combined(h(1:2,2))

  print*,'Initialize decomposition..'
  call d%init(1,1,MPI_COMM_WORLD)
  print*,'Add send..'
  call d%add_send(right,h(3,1))
  print*,'Add recv..'
  call d%add_recv(left,h(3,2))
  print*,'Create decomposition..'
  call d%create

  print*,'Check to see if validity check works'
  print*,'Result should be NOT valid'
  if (h(3,1)%is_valid_halo(c)) then
    print*,' Valid'
  else
    print*,'NOT valid'
  endif

  print*,'Update decomposition..'
  print*,'Result should be valid'
  if (h(3,1)%is_valid_halo(a)) then
    print*,' Valid'
  else
    print*,'NOT valid'
  endif

  print*,'Update decomposition..'
  call d%update(a,a)

  write(*,'(a,i3,a,12f13.3)')'AFTER  Rank',rank,' data=',a(:)
  deallocate(a)

  call MPI_Finalize(ierr)

end

