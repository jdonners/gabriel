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

  real,dimension(:,:,:),allocatable :: a
  integer ierr,rank,right,left,mpisize,i,j,k

  type(halo) :: h(2)
  type(decomposition) :: d
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

  write(*,'(a,i3,a,10f13.3)')'BEFORE Rank',rank,' data=',a(:,1,1)

  print*,'Define halos..'
  call h(1)%subarray(a,(/n,1,1/),(/n,5,6/))
  call h(2)%subarray(a,(/0,1,1/),(/0,5,6/))

  print*,'Define decomposition..'
  call d%init(1,1,MPI_COMM_WORLD)
  call d%add_send(right,h(1))
  call d%add_recv(left,h(2))
  call d%create

  print*,'Update decomposition..'
  call d%update(a,a)

  write(*,'(a,i3,a,10f13.3)')'AFTER  Rank',rank,' data=',a(:,1,1)
  deallocate(a)

  call MPI_Finalize(ierr)

end

