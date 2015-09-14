!***********************************************************************
module halos
!***********************************************************************
! contains info for the message passaging                               
!-----------------------------------------------------------------------
      use mpi, only : MPI_ADDRESS_KIND
      implicit none 

      private

!      logical, parameter :: debug=.true.
!      logical, parameter :: check=.true.
!      integer, parameter :: verbose=1
      logical, parameter :: debug=.false.
      logical, parameter :: check=.true.
      integer, parameter :: verbose=1

!> Class to define halos
      type, public :: halo
        integer, private :: m                         !< MPI type
        logical, private :: initialized = .false.              !< is m a valid MPI type?
        class(halo), pointer, private :: i            !< pointer to extended halo type for verification (boundary checking)
        contains
          procedure, public :: mpitype                        !< Return the MPI type
          procedure, private :: copy_halo                      !< Copy a halo
          generic :: assignment(=) => copy_halo       !< Define assignment of a halo
          procedure :: combined => create_combined                !< Combine different halos of the same variable
          procedure :: joined => create_joined                  !< Combine the same halo of different variables
          procedure :: subarray => create_subarray                !< Create a subarray type
          procedure :: is_valid_halo => check_halo                 !< Verify a halo
          final :: finalize_halo                      !< Finalize a halo
          procedure :: print => print_halo            !< print halo
      end type halo

!> Type to describe subarrays
      type, extends(halo) :: subarray
        integer                 :: ndim               !< dimension of array
        integer                 :: oldtype            !< MPI type of elements of array
        integer,dimension(:),allocatable :: lb                 !< lower bounds for array
        integer,dimension(:),allocatable :: ub                 !< upper bounds for array
        integer,dimension(:),allocatable :: sizes              !< sizes of array
        integer,dimension(:),allocatable :: subsizes           !< sizes of subarray
        integer,dimension(:),allocatable :: starts             !< starting indexes of subarray
        contains
          procedure :: is_valid_halo => check_subarray         !< Verify a subarray
          procedure :: print => print_subarray                 !< print a subarray
      end type

!> Type to join multiple halos of the same variable
      type, extends(halo) :: combined
        integer :: n                                           !< number of halos
        type(halo),dimension(:),allocatable :: halos           !< array of halos
        contains
          procedure :: is_valid_halo => check_combined         !< verify a combined halo
      end type

!> Type to join the same halo of multiple variables
      type, extends(halo) :: joined(length)
        integer,len :: length                                  !< max number of variables
        integer     :: n = 0                                   !< actual number of variables
        type(halo)  :: h                                       !< halo for each variable
        integer(MPI_ADDRESS_KIND),dimension(length) :: variables !< multiple variables
      end type

!> Type to define decomposition
      type, public :: decomposition
        private
        integer :: comm_parent             !< parent communicator
        integer :: comm                    !< communicator
        integer :: sends=0                 !< actual number of neighbors to send to
        integer :: recvs=0                 !< actual number of neighbors to receive from
        integer :: maxsends                !< max number of neighbors to send to
        integer :: maxrecvs                !< max number of neighbors to receive from
        integer,allocatable,dimension(:) :: sendranks !< ranks of neighbors to send to
        integer,allocatable,dimension(:) :: recvranks !< ranks of neighbors to send to
        integer,allocatable,dimension(:) :: sendcnts  !< type count to send
        integer,allocatable,dimension(:) :: recvcnts  !< type count to receive
        integer,allocatable,dimension(:) :: sendweights !< weight of sends 
        integer,allocatable,dimension(:) :: recvweights !< weight of recvs
        integer(MPI_ADDRESS_KIND),allocatable,dimension(:) :: senddispls !< displacement of send types
	integer(MPI_ADDRESS_KIND),allocatable,dimension(:) :: recvdispls !< displacement of receive types
        type(halo),allocatable,dimension(:) :: sendhalos  !< halo types to send
        type(halo),allocatable,dimension(:) :: recvhalos  !< halo types to receive
        contains
          procedure :: init => init_decomposition_        !< initialize decomposition
          procedure :: add_send => add_decomposition_send_  !< add a neighbor to send to
          procedure :: add_recv => add_decomposition_recv_  !< add a neighbor to recv from
          procedure :: create => create_decomposition_      !< finalize the decomposition
          procedure :: update => update_decomposition_      !< update the decomposition
      end type decomposition                                                                        

      contains 

!> Function to return the MPI type
      function mpitype(self)
        class(halo), intent(in) :: self
        integer :: mpitype

        mpitype=self%m
      end function mpitype

!> Create type to combine different halos.
!! This is usually used to communicate different parts of the same variable
      subroutine create_subarray(self,array,starts,stops,subsizes)
        use mpi
        integer, parameter           :: ndim=3

        class(halo)                                                :: self
        real, intent(in), dimension(:,:,:), allocatable            :: array
        integer, intent(in), dimension(ndim)                       :: starts
        integer, intent(in), dimension(ndim), optional             :: stops,subsizes

        integer                 :: mpierr

! check arguments
        if (present(stops).and.present(subsizes)) then
          print*,'Error: both stops and subsizes defined. Please define only one.'
          call MPI_Abort(MPI_COMM_WORLD,4,mpierr)
        elseif (.not.(present(stops).or.present(subsizes))) then
          print*,'Error: Please define stops or subsizes.'
          call MPI_Abort(MPI_COMM_WORLD,5,mpierr)
        endif


! allocate subarray halo type
        allocate(subarray :: self%i)

        select type(sub => self%i)
        type is (subarray)
! initialize type
          allocate(sub%lb(ndim))
          allocate(sub%ub(ndim))
          allocate(sub%sizes(ndim))
          allocate(sub%subsizes(ndim))
          allocate(sub%starts(ndim))
          sub%lb=lbound(array)
          sub%ub=ubound(array)
          sub%sizes=shape(array)

! check bounds
          if (any(starts.lt.sub%lb)) then
            print*,'Error: starts lower than lower bound of array'
            print*,'   Starts     :',starts
            print*,'   Lower bound:',sub%lb
            call MPI_Abort(MPI_COMM_WORLD,6,mpierr)
          endif
          sub%starts=starts-sub%lb

          if (present(stops)) then
            if (any(stops.gt.sub%ub)) then
              print*,'Error: stops higher than upper bound of array'
              print*,'   Stops      :',stops
              print*,'   Upper bound:',sub%ub
              call MPI_Abort(MPI_COMM_WORLD,7,mpierr)
            endif
            sub%subsizes=stops-starts+1
          endif

          if (present(subsizes)) then
            if (any(starts+subsizes-1.gt.sub%ub)) then
              print*,'Error: starts+subsizes higher than upper bound of array'
              print*,'   Starts+subsizes:',starts+subsizes
              print*,'   Upper bound:',sub%ub
              call MPI_Abort(MPI_COMM_WORLD,8,mpierr)
            endif
            sub%subsizes=subsizes
          endif

! create type
          call MPI_Type_create_subarray(ndim,sub%sizes,sub%subsizes,sub%starts,    &
     &       MPI_ORDER_FORTRAN,MPI_REAL,self%m,mpierr)
! note that the mpitype is stored in the halo type that points to the subarray type, NOT in the subarray type.
        end select ! sub => h%i
        call MPI_Type_commit(self%m,mpierr)
        self%initialized=.true.
      end subroutine create_subarray

!> Create type to combine different halos.
!! This is usually used to communicate different parts of the same variable
      subroutine create_combined(self,halos)
        use mpi

        class(halo)              :: self
        type(halo),dimension(:)  :: halos

        integer                  :: i,sz,mpierr
        integer,allocatable,dimension(:) :: hh,ones
        integer(kind=MPI_ADDRESS_KIND),allocatable,dimension(:) :: zeroes

        sz=size(halos)
! check arguments
        if (sz.le.1) then
          print*,'Error: Combined halo needs at least 2 halos as input.'
          call MPI_Abort(MPI_COMM_WORLD,9,mpierr)
        endif

! allocate subarray halo type
        allocate(combined :: self%i)
        select type (c=> self%i)
        type is (combined)
          c%n=sz
          allocate(c%halos(sz))
          allocate(hh(sz),ones(sz),zeroes(sz))
          do i=1,sz
            hh(i)=halos(i)%m
          enddo
          ones=1
          zeroes=0_MPI_ADDRESS_KIND
          c%halos=halos
        end select

! create type
          call MPI_Type_create_struct(sz,ones,zeroes,hh,self%m,mpierr)
! note that the mpitype is stored in the halo type that points to the subarray type, NOT in the subarray type.
        call MPI_Type_commit(self%m,mpierr)
        self%initialized=.true.
        deallocate(hh)
      end subroutine create_combined

!> Create type to join the same halo of multiple variables.
!! This is usually used to communicate the same part of different variables.
!! @param n number of variables to communicate the halos
! The resulting halo is based on absolute addresses, so it will be communicated
! with MPI_BOTTOM as both sending and receiving buffers.
      subroutine create_joined(self,n)
        use mpi

        class(halo)               :: self
        integer,intent(in)        :: n

        type(halo)                :: h
        integer                   :: i,mpierr

! check arguments
        if (n.lt.2) then
          print*,'Error: Joined halo needs at least room for 2 halo variables as input.'
          call MPI_Abort(MPI_COMM_WORLD,10,mpierr)
        endif

! allocate subarray halo type
        allocate(joined(n) :: h%i)

! create type
!          call MPI_Type_create_struct(sz,ones,zeroes,hh,h%m)
! note that the mpitype is stored in the halo type that points to the subarray type, NOT in the subarray type.
!        call MPI_Type_commit(h%m,mpierr)

        self=h
      end subroutine create_joined

!> check validity of a halo for a variable
      logical function check_halo(self,v)
        class(halo),intent(in) :: self
        real,dimension(:,:,:),allocatable,intent(in) :: v

        if(debug)print*,'check_halo'

        check_halo=.false.
        select type(h=>self%i)
        type is (subarray)
          check_halo=h%is_valid_halo(v)
        type is (combined)
          check_halo=h%is_valid_halo(v)
        class default
          if (h%initialized) check_halo=.true.
        end select

      end function check_halo

      subroutine print_halo(self)
        class(halo),intent(in) :: self

        if(debug)print*,'print_halo'
        select type(h=>self%i)
        type is (subarray)
          call h%print
        class default
          if (h%initialized) print*,'is initialized. mpitype=',h%m
        end select

      end subroutine print_halo

      subroutine print_subarray(self)
        class(subarray),intent(in) :: self

        if(debug)print*,'print_subarray'
        print*,'ndim=',self%ndim

      end subroutine print_subarray


!> add a variable to a joined halo type
      subroutine add_joined(self,v)
        use mpi

        class(halo), intent(inout) :: self
        real, dimension(:,:,:), allocatable, intent(in)    :: v
                                                                        
        integer mpierr,n 

      select type(j=>self%i)
      type is (joined)
        if (.not.j%is_valid_halo(v)) then
          print*,"Halo is not valid for variable" 
          call MPI_Abort(MPI_COMM_WORLD,11,mpierr) 
        endif        

        j%n=j%n+1
        if (j%n.gt.j%length) then 
          print*,"Too many halos for joined halo type" 
          call MPI_Abort(MPI_COMM_WORLD,12,mpierr) 
        endif 
        call MPI_Get_address(v,j%variables(n),mpierr) 
      class default
          print*,"This is not a joined halo type"
          call MPI_Abort(MPI_COMM_WORLD,13,mpierr) 
      end select
                                                                        
      end subroutine add_joined

      subroutine finalize_halo(h)
        type(halo) :: h

        integer mpierr

        print*,'Finalizing halo ',h%m
        if (associated(h%i)) then
          deallocate(h%i)
          nullify(h%i)
        endif
        if (h%initialized) then
          call MPI_Type_free(h%m,mpierr)
          h%initialized=.false.
        endif

      end subroutine finalize_halo

      subroutine copy_halo(hout,hin)
        class(halo),intent(inout) :: hout
        class(halo),intent(in)    :: hin

        integer mpierr

        if(debug)print*,'copy_halo ',hin%initialized
        if (hin%initialized) then
          call MPI_Type_dup(hin%m,hout%m,mpierr)
          allocate(hout%i,source=hin%i)
        endif
      end subroutine copy_halo

      subroutine init_decomposition_(d,sends,recvs,comm)
        class(decomposition), intent(out)                  :: d
        integer, intent(in)                                :: sends,recvs,comm

        d%comm_parent=comm
        d%maxsends=sends
        d%maxrecvs=recvs
        d%sends=0
        d%recvs=0
        allocate(d%sendranks(sends))
        allocate(d%recvranks(recvs))
        allocate(d%sendhalos(sends))
        allocate(d%recvhalos(recvs))
        allocate(d%sendcnts(sends))
        allocate(d%recvcnts(recvs))
        allocate(d%senddispls(sends))
        allocate(d%recvdispls(recvs))
        allocate(d%sendweights(sends))
        allocate(d%recvweights(recvs))

      end subroutine init_decomposition_

      subroutine add_decomposition_send_(d,rank,h)
        use mpi
        class(decomposition), intent(inout)                 :: d
        type(halo), intent(in)                             :: h
        integer, intent(in)                                :: rank

        integer :: n,mpierr

        n=d%sends+1
        if (n.gt.d%maxsends) then
          print*,"Exceeded maximum number of sending neighbours!"
          call MPI_Abort(MPI_COMM_WORLD,1,mpierr)
        endif

        d%sendranks(n)=rank
        d%sendhalos(n)=h
        d%senddispls(n)=0
        d%sendcnts(n)=1
        d%sendweights(n)=1
        d%sends=n

      end subroutine add_decomposition_send_

      subroutine add_decomposition_recv_(d,rank,h)
        use mpi
        class(decomposition), intent(inout)                 :: d
        type(halo), intent(in)                             :: h
        integer, intent(in)                                :: rank

        integer n,mpierr

        n=d%recvs+1
        if (n.gt.d%maxrecvs) then
          print*,"Exceeded maximum number of receiving neighbours!"
          call MPI_Abort(MPI_COMM_WORLD,2,mpierr)
        endif
        
        d%recvranks(n)=rank
        d%recvhalos(n)=h
        d%recvdispls(n)=0
        d%recvcnts(n)=1
        d%recvweights(n)=1
        d%recvs=n

      end subroutine add_decomposition_recv_

      subroutine create_decomposition_(d,i,reorder)
        use mpi
        class(decomposition), intent(inout)            :: d
        integer, optional, intent(in)                 :: i
        logical, optional, intent(in)                 :: reorder

        integer info,ierr
        logical mpi_reorder

        if (present(i)) then
          info=i
        else
          info=MPI_INFO_NULL
        endif

        if (present(reorder)) then
          mpi_reorder=reorder
        else
          mpi_reorder=.true.
        endif

        call MPI_Dist_graph_create_adjacent(d%comm_parent,d%recvs,d%recvranks,d%recvweights, &
     &           d%sends,d%sendranks,d%sendweights,info,mpi_reorder,d%comm,ierr)
        if (ierr.ne.MPI_SUCCESS) then
          print*,'MPI_Dist_graph_create_adjacent error: ',ierr
          call MPI_Abort(MPI_COMM_WORLD,3,ierr)
        endif

      end subroutine create_decomposition_

      subroutine update_decomposition_(self,vsend,vrecv)
      use mpi

      real, dimension(:,:,:), allocatable, intent(in)    :: vsend
      real, dimension(:,:,:), allocatable, intent(inout) :: vrecv
      class(decomposition), intent(in)                    :: self

      integer mpierr,status(MPI_STATUS_SIZE)
      integer i

      if(debug)print*,'update_decomposition'
      if (check) then
        do i=1,self%sends
          if(debug)print*,'i=',i
          if (.not.self%sendhalos(i)%is_valid_halo(vsend)) then
            print*,'Not valid '
            call MPI_Abort(MPI_COMM_WORLD,33,mpierr)
          endif 
        enddo
        do i=1,self%recvs
          if (.not.self%recvhalos(i)%is_valid_halo(vrecv)) then
            print*,'Not valid '
            call MPI_Abort(MPI_COMM_WORLD,34,mpierr)
          endif 
        enddo
      endif
      
      call MPI_Neighbor_alltoallw(vsend,self%sendcnts,self%senddispls,self%sendhalos%m, &
     &  vrecv,self%recvcnts,self%recvdispls,self%recvhalos%m,self%comm,mpierr)
                   
      end subroutine update_decomposition_
                                                                        
      logical function check_subarray(self,v)
        use mpi 
        integer, parameter           :: ndim=3
                                                                        
        class(subarray), intent(in)                         :: self
        real, dimension(:,:,:), allocatable, intent(in)    :: v 
                                                                        
      integer ierr,status(MPI_STATUS_SIZE) 
                                                          ! bounds for s
      integer lb(ndim),ub(ndim)
       
      if (debug) print*,'check_subarray'           
      ierr=.true.      
                                                           
      lb=lbound(v) 
      ub=ubound(v) 
                                                                        
!check bounds of array                                             
      if (.not.(all(lb.eq.self%lb))) then 
         if (verbose.gt.0)print*,"Lower bounds of array not as assumed.." 
         if (verbose.gt.0)print*,"Found:    lb=",lb 
         if (verbose.gt.0)print*,"Expected: lb=",self%lb 
         ierr=.false.
      endif 
      if (.not.(all(ub.eq.self%ub))) then 
         if (verbose.gt.0)print*,"Upper bounds of array not as assumed.."
         if (verbose.gt.0)print*,"Found:    ub=",ub 
         if (verbose.gt.0)print*,"Expected: ub=",self%ub 
         ierr=.false. 
      endif 
!      if (ierr.eq.1) then
!        flush(6)
!        call MPI_Abort(MPI_COMM_WORLD,10,ierr) 
!      endif

        check_subarray=ierr
      end function check_subarray

      logical function check_combined(self,v)
        use mpi 
        integer, parameter           :: ndim=3
                                                                        
        class(combined), intent(in)                         :: self
        real, dimension(:,:,:), allocatable, intent(in)    :: v 
                                                                        
      integer i,ierr,status(MPI_STATUS_SIZE) 
                                                          ! bounds for s
      if (debug) print*,'check_combined'
      ierr=.true.      

      if (debug) print*,'n=',self%n                                                          
      do i=1,self%n
        if(.not.self%halos(i)%is_valid_halo(v)) ierr=.false.
      enddo

        check_combined=ierr
      end function check_combined
                                                                        
      subroutine create_joined_halo(self) 
        use mpi 
                                                                        
        class(halo),intent(inout) :: self

        integer mpierr
                                                                        
        select type (j=>self%i)
        type is (joined)                                                                        
          call MPI_Type_create_hindexed_block(j%n,         & 
                  1,j%variables,j%m,self%m,mpierr) 
        end select
        call MPI_Type_commit(self%m,mpierr) 
                                                                        
      end subroutine create_joined_halo 
                                                                        
      END
