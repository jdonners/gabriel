!***********************************************************************
module halos
!***********************************************************************
! contains info for the message passing
!-----------------------------------------------------------------------
      use mpi, only : MPI_ADDRESS_KIND
      implicit none 

      private

!      logical, parameter :: debug=.true.
!      logical, parameter :: check=.true.
!      integer, parameter :: verbose=1
      logical, parameter :: check=.true.

! Verbosity levels
! 0 : no verbosity
! 1 : only errors
! 2 : add warnings
! 3 : add info
! 4 : add debug 
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
          procedure :: print => print_halo            !< print halo
          final :: finalize_halo                      !< Finalize a halo
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
      type, extends(halo) :: joined
        integer     :: length                                  !< max number of variables
        integer     :: n = 0                                   !< actual number of variables
        type(halo)  :: h                                       !< halo for each variable
        integer(MPI_ADDRESS_KIND),dimension(:),allocatable :: variables !< multiple variables
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
          procedure :: autocreate => create_decomposition_halo_      !< create autodecomposition
      end type decomposition          


      contains 

      logical function isdebug()

        isdebug=.false.
        if (verbose.eq.4) isdebug=.true.

      end function

      subroutine debug(s)
        character(len=*),intent(in) :: s

        if (isdebug()) print*,s
      end subroutine

!> Create type to combine different halos.
!! This is used to communicate a subarray of a larger array
!! @relates halos::halo
      subroutine create_subarray(self,array,starts,stops,subsizes)
        use mpi

        class(halo)                                             :: self
        real, intent(in), dimension(..), allocatable            :: array
        integer, intent(in), dimension(:)                       :: starts
        integer, intent(in), dimension(:), optional             :: stops,subsizes

        integer                 :: mpierr
        integer           :: ndim

! check arguments
        if (present(stops).and.present(subsizes)) then
          print*,'Error: both stops and subsizes defined. Please define only one.'
          call MPI_Abort(MPI_COMM_WORLD,4,mpierr)
        elseif (.not.(present(stops).or.present(subsizes))) then
          print*,'Error: Please define stops or subsizes.'
          call MPI_Abort(MPI_COMM_WORLD,5,mpierr)
        endif

        ndim=rank(array)

! allocate subarray halo type
        allocate(subarray :: self%i)

        select type(sub => self%i)
        type is (subarray)
! initialize type
          sub%ndim=ndim
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
!! @relates halos::halo
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
!! @param self halo object
!! @param n number of variables to communicate the halos
!! @relates halos::halo
      subroutine create_joined(self,n)
! The resulting halo is based on absolute addresses, so it will be communicated
! with MPI_BOTTOM as both sending and receiving buffers.
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
        allocate(joined :: h%i)
        select type (j=>h%i)
        type is (joined)
          allocate(j%variables(n))
        end select

! create type
!          call MPI_Type_create_struct(sz,ones,zeroes,hh,h%m)
! note that the mpitype is stored in the halo type that points to the subarray type, NOT in the subarray type.
!        call MPI_Type_commit(h%m,mpierr)

        self=h
      end subroutine create_joined

!> print a halo
!! @relates halos::halo
      subroutine print_halo(self)
        class(halo),intent(in) :: self

        call debug('print_halo')
        select type(h=>self%i)
        type is (subarray)
          call h%print
        class default
          if (h%initialized) print*,'is initialized. mpitype=',h%m
        end select

      end subroutine print_halo

!> print a subarray
!! @relates halos::subarray
!! @private
      subroutine print_subarray(self)
        class(subarray),intent(in) :: self

        call debug('print_subarray')
        print*,'ndim=',self%ndim

      end subroutine print_subarray


!> add a variable to a joined halo type
!! @relates halos::subarray
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

!> Initialize a decomposition
!! @relates halos::decomposition
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

!> This is meant to trick doxygen in documenting the next routine
!!@ private
      subroutine dummy1
        use mpi
      end subroutine dummy1

!> Add a receive to a decomposition
!! @relates halos::decomposition
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

!> Add a send to a decomposition
!! @relates halos::decomposition
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

!> Commit a decomposition
!! @relates halos::decomposition
      subroutine create_decomposition_(d,i,reorder)
        use mpi
        class(decomposition), intent(inout)           :: d            !> decomposition type
        integer, optional, intent(in)                 :: i            !> MPI_Info, default MPI_INFO_NULL
        logical, optional, intent(in)                 :: reorder      !> reorder, default .true.

        integer info,ierr
        integer commsize
        logical mpi_reorder
        integer cnt,n

        type(halo) h

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

        call MPI_Comm_size(d%comm_parent,commsize,ierr)
        do n=0,commsize-1
          
          if (count(d%recvranks.eq.n).gt.1) then
            if(isdebug())print*,'n,d%recvs,count(d%recvranks.eq.n)=',n,d%recvs,count(d%recvranks.eq.n)
!combine receive halos from the same rank into one combined type
            call h%combined(pack(d%recvhalos(1:d%recvs),d%recvranks(1:d%recvs).eq.n))
            cnt=count(d%recvranks(1:d%recvs).ne.n)
            d%recvhalos(1:cnt)=pack(d%recvhalos,d%recvranks(1:d%recvs).ne.n)
            d%recvdispls(1:cnt)=pack(d%recvdispls,d%recvranks(1:d%recvs).ne.n)
            d%recvcnts(1:cnt)=pack(d%recvcnts,d%recvranks(1:d%recvs).ne.n)
            d%recvweights(1:cnt)=pack(d%recvweights,d%recvranks(1:d%recvs).ne.n)
            d%recvranks(1:cnt)=pack(d%recvranks,d%recvranks(1:d%recvs).ne.n)
            d%recvs=cnt
!add combined type
            call d%add_recv(n,h)
          endif
          if (count(d%sendranks.eq.n).gt.1) then
            if(isdebug())print*,'n,d%sends,count(d%sendranks.eq.n)=',n,d%sends,count(d%sendranks.eq.n)
!combine send halos to the same rank into one combined type
            call h%combined(pack(d%sendhalos(1:d%sends),d%sendranks(1:d%sends).eq.n))
!remove separate types from send arrays
            cnt=count(d%sendranks(1:d%sends).ne.n)
            d%sendhalos(1:cnt)=pack(d%sendhalos,d%sendranks(1:d%sends).ne.n)
            d%senddispls(1:cnt)=pack(d%senddispls,d%sendranks(1:d%sends).ne.n)
            d%sendcnts(1:cnt)=pack(d%sendcnts,d%sendranks(1:d%sends).ne.n)
            d%sendweights(1:cnt)=pack(d%sendweights,d%sendranks(1:d%sends).ne.n)
            d%sendranks(1:cnt)=pack(d%sendranks,d%sendranks(1:d%sends).ne.n)
            d%sends=cnt
!add combined type
            call d%add_send(n,h)
          endif
        enddo

        call MPI_Dist_graph_create_adjacent(d%comm_parent,d%recvs,d%recvranks,d%recvweights, &
     &           d%sends,d%sendranks,d%sendweights,info,mpi_reorder,d%comm,ierr)
        if (ierr.ne.MPI_SUCCESS) then
          print*,'MPI_Dist_graph_create_adjacent error: ',ierr
          call MPI_Abort(MPI_COMM_WORLD,3,ierr)
        endif

      end subroutine create_decomposition_

!> Automatically create a decomposition with all halos.
!! This is a collective MPI call.
!! @relates halos::decomposition
!      subroutine create_decomposition_halo_(d,v,lower,upper,comm,offset,periodic,lower_global,upper_global)
      subroutine create_decomposition_halo_(d,v,lower,upper,comm,offset,periodic)
        use mpi
        class(decomposition), intent(inout)            :: d    !> Resulting decomposition
        real, dimension(..), allocatable, intent(in)   :: v    !> variable to create halos for
        integer, dimension(:), intent(in) :: lower             !> lower bound of active domain
        integer, dimension(:), intent(in) :: upper             !> upper bound of active domain
        integer, intent(in)               :: comm              !> communicator
        integer, dimension(:), intent(in), optional :: offset  !> offset of array indices
        logical, dimension(:), intent(in), optional :: periodic       !> periodicity of global domain
!        integer, dimension(:), intent(in), optional :: global_lower   !> lower bound of global domain
!        integer, dimension(:), intent(in), optional :: global_upper   !> upper bound of global domain

        integer, parameter            :: MAX_HALOS = 30
        integer :: sendcount=0
        integer :: recvcount=0
        integer,dimension(MAX_HALOS) :: sends
        integer,dimension(MAX_HALOS) :: recvs
        integer :: r
        integer,dimension(:),allocatable :: lb,ub,off
        integer,dimension(:),allocatable :: shf
        integer,dimension(:),allocatable :: low,up
        integer,dimension(:,:),allocatable :: lshalo,ushalo
        integer,dimension(:,:),allocatable :: lrhalo,urhalo
        integer,dimension(:,:),allocatable :: lowers,uppers
        integer,dimension(:,:),allocatable :: lbs,ubs
        integer,dimension(:),allocatable :: global_low,global_up
        logical,dimension(:),allocatable :: per
        integer :: commsize,mpierr,commrank
        integer :: i
        type(halo) :: h

        r=rank(v)
        if (size(lower).ne.r) call error("Size of lower bound array not equal to rank!",12)
        if (size(upper).ne.r) call error("Size of upper bound array not equal to rank!",13)

        allocate(off(r))
        off=0
        if (present(offset)) then
          if (size(offset).ne.r) call error("Size of offset array not equal to rank!",16)
          off=offset
        endif

        allocate(lb(r),ub(r))
        allocate(low(r),up(r))

        lb=lbound(v)+off
        ub=ubound(v)+off
        low=lower+off
        up=upper+off

        if (any(low.lt.lb)) call error("Lower bound array incorrect!",14)
        if (any(up.gt.ub)) call error("Upper bound array incorrect!",15)

! Possibly a check with a warning to see if the active domain equals the variable bounds, i.e. no halo regions

        call MPIcheck
        call MPI_Comm_Size(comm,commsize,mpierr)
        call MPI_Comm_Rank(comm,commrank,mpierr)

        allocate(lowers(r,commsize),uppers(r,commsize))
        allocate(lbs(r,commsize),ubs(r,commsize))

!Of course, this could be implemented more efficiently..
        call MPI_Allgather(low,r,MPI_INTEGER,lowers,r,MPI_INTEGER,comm,mpierr)
        call MPI_Allgather(up,r,MPI_INTEGER,uppers,r,MPI_INTEGER,comm,mpierr)
        call MPI_Allgather(lb,r,MPI_INTEGER,lbs,r,MPI_INTEGER,comm,mpierr)
        call MPI_Allgather(ub,r,MPI_INTEGER,ubs,r,MPI_INTEGER,comm,mpierr)

        allocate(lshalo(r,MAX_HALOS),ushalo(r,MAX_HALOS))
        allocate(lrhalo(r,MAX_HALOS),urhalo(r,MAX_HALOS))
        do i=1,commsize
        if (i-1.ne.commrank) then  
! check for overlap of active domains, if so, error!
          if (all(up.ge.lowers(:,i)).and.all(low.le.uppers(:,i))) then
             print*,'upper=',up
             print*,'uppers=',uppers(:,i)
             print*,'lower=',low
             print*,'lowers=',lowers(:,i)
             call error("Overlap of active domains!")
          endif
! check for overlap of my active domain with other domains
          if (all(up.ge.lbs(:,i)).and.all(low.le.ubs(:,i))) then
! send overlapping data from my active domain
            sendcount=sendcount+1
            if (sendcount.gt.MAX_HALOS) call error("Too many sends!")
            sends(sendcount)=i-1
            lshalo(:,sendcount)=max(lbs(:,i),low)-off
            ushalo(:,sendcount)=min(ubs(:,i),up)-off
          endif
! check for overlap of my full domain with other active domains
          if (all(ub.ge.lowers(:,i)).and.all(lb.le.uppers(:,i))) then
! receive overlapping data from other active domain
            recvcount=recvcount+1
            if (recvcount.gt.MAX_HALOS) call error("Too many receives!")
            recvs(recvcount)=i-1
            lrhalo(:,recvcount)=max(lb,lowers(:,i))-off
            urhalo(:,recvcount)=min(ub,uppers(:,i))-off
          endif
        endif
        enddo

! check for periodicity arguments
! if global bounds are used: check if defined global upper and lower bounds really are the maximum boundary values
!        if (present(periodic).and.present(global_lower).and.present(global_upper)) then
        if (present(periodic)) then
          allocate(global_low(r),global_up(r))
          global_low=minval(lowers,dim=2)
          global_up=maxval(uppers,dim=2)
          if (verbose.eq.1.and.commrank.eq.0) then
            print*,'Requested periodicity: ',periodic
            print*,'Global lower bounds: ',global_low
            print*,'Global upper bounds: ',global_up
          endif
          allocate(per(r))
          per=periodic
          allocate(shf(r))
          do
          shf=merge(global_up-global_low+1,0,per)
          do
          if (isdebug())print*,'shf=',shf
          do i=1,commsize
! check for overlap of my active domain with other domains
            if (all(up.ge.lbs(:,i)+shf).and.all(low.le.ubs(:,i)+shf)) then
! send overlapping data from my active domain
              sendcount=sendcount+1
              if (sendcount.gt.MAX_HALOS) call error("Too many sends!")
              sends(sendcount)=i-1
              lshalo(:,sendcount)=max(lbs(:,i)+shf,low)-off
              ushalo(:,sendcount)=min(ubs(:,i)+shf,up)-off
              if(isdebug())write(*,'(a,8i4)')'per1,lshalo,ushalo=',i-1,commrank,lshalo(:,sendcount),ushalo(:,sendcount)
            endif
! check for overlap of my full domain with other active domains
            if (all(ub.ge.lowers(:,i)-shf).and.all(lb.le.uppers(:,i)-shf)) then
! receive overlapping data from other active domain
              recvcount=recvcount+1
              if (recvcount.gt.MAX_HALOS) call error("Too many receives!")
              recvs(recvcount)=i-1
              lrhalo(:,recvcount)=max(lb,lowers(:,i)-shf)-off
              urhalo(:,recvcount)=min(ub,uppers(:,i)-shf)-off
              if(isdebug())write(*,'(a,8i4)')'per4,lrhalo,urhalo=',i-1,commrank,lrhalo(:,recvcount),urhalo(:,recvcount)
            endif
          enddo
          if (.not.signs(shf)) exit
          enddo ! do signed
          if (.not.combo(periodic,per)) exit
          enddo ! do combo
          deallocate(global_low,global_up,shf,per)
        endif

        call d%init(sendcount,recvcount,comm)
        do i=1,recvcount
          call h%subarray(v,lrhalo(:,i),urhalo(:,i))
          call d%add_recv(recvs(i),h)
        enddo
        do i=1,sendcount
          call h%subarray(v,lshalo(:,i),ushalo(:,i))
          call d%add_send(recvs(i),h)
        enddo

        deallocate(lshalo,ushalo,lrhalo,urhalo)
        deallocate(lb,ub,off)
        deallocate(low,up)
        deallocate(lbs,ubs)
        deallocate(lowers,uppers)
        call d%create

      end subroutine create_decomposition_halo_

logical recursive function combo(c,d,n) result(combor)
! this function makes a switch in the logical array.
! c is the original array
! d is the work array
! n is an internal argument used to remember the position in the array.
! returns .false. if no switch is made (to the calling routine this is the moment to stop)
! returns .true. if a switch is made
  logical, dimension(:),intent(in) :: c
  logical, dimension(:),intent(inout) :: d
  integer, optional :: n

  integer i

  if (present(n)) then
    i=n
  else
    i=1
  endif

  if (i.gt.size(c)) then
    combor=.false.
  elseif (combo(c,d,i+1)) then 
    combor=.true.
  elseif (d(i)) then
    d(i)=.false.
    d(i+1:size(d))=c(i+1:size(c))
    combor=.true.
  else
    combor=.false.
  endif

  if (count(d).eq.0) combor=.false.

  return
  end function

logical recursive function signs(d,n) result(signsr)
! this function creates all combinations of a signed array.
! d is the work array, initially with either 0's or 1's.
! n is an internal argument used to remember the position in the array.
! returns .false. if no switch is made (to the calling routine this is the moment to stop)
! returns .true. if a switch is made
  integer, dimension(:),intent(inout) :: d
  integer, optional :: n

  integer i
  logical firstsign

  if (present(n)) then
    i=n
  else
    i=1
  endif

  firstsign=all(d(1:i-1).eq.0)

  if (d(i).ge.1.and.firstsign) then 
    d=-d
    signsr=.true.
  elseif (d(i).le.-1.and.firstsign) then 
    d=-d
    signsr=signs(d,i+1)
  elseif (i.gt.size(d)) then
    signsr=.false.
  elseif (signs(d,i+1)) then 
    signsr=.true.
  elseif (d(i).gt.0) then
    d(i)=-d(i)
    where (d(i+1:size(d)).lt.0) d(i+1:size(d))=-d(i+1:size(d))
    signsr=.true.
  else
    signsr=.false.
  endif


  return
  end function

      subroutine MPIcheck
        use mpi
        logical initialized,finalized
        integer ierr

        call MPI_Initialized(initialized,ierr)
        if (.not.initialized) call error ('MPI library not yet initialized')
        call MPI_Finalized(finalized,ierr)
        if (finalized) call error ('MPI library already finalized')

      end subroutine
      
!> Update a decomposition
!! @relates halos::decomposition
      subroutine update_decomposition_(self,vsend,vrecv)
      use mpi
      use iso_c_binding, only : c_loc,c_f_pointer

      real, dimension(..), allocatable, target, intent(in)    :: vsend
      real, dimension(..), allocatable, target, intent(inout) :: vrecv
      class(decomposition), intent(in)                    :: self

      real,dimension(:),pointer :: psend,precv
      integer mpierr,status(MPI_STATUS_SIZE)
      integer i

      call debug('update_decomposition')
      if (check) then
        do i=1,self%sends
          if (isdebug()) print*,'i=',i
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
      
      call c_f_pointer(c_loc(vsend),psend,(/size(vsend)/))
      call c_f_pointer(c_loc(vrecv),precv,(/size(vrecv)/))
      call MPI_Neighbor_alltoallw(psend,self%sendcnts,self%senddispls,self%sendhalos%m, &
     &  precv,self%recvcnts,self%recvdispls,self%recvhalos%m,self%comm,mpierr)
                   
      end subroutine update_decomposition_

      subroutine error(s,e)
        use mpi 
        character(len=*),intent(in) :: s
        integer, intent(in), optional :: e

        integer mpierr

        print*,s
        if (present(e)) then        
          call MPI_Abort(MPI_COMM_WORLD,e,mpierr)
        else
          call MPI_Abort(MPI_COMM_WORLD,99,mpierr)
        endif

      end subroutine error

!> finalize halo
!! @private
      subroutine finalize_halo(h)
        type(halo) :: h

        integer mpierr

        if (verbose.ge.2) print*,'Finalizing halo ',h%m
        if (associated(h%i)) then
          deallocate(h%i)
          nullify(h%i)
        endif
        if (h%initialized) then
          call MPI_Type_free(h%m,mpierr)
          h%initialized=.false.
        endif

      end subroutine finalize_halo

!> copy halo
!! @private
      subroutine copy_halo(hout,hin)
        class(halo),intent(inout) :: hout
        class(halo),intent(in)    :: hin

        integer mpierr

        if (isdebug()) print*,'copy_halo ',hin%initialized
        if (hin%initialized) then
          call MPI_Type_dup(hin%m,hout%m,mpierr)
          allocate(hout%i,source=hin%i)
        endif
      end subroutine copy_halo

!> check subarray halo
!! @private
      function check_subarray(self,v)
        use mpi 
        logical :: check_subarray
                                                                        
        class(subarray), intent(in)                     :: self
        real, dimension(..), allocatable, intent(in)    :: v 
                                                                        
      integer           :: ndim
      integer status(MPI_STATUS_SIZE) 
      logical ierr

      integer,dimension(:),allocatable :: lb,ub
       
      call debug('check_subarray')
      ierr=.true.      
                                                           
      ndim=rank(v)

!check rank of array
      if (.not.(ndim.eq.self%ndim)) then
         if (verbose.gt.0)print*,"Rank of array not as assumed.." 
         if (verbose.gt.0)print*,"Found:    rank=",ndim
         if (verbose.gt.0)print*,"Expected: rank=",self%ndim
         ierr=.false.
         return
      endif 

!check contiguity of array
!      if (.not.is_contiguous(v)) then
!         if (verbose.gt.0)print*,"Array is not contiguous.."
!         ierr=.false.
!         return
!      endif

      allocate(lb(ndim),ub(ndim))
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

!> check subarray halo
!! @private
      function check_combined(self,v)
        use mpi 
        logical  :: check_combined
        integer, parameter           :: ndim=3
                                                                        
        class(combined), intent(in)                         :: self
        real, dimension(..), allocatable, intent(in)    :: v 
                                                                        
      integer i,status(MPI_STATUS_SIZE) 
      logical ierr
                                                          ! bounds for s
      call debug('check_combined')
      ierr=.true.      

      if (isdebug()) print*,'n=',self%n                                                          
      do i=1,self%n
        if(.not.self%halos(i)%is_valid_halo(v)) ierr=.false.
      enddo

        check_combined=ierr
      end function check_combined
                                                                        
!> create and commit a joined halo
!! @relates halos::halo
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

!> check validity of a halo for a variable
!! @relates halos::halo
      function check_halo(self,v)
        logical :: check_halo
        class(halo),intent(in) :: self
        real,dimension(..),allocatable,intent(in) :: v

        call debug('check_halo')

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

!> Function to return the MPI type
!! @relates halos::halo
!! @public
      function mpitype(self)
        class(halo), intent(in) :: self
        integer :: mpitype

        mpitype=self%m
      end function mpitype

end module halos
