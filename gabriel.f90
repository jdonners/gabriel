!***********************************************************************
module gabriel
!***********************************************************************
! contains info for the message passing
!-----------------------------------------------------------------------
      use mpi, only : MPI_ADDRESS_KIND,MPI_DATATYPE_NULL
      implicit none 

      private

!      logical, parameter :: debug=.true.
!      logical, parameter :: check=.true.
!      integer, parameter :: verbose=1
!      logical, parameter :: check=.true.

! Verbosity levels
! 0 : no verbosity
! 1 : only errors
! 2 : add warnings
! 3 : add info
! 4 : add debug 
!
!dox doxygen comments, not used now	 (this should be replaced with !!, should doxygen be used)

      public :: gabriel_init

      integer :: verbose=1
      public :: gabriel_set_verbosity

      logical :: check=.true.
      public :: gabriel_set_checking, gabriel_enable_checking, gabriel_disable_checking

      integer :: realtype = MPI_DATATYPE_NULL

!> Class to define halos
      type, public :: halo
        integer, private :: m=MPI_DATATYPE_NULL                   !< MPI type
        logical, private :: initialized = .false.              !< is m a valid MPI type?
        logical, private :: absolute = .false.        !< is m an absolute MPI type (i.e. should it be used with MPI_BOTTOM)?
        class(halo), pointer, private :: i            !< pointer to extended halo type for verification (boundary checking)
        contains
          procedure, public :: mpitype                        !< Return a duplicate of the MPI type
          procedure, private :: copy_halo                      !< Copy a halo
          generic :: assignment(=) => copy_halo       !< Define assignment of a halo
          procedure :: combined => create_combined                !< Combine different halos of the same variable
          procedure :: joined => halo_joined_init                  !< Join the same halo of different variables
          procedure :: add_joined
          procedure :: commit => halo_commit                      !< commit halo datatype and create joined datatype
          procedure :: subarray => create_subarray                !< Create a subarray type
          procedure :: is_valid_halo => check_halo                 !< Verify a halo
          procedure :: is_absolute => halo_is_absolute       !< Is a halo absolute?
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

!> Type to combine multiple halos of the same variable
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
          procedure :: init => init_decomposition_              !< initialize decomposition
          procedure :: add_send => add_decomposition_send_      !< add a neighbor to send to
          procedure :: add_recv => add_decomposition_recv_      !< add a neighbor to recv from
          procedure :: create => create_decomposition_          !< create the graph communicator of the decomposition
! Unfortunately, this generic triggers an Intel compiler bug because of assumed-rank arrays, see topic 595234 in the Intel Forum
!          generic   :: update => decomposition_update_single,decomposition_update_bottom,decomposition_update_sendrecv           !< update the decomposition
! This is a temporary fix to the compiler problem above. This fix resolves the right routine at runtime, not at compile-time,
! and will therefore be slower.
! As a workaround, the user could call the three underlying routines directly.
! These should be private if the compiler works correctly.
          procedure :: update => decomposition_update_alt           !< update the decomposition
          procedure :: autocreate => decomposition_autocreate   !< create autodecomposition
          procedure :: transform => create_reshuffle_           !< create transformation
          procedure :: halo => decomposition_halo               !< setup halos, but don't create
          procedure :: joined => decomposition_joined           !< setup joined decomposition
          procedure :: joined_add => decomposition_joined_add   !< add variable to joined decomposition
          procedure, public :: decomposition_update_single           !< update the decomposition
          procedure, public :: decomposition_update_bottom           !< update the decomposition
          procedure, public :: decomposition_update_sendrecv         !< update the decomposition
      end type decomposition          

      contains 

      subroutine gabriel_init
        character(LEN=256) :: verbosity
        integer :: v

        call get_environment_variable("GABRIEL_VERBOSE",verbosity)
        read(verbosity,'(i)')v
        call gabriel_set_verbosity(v)

      end subroutine

      subroutine gabriel_set_checking(flag)
        logical :: flag                        !< .true. or .false.

        check=flag
      end subroutine

      subroutine gabriel_disable_checking
        check=.false.
      end subroutine

      subroutine gabriel_enable_checking
        check=.true.
      end subroutine

      subroutine gabriel_set_verbosity(level)
        integer :: level

        if (level.ge.0 .and. level.le.4) verbose=level
        if (level.ge.4) print*,'gabriel verbosity: ',level
      end subroutine

      logical function isdebug()

        isdebug=.false.
        if (verbose.ge.4) isdebug=.true.

      end function

      logical function isinfo()

        isinfo=.false.
        if (verbose.ge.3) isinfo=.true.

      end function

      logical function iswarning()

        iswarning=.false.
        if (verbose.ge.2) iswarning=.true.

      end function

      logical function iserror()

        iserror=.false.
        if (verbose.ge.1) iserror=.true.

      end function

      subroutine debug(s)
        character(len=*),intent(in) :: s

        if (isdebug()) print*,s
      end subroutine

      subroutine error(errcode,s,err)
        use mpi
        use, intrinsic :: iso_fortran_env, only : stderr => error_unit

        integer, intent(in)            :: errcode
        character(len=*), optional     :: s
        integer, intent(out), optional :: err

        integer                        :: mpierr

        if (iserror().and.present(s)) write(stderr,*) 'Error: ',s
        if (.not.present(err)) then
          call MPI_Abort(MPI_COMM_WORLD,errcode,mpierr)
        else
          err=errcode
        endif
        
      end subroutine

      logical function isnonzero(err)
        integer,intent(in),optional :: err
        isnonzero=.false.
        if (present(err)) then
          if (err.ne.0) isnonzero=.true.
        endif

      end function

      subroutine create_subarray(self,array,starts,stops,subsizes,err)
!! Create type to combine different halos.
!! This is used to communicate a subarray of a larger array
        use mpi

        class(halo)                                             :: self           !< halo
        real, intent(in), dimension(..), allocatable            :: array          !< input array
        integer, intent(in), dimension(:)                       :: starts         !< starting indices of subarray
        integer, intent(in), dimension(:), optional             :: stops          !< stopping indices of subarray
        integer, intent(in), dimension(:), optional             :: subsizes       !< subsizes of subarray
        integer, intent(out), optional                          :: err            !< error code

        integer           :: mpierr
        integer           :: ndim
        real              :: dummyreal

        if (present(err)) err=0
! check arguments
        if (present(stops).and.present(subsizes)) then
          call error(4,"both stops and subsizes defined. Please define only one.",err)
          return
        elseif (.not.(present(stops).or.present(subsizes))) then
          call error(5,"please define stops or subsizes.",err)
          return
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
            if (iserror()) then
              print*,'Error: starting indices lower than lower bound of array'
              print*,'   Starts     :',starts
              print*,'   Lower bound:',sub%lb
            endif
            call error(6,"Starting indices lower than lower bound of array",err)
            return
          endif
          sub%starts=starts-sub%lb

          if (present(stops)) then
            if (any(stops.gt.sub%ub)) then
              if (iserror()) then
                print*,'Error: stopping indices higher than upper bound of array'
                print*,'   Stops      :',stops
                print*,'   Upper bound:',sub%ub
              endif
              call error(7,"Stopping indices higher than upper bound of array",err)
              return
            endif
            sub%subsizes=stops-starts+1
          endif

          if (present(subsizes)) then
            if (any(starts+subsizes-1.gt.sub%ub)) then
              if (iserror()) then
                print*,'Error: starts+subsizes higher than upper bound of array'
                print*,'   Starts+subsizes:',starts+subsizes
                print*,'   Upper bound:',sub%ub
              endif
              call error(8,"Starts+subsizes higher than upper bound of array",err)
              return
            endif
            sub%subsizes=subsizes
          endif

! create type
          if (realtype.eq.MPI_DATATYPE_NULL) then
            call MPI_Type_create_f90_real(precision(dummyreal),exponent(dummyreal),realtype,mpierr)
          endif
          call MPI_Type_create_subarray(ndim,sub%sizes,sub%subsizes,sub%starts,    &
     &       MPI_ORDER_FORTRAN,realtype,self%m,mpierr)
! note that the mpitype is stored in the halo type that points to the subarray type, NOT in the subarray type.
        end select ! sub => h%i
        call MPI_Type_commit(self%m,mpierr)
        self%initialized=.true.
      end subroutine create_subarray

!> Create type to combine different halos.
!> This is usually used to communicate different parts of the same variable
!dox @relates gabriel::halo
      subroutine create_combined(self,halos,err)
        use mpi

        class(halo),intent(inout)            :: self
        type(halo),dimension(:),intent(in)   :: halos
        integer,intent(out),optional         :: err

        integer                  :: i,sz,mpierr
        integer,allocatable,dimension(:) :: hh,ones
        integer(kind=MPI_ADDRESS_KIND),allocatable,dimension(:) :: zeroes
        integer nabs

        if (present(err)) err=0

        sz=size(halos)
! check arguments
        if (sz.le.1) then
          call error(9,"combined halo needs at least 2 halos as input.",err)
          return
        endif

! check for types with absolute addresses
        nabs=0
        do i=1,sz
          if (halos(i)%is_absolute()) nabs=nabs+1
        enddo
        if (nabs.gt.0 .and. nabs.ne.sz) then
          call error(37,"combined halo can not combine types with relative and absolute addresses.",err)
          return
        endif

! allocate subarray halo type
        allocate(combined :: self%i)
        if (nabs.eq.sz) self%absolute=.true.
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
!> This is usually used to communicate the same part of different variables.
!dox @param self halo object
!dox @param n number of variables to communicate the halos
!dox @relates gabriel::halo
      subroutine halo_joined_init(self,n,err)
! The resulting halo is based on absolute addresses, so it will be communicated
! with MPI_BOTTOM as both sending and receiving buffers.
        use mpi

        class(halo),intent(inout) :: self
        integer,intent(in)        :: n
        integer,intent(out),optional :: err

        class(halo),allocatable   :: h
        integer                   :: i,mpierr

        if (present(err)) err=0
! check arguments
        if (n.lt.2) then
          call error(10,'Error: Joined halo needs at least room for 2 halo variables as input.',err)
          return
        endif

! allocate joined halo type
        allocate(joined :: h)
        select type (h)
        type is (joined)
          allocate(h%variables(n))
          h%length=n
          allocate(h%i,source=self)
        end select
! create type
!          call MPI_Type_create_struct(sz,ones,zeroes,hh,h%m)
! note that the mpitype is stored in the halo type that points to the subarray type, NOT in the subarray type.
!        call MPI_Type_commit(h%m,mpierr)

        allocate(self%i,source=h)
        self%absolute=.true.
        if (isdebug()) then
        select type (h)
        type is (joined)
          print*,'original halo is joined'
        end select
        select type (j=>self%i)
        type is (joined)
          print*,'copied halo is joined'
        class default
          print*,'copied halo is not joined'
        end select
        endif

      end subroutine halo_joined_init

!> create a joined decomposition
!dox @ relates gabriel::decomposition
      subroutine decomposition_joined(self,n,err)
! The resulting halo is based on absolute addresses, so it will be communicated
! with MPI_BOTTOM as both sending and receiving buffers.
        use mpi

        class(decomposition),intent(inout) :: self
        integer,intent(in)        :: n
        integer,intent(out),optional :: err

        type(halo)                :: h
        integer                   :: i,mpierr

        if (present(err)) err=0
! check arguments
        if (n.lt.2) then
          call error(10,'Error: Joined halo needs at least room for 2 halo variables as input.',err)
          return
        endif

        do i=1,self%sends
          call self%sendhalos(i)%joined(n,err)
        enddo
        do i=1,self%recvs
! allocate subarray halo type
          call self%recvhalos(i)%joined(n,err)
        enddo

      end subroutine decomposition_joined

!> add a variable to a joined halo type
!dox @relates gabriel::subarray
      subroutine decomposition_joined_add(self,v,err)
        use mpi

        class(decomposition), intent(inout) :: self
        real, dimension(:,:,:), allocatable, intent(in)    :: v
        integer, intent(out), optional :: err
                                                                        
        integer i

        if (present(err)) err=0
  
        do i=1,self%sends
          call self%sendhalos(i)%add_joined(v,err)
        enddo
        do i=1,self%recvs
          call self%recvhalos(i)%add_joined(v,err)
        enddo

      end subroutine decomposition_joined_add

!> print a halo
!dox @relates gabriel::halo
      subroutine print_halo(self,err)
        class(halo),intent(in) :: self
        integer,intent(out),optional :: err

        if (present(err)) err=0
        call debug('print_halo')
        select type(h=>self%i)
        type is (subarray)
          call h%print(err=err)
        class default
          if (h%initialized) print*,'is initialized. mpitype=',h%m
        end select

      end subroutine print_halo

!> print a subarray
!dox @relates gabriel::subarray
!dox @private
      subroutine print_subarray(self,err)
        class(subarray),intent(in) :: self
        integer,intent(out),optional :: err

        if (present(err)) err=0

        call debug('print_subarray')
        print*,'ndim=',self%ndim

      end subroutine print_subarray


!> Add a variable to a joined halo type
!dox @relates gabriel::subarray
      subroutine add_joined(self,v,err)
        use mpi

        class(halo), intent(inout) :: self
        real, dimension(:,:,:), allocatable, intent(in)    :: v
        integer, intent(out), optional :: err
                                                                        
        integer mpierr,n 

      if (present(err)) err=0
      select type(j=>self%i)
      type is (joined)
        if (.not.j%i%is_valid_halo(v)) then
          call error(11,"Halo is not valid for variable",err)
          return
        endif        

        j%n=j%n+1
        if (isdebug()) print*,'n,length=',j%n,j%length
        if (j%n.gt.j%length) then 
          call error(12,"Too many halos for joined halo type",err)
          return
        endif 
        if (isdebug()) print*,'ln=',size(j%variables),lbound(j%variables)
        call MPI_Get_address(v,j%variables(j%n),mpierr) 
      class default
          call error(13,"This is not a joined halo type",err)
          return
      end select
                                                                        
      end subroutine add_joined

!> Initialize a decomposition
!dox @relates gabriel::decomposition
      subroutine init_decomposition_(d,sends,recvs,comm,err)
        class(decomposition), intent(out)                  :: d
        integer, intent(in)                                :: sends,recvs,comm
        integer, intent(out), optional                     :: err

        integer mpierr

        if (present(err)) err=0

        call MPI_Comm_dup(comm,d%comm_parent,mpierr)
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

!> Add a receive to a decomposition
!dox @relates gabriel::decomposition
      subroutine add_decomposition_recv_(d,rank,h,err)
        use mpi
        class(decomposition), intent(inout)                 :: d
        type(halo), intent(in)                             :: h
        integer, intent(in)                                :: rank
        integer, intent(out), optional                     :: err

        integer n,mpierr

        if (present(err)) err=0
        n=d%recvs+1
        if (n.gt.d%maxrecvs) then
          call error(2,"Exceeded maximum number of receiving neighbours!",err)
          return
        endif
        
        d%recvranks(n)=rank
        d%recvhalos(n)=h
        d%recvdispls(n)=0
        d%recvcnts(n)=1
        d%recvweights(n)=1
        d%recvs=n

      end subroutine add_decomposition_recv_

!> Add a send to a decomposition
!dox @relates gabriel::decomposition
      subroutine add_decomposition_send_(d,rank,h,err)
        use mpi
        class(decomposition), intent(inout)                 :: d
        type(halo), intent(in)                             :: h
        integer, intent(in)                                :: rank
        integer, intent(out), optional                     :: err

        integer :: n,mpierr

        if (present(err)) err=0
        n=d%sends+1
        if (n.gt.d%maxsends) then
          call error(1,"Exceeded maximum number of sending neighbours!",err)
          return
        endif

        d%sendranks(n)=rank
        d%sendhalos(n)=h
        d%senddispls(n)=0
        d%sendcnts(n)=1
        d%sendweights(n)=1
        d%sends=n

      end subroutine add_decomposition_send_

!> Commit a decomposition
!dox @relates gabriel::decomposition
      subroutine create_decomposition_(d,i,reorder,err)
        use mpi
        class(decomposition), intent(inout)           :: d            !< decomposition type
        integer, optional, intent(in)                 :: i            !< MPI_Info, default MPI_INFO_NULL
        logical, optional, intent(in)                 :: reorder      !< reorder, default .true.
        integer, intent(out), optional                 :: err         !< error indicator

        integer info,ierr
        integer commsize
        logical mpi_reorder
        integer cnt,n

        type(halo) h

        if (present(err)) err=0
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

        do n=1,d%recvs
          call d%recvhalos(n)%commit()
        enddo
        do n=1,d%sends
          call d%sendhalos(n)%commit()
        enddo
        do n=0,commsize-1
          if (count(d%recvranks(1:d%recvs).eq.n).gt.1) then
            if(isdebug())print*,'n,d%recvs,count(d%recvranks.eq.n)=',n,d%recvs,count(d%recvranks(1:d%recvs).eq.n)
!combine receive halos from the same rank into one combined type
            call h%combined(pack(d%recvhalos(1:d%recvs),d%recvranks(1:d%recvs).eq.n),err=err)
            if (isnonzero(err))return
            cnt=count(d%recvranks(1:d%recvs).ne.n)
            d%recvhalos(1:cnt)=pack(d%recvhalos(1:d%recvs),d%recvranks(1:d%recvs).ne.n)
            d%recvdispls(1:cnt)=pack(d%recvdispls(1:d%recvs),d%recvranks(1:d%recvs).ne.n)
            d%recvcnts(1:cnt)=pack(d%recvcnts(1:d%recvs),d%recvranks(1:d%recvs).ne.n)
            d%recvweights(1:cnt)=pack(d%recvweights(1:d%recvs),d%recvranks(1:d%recvs).ne.n)
            d%recvranks(1:cnt)=pack(d%recvranks(1:d%recvs),d%recvranks(1:d%recvs).ne.n)
            d%recvs=cnt
!add combined type
            call d%add_recv(n,h,err=err)
            if (isnonzero(err))return
          endif
          if (count(d%sendranks(1:d%sends).eq.n).gt.1) then
            if(isdebug())print*,'n,d%sends,count(d%sendranks.eq.n)=',n,d%sends,count(d%sendranks.eq.n)
!combine send halos to the same rank into one combined type
            call h%combined(pack(d%sendhalos(1:d%sends),d%sendranks(1:d%sends).eq.n),err=err)
            if (isnonzero(err))return
!remove separate types from send arrays
            cnt=count(d%sendranks(1:d%sends).ne.n)
            d%sendhalos(1:cnt)=pack(d%sendhalos(1:d%sends),d%sendranks(1:d%sends).ne.n)
            d%senddispls(1:cnt)=pack(d%senddispls(1:d%sends),d%sendranks(1:d%sends).ne.n)
            d%sendcnts(1:cnt)=pack(d%sendcnts(1:d%sends),d%sendranks(1:d%sends).ne.n)
            d%sendweights(1:cnt)=pack(d%sendweights(1:d%sends),d%sendranks(1:d%sends).ne.n)
            d%sendranks(1:cnt)=pack(d%sendranks(1:d%sends),d%sendranks(1:d%sends).ne.n)
            d%sends=cnt
!add combined type
            call d%add_send(n,h,err=err)
            if (isnonzero(err))return
          endif
        enddo

        call MPI_Dist_graph_create_adjacent(d%comm_parent,d%recvs,d%recvranks,d%recvweights, &
     &           d%sends,d%sendranks,d%sendweights,info,mpi_reorder,d%comm,ierr)
        if (ierr.ne.MPI_SUCCESS) then
          call error(3,"MPI_Dist_graph_create_adjacent error",err)
          return
        endif

      end subroutine create_decomposition_

!> Automatically create a decomposition with all halos.
!> This is a collective MPI call.
!dox @relates gabriel::decomposition
      subroutine decomposition_halo(d,v,lower,upper,comm,offset,periodic,err)
        use mpi
        class(decomposition), intent(inout)            :: d    !< Resulting decomposition
        real, dimension(..), allocatable, intent(in)   :: v    !< variable to create halos for
        integer, dimension(:), intent(in) :: lower             !< lower bound of active domain
        integer, dimension(:), intent(in) :: upper             !< upper bound of active domain
        integer, intent(in)               :: comm              !< communicator
        integer, dimension(:), intent(in), optional :: offset  !< offset of array indices
        logical, dimension(:), intent(in), optional :: periodic       !< periodicity of global domain
        integer, intent(out), optional :: err  !< error indicator
!        integer, dimension(:), intent(in), optional :: global_lower   !< lower bound of global domain
!        integer, dimension(:), intent(in), optional :: global_upper   !< upper bound of global domain

        integer, parameter            :: MAX_HALOS = 30
        integer :: sendcount
        integer :: recvcount
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

        if (present(err)) err=0
        sendcount=0
        recvcount=0

        r=rank(v)
        if (size(lower).ne.r) then
          call error(100,"Size of lower bound array not equal to rank!",err)
          return
        endif
        if (size(upper).ne.r) then
          call error(101,"Size of upper bound array not equal to rank!",err)
          return
        endif
        if (any(lower.gt.upper)) then
          if (isinfo())print*,'lower=',lower
          if (isinfo())print*,'upper=',upper
          call error(102,"Lower bound greater than upper bound!",err)
          return
        endif

        allocate(off(r))
        off=0
        if (present(offset)) then
          if (size(offset).ne.r) then
            call error(103,"Size of offset array not equal to rank!",err)
            return
          endif
          off=offset
        endif

        allocate(lb(r),ub(r))
        allocate(low(r),up(r))

        lb=lbound(v)+off
        ub=ubound(v)+off
        low=lower+off
        up=upper+off

        if (any(low.lt.lb)) then
          call error(104,"Lower bound array incorrect!",err)
          return
        endif
        if (any(up.gt.ub)) then
          call error(105,"Upper bound array incorrect!",err)
          return
        endif

! Possibly a check with a warning to see if the active domain equals the variable bounds, i.e. no halo regions

        call MPIcheck(err)
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
             if (iserror())print*,'upper=',up
             if (iserror())print*,'uppers=',uppers(:,i)
             if (iserror())print*,'lower=',low
             if (iserror())print*,'lowers=',lowers(:,i)
             call error(106,"Overlap of active domains!",err)
             return
          endif
! check for overlap of my active domain with other domains
          if (all(up.ge.lbs(:,i)).and.all(low.le.ubs(:,i))) then
! send overlapping data from my active domain
            sendcount=sendcount+1
             if (isinfo())print*,'rank to send to=',i-1
             if (isinfo())print*,'lower=',max(lbs(:,i),low)-off
             if (isinfo())print*,'upper=',min(ubs(:,i),up)-off
            if (sendcount.gt.MAX_HALOS) then
              call error(107,"Too many sends!",err)
              return
            endif
            sends(sendcount)=i-1
            lshalo(:,sendcount)=max(lbs(:,i),low)-off
            ushalo(:,sendcount)=min(ubs(:,i),up)-off
          endif
! check for overlap of my full domain with other active domains
          if (all(ub.ge.lowers(:,i)).and.all(lb.le.uppers(:,i))) then
! receive overlapping data from other active domain
            recvcount=recvcount+1
             if (isinfo())print*,'rank to receive from=',i-1
             if (isinfo())print*,'lower=',max(lb,lowers(:,i))-off
             if (isinfo())print*,'upper=',min(ub,uppers(:,i))-off
            if (recvcount.gt.MAX_HALOS) then
              call error(108,"Too many receives!",err)
              return
            endif
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
          if (isinfo().and.commrank.eq.0) then
            print*,'Requested periodicity: ',periodic
            print*,'Global lower bounds: ',global_low
            print*,'Global upper bounds: ',global_up
          endif
          allocate(per(r))
          per=periodic
          allocate(shf(r))
          if (any(periodic)) then
          do
          shf=merge(global_up-global_low+1,0,per)
          do
          if (isdebug())print*,'shf=',shf
          do i=1,commsize
! check for overlap of my active domain with other domains
            if (all(up.ge.lbs(:,i)+shf).and.all(low.le.ubs(:,i)+shf)) then
! send overlapping data from my active domain
              sendcount=sendcount+1
              if (sendcount.gt.MAX_HALOS) then
                call error(109,"Too many sends!",err)
                return
              endif
              sends(sendcount)=i-1
              lshalo(:,sendcount)=max(lbs(:,i)+shf,low)-off
              ushalo(:,sendcount)=min(ubs(:,i)+shf,up)-off
              if(isdebug())write(*,'(a,8i7)')'per1,lshalo,ushalo=',i-1,commrank,lshalo(:,sendcount),ushalo(:,sendcount)
            endif
! check for overlap of my full domain with other active domains
            if (all(ub.ge.lowers(:,i)-shf).and.all(lb.le.uppers(:,i)-shf)) then
! receive overlapping data from other active domain
              recvcount=recvcount+1
              if (recvcount.gt.MAX_HALOS) then
                call error(110,"Too many receives!",err)
                return
              endif
              recvs(recvcount)=i-1
              lrhalo(:,recvcount)=max(lb,lowers(:,i)-shf)-off
              urhalo(:,recvcount)=min(ub,uppers(:,i)-shf)-off
              if(isdebug())write(*,'(a,8i7)')'per4,lrhalo,urhalo=',i-1,commrank,lrhalo(:,recvcount),urhalo(:,recvcount)
            endif
          enddo
          if (.not.signs(shf)) exit
          enddo ! do signed
          if (.not.combo(periodic,per)) exit
          enddo ! do combo
          endif
          deallocate(global_low,global_up,shf,per)
        endif

        call d%init(sendcount,recvcount,comm,err=err)
        if (isnonzero(err)) return
        do i=1,recvcount
          call h%subarray(v,lrhalo(:,i),urhalo(:,i),err=err)
          if (isnonzero(err)) return
          call d%add_recv(recvs(i),h,err=err)
          if (isnonzero(err)) return
        enddo
        do i=1,sendcount
          call h%subarray(v,lshalo(:,i),ushalo(:,i),err=err)
          if (isnonzero(err)) return
          call d%add_send(sends(i),h,err=err)
          if (isnonzero(err)) return
        enddo

        deallocate(lshalo,ushalo,lrhalo,urhalo)
        deallocate(lb,ub,off)
        deallocate(low,up)
        deallocate(lbs,ubs)
        deallocate(lowers,uppers)

      end subroutine decomposition_halo

      subroutine decomposition_autocreate(d,v,lower,upper,comm,offset,periodic,err)
        use mpi
        class(decomposition), intent(inout)            :: d    !< Resulting decomposition
        real, dimension(..), allocatable, intent(in)   :: v    !< variable to create halos for
        integer, dimension(:), intent(in) :: lower             !< lower bound of active domain
        integer, dimension(:), intent(in) :: upper             !< upper bound of active domain
        integer, intent(in)               :: comm              !< communicator
        integer, dimension(:), intent(in), optional :: offset  !< offset of array indices
        logical, dimension(:), intent(in), optional :: periodic       !< periodicity of global domain
        integer, intent(out), optional :: err  !< error indicator
!        integer, dimension(:), intent(in), optional :: global_lower   !< lower bound of global domain
!        integer, dimension(:), intent(in), optional :: global_upper   !< upper bound of global domain

        call d%halo(v,lower,upper,comm,offset,periodic,err)
        call d%create(err=err)
      end subroutine decomposition_autocreate
   
!> Automatically create a transformation with all halos.
!> This is a collective MPI call.
!dox @relates gabriel::decomposition
      subroutine create_reshuffle_(d,vfrom,vto,lower,upper,comm,offset,err)
        use mpi
        class(decomposition), intent(inout)            :: d    !< Resulting decomposition
        real, dimension(..), allocatable, intent(in)   :: vfrom    !< variable to create halos for
        real, dimension(..), allocatable, intent(in)   :: vto    !< variable to create halos for
        integer, dimension(:), intent(in) :: lower             !< lower bound of active domain
        integer, dimension(:), intent(in) :: upper             !< upper bound of active domain
        integer, intent(in)               :: comm              !< communicator
        integer, dimension(:), intent(in), optional :: offset  !< offset of array indices
        integer, intent(out), optional :: err  !< error indicator
        integer, parameter            :: MAX_HALOS = 100
        integer :: sendcount
        integer :: recvcount
        integer,dimension(MAX_HALOS) :: sends
        integer,dimension(MAX_HALOS) :: recvs
        integer :: r
        integer,dimension(:),allocatable :: lb,ub,off
        integer,dimension(:),allocatable :: to_lb,to_ub
        integer,dimension(:),allocatable :: low,up
        integer,dimension(:,:),allocatable :: lshalo,ushalo
        integer,dimension(:,:),allocatable :: lrhalo,urhalo
        integer,dimension(:,:),allocatable :: lowers,uppers
        integer,dimension(:,:),allocatable :: lbs,ubs
        integer,dimension(:,:),allocatable :: to_lbs,to_ubs
        integer,dimension(:),allocatable :: global_low,global_up
        logical,dimension(:),allocatable :: per
        integer :: commsize,mpierr,commrank
        integer :: i
        type(halo) :: h

        if (present(err)) err=0
        sendcount=0
        recvcount=0

        call debug("Creating transform..")
        r=rank(vfrom)
        if (rank(vto).ne.r) then
          call error(21,"Ranks of from- and to-arrays do not match!",err)
          return
        endif
        if (size(lower).ne.r) then
          call error(22,"Size of lower bound array not equal to rank!",err)
          return
        endif
        if (size(upper).ne.r) then
          call error(23,"Size of upper bound array not equal to rank!",err)
          return
        endif

        allocate(lb(r),ub(r))
        allocate(low(r),up(r))
        allocate(to_lb(r),to_ub(r))

        lb=lbound(vfrom)
        ub=ubound(vfrom)
        low=lower
        up=upper
        to_lb=lbound(vto)
        to_ub=ubound(vto)

        if (any(low.lt.lb)) then
             if (iserror())print*,'lower=',low
             if (iserror())print*,'lbound=',lb
          call error(24,"Lower bound array incorrect!",err)
          return
        endif
        if (any(up.gt.ub)) then
             if (iserror())print*,'upper=',up
             if (iserror())print*,'ubound=',ub
             call error(25,"Upper bound array incorrect!",err)
          return
        endif

        allocate(off(r))
        off=0
        if (present(offset)) then
          call error(26,"Offset not yet implemented for transform!",err)
          return
          if (size(offset).ne.r) then
            call error(27,"Size of offset array not equal to rank!",err)
            return
          endif
          off=offset
        endif

        call MPI_Comm_Size(comm,commsize,mpierr)
        call MPI_Comm_Rank(comm,commrank,mpierr)

        allocate(lowers(r,commsize),uppers(r,commsize))
        allocate(lbs(r,commsize),ubs(r,commsize))
        allocate(to_lbs(r,commsize),to_ubs(r,commsize))

!Of course, this could be implemented more efficiently..
        call MPI_Allgather(low,r,MPI_INTEGER,lowers,r,MPI_INTEGER,comm,mpierr)
        call MPI_Allgather(up,r,MPI_INTEGER,uppers,r,MPI_INTEGER,comm,mpierr)
        call MPI_Allgather(lb,r,MPI_INTEGER,lbs,r,MPI_INTEGER,comm,mpierr)
        call MPI_Allgather(ub,r,MPI_INTEGER,ubs,r,MPI_INTEGER,comm,mpierr)
        call MPI_Allgather(to_lb,r,MPI_INTEGER,to_lbs,r,MPI_INTEGER,comm,mpierr)
        call MPI_Allgather(to_ub,r,MPI_INTEGER,to_ubs,r,MPI_INTEGER,comm,mpierr)

        allocate(lshalo(r,MAX_HALOS),ushalo(r,MAX_HALOS))
        allocate(lrhalo(r,MAX_HALOS),urhalo(r,MAX_HALOS))
        if(isdebug())print*,'commsize=',commsize
        do i=1,commsize
! check for overlap of active from-domains with other from-domains, if so, error!
          if (i-1.ne.commrank .and. all(up.ge.lowers(:,i)).and.all(low.le.uppers(:,i))) then
             if(isinfo())print*,'upper=',up
             if(isinfo())print*,'uppers=',uppers(:,i)
             if(isinfo())print*,'lower=',low
             if(isinfo())print*,'lowers=',lowers(:,i)
             call error(28,"Overlap of active domains!",err)
             return
          endif
! check for overlap of my active from-domain with other to-domains
          if (all(up.ge.to_lbs(:,i)).and.all(low.le.to_ubs(:,i))) then
! send overlapping data from my active domain
            sendcount=sendcount+1
            if (sendcount.gt.MAX_HALOS) then
              call error(29,"Too many sends!",err)
              return
            endif
            sends(sendcount)=i-1
            lshalo(:,sendcount)=max(to_lbs(:,i),low)-off
            ushalo(:,sendcount)=min(to_ubs(:,i),up)-off
          endif
! check for overlap of my full to-domain with active from-domains
          if (all(to_ub.ge.lowers(:,i)).and.all(to_lb.le.uppers(:,i))) then
! receive overlapping data from other active domain
            recvcount=recvcount+1
            if (recvcount.gt.MAX_HALOS) then
              call error(30,"Too many receives!",err)
              return
            endif
            recvs(recvcount)=i-1
            lrhalo(:,recvcount)=max(to_lb,lowers(:,i))-off
            urhalo(:,recvcount)=min(to_ub,uppers(:,i))-off
          endif
        enddo

        call d%init(sendcount,recvcount,comm,err=err)
        if (isnonzero(err)) return
        do i=1,recvcount
          call h%subarray(vto,lrhalo(:,i),urhalo(:,i),err=err)
          if (isnonzero(err)) return
          call d%add_recv(recvs(i),h,err=err)
          if (isnonzero(err)) return
        enddo
        do i=1,sendcount
          call h%subarray(vfrom,lshalo(:,i),ushalo(:,i),err=err)
          if (isnonzero(err)) return
          call d%add_send(sends(i),h,err=err)
          if (isnonzero(err)) return
        enddo

        deallocate(lshalo,ushalo,lrhalo,urhalo)
        deallocate(lb,ub)
        deallocate(low,up)
        deallocate(lbs,ubs)
        deallocate(to_lb,to_ub,to_lbs,to_ubs)
        deallocate(lowers,uppers)
        call d%create(err=err)

      end subroutine create_reshuffle_

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

      subroutine MPIcheck(err)
        use mpi
        integer, intent(out),optional :: err
        logical initialized,finalized
        integer ierr

        if (present(err)) err=0
        call MPI_Initialized(initialized,ierr)
        if (.not.initialized) then
          call error(31,"MPI library not yet initialized",err)
          return
        endif
        call MPI_Finalized(finalized,ierr)
        if (finalized) then
          call error(32,"MPI library already finalized",err)
          return
        endif

      end subroutine
      
!> Update a decomposition with one array as source and destination
!dox @relates gabriel::decomposition
      subroutine decomposition_update_single(self,v,err)
        use mpi
        use iso_c_binding, only : c_loc,c_f_pointer
 
        real, dimension(..), allocatable, target, intent(inout) :: v
        class(decomposition), intent(in)                    :: self
        integer, intent(out), optional                      :: err

        real,dimension(:),pointer :: psend,precv
        integer mpierr,status(MPI_STATUS_SIZE)
        integer i

        call self%update(v,v,err)

      end subroutine decomposition_update_single

!> Update a decomposition with types that use absolute addressing
!dox @relates gabriel::decomposition
      subroutine decomposition_update_bottom(self,err)
      use mpi
      use iso_c_binding, only : c_loc,c_f_pointer

      class(decomposition), intent(in)                    :: self
      integer, intent(out), optional                      :: err

      integer mpierr,status(MPI_STATUS_SIZE)
      integer i

      if (present(err)) err=0

      call debug('decomposition_update_bottom')
      if (check) then
        do i=1,self%sends
          if (isdebug()) print*,'i=',i
          if (.not.self%sendhalos(i)%is_absolute()) then
            call error(35,"Send halo not absolute",err)
            return
          endif 
        enddo
        do i=1,self%recvs
          if (.not.self%recvhalos(i)%is_absolute()) then
            call error(36,"Receive halo not absolute",err)
            return
          endif 
        enddo
      endif
      
      call MPI_Neighbor_alltoallw(MPI_BOTTOM,self%sendcnts,self%senddispls,self%sendhalos%m, &
     &  MPI_BOTTOM,self%recvcnts,self%recvdispls,self%recvhalos%m,self%comm,mpierr)

      end subroutine decomposition_update_bottom

!> Update a decomposition with separate arrays as source and destination
      subroutine decomposition_update_sendrecv(self,vsend,vrecv,err)
      use mpi
      use iso_c_binding, only : c_loc,c_f_pointer

      class(decomposition), intent(in)                    :: self
      real, dimension(..), allocatable, target, intent(in)    :: vsend
      real, dimension(..), allocatable, target, intent(inout) :: vrecv
      integer, intent(out), optional                      :: err

      real,dimension(:),pointer :: psend,precv
      integer mpierr,status(MPI_STATUS_SIZE)
      integer i

      if (present(err)) err=0

      call debug('update_decomposition')
      if (check) then
        do i=1,self%sends
          if (isdebug()) print*,'i=',i
          if (.not.self%sendhalos(i)%is_valid_halo(vsend)) then
            call error(33,"Send halo not valid",err)
            return
          endif 
        enddo
        do i=1,self%recvs
          if (.not.self%recvhalos(i)%is_valid_halo(vrecv)) then
            call error(34,"Receive halo not valid",err)
            return
          endif 
        enddo
      endif
      
      call c_f_pointer(c_loc(vsend),psend,(/size(vsend)/))
      call c_f_pointer(c_loc(vrecv),precv,(/size(vrecv)/))
      call MPI_Neighbor_alltoallw(psend,self%sendcnts,self%senddispls,self%sendhalos%m, &
     &  precv,self%recvcnts,self%recvdispls,self%recvhalos%m,self%comm,mpierr)
                   
      end subroutine decomposition_update_sendrecv

      subroutine decomposition_update_alt(self,vsend,vrecv,err)
      use mpi
      use iso_c_binding, only : c_loc,c_f_pointer

      class(decomposition), intent(in)                    :: self
      real, dimension(..), allocatable, target, optional, intent(inout) :: vsend
      real, dimension(..), allocatable, target, optional, intent(inout) :: vrecv
      integer, intent(out), optional                      :: err

      if (present(vsend).and.present(vrecv)) then
        call self%decomposition_update_sendrecv(vsend,vrecv,err)
      elseif (present(vsend)) then
        call self%decomposition_update_single(vsend,err)
      else
        call self%decomposition_update_bottom(err)
      endif

      end subroutine decomposition_update_alt

!> finalize halo
!dox @private
      subroutine finalize_halo(h)
        type(halo) :: h

        integer mpierr

        if (isdebug()) print*,'Finalizing halo ',h%m
        if (associated(h%i)) then
!          deallocate(h%i)
          nullify(h%i)
        endif
        if (h%initialized) then
          call MPI_Type_free(h%m,mpierr)
          h%initialized=.false.
        endif

      end subroutine finalize_halo

!> copy halo
      subroutine copy_halo(hout,hin)
        class(halo),intent(inout) :: hout
        class(halo),intent(in)    :: hin

        integer mpierr

        if (isdebug()) print*,'copy_halo ',hin%initialized,hin%m
        hout%absolute=hin%absolute
        if (hin%m.ne.MPI_DATATYPE_NULL) then
          call MPI_Type_dup(hin%m,hout%m,mpierr)
        endif
        if(associated(hin%i))allocate(hout%i,source=hin%i)
      end subroutine copy_halo

!> check subarray halo
!dox @private
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
!dox @private
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

!> is halo absolute?
!dox @private
      function halo_is_absolute(self)
        use mpi 
        logical  :: halo_is_absolute

        class(halo), intent(in)                         :: self

        call debug('halo_is_absolute')

        halo_is_absolute=self%absolute

      end function halo_is_absolute
                                                                        
!> Commit a halo. If the object is a joined halo, it first creates the derived datatype.
!dox @relates gabriel::halo
      subroutine halo_commit(self,err) 
        use mpi 
                                                                        
        class(halo),intent(inout) :: self
        integer, intent(out), optional :: err

        integer mpierr
                                                                        
        select type (j=>self%i)
        type is (joined)                                                                        
          if (.not.j%initialized) then
            call MPI_Type_create_hindexed_block(j%n,1,j%variables,j%i%m,self%m,mpierr) 
            self%initialized=.true.
            if(.not.self%is_absolute()) call error(999,"Sanity check: type is joined, but not absolute",err)
          endif
        end select

        call MPI_Type_commit(self%m,mpierr) 
                                                                        
      end subroutine halo_commit

!> check validity of a halo for a variable
!dox @private
      function check_halo(self,v)
        logical :: check_halo
        class(halo),intent(in) :: self
        real,dimension(..),allocatable,intent(in) :: v

        call debug('check_halo')

        check_halo=.false.
! if the halo is of absolute type, return .false.
        if (self%is_absolute()) return
        select type(h=>self%i)
        type is (subarray)
          check_halo=h%is_valid_halo(v)
        type is (combined)
          check_halo=h%is_valid_halo(v)
        class default
          if (h%initialized) check_halo=.true.
        end select

      end function check_halo

!> Function to return a duplicate of the MPI type
!dox @relates gabriel::halo
!dox @public
      function mpitype(self)
        class(halo), intent(in) :: self
        integer :: mpitype
        integer :: mpierr

       call MPI_Type_dup(self%m,mpitype,mpierr)

      end function mpitype

end module
