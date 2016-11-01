program ex

  use gabriel
  use MPI
  use iso_fortran_env
  implicit none

  namelist/halo_exchange/array_size,nproc_row,halo_width,warmloop,timeloop,vars

  integer :: array_size=50
  integer :: nproc_row=8
  integer :: halo_width=2
  integer :: warmloop=100
  integer :: timeloop=1000
  integer :: vars=4
  real*8  :: gbpers

  real,dimension(:,:,:),allocatable :: a,b,c,e
  real :: dummyreal
  integer ierr,rank,right,left,mpisize,i,j,k
  integer hor,ver,pow2,s,realsize,realtype
  real*8 t0, t1,t,datasize

  type(box) :: bo
  type(distribution) :: d
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mpisize,ierr)
  right=rank+1
  left=rank-1
  if (right.ge.mpisize)right=0
  if (left.lt.0)left=mpisize-1
!  print*,'left,rank,right=',left,rank,right

  open(11, file='input')
  read(11, nml = halo_exchange)
  s=array_size
  if(rank.eq.0)then
    write(output_unit,nml=halo_exchange)
    print*,'MPI size: ',mpisize
  endif
  if (vars.lt.2.or.vars.gt.4)then
    print*,'Vars should be between 2 and 4'
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  endif

  hor=mod(rank,nproc_row)
  ver=rank/nproc_row
  allocate(a(hor*s-halo_width:hor*s+s+halo_width-1,ver*s-halo_width+1:ver*s+s+halo_width,s))
  allocate(b(hor*s-halo_width:hor*s+s+halo_width-1,ver*s-halo_width+1:ver*s+s+halo_width,s))
  if(vars.ge.3)allocate(c(hor*s-halo_width:hor*s+s+halo_width-1,ver*s-halo_width+1:ver*s+s+halo_width,s))
  if(vars.ge.4)allocate(e(hor*s-halo_width:hor*s+s+halo_width-1,ver*s-halo_width+1:ver*s+s+halo_width,s))
  a=rank
  b=10+rank
  call zet_array(a,halo_width)
  call zet_array(b,halo_width)
  if(vars.ge.3)call zet_array(c,halo_width)
  if(vars.ge.4)call zet_array(e,halo_width)

!  write(*,'(a,4i3)')'Size,Rank,hor,ver=',mpisize,rank,hor,ver

  call gabriel_init
  
  call bo%init(a,(/hor*s,ver*s+1,1/),(/hor*s+s-1,ver*s+s,s/),MPI_COMM_WORLD,periodic=(/.true.,.true.,.true./))
  call d%halo(bo)
!  call d%autocreate(a,(/hor*s,ver*s+1,1/),(/hor*s+s-1,ver*s+s,6/),MPI_COMM_WORLD,periodic=(/.true.,.false.,.false./))
  call d%joined(vars)
  call d%joined_add(a)
  call d%joined_add(b)
  if(vars.ge.3)call d%joined_add(c)
  if(vars.ge.4)call d%joined_add(e)
  call d%create()

  if(rank.eq.0)print*,'Warming up..'
  do i=1,warmloop
    call d%update()
  enddo
  if(rank.eq.0)print*,'Apply distribution..'
  call gabriel_disable_checking

  t=0.0
  pow2=1
  do i=1,timeloop
    if (rank.eq.0.and.i.eq.pow2) then
      write(*,'(i5,a1)')i,' '
      pow2=pow2*2
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    t0=MPI_Wtime()
    call d%update()

    t1=MPI_Wtime()
    t=t+(t1-t0)
!    a=0.0
    call zet_array(a,halo_width)
    call zet_array(b,halo_width)
    if(vars.ge.3)call zet_array(c,halo_width)
    if(vars.ge.4)call zet_array(e,halo_width)
  enddo

  call MPI_Type_create_f90_real(precision(dummyreal),exponent(dummyreal),realtype,ierr)
  call MPI_Type_size(realtype,realsize,ierr)
!          vars                                      directions*type_size)/kb
  datasize=vars*(halo_width*s*s+halo_width*halo_width*s)*4.d0*realsize/1024.d0
  if(rank.eq.0)then
    call MPI_Reduce(MPI_IN_PLACE,t,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  else
    call MPI_Reduce(t,t,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  endif
  if(rank.eq.0)then
    print*,'Bytes per value         :',realsize
    print*,'Data size per comm (kB) :',datasize
    print*,'Comm. time (s)          :',t
    gbpers=datasize*timeloop/(1024.d0*1024.d0*t)
    write(*,'(a,2f10.2)')'Comm. speed (GB/s, GBit/s)      :',gbpers,gbpers*8
  endif

!  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s,a(:,ver*s,1)
!  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s+1,a(:,ver*s+1,1)
!  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s+2,a(:,ver*s+2,1)
!  write(*,'(a,i3,a,i4,4f13.3)')' AFTER Rank',rank,' data=',ver*s+3,a(:,ver*s+3,1)
  deallocate(a,b)
  if(vars.ge.3)deallocate(c)
  if(vars.ge.4)deallocate(e)
  call MPI_Finalize(ierr)

contains

  subroutine zet_array(a,halo_width)
    real,dimension(:,:,:),allocatable :: a
    integer halo_width
    integer ii,j,k

    a=0.0
    do k=lbound(a,3)+halo_width,ubound(a,3)-halo_width
    do j=lbound(a,2)+halo_width,ubound(a,2)-halo_width
    do ii=lbound(a,1)+halo_width,ubound(a,1)-halo_width
      a(ii,j,k)=a(ii,j,k)+ii*0.1+j*0.01+k*0.001
    enddo
    enddo
    enddo
  end subroutine zet_array    

  logical function check_array_halo(a,halo_width)
    real,dimension(:,:,:),allocatable :: a,b
    integer :: nprocs,nproc_row,row,col,halo_width
    logical :: check
    integer ii,j,k

    check=.true.
    do k=lbound(a,3),ubound(a,3)
    do j=lbound(a,2),ubound(a,2)
    do ii=lbound(a,1),ubound(a,1)
      if (ii.lt.lbound(a,1)+halo_width .or.ii.gt.lbound(a,1)-halo_width.or. &
     &    j.lt.lbound(a,2)+halo_width .or.j.gt.lbound(a,2)-halo_width.or.&
     &    k.lt.lbound(a,3)+halo_width .or.k.gt.lbound(a,3)-halo_width)then
!         point is part of the halo
        check=check.and.(a(ii,j,k).eq.ii*0.1+j*0.01+k*0.001)
      endif
    enddo
    enddo
    enddo
    check_array_halo=check
  end function check_array_halo
end

