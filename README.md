Gabriel library
===============
Gabriel is a Fortran library to simplify the use of MPI.
It offers verification, ease of use and a high performance, 
especially for models on regular grids.

Compilation
-----------
The requirements to compile the gabriel library are

* a recent Fortran compiler supporting Fortran 2003 and TS 29113:
  - Intel 16.0.0 or newer
  - GNU gfortran 4.9 or newer
  - recent Cray Fortran
  - recent IBM XL Fortran, can anyone test this?
* a recent MPI library supporting the MPI-3.0 standard:
  - Intel MPI 5.0.0 or newer
  - MPICH2 3.0.4 or newer (and maybe older)
  - recent Cray MPI
  - recent IBM POE, can anyone test this?

It's recommended to download a release tarball of gabriel. The latest
release tarball is [gabriel-1.2.tar.gz](https://github.com/jdonners/gabriel/releases/download/1.2/gabriel-1.2.tar.gz).
An example of a build of the library is 

    ./configure                    # use --enable-real8 for REAL*8 arrays
    make
    make install prefix=[prefix]

where [prefix] is the directory where you want to install the library
and the Fortran module.

Should you want to build from the repository, you'll need recent versions
of the GNU autotools as well.

Usage
-----
There are several examples included in the distribution.

The library exposes two derived types to the user: HALO and DECOMPOSITION.
Each type has several type-bound procedures.

I'll explain the most important procedure in some more detail: autocreate.

In a regular 2D decomposition, each rank has 8 neighboring ranks:

    .-------.-------.-------.
    |       |       |       |
    |   1   |   2   |   3   |
    |      X|XXXXXXX|X      |
    .-------.-------.-------.
    |      X|OOOOOOO|X      |
    |   4  X|OOO5OOO|X  6   |
    |      X|OOOOOOO|X      |
    .-------.-------.-------.
    |      X|XXXXXXX|X      |
    |   7   |   8   |   9   |
    |       |       |       |
    .-------.-------.-------.

    [1-9] - rank
        O - local domain for rank 5
        X - halo regions for rank 5

The routine autocreate sets up all communication to update all halo regions of a 
variable. Note that this includes the elements located on diagonal ranks. The
gabriel library support arrays of any dimension and rank.

The input arguments are the variable to be updated, the lower and upper indices of the
array that are calculated locally and the communicator of all processes that
will calculate part of the full domain. E.g. for the simple case of a 1D array, 
it would be something like:

    MPI rank                0          1
    Array A              0123456    3456789
    Locally computed     .....        .....
    Halo region               ..    ..

So rank 0 calculates indices 0-4 and rank 1 calculates indices 5-9.
Both ranks need some data from the neighbor in the halo region
to be able to calculate its local indices. This would look something like:

    program one_two
      use mpi
      use gabriel

      real,dimension(:),allocatable :: a
      type(decomposition) :: dec
      integer ierr,rank

      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

      if (rank.eq.0) then
        allocate(a(0:6))
        a=0.0
        call dec%autocreate(a,(/0/),(/4/),MPI_COMM_WORLD)
      endif
      if (rank.eq.1) then
        allocate(a(3:9))
        a=1.0
        call dec%autocreate(a,(/5/),(/9/),MPI_COMM_WORLD)
      endif

      call dec%update(a,a)
      write(*,'(a,i2,a,5f5.1)')'Rank=',rank,'   My values:=',a
    end program

and the output of the program will look like:

    Rank= 0   My values:=  0.0  0.0  0.0  0.0  0.0  1.0  1.0
    Rank= 1   My values:=  0.0  0.0  1.0  1.0  1.0  1.0  1.0

And this works analogously for multi-dimensional arrays. 

The routine has an optional argument 'offset' to set the offset for each rank. E.g. if 
your arrays are defined to start at index 1 on each process, the offset will 
be added to the indices of the array to determine its position in the 
global domain across all processes. 

The optional argument 'periodic' can be used to indicate periodic boundary conditions.

