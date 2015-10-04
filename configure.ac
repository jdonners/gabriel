#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([gabriel], [1.0], [john.donners@surfsara.nl])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
AM_SILENT_RULES([yes])

AC_LANG([Fortran])
#set the extension to f90, to compile the tests
AC_FC_SRCEXT(f90)

LT_INIT

# AC_PROG_FC
AX_PROG_FC_MPI([true],
  [ use_mpi=yes],
  [AC_MSG_FAILURE([MPI compiler needed, but could not find MPI.],1)
  ])

    # test if MPI-3 is available
    # We do not use AC_SEARCH_LIBS here, as it caches its outcome and
    # thus disallows corresponding calls in the other AX_PROG_*_MPI
    # macros.
    AC_MSG_CHECKING([for function MPI_Neighbor_alltoallw])
    AC_LINK_IFELSE([AC_LANG_CALL([],[MPI_NEIGHBOR_ALLTOALLW])],
      [ _ax_mpi_neighbor_found=yes ],
      [ _ax_mpi_neighbor_found=no ])
    AC_MSG_RESULT($_ax_mpi_neighbor_found)

    AC_MSG_CHECKING([for function MPI_Type_create_hindexed_block])
    AC_LINK_IFELSE([AC_LANG_CALL([],[MPI_TYPE_CREATE_HINDEXED_BLOCK])],
      [ _ax_mpi_typehblock_found=yes ],
      [ _ax_mpi_typehblock_found=no ])
    AC_MSG_RESULT($_ax_mpi_typehblock_found)

if test "x$_ax_mpi_typehblock_found" == "xno" -o "x$_ax_mpi_neighbor_found" == "xno" ; then
  AC_MSG_ERROR([Unable to find MPI-3 compatible library.])
fi

    AC_MSG_CHECKING([for module mpi_f08])
    AC_LINK_IFELSE(
     [AC_LANG_SOURCE([[
      program conftest
      use mpi_f08
      implicit none

      end program]])],
    [_ax_mpi_mpi_f08_found=yes],[_ax_mpi_mpi_f08_found=no])
    AC_MSG_RESULT($_ax_mpi_mpi_f08_found)

AC_FC_MODULE_FLAG
# Lookup extension for modules, used for 'make clean'
AC_FC_MODULE_EXTENSION
AC_SUBST(FC_MODEXT)

# Lookup flag for passing preprocessor flags to the Fortran compiler. Used to define HAVE_CONFIG_H
AC_FC_PP_DEFINE
AC_SUBST(FC_DEFINE)

AC_MSG_CHECKING([for flag to promote REAL to DOUBLE PRECISION])
save_FCFLAGS="$FCFLAGS"
for ax_option in -r8 -qautodbl=dbl4 "-fdefault-real-8 -fdefault-double-8" "-s real64" -Mr8; do
  FCFLAGS="$save_FCFLAGS $ax_option $dummy_FCFLAGS"
  AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
      program conftest
      implicit none

      print*,dsin(0.0)

      end program]])],
    [real_to_double=$ax_option])
  if test "$real_to_double" != ""; then
    AC_MSG_RESULT([$real_to_double])
    break
  fi
done

AC_MSG_CHECKING([for assumed rank support])
AC_LINK_IFELSE(
  [AC_LANG_SOURCE([[
      program conftest
      contains

      subroutine arr(a)
        real,dimension(..),allocatable :: a 

        print*,rank(a),lbound(a),ubound(a)
      end subroutine
      end]])],
  [assumedrank_works=yes],
  [assumedrank_works=no AC_MSG_ERROR([Compiler doesn''t support assumed rank arrays. Please use gcc 4.9, Intel 2016 or higher])])
AC_MSG_RESULT([$assumedrank_works])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT