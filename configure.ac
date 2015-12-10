#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([gabriel], [1.2], [john.donners@surfsara.nl])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
AM_SILENT_RULES([yes])

AC_LANG([Fortran])

# AC_PROG_FC
AX_PROG_FC_MPI([true],
  [ use_mpi=yes],
  [AC_MSG_FAILURE([MPI compiler needed, but could not find MPI.],1)
  ])

#set the extension to f90, to compile the tests
AC_FC_SRCEXT(f90)

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

r8_wanted=yes
AC_ARG_ENABLE([real8], AS_HELP_STRING([--enable-real8],
  [enable 8-byte real arrays]),
  [if test "x$enableval" == "xno"; then r8_wanted=no; fi],[])
if test "x$enableval" == "xyes"; then
AC_MSG_CHECKING([for flag to promote REAL to REAL*8])
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
fi

    AC_MSG_CHECKING([if ${FC} supports finalization])
    AC_LINK_IFELSE(
     [AC_LANG_SOURCE([[
      module conftest
      implicit none

      type :: a
        integer b
        contains
        final :: finalize_a
      end type
      contains
      subroutine finalize_a(x)
        type(a) :: x
          x%b=0
      end subroutine
      end module

      program t
        use conftest
      end program]])],
    [_ax_f03_finalization_found=yes],[_ax_f03_finalization_found=no])
    AC_MSG_RESULT($_ax_f03_finalization_found)
if test "x$_ax_f03_finalization_found" == "xno" ; then
  AC_MSG_ERROR([${FC} does not support finalization.])
fi

AC_MSG_CHECKING([if ${FC} supports assumed rank arrays])
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
  [assumedrank_works=no AC_MSG_ERROR([${FC} does not support assumed rank arrays.])])
AC_MSG_RESULT([$assumedrank_works])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
LT_INIT

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
