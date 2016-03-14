lib_LTLIBRARIES = libgabriel.la
libgabriel_la_LDFLAGS = -version-info 2:0:0
libgabriel_la_SOURCES = gabriel.f90

nodist_include_HEADERS = gabriel.$(FC_MODEXT)

TESTS = example1d example_auto example_auto_offset example_transform example_err example_joined example_auto_joined bandwidth
check_PROGRAMS = $(TESTS)
AM_DEFAULT_SOURCE_EXT = .f90
LDADD = -lgabriel

EXTRA_DIST =  gabriel.docset

# Change preprocessor flag. This is specific to some Fortran compilers.
# This hack shouldn't be used if C/C++ code is included!
dist-hook:
	chmod u+w $(distdir)/configure
	sed -i 's/^DEFS=-DHAVE_CONFIG_H$$/DEFS=$${FC_DEFINE}HAVE_CONFIG_H/' $(distdir)/configure	

# Clean module files
mostlyclean-local:
	-test -z "$(FC_MODEXT)" || rm -f gabriel.$(FC_MODEXT)


