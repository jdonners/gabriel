lib_LTLIBRARIES = libgabriel.la
libgabriel_la_LDFLAGS = -version-info 0:0:0
libgabriel_la_SOURCES = halos.f90

include_HEADERS = gabriel.$(FC_MODEXT)

EXTRA_DIST =  halos.docset

# Change preprocessor flag. This is specific to some Fortran compilers.
# This hack shouldn't be used if C/C++ code is included!
dist-hook:
	chmod u+w $(distdir)/configure
	sed -i 's/^DEFS=-DHAVE_CONFIG_H$$/DEFS=$${FC_DEFINE}HAVE_CONFIG_H/' $(distdir)/configure	

# Clean module files
mostlyclean-local:
	-test -z "$(FC_MODEXT)" || rm -f halos.$(FC_MODEXT)


