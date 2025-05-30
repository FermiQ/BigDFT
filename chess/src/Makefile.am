# To be installed:
dynlibdir        = $(libdir)
lib_LIBRARIES    = libCheSS-1.a
module_DATA      = $(modules)
dynlib_DATA      = $(chess_dynlib)
noinst_DATA      = $(chess_linklib)
bin_PROGRAMS     = $(binaries_sources)

# Static or dynamic library
chess_dynamic_library = libCheSS-1.so.@CHESS_MINOR_VERSION@.0.@CHESS_MICRO_VERSION@
chess_library_soname  = libCheSS-1.so.@CHESS_MINOR_VERSION@

if BUILD_DYNAMIC_LIBS
chess_dynlib = $(chess_dynamic_library)
chess_linklib = $(chess_library_soname)
chess_library = $(chess_library_soname)
else
chess_dynlib =
chess_linklib =
chess_library = libCheSS-1.a @LIB_CHESS_DEPS@
endif

if BUILD_BINARIES 
binaries_sources = chess_toolbox
else
binaries_sources =
endif

SUBDIRS=python

# C objects
c_level = chess_input_variables_fetch.c

chess_input_variables_fetch.o: chess_input_variables_definition-inc.h
#build the fortran include file with the stream of input variables into a local string
#the sed command below will transform the entire yaml document into a c string, to be passed to fortran
chess_input_variables_definition-inc.h: chess_input_variables_definition.yaml
	$(SED) -e "s/^/\"/;s/$$/\\\n\"/" $< > $@

libCheSS_1_a_SOURCES = chebyshev.f90 \
	  foe_base.f90 \
	  foe_common.f90 \
	  fermi_level.f90 \
	  foe.f90 \
	  ice.f90 \
	  coeffs.f90 \
	  matrix_operations.f90 \
	  sparsematrix_highlevel.f90 \
	  highlevel_wrappers.f90 \
	  chess_base.f90 \
	  $(c_level)

# Module files.
if CAPITALIZE
modules = CHEBYSHEV.@MODULE_EXT@ \
	  FOE_BASE.@MODULE_EXT@ \
	  FOE_COMMON.@MODULE_EXT@ \
	  FERMI_LEVEL.@MODULE_EXT@ \
	  FOE.@MODULE_EXT@ \
      MODULE_FUNC.@MODULE_EXT@ \
	  ICE.@MODULE_EXT@ \
	  COEFFS.@MODULE_EXT@ \
	  MATRIX_OPERATIONS.@MODULE_EXT@ \
	  SPARSEMATRIX_HIGHLEVEL.@MODULE_EXT@ \
	  HIGHLEVEL_WRAPPERS.@MODULE_EXT@ \
	  CHESS_BASE.@MODULE_EXT@
else
modules =  chebyshev.@MODULE_EXT@ \
	  foe_base.@MODULE_EXT@ \
	  foe_common.@MODULE_EXT@ \
	  fermi_level.@MODULE_EXT@ \
	  foe.@MODULE_EXT@ \
      module_func.@MODULE_EXT@ \
	  ice.@MODULE_EXT@ \
	  coeffs.@MODULE_EXT@ \
	  matrix_operations.@MODULE_EXT@ \
	  sparsematrix_highlevel.@MODULE_EXT@ \
	  highlevel_wrappers.@MODULE_EXT@ \
	  chess_base.@MODULE_EXT@
endif

EXTRA_DIST = \
    chess_input_variables_definition.yaml \
	configure.inc.in

subs_static_library = $(top_builddir)/sparsematrix/libsparsematrix.a \
                      $(top_builddir)/chesspexsi/libchesspexsi.a


chess_toolbox_SOURCES = chess_toolbox.f90
chess_toolbox_LDADD = $(chess_library) @LIB_FUTILE_LIBS@ @LINALG_LIBS@

CLEANFILES = *.@MODULE_EXT@ $(chess_dynamic_library) $(chess_library_soname) libCheSS-1.so
clean-local:
	rm -rf tmp-libCheSS

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@

if BUILD_DYNAMIC_LIBS
LD_RPATH=-Wl,-rpath=$(DESTDIR)$(dynlibdir)
AM_LDFLAGS = $(LD_RPATH)
endif

$(modules): $(libCheSS_1_a_OBJECTS)

# Dynamic library building
$(chess_dynamic_library): $(libCheSS_1_a_OBJECTS) tmp-libCheSS/extract.stamp
	$(FC) $(FCFLAGS) -shared $(LDFLAGS) $(AM_LDFLAGS) -Wl,-soname=$(chess_library_soname) $(LD_RPATH) -o $@ $(libCheSS_1_a_OBJECTS) tmp-libCheSS/*.o @LIB_CHESS_DEPS@ $(LIBS)
	@chmod a+x $@

$(chess_library_soname): $(chess_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libCheSS-1.so

install-data-hook:
	if test -n "$(chess_dynlib)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(chess_dynamic_library) ; \
	  ln -fs $(chess_dynamic_library) $(chess_library_soname) ; \
	  ln -fs $(chess_library_soname) libCheSS-1.so ; \
	fi

tmp-libCheSS/extract.stamp: $(subs_static_library)
	test -e "tmp-libCheSS" || $(INSTALL) -d -m 755 tmp-libCheSS
	cd tmp-libCheSS ; \
	rm -f *; \
	for lib in $^; do \
	  $(AR) x ../$$lib ; \
	done ; \
	touch extract.stamp

libCheSS-1.a: $(libCheSS_1_a_OBJECTS) tmp-libCheSS/extract.stamp insmod
	$(AR) $(ARFLAGS) $@ $(libCheSS_1_a_OBJECTS) tmp-libCheSS/*.o
	$(RANLIB) $@

#dependencies
$(modules): $(libCheSS_1_a_OBJECTS)
foe_common.o: foe_base.o chebyshev.o fermi_level.o chess_base.o
ice.o: chebyshev.o foe_base.o foe_common.o fermi_level.o
matrix_operations.o: ice.o foe_base.o
foe.o: chebyshev.o foe_base.o foe_common.o fermi_level.o matrix_operations.o
sparsematrix_highlevel.o: ice.o foe.o
highlevel_wrappers.o: sparsematrix_highlevel.o
coeffs.o: sparsematrix_highlevel.o
chess_toolbox.o: $(libCheSS_1_a_OBJECTS)

# Copy modules in a common include directory.
insmod: $(modules)
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
