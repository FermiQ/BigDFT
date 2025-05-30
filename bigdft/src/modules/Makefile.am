# To be installed:
noinst_LIBRARIES = libmodules.a
module_DATA = $(modules)
dynlibdir = $(libdir)
dynlib_DATA = $(openbabel_dynamic_library)
noinst_DATA = $(openbabel_soname_library)

#if USE_MPI
#mpi_include =
#else
#mpi_include = mpif.h
#endif

if HAVE_CP2K
cp2k_interfacemod = module_cp2k.f90
cp2k_interfacemodCAP = MODULE_CP2k.F90
cp2k_interfacemodmod = module_cp2k
cp2k_interfacemodCAPmod = MODULE_CP2k
else
cp2k_interfacemod = cp2k_fake.f90 module_cp2k.f90
cp2k_interfacemodCAP = CP2K_FAKE.F90 MODULE_CP2K.F90
cp2k_interfacemodmod = cp2k_fake module_cp2k
cp2k_interfacemodCAPmod = CP2K_FAKE MODULE_CP2K
endif

if HAVE_ALBORZ
alborz_interfacemod = module_alborz.f90
alborz_interfacemodCAP = MODULE_ALBORZ.F90
alborz_interfacemodmod = module_alborz
alborz_interfacemodCAPmod = MODULE_ALBORZ
else
alborz_interfacemod = alborz_fake.f90 module_alborz.f90
alborz_interfacemodCAP = alborz_FAKE.F90 MODULE_ALBORZ.F90
alborz_interfacemodmod = alborz_fake module_alborz
alborz_interfacemodCAPmod = alborz_FAKE MODULE_alborz
endif

#if HAVE_PEXSI
#PEXSI_SRC = pexsi_base.f90 pexsi_interfaces.f90 pexsi.f90
#else
#PEXSI_SRC = pexsi_base.f90 pexsi_interfaces_fake.f90 pexsi.f90
#endif

libmodules_a_SOURCES = ab7_mixing.f90 \
	dpbox.f90 \
	psp_projectors.f90 \
	psp_projectors_base.f90 \
	input_keys.f90 \
	gaussdaub.f90 \
	input_dicts.f90 \
	types.f90 \
	input.f90 \
	interfaces.f90 \
	private_api.f90 \
	module_lj.f90 \
	module_morse_bulk.f90 \
	module_lenosky_si.f90 \
	module_tersoff.f90 \
	module_alborz.f90 \
	module_bornmayerhugginstosifumi.f90 \
	module_tdpot.f90 \
	module_coulomb.f90 \
	module_bazant.f90 \
	BigDFT_API.f90 \
	vdwcorrection.f90 \
	fragments.f90 \
	fragment_base.f90 \
	constrained_dft.f90 \
	diis_module.f90 \
	reformatting.f90 \
	atoms_data.f90 \
	communications_init.f90 \
	communications.f90 \
	sparsematrix_wrappers.f90 \
	internal_coordinates.f90 \
	module_sqn.f90 \
	transposed_operations.f90 \
	module_forces.f90 \
	unitary_tests.f90 \
	io.f90 \
	postprocessing_linear.f90 \
	orbital_basis.f90 \
	bigdft_run.f90 \
	sockets.c \
	fsockets.f90 \
	rhopotential.f90 \
	multipole_base.f90 \
	multipole.f90 \
	orthonormalization.f90 \
	locregs_init.f90 \
	sort.f90 \
	forces_linear.f90 \
	bigdft_matrices.f90 \
	configure.c \
	$(alborz_interfacemod) \
	$(cp2k_interfacemod) \
	$(openbabel_static_sources) \
	module_dftbp.f90 \
	constrained_magnetic_field.f90 \
	sw_potential.f90 \
	get_basis.f90 \
	get_kernel.f90

# Modules
if CAPITALIZE
modules = MODULE_TYPES.@MODULE_EXT@ \
	MODULE_MIXING.@MODULE_EXT@ \
	MODULE_INTERFACES.@MODULE_EXT@ \
	MODULE_PRIVATE_API.@MODULE_EXT@ \
	MODULE_INPUT.@MODULE_EXT@ \
	GAUSSDAUB.@MODULE_EXT@ \
	BIGDFT_API.@MODULE_EXT@ \
	VDWCORRECTION.@MODULE_EXT@ \
	MODULE_FRAGMENTS.@MODULE_EXT@ \
	CONSTRAINED_DFT.@MODULE_EXT@ \
	DIIS_SD_OPTIMIZATION.@MODULE_EXT@ \
	INPUT_OLD_TEXT_FORMAT.@MODULE_EXT@ \
	MODULE_INPUT_KEYS.@MODULE_EXT@ \
	MODULE_INPUT_DICTS.@MODULE_EXT@ \
	PSP_PROJECTORS.@MODULE_EXT@ \
	PSP_PROJECTORS_BASE.@MODULE_EXT@ \
	MODULE_DPBOX.@MODULE_EXT@ \
	REFORMATTING.@MODULE_EXT@ \
	MODULE_ATOMS.@MODULE_EXT@ \
	MODULE_CFD.@MODULE_EXT@ \
	COMMUNICATIONS_INIT.@MODULE_EXT@ \
	COMMUNICATIONS.@MODULE_EXT@ \
	SPARSEMATRIX_WRAPPERS.@MODULE_EXT@ \
	INTERNAL_COORDINATES.@MODULE_EXT@ \
	MODULE_SQN.@MODULE_EXT@ \
	MODULE_FORCES.@MODULE_EXT@ \
	TRANSPOSED_OPERATIONS.@MODULE_EXT@ \
	UNITARY_TESTS.@MODULE_EXT@ \
	IO.@MODULE_EXT@ \
	POSTPROCESSING_LINEAR.@MODULE_EXT@ \
	BIGDFT_RUN.@MODULE_EXT@ \
	RHOPOTENTIAL.@MODULE_EXT@ \
	MULTIPOLE_BASE.@MODULE_EXT@ \
	FRAGMENT_BASE.@MODULE_EXT@ \
	MULTIPOLE.@MODULE_EXT@ \
	ORTHONORMALIZATION.@MODULE_EXT@ \
	LOCREGS_INIT.@MODULE_EXT@ \
	SORT.@MODULE_EXT@ \
	FORCES_LINEAR.@MODULE_EXT@ \
	MODULE_LJ.@MODULE_EXT@ \
	F90SOCKETS.@MODULE_EXT@ \
	MODULE_LENOSKY_SI.@MODULE_EXT@\
	MODULE_MORSE_BULK.@MODULE_EXT@\
	MODULE_TERSOFF.@MODULE_EXT@\
	BIGDFT_MATRICES.@MODULE_EXT@ \
	MODULE_BORNMAYERHUGGINSTOSIFUMI.@MODULE_EXT@\
	MODULE_DFTBP.@MODULE_EXT@ \
	MODULE_COULOMB.@MODULE_EXT@ \
	MODULE_BAZANT.@MODULE_EXT@ \
	ORBITALBASIS.@MODULE_EXT@ \
	MODULE_TDPOT.@MODULE_EXT@ \
	MODULE_CP2K.@MODULE_EXT@ \
	SWPOTENTIAL.@MODULE_EXT@ \
	GET_BASIS.@MODULE_EXT@ \
	GET_KERNEL.@MODULE_EXT@
else
modules = module_types.@MODULE_EXT@ \
	module_mixing.@MODULE_EXT@ \
	module_interfaces.@MODULE_EXT@ \
	module_private_api.@MODULE_EXT@ \
	module_input.@MODULE_EXT@ \
	gaussdaub.@MODULE_EXT@ \
	bigdft_api.@MODULE_EXT@ \
	vdwcorrection.@MODULE_EXT@ \
	module_fragments.@MODULE_EXT@ \
	constrained_dft.@MODULE_EXT@ \
	diis_sd_optimization.@MODULE_EXT@ \
	input_old_text_format.@MODULE_EXT@ \
	module_input_keys.@MODULE_EXT@ \
	module_input_dicts.@MODULE_EXT@ \
	psp_projectors.@MODULE_EXT@ \
	psp_projectors_base.@MODULE_EXT@ \
	module_dpbox.@MODULE_EXT@ \
	reformatting.@MODULE_EXT@ \
	module_atoms.@MODULE_EXT@ \
	module_cfd.@MODULE_EXT@ \
	communications_init.@MODULE_EXT@ \
	communications.@MODULE_EXT@ \
	sparsematrix_wrappers.@MODULE_EXT@ \
	internal_coordinates.@MODULE_EXT@ \
	module_sqn.@MODULE_EXT@ \
	module_forces.@MODULE_EXT@ \
	transposed_operations.@MODULE_EXT@ \
	unitary_tests.@MODULE_EXT@ \
	io.@MODULE_EXT@ \
	postprocessing_linear.@MODULE_EXT@ \
	bigdft_run.@MODULE_EXT@ \
	rhopotential.@MODULE_EXT@ \
	multipole_base.@MODULE_EXT@ \
	fragment_base.@MODULE_EXT@ \
	multipole.@MODULE_EXT@ \
	orthonormalization.@MODULE_EXT@ \
	locregs_init.@MODULE_EXT@ \
	forces_linear.@MODULE_EXT@ \
	bigdft_matrices.@MODULE_EXT@ \
	sort.@MODULE_EXT@ \
	module_lj.@MODULE_EXT@ \
	f90sockets.@MODULE_EXT@ \
	module_morse_bulk.@MODULE_EXT@ \
	module_tersoff.@MODULE_EXT@ \
	module_alborz.@MODULE_EXT@ \
	module_bornmayerhugginstosifumi.@MODULE_EXT@ \
	module_lenosky_si.@MODULE_EXT@ \
	orbitalbasis.@MODULE_EXT@ \
	module_dftbp.@MODULE_EXT@ \
	module_tdpot.@MODULE_EXT@ \
	module_cp2k.@MODULE_EXT@ \
	module_coulomb.@MODULE_EXT@ \
	module_bazant.@MODULE_EXT@ \
	swpotential.@MODULE_EXT@ \
	get_basis.@MODULE_EXT@ \
	get_kernel.@MODULE_EXT@
endif

ATOMS_INCLUDES = astruct-inc.f90
COMMS_INCLUDES = check_array-inc.f90

if HAVE_OPENBABEL
openbabel_dynamic_sources = openbabel.cpp
openbabel_dynamic_library = libbabel.so.@BIGDFT_MINOR_VERSION@.0.@BIGDFT_MICRO_VERSION@
openbabel_soname_library = libbabel.so.@BIGDFT_MINOR_VERSION@
openbabel_library = $(openbabel_soname_library)
AM_LDFLAGS = -Wl,-rpath=$(DESTDIR)$(dynlibdir)
else
openbabel_static_sources = openbabel_fake.cpp
endif

EXTRA_DIST = $(ATOMS_INCLUDES) $(COMMS_INCLUDES) openbabel.cpp
CLEANFILES = *.mod libbabel.so*

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@ @LIB_FUTILE_CFLAGS@
AM_CXXFLAGS = @LIB_FUTILE_CFLAGS@ @LIB_OPENBABEL_CFLAGS@

# Hack for XLF strange DEFS option.
PPFCCOMPILE = $(FC) $(DEFS:-D%=@FCDEFS@%) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)

#dependencies
$(modules): $(libmodules_a_OBJECTS)
#multipole_base.o: base.o
#fragment_base.o: public_keys.o

input.o: fragment_base.o
input_dicts.o: fragment_base.o input.o
atoms_data.o: $(ATOMS_INCLUDES) internal_coordinates.o input_dicts.o
input_keys.o: multipole_base.o atoms_data.o input.o fragment_base.o vdwcorrection.o
types.o: psp_projectors_base.o atoms_data.o input_keys.o \
	ab7_mixing.o multipole_base.o dpbox.o
interfaces.o: types.o fragments.o diis_module.o constrained_dft.o input_keys.o

BigDFT_API.o: types.o interfaces.o atoms_data.o input_dicts.o psp_projectors.o communications_init.o bigdft_run.o
vdwcorrection.o: atoms_data.o
private_api.o: types.o interfaces.o
diis_module.o: types.o
fragments.o: input_keys.o types.o fragment_base.o
constrained_dft.o: input_keys.o
#psp_projectors_base.o:  gaussians.o
psp_projectors.o: locregs_init.o psp_projectors_base.o
#reformatting.o: base.o
#communications_base.o: base.o
#dpbox.o: base.o bounds.o
communications_init.o: types.o
communications.o: types.o
sparsematrix_wrappers.o: types.o interfaces.o
gaussdaub.o: psp_projectors_base.o
#module_sqn.o: base.o
bigdft_run.o: types.o input_dicts.o interfaces.o module_forces.o module_lenosky_si.o module_lj.o module_morse_bulk.o \
		module_tersoff.o module_alborz.o module_bornmayerhugginstosifumi.o module_cp2k.o module_dftbp.o sw_potential.o module_tdpot.o \
		module_coulomb.o module_bazant.o
#ab7_mixing.o: defs.o
#module_lj.o: base.o defs.o
#module_tdpot.o: base.o defs.o
#module_cp2k.o: base.o defs.o
#module_morse_bulk.o: base.o
module_coulomb.o: types.o
module_bazant.o: types.o
module_dftbp.o: types.o 
module_tersoff.o: atoms_data.o
module_alborz.o: atoms_data.o
module_bornmayerhugginstosifumi.o: atoms_data.o
#module_lenosky_si.o: base.o defs.o
transposed_operations.o: types.o bigdft_matrices.o
module_forces.o: types.o 
unitary_tests.o: types.o interfaces.o communications_init.o communications.o transposed_operations.o rhopotential.o
io.o: types.o interfaces.o input_keys.o
postprocessing_linear.o: types.o interfaces.o communications.o transposed_operations.o rhopotential.o sort.o communications_init.o locregs_init.o multipole_base.o multipole.o io.o get_kernel.o
rhopotential.o: types.o interfaces.o communications.o bigdft_matrices.o constrained_magnetic_field.o atoms_data.o
multipole.o: multipole_base.o interfaces.o orthonormalization.o io.o orbital_basis.o
orbital_basis.o:  types.o interfaces.o communications.o
orthonormalization.o: types.o communications.o transposed_operations.o atoms_data.o
locregs_init.o: types.o communications.o sort.o interfaces.o
forces_linear.o: types.o psp_projectors.o transposed_operations.o
sw_potential.o: atoms_data.o
bigdft_matrices.o: types.o sparsematrix_wrappers.o
get_basis.o: types.o interfaces.o communications.o rhopotential.o constrained_dft.o fragment_base.o fragments.o transposed_operations.o locregs_init.o get_kernel.o
get_kernel.o: types.o interfaces.o diis_module.o communications.o constrained_dft.o rhopotential.o transposed_operations.o locregs_init.o orthonormalization.o
##pexsi.o also depends on pexsi_interfaces, but it is not known which one...
#if HAVE_PEXSI
#pexsi_interfaces.o: pexsi_base.o
#pexsi.o: pexsi_base.o pexsi_interfaces.o
#else
#pexsi_interfaces_fake.o: pexsi_base.o
#pexsi.o: pexsi_base.o pexsi_interfaces_fake.o
#endif


# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: $(modules) libmodules.a
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
#	test -z "$(mpi_include)" || $(INSTALL) -m 644 mpif.h $(top_builddir)/includes

$(openbabel_dynamic_library): $(openbabel_dynamic_sources)
	$(CXX) $(CXXFLAGS) @CC_PIC_FLAG@ @LIB_OPENBABEL_CFLAGS@ -shared -I$(top_builddir) @GLIB_CFLAGS@ @LIB_FUTILE_CFLAGS@ $(LDFLAGS) $(AM_LDFLAGS) -Wl,-soname=$(openbabel_soname_library) -o $@ $< @LIB_OPENBABEL_LIBS@ @LIB_LIBDICTS_LIBS@ 
	@chmod a+x $@

$(openbabel_soname_library): $(openbabel_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libbabel.so

install-data-hook:
	if test -n "$(openbabel_dynamic_library)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(openbabel_dynamic_library) ; \
	  ln -fs $(openbabel_dynamic_library) $(openbabel_soname_library) ; \
	  ln -fs $(openbabel_soname_library) libbabel.so ; \
	fi
