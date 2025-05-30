## Process this file with automake to produce Makefile.in

CHECKS = test-ccs.report.yaml test-foeccs.report.yaml test-random.report.yaml test-eigenvalues.report.yaml

#include $(srcdir)/../check.mk

check: $(CHECKS) report

EXTRA_DIST = test-ccs.ref.yaml \
	matrix1-base.dat matrix2-base.dat \
	test-foeccs.ref.yaml \
    overlap_ccs-base.txt hamiltonian_ccs-base.txt \
	density_kernel_ccs-base.txt density_kernel_matmul_ccs-base.txt \
	data1_base.dat data2_base.dat

check_PROGRAMS = driver_ccs driver_foe driver_single driver_random driver_eigenvalues
EXTRA_PROGRAMS = driver

AM_FCFLAGS = -I$(top_builddir)/includes -I. @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
chess_library = $(top_builddir)/src/libCheSS-1.so.@CHESS_MINOR_VERSION@
else
chess_library = $(top_builddir)/src/libCheSS-1.a @LIB_CHESS_DEPS@
endif

driver_ccs_SOURCES = utilities.f90 driver_ccs.f90
driver_ccs_LDADD = $(chess_library) @LINALG_LIBS@

driver_SOURCES = driver.f90

driver_foe_SOURCES = utilities.f90 driver_foe.f90
driver_foe_LDADD = $(chess_library) @LINALG_LIBS@

driver_single_SOURCES = utilities.f90 driver_single.f90
driver_single_LDADD = $(chess_library) @LIB_FUTILE_LIBS@ @LINALG_LIBS@

driver_random_SOURCES = utilities.f90 driver_random.f90
driver_random_LDADD = $(chess_library) @LINALG_LIBS@

driver_eigenvalues_SOURCES = utilities.f90 driver_eigenvalues.f90
driver_eigenvalues_LDADD = $(chess_library) @LINALG_LIBS@

CLEANFILES = test-ccs.out.yaml test-ccs.report.yaml test-ccs.report \
             test-foeccs.out.yaml test-foeccs.report.yaml test-foeccs.report \
			 test-random.out.yaml test-random.report.yaml test-random.report \
			 test-eigenvalues.out.yaml test-eigenvalues.report.yaml test-eigenvalues.report *.@MODULE_EXT@

test-ccs.report.yaml: $(srcdir)/test-ccs.ref.yaml test-ccs.out.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $< -d test-ccs.out.yaml -t $(abs_top_srcdir)/tests/tols.yaml --label=test-ccs -o $@

test-foeccs.report.yaml: $(srcdir)/test-foeccs.ref.yaml test-foeccs.out.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $< -d test-foeccs.out.yaml -t $(abs_top_srcdir)/tests/tols.yaml --label=test-foeccs -o $@ 

test-random.report.yaml: $(srcdir)/test-random.ref.yaml test-random.out.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $< -d test-random.out.yaml -t $(abs_top_srcdir)/tests/tols.yaml --label=test-random -o $@ 

test-eigenvalues.report.yaml: $(srcdir)/test-eigenvalues.ref.yaml test-eigenvalues.out.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $< -d test-eigenvalues.out.yaml -t $(abs_top_srcdir)/tests/tols.yaml --label=test-eigenvalues -o $@ 

test-ccs.out.yaml: driver_ccs
	cp -f $(srcdir)/matrix1-base.dat $(PWD)/matrix1.dat
	cp -f $(srcdir)/matrix2-base.dat $(PWD)/matrix2.dat
	$(run_parallel) ./driver_ccs > test-ccs.out.yaml

test-foeccs.out.yaml: driver_foe
	cp -f $(srcdir)/overlap_ccs-base.txt $(PWD)/overlap_ccs.txt
	cp -f $(srcdir)/hamiltonian_ccs-base.txt $(PWD)/hamiltonian_ccs.txt
	cp -f $(srcdir)/density_kernel_ccs-base.txt $(PWD)/density_kernel_ccs.txt
	cp -f $(srcdir)/density_kernel_matmul_ccs-base.txt $(PWD)/density_kernel_matmul_ccs.txt
	cp -f $(srcdir)/sparsematrix_metadata-base.dat $(PWD)/sparsematrix_metadata.dat
	$(run_parallel) ./driver_foe \
        --sparsity_format=ccs \
        --matrix_format=serial_text \
        --metadata_file=sparsematrix_metadata.dat \
        --overlap_file=overlap_ccs.txt \
        --hamiltonian_file=hamiltonian_ccs.txt \
        --kernel_file=density_kernel_ccs.txt \
        --kernel_matmul_file=density_kernel_matmul_ccs.txt \
        --check_spectrum=yes \
        > test-foeccs.out.yaml

test-random.out.yaml: driver_random
	cp -f $(srcdir)/data1_base.dat $(PWD)/data1.dat
	cp -f $(srcdir)/data2_base.dat $(PWD)/data2.dat
	$(run_parallel) ./driver_random \
       --sparsegen_method=random \
       --matgen_method=random \
       --nfvctr=100 \
       --nvctr=10000 \
       --condition_number=50 \
       --expo=-1 \
       > test-random.out.yaml
	$(run_parallel) ./driver_random \
       --sparsegen_method=random \
       --matgen_method=random \
       --nfvctr=2000 \
       --nvctr=100000 \
       --condition_number=100 \
       --expo=-0.5 \
       >> test-random.out.yaml
