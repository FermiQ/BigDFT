## Process this file with automake to produce Makefile.in

EXTRA_DIST = \
	README \
	styles \
	images \
	tols-BigDFT.yaml

SUBDIRS = \
	  DFT \
	  overDFT \
	  tutorials \
	  libs \
	  benchs \
	  plugins

python_PYTHON = parse_log.py

if USE_OCL
oclrun_message = oclrun
else
oclrun_message =
endif

if USE_MPI
mpirun_message = mpirun
else
mpirun_message =
endif

head_message: all-recursive
	@echo "========================================================"
	@echo "To run the tests, type 'make check' (might take a while)"

foot_message: $(mpirun_message) $(oclrun_message) head_message
	@echo "========================================================"

mpirun: head_message
	@echo "     Use the environment variable run_parallel"
	@echo "         ex: export run_parallel='mpirun -np 2'  "

oclrun: head_message $(mpirun_message)
	@echo "     Use the environment variable run_ocl"
	@echo "         ex: export run_ocl='on' to use OpenCL acceleration"
	@echo "     Use run_ocl='CPU' or 'ACC' to force use of hardware"
	@echo "         different than GPU"
	@echo "     and the environment variable ocl_platform"
	@echo "         ex: export ocl_platform='NVIDIA'"
	@echo "     and the environment variable ocl_devices"
	@echo "         ex: export ocl_devices='K20'"

CLEANFILES =

all: head_message foot_message

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(pythondir)/report.py
