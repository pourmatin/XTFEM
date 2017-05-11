################################################################################
# Makefile for
#   XTFEM: Extended Space-Time Finite Element Method 3D (MPI)
#
# Developed by Rui Zhang, Jun. 2015, at UT Dallas, Richardson, TX.
#
# E-mail: rui.zhang4@utdallas.edu or ruizhang@mail.nwpu.edu.cn
#
# Supervised by Dr. Dong Qian,
#   Department of Mechanical Engineering, 
#   The University of Texas at Dallas, Richardson, Texas, USA
# and Dr. Lihua Wen,
#   school of Astronautics,
#   Northwestern Polytechnical University, Xi'an, Shaanxi, China
#
# All rights reserved (c) 2016
#_______________________________________________________________________________
#
# Makefile project supported on Linux and Mac OS X Platforms
#_______________________________________________________________________________
#
# Options:
#   all             build
#   run             run program in work directory
#   nohup_run       run program in background
#   clean           clean source directory
#   veryclean       clean source and work directory
#   archive_src     archive source files and copy to backup directory
#   archive_rst     archive result files and copy to backup directory
#
################################################################################

# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# These flags will override any settings
ifeq ($(i386),1)
	OS_SIZE = 32
	OS_ARCH = i686
endif

ifeq ($(x86_64),1)
	OS_SIZE = 64
	OS_ARCH = x86_64
endif

# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))

# Directory settings
WORK = /work/04161/hoss3in/lonestar
TOPDIR = $(WORK)
WORKDIR = $(TOPDIR)/workspace
LIBDIR = $(TOPDIR)/libs
METISDIR = $(LIBDIR)/metis-5.1.0-gcc
MUMPSDIR = $(LIBDIR)/MUMPS_5.0.1-gcc
SCADIR = $(LIBDIR)/scalapack-2.0.2
BLASDIR = $(LIBDIR)/BLAS-3.6.0
LAPDIR = $(LIBDIR)/lapack-3.6.1
MPIDIR = /opt/apps/gcc4_9/mvapich2/2.1/lib
FEMDIR = $(WORKDIR)/Thermal/XTFEM
SRCDIR = $(FEMDIR)/src
MESHDIR = $(FEMDIR)/inp
BACKUPDIR = $(FEMDIR)/backup
OOCDIR = $(TEMPDIR)/ooc

# Mesh File
# MESHFILE = prismaticrod-325-192.inp
# MESHFILE = prismaticrod-2025-1536.inp
# MESHFILE = prismaticrod-3751-3000.inp
# MESHFILE = prismaticrod-26901-24000.inp
 MESHFILE = sent_c3d8_1824_2435.inp
# MESHFILE = sent_c3d8_9380_10934.inp
# MESHFILE = sent_c3d8_20658_18160.inp
# MESHFILE = sent_c3d8_48516_45100.inp
# MESHFILE = sent_c3d8_133671_126650.inp
# MESHFILE = sent_c3d8_948693_924600.inp
# MESHFILE = panel_4424_2144.inp
# MESHFILE = panel_7634_3749.inp
# MESHFILE = panel_29286_14513.inp
# MESHFILE = panel_38318_18903.inp
# MESHFILE = panel_171015_113490.inp
# MESHFILE = panel_230685_152766.inp
# MESHFILE = panel_519324_343684.inp
# MESHFILE = panel_1689895_1346820.inp
# MESHFILE = panel98114_47923.inp
#MESHFILE = rod6171_5000.inp
#MESHFILE = rod4851_4000.inp
#MESHFILE = rod525_320.inp
#MESHFILE = rod99_40.inp
#MESHFILE = rod1.inp
#MESHFILE = rod5025_3200.inp
MESHINP = sent.inp

# file operations
CP = cp -f
MV = mv -f
RM = rm -fr
MKDIR = mkdir -p
ZIP = tar -czf
ZIPEXT = tar.gz

ifneq (,${findstring Cuda,${KOKKOS_DEVICES}})
   CXX = $(KOKKOSDIR)/bin/nvcc_wrapper
   LINK = ${CXX}
   LINKFLAGS = --std=c++11
   KOKKOS_FLAGS = -mavx
   EXE = cuda
   KOKKOS_DEVICES = "Cuda,OpenMP"
   CCFLAGS :=
   FCFLAGS := -fopenmp
   LDFLAGS := -fopenmp
   KOKKOS_ARCH = -arch=sm_35
   KOKKOSDIR = $(LIBDIR)/kokkos-lib-cuda
else
   CXX = mpicxx
   LINK = ${CXX}
   LINKFLAGS = --std=c++11
   KOKKOS_FLAGS = -mavx
   EXE = host
   KOKKOS_DEVICES = "OpenMP"
   CCFLAGS := -fopenmp
   FCFLAGS := -fopenmp
   LDFLAGS := -fopenmp
   KOKKOS_ARCH = -fopenmp
   KOKKOSDIR = $(LIBDIR)/kokkos-lib-openmp
endif

TEMPDIR = $(FEMDIR)/temp-${EXE}-${MESHFILE}

# The compiler and linker
FC = mpif90
CC = mpicc
LD = mpif90
NVCC = nvcc

# timestamp
TIMESTAMP = $(shell date +'%Y_%m_%d-%H_%M_%S')

# backup setting
BACKUPFILESRC = $(PROG)-src-$(TIMESTAMP).$(ZIPEXT)
BACKUPFILERST = $(PROG)-rst-$(TIMESTAMP).$(ZIPEXT)

# files that will be cleaned
CLEAN = *.o *~ *.mod
VERYCLEAN = *

# OS-specific build flags
ifneq ($(DARWIN),) 
      LDFLAGS   += -Xlinker -rpath
      CCFLAGS   += -arch $(OS_ARCH) 
      FCFLAGS   += -arch $(OS_ARCH)
else
  ifeq ($(OS_SIZE),32)
      LDFLAGS   += 
      CCFLAGS   += -m32
      FCFLAGS   += -m32
      NVCCFLAGS += -m32 
  else
      LDFLAGS   += 
      CCFLAGS   += -m64
      FCFLAGS   += -m64
      NVCCFLAGS += -m64
  endif
endif

# flags for compiler and linker options
ifneq (,${findstring debug,${Mode}})
   CCFLAGS += -g -O0
   FCFLAGS += -g -O0
   LDFLAGS += -g
else
   LDFLAGS +=
   FCFLAGS += -O3
   CCFLAGS += -O3
endif

# libraries

# SCALAPACK, BLACS, LAPACK, BLAS
#LIBS = -L$(SCADIR)
#LIBS += -L$(LAPDIR) -llapack
#LIBS += -L$(BLASDIR) -lblas

#MPI
#LIBS += $(MPIDIR)/libmpi.so $(MPIDIR)/libmpich.so -lpthread

#SCALAP should define the SCALAPACK and  BLACS libraries.
SCADIR = /work/04161/hoss3in/lonestar/libs/scalapack-2.0.2
LIBS  = -L$(SCADIR) -lscalapack
LIBS += -L$(LAPDIR) -llapack
# INCLUDE DIRECTORY FOR MPI
#MPIDIR = /opt/apps/gcc4_9/mvapich2/2.1/lib
#INCPAR = -I$(MPIDIR)

# LIBRARIES USED BY THE PARALLEL VERSION OF MUMPS: $(SCALAP) and MPI
#LIBS += $(SCALAP) -L$(MPIDIR) -lmpi
#LIBPAR = $(SCALAP)   $(MPIdir)/libmpi.so $(MPIdir)/libmpich.so

BLASDIR = /work/04161/hoss3in/lonestar/libs/BLAS-3.6.0
#LIBLAP = -L$(LAPDIR) -llapack
LIBS += -L$(BLASDIR) -lblas

# DEFINE YOUR PTHREAD LIBRARY
LIBS += -lpthread

# METIS
LMETIS    = $(METISDIR)/build/Linux-x86_64/libmetis/libmetis.a
IMETIS    = -I$(METISDIR)/include

# MUMPS
LIBMUMPS_COMMON = $(MUMPSDIR)/lib/libmumps_common.a
LIBDMUMPS = $(MUMPSDIR)/lib/libdmumps.a $(LIBMUMPS_COMMON)
INCLUDES += -I$(MUMPSDIR)/include

LORDERINGS = -L$(MUMPSDIR)/PORD/lib -lpord $(LMETIS)
INCLUDES += $(IMETIS)
INCLUDES += -I/usr/include/openmpi-x86_64/ $(IMETIS)

# KOKKOS
LKOKKOS	= -L$(KOKKOSDIR)/lib -L/opt/apps/cuda/6.5/lib64 -L/opt/apps/cuda/6.5/lib64/stubs  -lkokkos -lcudart -lcuda
IKOKKOS	= -I$(KOKKOSDIR)/include 
LIBS += -lrt -lm -lstdc++

# Program name
PROG = xtfem3d
BATCHFILE = $(EXE).slurm

# List of generated object files
# general
OBJS = kinds.o interface_definitions.o procs.o thermoelastic.o
# sparse matrix
OBJS += spmtx.o
# debugging
OBJS += debug.o
# MUMPS-based solver
OBJS += ilut.f kronecker.o kronecker_mumps.o mumps.o
# MPI
OBJS += mpisubs.o
# post-process and damage
OBJS += postprocess.o kokkos_damage.o kokkos_host.o
# element deletion
OBJS += delelement.o
# adddiagonal.o
# main and others
OBJS += stcoeff.o getfextc.o getdfinal.o elesubs.o main3d.o  

#include ${KOKKOS_PATH}/Makefile.kokkos

# all
all: build
build: ${PROG}

# Rue to link the program
${PROG}: $(OBJS)
	$(LD) $(LDFLAGS) $(KOKKOS_LDFLAGS) $(OBJS) $(KOKKOS_FLAGS) $(LKOKKOS) $(LIBDMUMPS) $(LORDERINGS) ${LIBS} -o ${PROG}.${EXE}
	$(MKDIR) $(TEMPDIR)
	cd $(TEMPDIR) && $(RM) $(VERYCLEAN)
	$(MKDIR) $(OOCDIR)
	$(CP) ${PROG}.${EXE} $(TEMPDIR)
	$(CP) $(BATCHFILE) $(TEMPDIR)/$(BATCHFILE)
	$(CP) $(MESHDIR)/$(MESHFILE) $(TEMPDIR)/$(MESHINP)


# General rules 
%.o: %.f
	$(FC) $(FCFLAGS) $(INCLUDES) -o $@ -c $<
%.o: %.f90
	$(FC) $(FCFLAGS) $(INCLUDES) -o $@ -c $<
%.o: %.c
	$(CC) $(CCFLAGS) ${INCLUDES} -o $@ -c $<
%.o: %.cpp
	$(CXX) $(CCFLAGS) $(IKOKKOS) $(LINKFLAGS) $(KOKKOS_FLAGS) $(KOKKOS_ARCH) $(CXXFLAGS) $(EXTRA_INC) -o $@ -c $<


# Utility targets	
archive_src: clean
	$(ZIP) $(BACKUPFILESRC) *
	$(MKDIR) $(BACKUPDIR)
	$(MV) $(BACKUPFILESRC) $(BACKUPDIR)/$(BACKUPFILESRC)
	
archive_rst:
	cd $(TEMPDIR) && $(ZIP)  $(BACKUPFILERST) *
	$(MKDIR) $(BACKUPDIR)
	$(MV) $(TEMPDIR)/$(BACKUPFILERST) $(BACKUPDIR)/$(BACKUPFILERST)

cleantemp:
	cd $(TEMPDIR) && $(RM) $(VERYCLEAN) *.cuda *.host
      
buildtemp: cleantemp
	$(CP) $(PROG) $(TEMPDIR)
	$(CP) $(BATCHFILE) $(TEMPDIR)/$(BATCHFILE)
	$(CP) $(MESHDIR)/$(MESHFILE) $(TEMPDIR)/$(MESHINP)
	$(MKDIR) $(OOCDIR)

.PHONY: clean veryclean
clean::
	$(RM) $(CLEAN) $(PROG)
veryclean: clean
	cd $(TEMPDIR) && $(RM) $(VERYCLEAN)
