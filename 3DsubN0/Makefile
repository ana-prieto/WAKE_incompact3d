#=======================================================================
# Makefile for Imcompact3D
#=======================================================================

# Choose pre-processing options
#   -DSHM	   - enable shared-memory implementation
#   -DDOUBLE_PREC  - use double-precision
OPTIONS = -DDOUBLE_PREC

# Choose an FFT engine, available options are:
#   fftw3      - FFTW version 3.x
#   generic    - A general FFT algorithm (no 3rd-party library needed)
FFT= generic

# Paths to FFTW 3
FFTW3_PATH=   # full path of FFTW installation if using fftw3 engine above
FFTW3_INCLUDE = -I$(FFTW3_PATH)/include
FFTW3_LIB = -L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f

# Specify Fortran and C compiler names and flags here
# Normally, use MPI wrappers rather than compilers themselves 
# Supply a Fortran pre-processing flag together with optimisation level flags
# Some examples are given below:

#FC =  
#OPTFC = 
#CC = 
#CFLAGS = 

# PGI
#FC = ftn
#OPTFC = -fast -O3 -Mpreprocess
#CC = cc
#CFLAGS = -O3

# PathScale
#FC = ftn
#OPTFC = -Ofast -cpp
#CC = cc
#CFLAGS = -O3

# GNU
#FC = mpif90
#OPTFC =  -fpp -# xAVX -axSSE4.1 -xHost
#CC = cc
#CFLAGS = -O3
#fpp:  Fortran preprocessor
#xAVX: In Supermuc, generate Intel Advanced Vector Extension in SandyBridge
#-xSSE4.1 : turn on the auto-vectorizing
#cflag -03: Optimize at level 3. -0 is the optimization level.


#IBM
#FC=mpixlf95_r
#OPTFC = -O3 -qextname=flush -qsuffix=cpp=f90 -qarch=qp -qtune=qp -WF,-DDOUBLE_PREC
#CC = mpixlc_r
#CFLAGS = -O3
#Location: MPI compile scripts: /bgsys/drivers/ppcfloor/comm/xl/bin 
# Serial Fortran: /opt/ibmcmp/xlf/bg/14.1/bin


# Cray
#FC = ftn
#OPTFC = -e Fm
#CC = cc
#CFLAGS = 

# Archer
FC = ftn
OPTFC = -cpp -xAVX -axSSE4.1
CC = cc
CFLAGS = -O3

#-----------------------------------------------------------------------
# Normally no need to change anything below

# include PATH 
ifeq ($(FFT),generic)
  INC=
else ifeq ($(FFT),fftw3)
  INC=$(FFTW3_INCLUDE)
endif

# library path
ifeq ($(FFT),generic)
   LIBFFT=
else ifeq ($(FFT),fftw3)
   LIBFFT=$(FFTW3_LIB)
endif

# List of source files
SRC = decomp_2d.f90 glassman.f90 fft_$(FFT).f90 module_param.f90 io.f90 variables.f90 poisson.f90 schemes.f90 convdiff.f90 incompact3d.f90 navier.f90 filter.f90 derive.f90 parameters.f90 tools.f90 visu.f90

#-----------------------------------------------------------------------
# Normally no need to change anything below

ifneq (,$(findstring DSHM,$(OPTIONS)))
SRC := FreeIPC.f90 $(SRC)  
OBJ =	$(SRC:.f90=.o) alloc_shm.o FreeIPC_c.o
else
OBJ =	$(SRC:.f90=.o)
endif	

all: incompact3d

alloc_shm.o: alloc_shm.c
	$(CC) $(CFLAGS) -c $<

FreeIPC_c.o: FreeIPC_c.c
	$(CC) $(CFLAGS) -c $<

incompact3d : $(OBJ)
	$(FC) -O3 -o $@ $(OBJ) $(LIBFFT)

%.o : %.f90
	$(FC) $(OPTFC) $(OPTIONS) $(INC) -c $<

.PHONY: clean 
clean:
	rm -f *.o *.mod incompact3d

.PHONY: realclean
realclean: clean
	rm -f *~ \#*\#
