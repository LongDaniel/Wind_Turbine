########################################################################
# Makefile for building the 2DECOMP&FFT library
########################################################################

# Normally, do not change anything here. Modify the platform-dependent
# Makefile.inc file instead

include ../Makefile.inc

SRCS = decomp.f90 fft_fftw3.f90 discontinuity_smooth.f90 spectral.f90 io.f90 utils.f90 

FLL_INCLUDE = ../include

OBJS =	$(SRCS:.f90=.o)
LIB = lib2decomp_fft.a

.PHONY: all
all: $(LIB)
	$(MAKE) install

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -I$(FLL_INCLUDE) -I$(HDF5_INC)

.PHONY: objects
objects: $(OBJS)

## special treatment for Intel MKL (need to build MKL mod files)
#MKL_MOD=
#ifeq ($(FFT),mkl)
#   MKL_MOD=mkl_mod
#endif

$(LIB): $(OBJS)
	ar qc $@ $(OBJS)

.PHONY: lib
lib: $(LIB)

#lib2decomp_fft.a: $(MKL_MOD) $(OBJS)
#	ar qc $@ $(OBJS)

.PHONY: install
install:  $(LIB)
	mkdir -p ../include
	mkdir -p ../lib	
	cp -u *.mod ../include
	cp -u lib*.a ../lib

#alloc_shm.o: alloc_shm.c
#	$(CC) $(CFLAGS) -c $<

#FreeIPC_c.o: FreeIPC_c.c
#	$(CC) $(CFLAGS) -c $<

#.PHONY: mkl_mod
#mkl_mod:
#	echo Building Intel MKL mod files...
#	$(IFORT) -c $(MKL_PATH)/include/mkl_dfti.f90

#.PHONY: examples
#examples: lib

#.PHONY: test
#test: examples

fft_fftw3.o: decomp.o

discontinuity_smooth.o: decomp.o

spectral.o: decomp.o fft_fftw3.o

io.o: decomp.o

utils.o: decomp.o

.PHONY: clean
clean:
	rm -f *.o *.mod lib*.a

.PHONY: realclean
realclean: clean
	rm -f *~ \#*\#

