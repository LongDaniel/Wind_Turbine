include Makefile.inc

ICFLAGS= -I./include -I$(FFTW_INC) -I$(HDF5_INC)
LDFLAGS= -L$(FFTW_LIB) -lfftw3 -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

SRCS =  fft_hos.f90 hdf5_hos.f90 hos.f90 hos_param.f90 io_hos.f90 \
        post.f90 smooth.f90 spec_src.f90 spectral_hos.f90 \
        utils_hos.f90

MAIN1 = main.f90

OBJS := $(patsubst %.f90,%.o,$(SRCS))

#Define the name of the executable file
#EXECUTABLE=main.exe
EXECUTABLE=main_hos

.PHONY: default
default:
	$(MAKE) main_hos 

main_hos: $(OBJS) main.f90
	$(FC) $(FFLAGS) -o $@ main.f90 $(OBJS) $(LDFLAGS) $(ICFLAGS)
	@echo "Compilation finished."

.PHONY: all
all:  $(EXECUTABLE)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< $(ICFLAGS)

#.PHONY: objects
#objects: $(OBJS)


#	-mv *.mod include/
#	-mv *.o lib/

#hos_param.o param.o poisson.o fft_hos.o grid.f90 wavecontrol.o common.o navier.o smooth_hos.f90: 2decomp
fft_hos.o hdf5_hos.f90 spectral_hos.o utils_hos.f90 spec_src.o: hos_param.o
spec_src.o: spectral_hos.o
spectral_hos.o utils_hos.o: fft_hos.o
hos.o: smooth.o
smooth.o: spectral_hos.o
io_hos.o: utils_hos.o
post.o: spec_src.o hos.o

.PHONY: clean
clean:
	rm -f *.o *.mod *.MOD main_hos 

.PHONY: lightclean
lightclean:
	rm -f *.o *.mod *.MOD 
	
