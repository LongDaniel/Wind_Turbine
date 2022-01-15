include Makefile.inc

ICFLAGS= -I./include -I$(FFTW_INC) -I$(HDF5_INC)
LDFLAGS= -L$(FFTW_LIB) -lfftw3 -L$(HDF5_LIB) -lhdf5_fortran -lhdf5 ./lib/lib2decomp_fft.a ./lib/libfll.a

#Define the fixed form source files .f
#INDSRC=hos_param.f90 fft_hos.f90
#DPDSRC=spectral_hos.f90 hos.f90 main.f90 test.f90 
#F77SRC=

SRCS = hos_param.f90 param.f90 poisson.f90 lib_array.f90 \
       hdf5_hos.f90 grid.f90 fft_hos.f90 common.f90 wavecontrol.f90 \
       spectral_hos.f90 turbine_model.f90 navier.f90 io.f90 \
       smooth_hos.f90 hos.f90 io_hos.f90 main_hos.f90 post.f90

MAIN1 = main.f90

FLL = lib/libfll.a
DECOMP = lib/lib2decomp_fft.a

OBJS := $(patsubst %.f90,%.o,$(SRCS))

#Define the name of the executable file
#EXECUTABLE=main.exe
EXECUTABLE=fst

.PHONY: default
default:
	cd fll_simplified && $(MAKE)  
	cd 2decomp && $(MAKE) && $(MAKE)  
	$(MAKE) fst

fst: $(FLL) $(DECOMP) $(OBJS) main.f90
	$(FC) $(FFLAGS) -o $@ main.f90 $(OBJS) $(LDFLAGS) $(ICFLAGS)
	@echo "Compilation finished."

.PHONY: all
all:  $(FLL) $(DECOMP) $(EXECUTABLE)

$(OBJS): $(DECOMP)

$(DECOMP): $(FLL) 
	cd 2decomp && $(MAKE) && $(MAKE) install 

$(fll):
	cd fll_simplified && $(MAKE) && $(MAKE) install 

%.o : %.f90
	$(FC) $(FFLAGS) -c $< $(ICFLAGS)

.PHONY: objects
objects: $(DECOMP) $(OBJS)


#	-mv *.mod include/
#	-mv *.o lib/

#hos_param.o param.o poisson.o fft_hos.o grid.f90 wavecontrol.o common.o navier.o smooth_hos.f90: 2decomp
hdf5_hos.o: hos_param.o
spectral_hos.o: hos_param.o fft_hos.o
io_hos.o:  hos_param.o spectral_hos.o hdf5_hos.o hos.o
main_hos.o : io_hos.o hos.o
hos.o: param.o hos_param.o fft_hos.o smooth_hos.o spectral_hos.o

smooth_hos.o: param.o hos_param.o spectral_hos.o fft_hos.o
turbine_model.o: param.o spectral_hos.o
navier.o: param.o common.o hos_param.o turbine_model.o
io.o: navier.o
post.o: io.o
main.o: post.o


.PHONY: clean
clean:
	rm -f *.o *.mod *.MOD lib/* include/* fst
	$(MAKE) -C 2decomp clean
	$(MAKE) -C fll_simplified clean

.PHONY: lightclean
lightclean:
	rm -f *.o *.mod *.MOD lib/* include/*
	
