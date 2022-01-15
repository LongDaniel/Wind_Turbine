# find_package(MPI REQUIRED)
# include_directories(${MPI_Fortran_INCLUDE_PATH})

find_path(HDF5_INCLUDE_DIRS hdf5.mod PATHS ENV HDF5_INCLUDE_OPTS)
include_directories(${HDF5_INCLUDE_DIRS})

#set(FFTW_INCLUDE_DIRS /home/shenl/plyu/usr/fftw-3.3.5/build/include)

find_path (FFTW_INCLUDES NAMES fftw3.h PATHS ENV FFTW_INCLUDES)
find_library (FFTW_LIBRARIES NAMES fftw3 PATHS ENV FFTW_LIBRARIES)
#find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})
message(STATUS ${FFTW_INCLUDE_DIRS})

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(decomp2d_SRCS decomp.F90 fft_fftw3.f90 io.f90 spectral.f90 utils.f90
discontinuity_smooth.f90)

add_library(2decomp STATIC ${decomp2d_SRCS})
target_link_libraries(2decomp fll)
set_target_properties(2decomp PROPERTIES Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")

