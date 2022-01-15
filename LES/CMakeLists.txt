add_subdirectory(fll_simplified)
add_subdirectory(2decomp)

set(fst_SRCS grid.f90 main.f90 io.f90 param.f90 navier.f90 common.f90
poisson.f90 wavecontrol.f90 fft_hos.f90 hos_param.f90 io_hos.f90
smooth_hos.f90 spectral_hos.f90 hos.f90 turbine_model.f90  main_hos.f90 hdf5_hos.f90 
lib_array.f90 post.f90)

# find_package(MPI REQUIRED)
# include_directories(${MPI_Fortran_INCLUDE_PATH})

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories("${PROJECT_BINARY_DIR}/mod")

add_executable(fst ${fst_SRCS})
add_dependencies(2decomp fll)
add_dependencies(fst 2decomp)

target_link_libraries(fst fll 2decomp ${FFTW_LIBRARIES} ${HDF5_Fortran_LIBRARIES})

#add_executable(post ${fst_SRCS})
#add_dependencies(post 2decomp)
#target_compile_definitions(post PRIVATE POST=1)
#
#target_link_libraries(post 2decomp ${FFTW_LIBRARIES} ${HDF5_Fortran_LIBRARIES})
