# find_package(MPI REQUIRED)
# include_directories(${MPI_Fortran_INCLUDE_PATH})

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(fll_SRCS fll_type.f90 fll_out.f90 fll_cat.f90 fll_stich.f90 fll_rm.f90 
    fll_mv.f90 fll_mk.f90 fll_duplicate.f90 fll_cp.f90 fll_deattach.f90 
    fll_funct_prt.f90 fll_getnbytes.f90 fll_locate.f90 fll_getndata.f90 
    fll_match_pattern.f90 fll_mkdir.f90 fll_nnodes.f90 fll_read.f90 
    fll_read_record.f90 fll_read_ffa.f90 fll_sweep.f90 fll_write.f90 
    fll_write_ffa.f90 fll_scan_file.f90 fll_read_ucd.f90 fll_rename.f90 
    fll_mods.f90)

add_library(fll STATIC ${fll_SRCS})
set_target_properties(fll PROPERTIES Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")

