cmake_minimum_required(VERSION 3.12)
project(PDAF VERSION 2.1.0 LANGUAGES Fortran)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Find MPI package
# use customised FindMPI.cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
find_package(MPI REQUIRED)
# Add the library target
if(WIN32)
  add_library(mpimod STATIC ${MPI_Fortran_MODULE_DIR}/mpi.f90)
  set_source_files_properties(${MPI_Fortran_MODULE_DIR}/mpi.f90 PROPERTIES COMPILE_OPTIONS "")
  target_include_directories(mpimod PUBLIC ${MPI_Fortran_INCLUDE_PATH})
endif()

# # Customize compiler flags
# set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
# if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
#    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
#    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
# elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
#    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -N 1023")
# elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
#    if(WIN32)
#       set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /4R8")
#       set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
#    else()
#       set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
#       set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -mkl")
#       #KB For newer versions of CMake - this is the way to go
#       enable_language(C)
#       set(BLA_VENDOR Intel10_64lp_seq)
#       #KB https://cmake.org/cmake/help/latest/module/FindBLAS.html
#   endif()
# endif()

set(MOD_PDAF 
    PDAF_timer.F90
    PDAF_memcount.F90
    PDAF_mod_filtermpi.F90
    PDAF_mod_filter.F90
)

# Generic routines in PDAF
set(OBJ_PDAF_GEN 
    PDAF_init.F90
    PDAF_init_filters.F90
    PDAF_alloc_filters.F90
    PDAF_print_info.F90
    PDAF_print_version.F90
    PDAF_communicate_ens.F90
    PDAF_set_comm_pdaf.F90
    PDAF_options_filters.F90
    PDAF_get_state.F90
    PDAF_incremental.F90
    PDAF_set_forget.F90
    PDAF_set_forget_local.F90
    PDAF_add_increment.F90
    PDAF_generate_rndmat.F90
    PDAF_local_weights.F90
    PDAF_local_weight.F90
    PDAF_force_analysis.F90
    PDAF_set_memberid.F90
    PDAF_get_memberid.F90
    PDAF_get_obsmemberid.F90
    PDAF_smoother_shift.F90
    PDAF_smoother.F90
    PDAF_smoother_local.F90
    PDAF_set_smootherens.F90
    PDAF_get_smootherens.F90
    PDAF_set_ens_pointer.F90
    PDAF_put_state_prepost.F90
    PDAF_assimilate_prepost.F90
    PDAF_prepost.F90
    PDAF_inflate_ens.F90
    PDAF_sampleens.F90
    PDAF_mvnormalize.F90
    PDAF_eofcovar.F90
    PDAF_diag_histogram.F90
    PDAF_diag_ensstats.F90
    PDAF_diag_effsample.F90
    PDAF_diag_crps.F90
    PDAF_gather_dim_obs_f.F90
    PDAF_gather_obs_f.F90
    PDAF_gather_obs_f2.F90
    PDAF_gather_obs_f_flex.F90
    PDAF_gather_obs_f2_flex.F90
    PDAF_allreduce.F90
    PDAF_deallocate.F90
    PDAF_get_assim_flag.F90
    PDAF_get_localfilter.F90
    PDAF_get_globalobs.F90
    PDAF_inflate_weights.F90
    PDAFomi_put_state_global.F90
    PDAFomi_put_state_local.F90
    PDAFomi_assimilate_global.F90
    PDAFomi_assimilate_local.F90
    PDAF_reset_forget.F90
    PDAF_get_ensstats.F90
    PDAF_set_debug_flag.F90
)

# Specific PDAF-routines for SEIK
set(OBJ_SEIK 
    PDAF_seik_init.F90
    PDAF_seik_alloc.F90
    PDAF_seik_options.F90
    PDAF_seik_memtime.F90
    PDAF_put_state_seik.F90
    PDAF_assimilate_seik.F90
    PDAF_seik_update.F90
    PDAF_seik_analysis.F90
    PDAF_seik_resample.F90
    PDAF_seik_analysis_newT.F90
    PDAF_seik_resample_newT.F90
    PDAF_seik_analysis_trans.F90
    PDAF_seik_matrixT.F90
    PDAF_seik_uinv.F90
    PDAF_seik_omega.F90
    PDAF_seik_TtimesA.F90
)

# Specific PDAF-routines for local SEIK
set(OBJ_LSEIK 
    PDAF_lseik_init.F90
    PDAF_lseik_alloc.F90
    PDAF_lseik_options.F90
    PDAF_lseik_memtime.F90
    PDAF_put_state_lseik.F90
    PDAF_assimilate_lseik.F90
    PDAF_lseik_update.F90
    PDAF_lseik_analysis.F90
    PDAF_lseik_resample.F90
    PDAF_lseik_analysis_trans.F90
)

# Specific PDAF-routines for SEEK
set(OBJ_SEEK 
    PDAF_seek_init.F90
    PDAF_seek_alloc.F90
    PDAF_seek_options.F90
    PDAF_seek_memtime.F90
    PDAF_put_state_seek.F90
    PDAF_assimilate_seek.F90
    PDAF_seek_update.F90
    PDAF_seek_analysis.F90
    PDAF_seek_rediag.F90
)

# Specific PDAF-routines for EnKF
set(OBJ_ENKF 
    PDAF_enkf_init.F90
    PDAF_enkf_alloc.F90
    PDAF_enkf_options.F90
    PDAF_enkf_memtime.F90
    PDAF_put_state_enkf.F90
    PDAF_assimilate_enkf.F90
    PDAF_enkf_update.F90
    PDAF_enkf_obs_ensemble.F90
    PDAF_enkf_analysis_rlm.F90
    PDAF_enkf_analysis_rsm.F90
    PDAF_enkf_omega.F90
    PDAF_enkf_Tleft.F90
    PDAF_smoother_enkf.F90
)

# Specific PDAF-routines for ETKF
set(OBJ_ETKF 
    PDAF_etkf_init.F90
    PDAF_etkf_alloc.F90
    PDAF_etkf_options.F90
    PDAF_etkf_memtime.F90
    PDAF_put_state_etkf.F90
    PDAF_assimilate_etkf.F90
    PDAF_etkf_update.F90
    PDAF_etkf_analysis.F90
    PDAF_etkf_analysis_T.F90
    PDAF_etkf_analysis_fixed.F90
    PDAF_etkf_Tright.F90
    PDAF_etkf_Tleft.F90
)

# Specific PDAF-routines for LETKF
set(OBJ_LETKF 
    PDAF_letkf_init.F90
    PDAF_letkf_alloc.F90
    PDAF_letkf_options.F90
    PDAF_letkf_memtime.F90
    PDAF_put_state_letkf.F90
    PDAF_assimilate_letkf.F90
    PDAF_letkf_update.F90
    PDAF_letkf_analysis.F90
    PDAF_letkf_analysis_T.F90
    PDAF_letkf_analysis_fixed.F90
)

# Specific PDAF-routines for ESTKF
set(OBJ_ESTKF 
    PDAF_estkf_init.F90
    PDAF_estkf_alloc.F90
    PDAF_estkf_options.F90
    PDAF_estkf_memtime.F90
    PDAF_put_state_estkf.F90
    PDAF_assimilate_estkf.F90
    PDAF_estkf_update.F90
    PDAF_estkf_analysis.F90
    PDAF_estkf_analysis_fixed.F90
    PDAF_estkf_AOmega.F90
    PDAF_estkf_OmegaA.F90
)

# Specific PDAF-routines for LESTKF
set(OBJ_LESTKF 
    PDAF_lestkf_init.F90
    PDAF_lestkf_alloc.F90
    PDAF_lestkf_options.F90
    PDAF_lestkf_memtime.F90
    PDAF_put_state_lestkf.F90
    PDAF_assimilate_lestkf.F90
    PDAF_lestkf_update.F90
    PDAF_lestkf_analysis.F90
    PDAF_lestkf_analysis_fixed.F90
)

# Specific PDAF-routines for EnKF
set(OBJ_LENKF 
    PDAF_lenkf_init.F90
    PDAF_lenkf_alloc.F90
    PDAF_lenkf_options.F90
    PDAF_lenkf_memtime.F90
    PDAF_put_state_lenkf.F90
    PDAFomi_put_state_lenkf.F90
    PDAF_assimilate_lenkf.F90
    PDAFomi_assimilate_lenkf.F90
    PDAF_lenkf_update.F90
    PDAF_enkf_gather_resid.F90
    PDAF_lenkf_analysis_rsm.F90
)

# Specific PDAF-routines for NETF
set(OBJ_NETF 
    PDAF_netf_init.F90
    PDAF_netf_alloc.F90
    PDAF_netf_options.F90
    PDAF_netf_memtime.F90
    PDAF_put_state_netf.F90
    PDAF_assimilate_netf.F90
    PDAF_netf_update.F90
    PDAF_netf_smootherT.F90
    PDAF_netf_analysis.F90
    PDAF_smoother_netf.F90
)

# Specific PDAF-routines for LNETF
set(OBJ_LNETF 
    PDAF_lnetf_init.F90
    PDAF_lnetf_alloc.F90
    PDAF_lnetf_options.F90
    PDAF_lnetf_memtime.F90
    PDAF_put_state_lnetf.F90
    PDAF_assimilate_lnetf.F90
    PDAF_lnetf_update.F90
    PDAF_lnetf_analysis.F90
    PDAF_lnetf_smootherT.F90
    PDAF_smoother_lnetf.F90
)

# Specific PDAF-routines for PF
set(OBJ_PF 
    PDAF_pf_init.F90
    PDAF_pf_alloc.F90
    PDAF_pf_options.F90
    PDAF_pf_memtime.F90
    PDAF_put_state_pf.F90
    PDAF_assimilate_pf.F90
    PDAF_pf_update.F90
    PDAF_pf_analysis.F90
    PDAF_pf_resampling.F90
    PDAF_pf_add_noise.F90
)

# Specific PDAF-routines for LKNETF
set(OBJ_LKNETF 
    PDAF_lknetf_init.F90
    PDAF_lknetf_alloc.F90
    PDAF_lknetf_options.F90
    PDAF_lknetf_memtime.F90
    PDAF_put_state_lknetf.F90
    PDAF_assimilate_lknetf.F90
    PDAF_lknetf_update.F90
    PDAF_lknetf_analysis_T.F90
    PDAF_lknetf_step_update.F90
    PDAF_lknetf_ana_lnetf.F90
    PDAF_lknetf_ana_letkfT.F90
    PDAF_lknetf_compute_gamma.F90
    PDAF_lknetf_set_gamma.F90
    PDAF_lknetf_alpha_neff.F90
    PDAF_lknetf_reset_gamma.F90
)

# Specific PDAF-routines for generating observations
set(OBJ_OBSGEN 
    PDAF_genobs_init.F90
    PDAF_genobs_alloc.F90
    PDAF_genobs_options.F90
    PDAF_put_state_generate_obs.F90
    PDAFomi_put_state_generate_obs.F90
    PDAF_generate_obs.F90
    PDAFomi_generate_obs.F90
    PDAF_gen_obs.F90
)

# Specific PDAF-routines for 3DVAR initialization part
set(OBJ_3DVAR_INI 
    PDAF_3dvar_init.F90
    PDAF_3dvar_alloc.F90
    PDAF_3dvar_options.F90
)

# Specific PDAF-routines for 3DVAR
set(OBJ_3DVAR 
    PDAF_3dvar_memtime.F90
    PDAF_put_state_3dvar.F90
    PDAF_assimilate_3dvar.F90
    PDAF_3dvar_update.F90
    PDAF_3dvar_analysis_cvt.F90
    PDAF_3dvar_optim_lbfgs.F90
    PDAF_3dvar_optim_cgplus.F90
    PDAF_3dvar_costf_cvt.F90
    PDAF_3dvar_costf_cg_cvt.F90
    PDAF_3dvar_optim_cg.F90
    PDAF_put_state_en3dvar_lestkf.F90
    PDAF_assimilate_en3dvar_lestkf.F90
    PDAF_en3dvar_update_lestkf.F90
    PDAF_put_state_en3dvar_estkf.F90
    PDAF_assimilate_en3dvar_estkf.F90
    PDAF_en3dvar_update_estkf.F90
    PDAF_en3dvar_analysis_cvt.F90
    PDAF_en3dvar_optim_lbfgs.F90
    PDAF_en3dvar_optim_cgplus.F90
    PDAF_en3dvar_optim_cg.F90
    PDAF_en3dvar_costf_cvt.F90
    PDAF_en3dvar_costf_cg_cvt.F90
    PDAF_put_state_hyb3dvar_lestkf.F90
    PDAF_assimilate_hyb3dvar_lestkf.F90
    PDAF_hyb3dvar_update_lestkf.F90
    PDAF_hyb3dvar_analysis_cvt.F90
    PDAF_put_state_hyb3dvar_estkf.F90
    PDAF_assimilate_hyb3dvar_estkf.F90
    PDAF_hyb3dvar_update_estkf.F90
    PDAF_hyb3dvar_optim_lbfgs.F90
    PDAF_hyb3dvar_optim_cgplus.F90
    PDAF_hyb3dvar_optim_cg.F90
    PDAF_hyb3dvar_costf_cvt.F90
    PDAF_hyb3dvar_costf_cg_cvt.F90
    PDAFomi_assimilate_3dvar.F90
    PDAFomi_assimilate_en3dvar_estkf.F90
    PDAFomi_assimilate_en3dvar_lestkf.F90
    PDAFomi_assimilate_hyb3dvar_estkf.F90
    PDAFomi_assimilate_hyb3dvar_lestkf.F90
    PDAFomi_put_state_3dvar.F90
    PDAFomi_put_state_en3dvar_estkf.F90
    PDAFomi_put_state_en3dvar_lestkf.F90
    PDAFomi_put_state_hyb3dvar_estkf.F90
    PDAFomi_put_state_hyb3dvar_lestkf.F90
)

# Routines for PDAF-OMI
set(OBJ_PDAFOMI
    PDAFomi_obs_f.F90
    PDAFomi_obs_l.F90
    PDAFomi_obs_op.F90
    PDAFomi.F90
    PDAFomi_callback.F90
)

# External optimizer libraries
set(OBJ_OPTIM
    ../external/CG+_mpi/cgfam.f
    ../external/CG+_mpi/cgsearch.f
    ../external/CG+/cgfam.f
    ../external/CG+/cgsearch.f
    ../external/LBFGS/lbfgsb.f
    ../external/LBFGS/linpack.f
    ../external/LBFGS/timer.f
)

# List all source files
set(SRC_FILES
    ${MOD_PDAF}
    ${OBJ_PDAF_GEN}
    ${OBJ_SEIK}
    ${OBJ_LSEIK}
    ${OBJ_SEEK}
    ${OBJ_ENKF}
    ${OBJ_ETKF}
    ${OBJ_LETKF}
    ${OBJ_ESTKF}
    ${OBJ_LESTKF}
    ${OBJ_LENKF}
    ${OBJ_NETF}
    ${OBJ_LNETF}
    ${OBJ_PF}
    ${OBJ_LKNETF}
    ${OBJ_OBSGEN}
    ${OBJ_3DVAR_INI}
    ${OBJ_3DVAR}
    ${OBJ_PDAFOMI}
    ${OBJ_OPTIM}
    ${OBJ_OPTIM}
)

# Customize compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    set_source_files_properties (${SRC_FILES} PROPERTIES COMPILE_OPTIONS "-O3;-ffree-line-length-none -fdefault-real-8")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
    set_source_files_properties (${SRC_FILES} PROPERTIES COMPILE_OPTIONS "-O3;-N;1023")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   if(WIN32)
      set_source_files_properties (${SRC_FILES} PROPERTIES COMPILE_OPTIONS "/O3;/4R8;/Qdiag-disable:10448")
      set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
   else()
      # set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
      set_source_files_properties (${SRC_FILES} PROPERTIES COMPILE_OPTIONS "-O3;-r8")
      set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -mkl")
      #KB For newer versions of CMake - this is the way to go
      enable_language(C)
      set(BLA_VENDOR Intel10_64lp_seq)
      #KB https://cmake.org/cmake/help/latest/module/FindBLAS.html
  endif()
endif()

# Define a target for pdaf-var
add_library(pdaf-var STATIC ${SRC_FILES})

# Set properties for pdaf-var
target_include_directories(pdaf-var PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pdaf-var PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(pdaf-var PUBLIC ${MPI_Fortran_INCLUDE_PATH})

set_target_properties(pdaf-var PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib
)
set_target_properties(pdaf-var PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include)