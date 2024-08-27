cmake_minimum_required(VERSION 3.12)
project(PDAF VERSION 2.2.1 LANGUAGES Fortran)

# Set the PDAF library version
set(PDAF_VERSION "2.2.1")

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/..)

set(MOD_PDAF
    PDAF_timer.F90
    PDAF_memcount.F90
    PDAF_mod_filtermpi.F90
    PDAF_mod_filter.F90
)

# Generic routines in PDAF
set(OBJ_PDAF_GEN
    PDAF_analysis_utils.F90
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
    PDAFomi_put_state_global_nondiagR.F90
    PDAFomi_put_state_nonlin_nondiagR.F90
    PDAFomi_put_state_local.F90
    PDAFomi_put_state_local_nondiagR.F90
    PDAFomi_assimilate_global.F90
    PDAFomi_assimilate_global_nondiagR.F90
    PDAFomi_assimilate_nonlin_nondiagR.F90
    PDAFomi_assimilate_local.F90
    PDAFomi_assimilate_local_nondiagR.F90
    PDAF_reset_forget.F90
    PDAF_get_ensstats.F90
    PDAF_set_debug_flag.F90
    PDAF_reset_dim_p.F90
    PDAF_set_offline_mode.F90
    PDAFlocal.F90
    PDAFlocal_set_indices.F90
    PDAFlocal_set_increment_weights.F90
    PDAFlocal_clear_increment_weights.F90
    PDAFlocal_g2l_cb.F90
    PDAFlocal_l2g_cb.F90
    PDAFlocalomi_assimilate.F90
    PDAFlocalomi_assimilate_nondiagR.F90
    PDAFlocalomi_put_state.F90
    PDAFlocalomi_put_state_nondiagR.F90
    PDAF_correlation_function.F90
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
    PDAFlocal_put_state_lseik.F90
    PDAFlocal_assimilate_lseik.F90
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
    PDAF_enkf_gather_resid.F90
    PDAF_enkf_analysis_rlm.F90
    PDAF_enkf_analysis_rsm.F90
    PDAF_enkf_omega.F90
    PDAF_enkf_Tleft.F90
    PDAF_smoother_enkf.F90
    PDAFomi_put_state_enkf_nondiagR.F90
    PDAFomi_assimilate_enkf_nondiagR.F90
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
    PDAFlocal_put_state_letkf.F90
    PDAFlocal_assimilate_letkf.F90
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
    PDAFlocal_put_state_lestkf.F90
    PDAFlocal_assimilate_lestkf.F90
)

# Specific PDAF-routines for EnKF
set(OBJ_LENKF
    PDAF_lenkf_init.F90
    PDAF_lenkf_alloc.F90
    PDAF_lenkf_options.F90
    PDAF_lenkf_memtime.F90
    PDAF_put_state_lenkf.F90
    PDAFomi_put_state_lenkf.F90
    PDAFomi_put_state_lenkf_nondiagR.F90
    PDAF_assimilate_lenkf.F90
    PDAFomi_assimilate_lenkf.F90
    PDAFomi_assimilate_lenkf_nondiagR.F90
    PDAF_lenkf_update.F90
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
    PDAF_mod_lnetf.F90
    PDAF_smoother_lnetf.F90
    PDAFomi_put_state_lnetf_nondiagR.F90
    PDAFomi_assimilate_lnetf_nondiagR.F90
    PDAFlocal_put_state_lnetf.F90
    PDAFlocal_assimilate_lnetf.F90
    PDAFlocalomi_assimilate_lnetf_nondiagR.F90
    PDAFlocalomi_put_state_lnetf_nondiagR.F90
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
    PDAFomi_put_state_lknetf_nondiagR.F90
    PDAFomi_assimilate_lknetf_nondiagR.F90
    PDAFlocal_put_state_lknetf.F90
    PDAFlocal_assimilate_lknetf.F90
    PDAFlocalomi_assimilate_lknetf_nondiagR.F90
    PDAFlocalomi_put_state_lknetf_nondiagR.F90
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
   PDAF_3dvar_memtime.F90
)

# Specific PDAF-routines for 3DVAR
set(OBJ_3DVAR
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
    PDAFlocal_put_state_en3dvar_lestkf.F90
    PDAFlocal_put_state_hyb3dvar_lestkf.F90
    PDAFlocal_assimilate_en3dvar_lestkf.F90
    PDAFlocal_assimilate_hyb3dvar_lestkf.F90
    PDAFomi_assimilate_3dvar.F90
    PDAFomi_assimilate_en3dvar_estkf.F90
    PDAFomi_assimilate_en3dvar_lestkf.F90
    PDAFomi_assimilate_hyb3dvar_estkf.F90
    PDAFomi_assimilate_hyb3dvar_lestkf.F90
    PDAFomi_assimilate_3dvar_nondiagR.F90
    PDAFomi_assimilate_en3dvar_estkf_nondiagR.F90
    PDAFomi_assimilate_en3dvar_lestkf_nondiagR.F90
    PDAFomi_assimilate_hyb3dvar_estkf_nondiagR.F90
    PDAFomi_assimilate_hyb3dvar_lestkf_nondiagR.F90
    PDAFlocalomi_assimilate_en3dvar_lestkf.F90
    PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR.F90
    PDAFlocalomi_assimilate_hyb3dvar_lestkf.F90
    PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR.F90
    PDAFomi_put_state_3dvar.F90
    PDAFomi_put_state_en3dvar_estkf.F90
    PDAFomi_put_state_en3dvar_lestkf.F90
    PDAFomi_put_state_hyb3dvar_estkf.F90
    PDAFomi_put_state_hyb3dvar_lestkf.F90
    PDAFomi_put_state_3dvar_nondiagR.F90
    PDAFomi_put_state_en3dvar_estkf_nondiagR.F90
    PDAFomi_put_state_en3dvar_lestkf_nondiagR.F90
    PDAFomi_put_state_hyb3dvar_estkf_nondiagR.F90
    PDAFomi_put_state_hyb3dvar_lestkf_nondiagR.F90
    PDAFlocalomi_put_state_en3dvar_lestkf.F90
    PDAFlocalomi_put_state_en3dvar_lestkf_nondiagR.F90
    PDAFlocalomi_put_state_hyb3dvar_lestkf.F90
    PDAFlocalomi_put_state_hyb3dvar_lestkf_nondiagR.F90
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

set(OBJ_SANGOMA
    ../external/SANGOMA/SANGOMA_quicksort.F90
)

set (OBJ_INTERFACE
    PDAFlocal_interfaces.F90
    PDAF_interfaces_module.F90
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
    ${OBJ_PDAFOMI}
    ${OBJ_INTERFACE}
    ${OBJ_3DVAR_INI}
    ${OBJ_SANGOMA}
)

# check if Config Path is provided
if(NOT Config_PATH)
    message(FATAL_ERROR "Config_PATH is not given. Providing a Config_PATH in commandline by specify -DConfig_PATH=PATH_TO_CONFIG_FILE")
endif()

# import configuration file
include(${Config_PATH})

if (${PDAF_NAME} STREQUAL "pdaf-var")
   set(SRC_FILES ${SRC_FILES} ${OBJ_3DVAR} ${OBJ_OPTIM})
endif()

if (WIN32)
    # This will just compile the MPI library modules
    # This is necessary for MS-MPI on windows and cannot be included in the os-dependent config file
    add_library(mpimod OBJECT ${MPI_Fortran_MODULE_SRC_FILE})
    target_include_directories(mpimod PUBLIC ${MPI_Fortran_MODULE_INCLUDE_PATH})
    set_target_properties(mpimod PROPERTIES
        Fortran_MODULE_DIRECTORY ${MPI_Fortran_MODULE_DIR})
endif()

# Define a target for pdaf library
add_library(${PDAF_NAME} STATIC ${SRC_FILES})

# Set include directory for MPI
target_include_directories(${PDAF_NAME} PUBLIC ${MPI_Fortran_INCLUDE_PATH})
target_include_directories(${PDAF_NAME} PUBLIC ${MPI_Fortran_MODULE_INCLUDE_PATH})
target_include_directories(${PDAF_NAME} PUBLIC ${MPI_Fortran_MODULE_DIR})
set_target_properties(${PDAF_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include)

# Move the static library and include files to given directories.
install(TARGETS ${PDAF_NAME}
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../lib
    ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../lib
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../bin  # For Windows DLLs
)