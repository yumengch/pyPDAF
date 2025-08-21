project('pyPDAF', 'cython', 'fortran', 'c', version: '0.1.0',
  meson_version: '>=1.5.0', default_options: ['b_asneeded=false'])

fs = import('fs')

py_mod = import('python')
python = py_mod.find_installation(pure: false)
dep_py = python.dependency()

cc = meson.get_compiler('c')
ff = meson.get_compiler('fortran')

mpimod = get_option('mpi_mod')

if host_machine.system() == 'windows'
  mpi = dependency('mpi', language : 'fortran', required: true)
else
  mpi = dependency('', language : 'fortran', required: false)
endif

blas_lib = get_option('blas_lib')
incdirs = get_option('incdirs')
libdirs = get_option('libdirs')

blas_lib_dep = []
foreach lib : blas_lib
  blas_lib_dep += ff.find_library(lib, dirs:libdirs, required : true)
endforeach

blas_dep = declare_dependency(
  dependencies: blas_lib_dep,
  include_directories: incdirs,
)

c_cython_args = cc.get_supported_arguments(['/O2', '/GL', '-O3'], checked: 'warn')
cython_args = ['-Xboundscheck=False', '-Xbinding=False', '-Xwraparound=False', '-Xinitializedcheck=False']

is_flang = ff.get_id() in ['llvm-flang', 'flang', 'flang-new']
fruntime_dep = dependency('', required : false)
if is_flang
  fortran_ext_args = ['-O0', '-g', ]
  fortran_args = ['-O0', '-g', '-fdefault-real-8', '-DUSE_PDAF']  # no probe
  fortran_runtime_lib = ff.find_library('FortranRuntime', dirs: libdirs, required : true)
  fortran_decimal_lib = ff.find_library('FortranDecimal', dirs: libdirs, required : true)
  fruntime_dep = declare_dependency(
    dependencies: [fortran_runtime_lib, fortran_decimal_lib],
  )
else
  fortran_args = ff.get_supported_arguments(
    ['-O3', '-fdefault-real-8', '-DUSE_PDAF'], checked: 'warn')
  fortran_ext_args = ff.get_supported_arguments(
    ['-O3', ], checked: 'warn')
  # fortran_args = ff.get_supported_arguments([
  #   '-O0','-g','-fcheck=all','-fbacktrace','-Wall','-Wextra',
  #   '-fdefault-real-8','-DUSE_PDAF'
  # ], checked: 'warn')
endif

incdir_numpy = run_command(python,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()


fortran_sources = files('PDAF/src/PDAF3_assimilate_3dvars.F90',
                        'PDAF/src/PDAF3_assimilate_3dvars_nondiagR.F90',
                        'PDAF/src/PDAF3_assimilate_ens.F90',
                        'PDAF/src/PDAF3_assimilate_ens_nondiagR.F90',
                        'PDAF/src/PDAF3_assim_offline_3dvars.F90',
                        'PDAF/src/PDAF3_assim_offline_3dvars_nondiagR.F90',
                        'PDAF/src/PDAF3_assim_offline_ens.F90',
                        'PDAF/src/PDAF3_assim_offline_ens_nondiagR.F90',
                        'PDAF/src/PDAF_3dvar_analysis_cvt.F90',
                        'PDAF/src/PDAF_3dvar.F90',
                        'PDAF/src/PDAF_3dvar_optim.F90',
                        'PDAF/src/PDAF_3dvar_update.F90',
                        'PDAF/src/PDAF3_put_state_3dvars.F90',
                        'PDAF/src/PDAF3_put_state_3dvars_nondiagR.F90',
                        'PDAF/src/PDAF3_put_state_ens.F90',
                        'PDAF/src/PDAF3_put_state_ens_nondiagR.F90',
                        'PDAF/src/PDAF_analysis_utils.F90',
                        'PDAF/src/PDAF_assimilate_3dvar.F90',
                        'PDAF/src/PDAF_assimilate_en3dvar_estkf.F90',
                        'PDAF/src/PDAF_assimilate_en3dvar_lestkf.F90',
                        'PDAF/src/PDAF_assimilate_enkf.F90',
                        'PDAF/src/PDAF_assimilate_ensrf.F90',
                        'PDAF/src/PDAF_assimilate_estkf.F90',
                        'PDAF/src/PDAF_assimilate_etkf.F90',
                        'PDAF/src/PDAF_assimilate_hyb3dvar_estkf.F90',
                        'PDAF/src/PDAF_assimilate_hyb3dvar_lestkf.F90',
                        'PDAF/src/PDAF_assimilate_lenkf.F90',
                        'PDAF/src/PDAF_assimilate_lestkf.F90',
                        'PDAF/src/PDAF_assimilate_letkf.F90',
                        'PDAF/src/PDAF_assimilate_lknetf.F90',
                        'PDAF/src/PDAF_assimilate_lnetf.F90',
                        'PDAF/src/PDAF_assimilate_lseik.F90',
                        'PDAF/src/PDAF_assimilate_netf.F90',
                        'PDAF/src/PDAF_assimilate_pf.F90',
                        'PDAF/src/PDAF_assimilate_prepost.F90',
                        'PDAF/src/PDAF_assimilate_seik.F90',
                        'PDAF/src/PDAF_assim_interfaces.F90',
                        'PDAF/src/PDAF_cb_procedures.F90',
                        'PDAF/src/PDAF_comm_obs.F90',
                        'PDAF/src/PDAF_communicate_ens.F90',
                        'PDAF/src/PDAF_da.F90',
                        'PDAF/src/PDAF_diag.F90',
                        'PDAF/src/PDAF_en3dvar_analysis_cvt.F90',
                        'PDAF/src/PDAF_en3dvar_optim.F90',
                        'PDAF/src/PDAF_en3dvar_update.F90',
                        'PDAF/src/PDAF_enkf_analysis_rlm.F90',
                        'PDAF/src/PDAF_enkf_analysis_rsm.F90',
                        'PDAF/src/PDAF_enkf.F90',
                        'PDAF/src/PDAF_enkf_update.F90',
                        'PDAF/src/PDAF_ensrf_analysis.F90',
                        'PDAF/src/PDAF_ensrf.F90',
                        'PDAF/src/PDAF_ensrf_update.F90',
                        'PDAF/src/PDAF_estkf_analysis.F90',
                        'PDAF/src/PDAF_estkf_analysis_fixed.F90',
                        'PDAF/src/PDAF_estkf.F90',
                        'PDAF/src/PDAF_estkf_update.F90',
                        'PDAF/src/PDAF_etkf_analysis.F90',
                        'PDAF/src/PDAF_etkf_analysis_fixed.F90',
                        'PDAF/src/PDAF_etkf_analysis_T.F90',
                        'PDAF/src/PDAF_etkf.F90',
                        'PDAF/src/PDAF_etkf_update.F90',
                        'PDAF/src/PDAF.F90',
                        'PDAF/src/PDAF_forecast.F90',
                        'PDAF/src/PDAF_generate_obs.F90',
                        'PDAF/src/PDAF_generate_obs_update.F90',
                        'PDAF/src/PDAF_genobs.F90',
                        'PDAF/src/PDAF_get.F90',
                        'PDAF/src/PDAF_get_state.F90',
                        'PDAF/src/PDAF_hyb3dvar_analysis_cvt.F90',
                        'PDAF/src/PDAF_hyb3dvar_optim.F90',
                        'PDAF/src/PDAF_hyb3dvar_update.F90',
                        'PDAF/src/PDAF_iau.F90',
                        'PDAF/src/PDAF_info.F90',
                        'PDAF/src/PDAF_init.F90',
                        'PDAF/src/PDAF_lenkf_analysis_rsm.F90',
                        'PDAF/src/PDAF_lenkf.F90',
                        'PDAF/src/PDAF_lenkf_update.F90',
                        'PDAF/src/PDAF_lestkf_analysis.F90',
                        'PDAF/src/PDAF_lestkf_analysis_fixed.F90',
                        'PDAF/src/PDAF_lestkf.F90',
                        'PDAF/src/PDAF_lestkf_update.F90',
                        'PDAF/src/PDAF_letkf_analysis.F90',
                        'PDAF/src/PDAF_letkf_analysis_fixed.F90',
                        'PDAF/src/PDAF_letkf_analysis_T.F90',
                        'PDAF/src/PDAF_letkf.F90',
                        'PDAF/src/PDAF_letkf_update.F90',
                        'PDAF/src/PDAF_lknetf_analysis_step.F90',
                        'PDAF/src/PDAF_lknetf_analysis_sync.F90',
                        'PDAF/src/PDAF_lknetf.F90',
                        'PDAF/src/PDAF_lknetf_update_step.F90',
                        'PDAF/src/PDAF_lknetf_update_sync.F90',
                        'PDAF/src/PDAF_lnetf_analysis.F90',
                        'PDAF/src/PDAF_lnetf.F90',
                        'PDAF/src/PDAF_lnetf_update.F90',
                        'PDAF/src/PDAFlocal_assimilate_3dvars.F90',
                        'PDAF/src/PDAFlocal_assimilate_ens.F90',
                        'PDAF/src/PDAFlocal_callback.F90',
                        'PDAF/src/PDAFlocal.F90',
                        'PDAF/src/PDAFlocalomi_assimilate_3dvars.F90',
                        'PDAF/src/PDAFlocalomi_assimilate_ens.F90',
                        'PDAF/src/PDAFlocalomi_put_state_3dvars.F90',
                        'PDAF/src/PDAFlocalomi_put_state_ens.F90',
                        'PDAF/src/PDAFlocal_put_state_3dvars.F90',
                        'PDAF/src/PDAFlocal_put_state_ens.F90',
                        'PDAF/src/PDAF_lseik_analysis.F90',
                        'PDAF/src/PDAF_lseik_analysis_trans.F90',
                        'PDAF/src/PDAF_lseik.F90',
                        'PDAF/src/PDAF_lseik_update.F90',
                        'PDAF/src/PDAF_memcount.F90',
                        'PDAF/src/PDAF_mod_core.F90',
                        'PDAF/src/PDAF_mod_parallel.F90',
                        'PDAF/src/PDAF_netf_analysis.F90',
                        'PDAF/src/PDAF_netf.F90',
                        'PDAF/src/PDAF_netf_update.F90',
                        'PDAF/src/PDAFobs.F90',
                        'PDAF/src/PDAFomi_assimilate_3dvars.F90',
                        'PDAF/src/PDAFomi_assimilate_ens.F90',
                        'PDAF/src/PDAFomi_assimilate_ens_nondiagR.F90',
                        'PDAF/src/PDAFomi_callback.F90',
                        'PDAF/src/PDAFomi_dim_obs_l.F90',
                        'PDAF/src/PDAFomi.F90',
                        'PDAF/src/PDAFomi_obs_diag.F90',
                        'PDAF/src/PDAFomi_obs_f.F90',
                        'PDAF/src/PDAFomi_obs_l.F90',
                        'PDAF/src/PDAFomi_obs_op.F90',
                        'PDAF/src/PDAFomi_put_state_3dvars.F90',
                        'PDAF/src/PDAFomi_put_state_ens.F90',
                        'PDAF/src/PDAFomi_put_state_ens_nondiagR.F90',
                        'PDAF/src/PDAF_pf_analysis.F90',
                        'PDAF/src/PDAF_pf.F90',
                        'PDAF/src/PDAF_pf_update.F90',
                        'PDAF/src/PDAF_prepost.F90',
                        'PDAF/src/PDAF_put_state_3dvar.F90',
                        'PDAF/src/PDAF_put_state_en3dvar_estkf.F90',
                        'PDAF/src/PDAF_put_state_en3dvar_lestkf.F90',
                        'PDAF/src/PDAF_put_state_enkf.F90',
                        'PDAF/src/PDAF_put_state_ensrf.F90',
                        'PDAF/src/PDAF_put_state_estkf.F90',
                        'PDAF/src/PDAF_put_state_etkf.F90',
                        'PDAF/src/PDAF_put_state_generate_obs.F90',
                        'PDAF/src/PDAF_put_state_hyb3dvar_estkf.F90',
                        'PDAF/src/PDAF_put_state_hyb3dvar_lestkf.F90',
                        'PDAF/src/PDAF_put_state_lenkf.F90',
                        'PDAF/src/PDAF_put_state_lestkf.F90',
                        'PDAF/src/PDAF_put_state_letkf.F90',
                        'PDAF/src/PDAF_put_state_lknetf.F90',
                        'PDAF/src/PDAF_put_state_lnetf.F90',
                        'PDAF/src/PDAF_put_state_lseik.F90',
                        'PDAF/src/PDAF_put_state_netf.F90',
                        'PDAF/src/PDAF_put_state_pf.F90',
                        'PDAF/src/PDAF_put_state_prepost.F90',
                        'PDAF/src/PDAF_put_state_seik.F90',
                        'PDAF/src/PDAF_sample.F90',
                        'PDAF/src/PDAF_seik_analysis.F90',
                        'PDAF/src/PDAF_seik_analysis_newT.F90',
                        'PDAF/src/PDAF_seik_analysis_trans.F90',
                        'PDAF/src/PDAF_seik.F90',
                        'PDAF/src/PDAF_seik_update.F90',
                        'PDAF/src/PDAF_set.F90',
                        'PDAF/src/PDAF_smoother.F90',
                        # 'PDAF/src/PDAF_timer.F90',
                        'PDAF/src/PDAF_timer_mpi.F90',
                        'PDAF/src/PDAF_utils.F90',
                        'PDAF/src/PDAF_utils_filters.F90',
                        'PDAF/external/SANGOMA/SANGOMA_quicksort.F90'
                        )
pdaf_ext_sources = files('PDAF/external/CG+/cgfam.f',
                         'PDAF/external/CG+/cgsearch.f',
                         'PDAF/external/CG+_mpi/cgfam.f',
                         'PDAF/external/CG+_mpi/cgsearch.f',
                         'PDAF/external/LBFGS/lbfgsb.f',
                         'PDAF/external/LBFGS/linpack.f',
                         'PDAF/external/LBFGS/timer.f'
                        )

pdafc_sources = files('src/fortran/pdaf3_c_assim.f90',
                      'src/fortran/pdaf3_c_put.f90',
                      'src/fortran/pdaf_c_assim.f90',
                      'src/fortran/pdaf_c_callback.f90',
                      'src/fortran/pdaf_c_cb_interface.f90',
                      'src/fortran/pdaf_c_diag.f90',
                      'src/fortran/pdaf_c.f90',
                      'src/fortran/pdaf_c_get.f90',
                      'src/fortran/pdaf_c_iau.f90',
                      'src/fortran/pdaf_c_iau_internal.f90',
                      'src/fortran/pdaf_c_internal.f90',
                      'src/fortran/pdaf_c_put.f90',
                      'src/fortran/pdaf_c_setter.f90',
                      'src/fortran/pdaflocal_c_assim.f90',
                      'src/fortran/pdaflocal_c.f90',
                      'src/fortran/pdaflocal_c_put.f90',
                      'src/fortran/pdaflocalomi_c_assim.f90',
                      'src/fortran/pdaflocalomi_c_put.f90',
                      'src/fortran/pdafomi_c_assim.f90',
                      'src/fortran/pdafomi_c_diag.f90',
                      'src/fortran/pdafomi_c.f90',
                      'src/fortran/pdafomi_c_internal.f90',
                      'src/fortran/pdafomi_c_legacy.f90',
                      'src/fortran/pdafomi_c_put.f90',
                      'src/fortran/pdafomi_c_setter.f90',
                      'src/fortran/pdaf_c_f_interface.f90'
                      )

if host_machine.system() == 'windows'
  pdaf_ext_sources += files(mpimod)
endif



pdaf_ext = static_library(
  'pdaf_ext',
  pdaf_ext_sources,
  dependencies : [mpi, blas_dep],
  fortran_args: fortran_ext_args
)

pdafc_lib = shared_library(
  'PDAFc',
   pdafc_sources + fortran_sources,
  dependencies : [mpi, blas_dep],
  include_directories: include_directories([incdir_numpy]),
  c_args: c_cython_args,
  fortran_args: fortran_args,
  cython_args : cython_args,
  link_language: 'fortran',
  link_whole: [pdaf_ext],
  install_dir: python.get_install_dir() / 'pyPDAF',
  install: true
)

cython_ext = python.extension_module(
  'pdaf_c_cb_interface',
  'src/pyPDAF/pdaf_c_cb_interface.pyx',
  link_with: [pdafc_lib],
  dependencies : [mpi, dep_py, blas_dep, fruntime_dep],
  include_directories: include_directories([incdir_numpy]),
  c_args: c_cython_args,
  fortran_args: fortran_args,
  cython_args : cython_args,
  link_language: 'fortran',
  subdir: 'pyPDAF',
  install: true
)

subdir('src/pyPDAF/PDAF')
subdir('src/pyPDAF/PDAF3')
subdir('src/pyPDAF/PDAFlocal')
subdir('src/pyPDAF/PDAFlocalomi')
subdir('src/pyPDAF/PDAFomi')

python.install_sources(['src/pyPDAF/__init__.py'], subdir: 'pyPDAF')
python.install_sources(['src/pyPDAF/PDAF/__init__.py'], subdir: 'pyPDAF/PDAF/')
python.install_sources(['src/pyPDAF/PDAFlocal/__init__.py'], subdir: 'pyPDAF/PDAFlocal/')
python.install_sources(['src/pyPDAF/PDAFlocalomi/__init__.py'], subdir: 'pyPDAF/PDAFlocalomi/')
python.install_sources(['src/pyPDAF/PDAFomi/__init__.py'], subdir: 'pyPDAF/PDAFomi/')
python.install_sources(['src/pyPDAF/PDAF3/__init__.py'], subdir: 'pyPDAF/PDAF3/')
