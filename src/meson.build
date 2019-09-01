
sources = [
    'main.f90',
    'atom.f90',
    'cla.f90',
    'desp.f90',
    'group.f90',
    'helanal.f90',
    'interp_hel.f90',
    'io.f90',
    'kinds.f90',
    'linalg.f90',
    'mem_mod.f90',
    'precision.f90',
    'stringmod.f90',
    'utilities.f90']

# openmp = dependency('openmp', language : 'fortran', required : false)

# if not openmp.found()
#  subdir_done()
# endif

framework = ''
biHelix_deps = []

if get_option('framework') != ''
    framework = get_option('framework')
endif

if get_option('blas') == 'openblas'
    biHelix_deps += fc.find_library('blas')
endif

if get_option('lapack') == 'netlib'
    lapackdep = dependency('lapack-netlib', required: true)
    lapacklib = fc.find_library('lapack')
    if lapackdep.found()
        biHelix_deps += lapackdep
    else
        biHelix_deps += lapacklib
    endif
    biHelix_deps += fc.find_library('lapacke')
    biHelix_deps += fc.find_library('gfortran')
endif

if get_option('blas') == 'mkl' or get_option('lapack') == 'mkl'
    if (fc.get_id() == 'intel')
        if framework == 'openmp'
            libmkl = declare_dependency(link_args: ['-mkl=parallel'])
        else
            libmkl = declare_dependency(link_args: ['-mkl=sequential'])
        endif
    else
        libmkl = declare_dependency(
            link_args: [
                '-lmkl_intel_lp64',
				'-lmkl_intel_thread',
				'-lmkl_core',
				'-liomp5',
				'-lpthread',
				'-lm',
				'-ldl',
            ]
        )
    endif
    biHelix_deps += libmkl
endif

if get_option('blas') == '' and get_option('lapack') == ''
    blas = fc.find_library('blas', required: true)
    if not blas.found()
        error('MESON_SKIP_TEST: BLAS Library not available')
    endif
    biHelix_deps += blas
    lapack = fc.find_library('lapack', required: true)
    if not lapack.found()
        error('MESON_SKIP_TEST: LAPACK Library not available')
    endif
    biHelix_deps += lapack
endif

bindir = get_option('bindir')
if bindir == ''
   bindir = join_paths(meson.source_root(), 'bin')
endif

exe_name = join_paths(bindir, meson.project_name())
inc = include_directories(join_paths('..','include'))



exe = executable(
    exe_name,
    sources,
    include_directories: inc,
    dependencies: [biHelix_deps])
test('biHelix',exe)
