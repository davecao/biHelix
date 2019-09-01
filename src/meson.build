
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

blas = fc.find_library('blas', required: true)
if not blas.found()
    error('MESON_SKIP_TEST: BLAS Library not available')
endif

lapack = fc.find_library('lapack', required: true)
if not lapack.found()
    error('MESON_SKIP_TEST: LAPACK Library not available')
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
    dependencies: [blas, lapack])
test('biHelix',exe)
