
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

exe = executable('biHelix', sources, dependencies: [blas, lapack])
test('biHelix',exe)
