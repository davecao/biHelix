project(
    'biHelix', 'fortran',
    version:'1.0.0',
    default_options:['default_library=static', 'buildtype=release'])

realbits = get_option('realbits')

# %% -------------
# Compiler flags
#-----------------

# find the fortran compiler
fc = meson.get_compiler('fortran')

if fc.get_id() == 'gcc'
    # add_global_arguments('-fimplicit-none', '-cpp', '-fbounds-check', '-Wextra', '-Wpedantic', language : 'fortran')
    add_global_arguments('-cpp', '-fbounds-check', '-Wextra', '-Wpedantic', language : 'fortran')
    if fc.version().version_compare('>=8.0')
        add_global_arguments('-std=f2018', language: 'fortran')
    endif
    oldargs = ['-std=legacy', '-w']

    # dialect = [ '-ffree-form','-std=f2003','-fimplicit-none']
    # bounds = '-fbounds-check'
elif fc.get_id() == 'intel'
    add_global_arguments('-stand', 'f18', '-warn declarations', language: 'fortran')
    oldargs = '-nowarn'
    # dialect = ['-stand f03','-free','-implicitnone']
    # bounds = '-check bounds'
elif fc.get_id() == 'intel-cl'
    # /fpp allows #include etc preprocessor lines
    add_global_arguments('/stand:f18', '/warn:declarations', '/fpp', language: 'fortran')
    oldargs = '/nowarn'
elif fc.get_id() == 'pgi'
    add_global_arguments('-C', '-Mdclchk', language: 'fortran')
    oldargs = '-w'
    # dialect = [ '-Mfreeform','-Mdclchk',
    #            '-Mstandard','-Mallocatable=03']
    # bounds = '-C'
elif fc.get_id() == 'flang'
    add_global_arguments('-W', language: 'fortran')
    oldargs = '-w'
endif

os = host_machine.system()

f18errorstop = fc.links('character :: b; error stop b; end', name: 'F2018 error stop')

includes = include_directories('include')
subdir('src')
# bindir = meson.source_root() # do not pass to sub meson.build

#program(
#    'biHelix',
#    ['src/main.f90'],
#    include_directories: includes,
#    link_with:
#)


#executable('biHelix', sources)
