project('LATfield2','cpp','c',
        default_options : ['cpp_std=c++17',
        'warning_level=3'],
        version: '1.1')
add_global_arguments(['-Ofast','-march=native' ],language : 'cpp')

cpp=meson.get_compiler('cpp')
mpi=dependency('mpi',language: 'cpp')
boost = dependency('boost')
boost_mpi = dependency('boost',modules: ['mpi','serialization'])
boost_ut = dependency('boost',modules: ['unit_test_framework'])
mpirun = find_program('mpirun')
fftw3 = dependency('fftw3')
hdf5 = dependency('hdf5')

include=include_directories('.')

headers = files(
    'Imag.hpp',
    'int2string.hpp',
    'LATfield2_Catalyst.hpp',
    'LATfield2_Field_decl.hpp',
    'LATfield2_Field.hpp',
    'LATfield2.hpp',
    'LATfield2_IO_server.hpp',
    'LATfield2_Lattice_decl.hpp',
    'LATfield2_Lattice.hpp',
    'LATfield2_parallel2d_decl.hpp',
    'LATfield2_parallel2d.hpp',
    'LATfield2_PlanFFT_decl.hpp',
    'LATfield2_PlanFFT.hpp',
    'LATfield2_save_hdf5.h',
    'LATfield2_save_hdf5_pixie.h',
    'LATfield2_SettingsFile.hpp',
    'LATfield2_Site_decl.hpp',
    'LATfield2_Site.hpp',
    'looping_macro.hpp',
    'timer.hpp')

subdir('particles')
subdir('benchmarks_tests')
subdir('examples')

subdir('tests')
