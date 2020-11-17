bench_files =[ 
#    'benchmarks.cpp',
#    'benchmarks_server.cpp',
    'benchmarks_particles.cpp',
    'wave_test_back.cpp',
    'wave_test.cpp']

foreach f : bench_files
    executable(
        f.split('.')[0],
        f,
        include_directories: [include],
        dependencies: [mpi,hdf5,fftw3])
endforeach

#subdir('particles')
