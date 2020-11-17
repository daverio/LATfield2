examples_files =[ 
    'boundary_test.cpp',
    'fft.cpp',
    'gettingStarted.cpp',
    'particles.cpp',
    'particles_IOissue.cpp',
    'poissonSolver.cpp',
    'test_sum.cpp']

foreach f : examples_files
    executable(
        f.split('.')[0]+'.exe',
        f,
        include_directories: [include],
        dependencies: [mpi,hdf5,fftw3])
endforeach


