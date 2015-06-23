Some short explanation to compile and run examples:

=========================================
gaussian.cpp

A small example, which fill a gaussian and compute derivative of it (Grad and Laplacian).
It also test the output functionality.

To compile with HDF5
mpic++ -o latgauss gaussian.cpp -lhdf5 -I../ 

To compile without HDF5
mpic++ -o latgauss gaussian.cpp -I../  -DWITHOUT_HDF5

To run:
mpirun -np 4 ./latgauss -n 2 -m 2 -o test

=========================================
poissonSolver.cpp

A small example which solve the poisson equation. This example need FFTW3.xx installed. (FFTW2.xx should work too). The executable will output the error of the solver.

To compile:
mpic++ -o poisson poissonSolver.cpp -lfftw3 -DFFT3D  -I../  -DWITHOUT_HDF5

To run:
mpirun -np 4 ./poisson -n 2 -m 2 -b 128

=========================================
IOserver.cpp

A small example of the usage of the IOserver. The executable will create some file named testFileXXX.dat (which are for this example plain text file)

To compile:
mpic++ -o ioserver IOserver.cpp -DEXTERNAL_IO -DWITHOUT_HDF5 -I../ 

To run:
mpirun -np 24 ./ioserver -n 4 -m 4 -i 8 -g 2

