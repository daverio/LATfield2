Some short explanation to compile and run examples:

=========================================
gettingStarted.cpp

Compile this example with e.g. mpic++:
mpic++ -o getStart gettingStarted.cpp -I../ -DHDF5 -lhdf5

This will compile using HDF5 for the Field I/O. If you have not HDF5 installed then you should compile with:

mpic++ -o getStart gettingStarted.cpp -I../

It can be executed using (here using "mpirun -np 4" to run with 4 process):

mpirun -np 4 ./getStart -n 2 -m 2

The executable will prompt the following text:

Parallel grid size: (2,2). 
Lattice size: (25,57,32);
Process ranks: 0,(0,0); Local lattice size: (25,28,16); First local point coordinate: (0,0,0).
Process ranks: 1,(1,0); Local lattice size: (25,28,16); First local point coordinate: (0,0,16).
Process ranks: 2,(0,1); Local lattice size: (25,29,16); First local point coordinate: (0,28,0).
Process ranks: 3,(1,1); Local lattice size: (25,29,16); First local point coordinate: (0,28,16).

The first line gives the size of the 2 dimension of the parallel grid. The second line gives the size of a Lattice. Then each process output its ranks and the description of the local part of the lattice (the part of the lattice which is stored on the given process).

=========================================
fft.cpp

Compile this example with e.g. mpic++:

Double precision:
mpic++ -o fft_exec fft.cpp -I../ -DFFT3D -lfftw 

Single precision:
mpic++ -o fft_exec fft.cpp -I../ -DFFT3D -DSINGLE -lfftwf

It can be executed using (here using "mpirun -np 4" to run with 4 processes):
mpirun -np 4 ./fft_exec -n 2 -m 2

The executable will not return anything, and the FFT are perform on fields which have not been assigned. It exist just to show the usage of the FFT wrapper of LATfield2.


=========================================
IOserver.cpp

Compile this example with e.g. mpic++:
mpic++ -o ioserver_exec IOserver.cpp -I../ -DEXTERNAL_IO

It can be executed using (here using "mpirun -np 4" to run with 4 process):
mpirun -np 24 ./ioserver_exec -n 4 -m 4 -i 8 -g 4

The n and m parameter are as usual parameters to initialize the parallel object. Then 2 additional parameters are passed. First -i which is the total number of MPI processes of the IO server, then -g is the number of IO process which write data in a single file. Therefor n*m+i need to be equal to the total number of MPI processes used by the job and i/g must be an integer.

