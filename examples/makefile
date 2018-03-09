all: ta2a

mpi_alltoallv_general.o: ../mpi_alltoallv_general.c
	cc -c ../mpi_alltoallv_general.c -O3

test_alltoall.o: test_alltoall.cpp $(wildcard ../*.hpp)
	CC -c test_alltoall.cpp -I../ -DFFT3D -DHDF5 -DH5_HAVE_PARALLEL -std=c++11 -DSINGLE -O3

ta2a: mpi_alltoallv_general.o test_alltoall.o
	CC -o ta2a mpi_alltoallv_general.o test_alltoall.o
