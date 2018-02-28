#ifndef MPI_ALLTOALLV_GENERAL_H_

#define MPI_ALLTOALLV_GENERAL_H_

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

int MY_Alltoallv_init(void *sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void *recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm, int my_cores_per_node_row, MPI_Comm comm2, int my_cores_per_node_column, int max_message_size_node);
int MY_Alltoall_done(int handle);
int MY_Alltoall(int handle);

#ifdef __cplusplus
}
#endif

#endif
