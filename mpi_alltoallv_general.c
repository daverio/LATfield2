#include <stdio.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <string.h>
#include <mpi.h>

#define NUM_BARRIERS 4
#define MY_REQUEST_MAX 100000
#define MY_HANDLE_MAX 10

int shmemid;
char volatile *shmem=NULL;

void **comm_code[MY_HANDLE_MAX+1];

int setup_shared_memory(MPI_Comm comm, int my_cores_per_node_row, MPI_Comm comm2, int my_cores_per_node_column, int size_shared, int *shmemid, char volatile **shmem){
    MPI_Comm my_comm_node;
    int my_mpi_rank_row, my_mpi_size_row, my_lrank_row, my_node, type_size, my_mpi_rank_column, my_mpi_size_column, my_lrank_column, my_lrank_node, my_mpi_size_global, my_mpi_rank_global;
    if ((*shmem)!=NULL){
        return 1;
    }
    MPI_Comm_size(comm, &my_mpi_size_row);
    MPI_Comm_rank(comm, &my_mpi_rank_row);
    if (comm2!=MPI_COMM_NULL){
        MPI_Comm_size(comm2, &my_mpi_size_column);
        MPI_Comm_rank(comm2, &my_mpi_rank_column);
    }else{
        my_mpi_size_column=1;
        my_mpi_rank_column=0;
    }
    MPI_Comm_split(comm, my_mpi_rank_row/my_cores_per_node_row, my_mpi_rank_row%my_cores_per_node_row, &my_comm_node);
    if ((my_mpi_rank_row%my_cores_per_node_row==0)&&(my_mpi_rank_column%my_cores_per_node_column==0)){
        (*shmemid)=shmget(IPC_PRIVATE, size_shared, IPC_CREAT | 0666);
    }
    MPI_Bcast(shmemid, 1, MPI_INT, 0, my_comm_node);
    MPI_Barrier(my_comm_node);
    MPI_Comm_free(&my_comm_node);
    if (comm2!=MPI_COMM_NULL){
        MPI_Comm_split(comm2, my_mpi_rank_column/my_cores_per_node_column, my_mpi_rank_column%my_cores_per_node_column, &my_comm_node);
        MPI_Bcast(shmemid, 1, MPI_INT, 0, my_comm_node);
        MPI_Barrier(my_comm_node);
        MPI_Comm_free(&my_comm_node);
    }
    (*shmem)=(char*)shmat(*shmemid, NULL, 0);
    if ((*shmem)==NULL) exit(2);
    if (!((my_mpi_rank_row%my_cores_per_node_row==0))&&(my_mpi_rank_column%my_cores_per_node_column==0)){
        (*shmemid)=-1;
    }
    MPI_Barrier(comm);
    if (comm2!=MPI_COMM_NULL){
        MPI_Barrier(comm2);
        MPI_Barrier(comm);
    }
    return 0;
}

void setup_rank_translation(MPI_Comm comm, int my_cores_per_node_row, MPI_Comm comm2, int my_cores_per_node_column, int *global_ranks){
    MPI_Comm my_comm_node;
    int my_mpi_size_row, grank, my_mpi_size_column, my_mpi_rank_column, *lglobal_ranks;
    MPI_Comm_size(comm, &my_mpi_size_row);
    if (comm2!=MPI_COMM_NULL){
        MPI_Comm_size(comm2, &my_mpi_size_column);
        MPI_Comm_rank(comm2, &my_mpi_rank_column);
        MPI_Comm_split(comm2, my_mpi_rank_column/my_cores_per_node_column, my_mpi_rank_column%my_cores_per_node_column, &my_comm_node);
        MPI_Comm_rank(MPI_COMM_WORLD, &grank);
        lglobal_ranks=(int*)malloc(sizeof(int)*my_cores_per_node_column);
        MPI_Gather(&grank, 1, MPI_INT, lglobal_ranks, 1, MPI_INT, 0, my_comm_node);
        MPI_Bcast(lglobal_ranks, my_cores_per_node_column, MPI_INT, 0, my_comm_node);
        MPI_Barrier(my_comm_node);
        MPI_Comm_free(&my_comm_node);
        MPI_Gather(lglobal_ranks, my_cores_per_node_column, MPI_INT, global_ranks, my_cores_per_node_column, MPI_INT, 0, comm);
        free(lglobal_ranks);
    }else{
        MPI_Comm_rank(MPI_COMM_WORLD, &grank);
        MPI_Gather(&grank, 1, MPI_INT, global_ranks, 1, MPI_INT, 0, comm);
    }
    MPI_Bcast(global_ranks, my_mpi_size_row*my_cores_per_node_column, MPI_INT, 0, comm);
}

int MY_Alltoallv_init(void *sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void *recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm, int my_cores_per_node_row, MPI_Comm comm2, int my_cores_per_node_column, int max_message_size_node){
    int my_mpi_rank_row, my_mpi_size_row, my_lrank_row, my_node, type_size, my_mpi_rank_column, my_mpi_size_column, my_lrank_column, my_lrank_node, my_mpi_size_global, my_mpi_rank_global;
    int ip, handle;
    char volatile *my_shared_sendbuf, *my_shared_recvbuf;
    int *global_ranks, i, j, k, l, m, port, new_counts_displs, add, isize, my_size_shared_sendbuf, my_size_shared_recvbuf;
    int lshmemid;
    int volatile *lshmem_sendcounts, *lshmem_recvcounts, *lshmem=NULL;
    MPI_Type_size(sendtype, &type_size);
    if (shmem==NULL){
        for (i=0; i<MY_HANDLE_MAX; i++){
            comm_code[i]=NULL;
        }
    }
    handle=0;
    while ((comm_code[handle]!=NULL)&&handle<MY_HANDLE_MAX){
        handle++;
    }
    if (handle>=MY_HANDLE_MAX){
        return(-1);
    }
    comm_code[handle]=(void**)malloc(sizeof(void *)*1000000);
    MPI_Comm_size(comm, &my_mpi_size_row);
    MPI_Comm_rank(comm, &my_mpi_rank_row);
    if (comm2!=MPI_COMM_NULL){
        MPI_Comm_size(comm2, &my_mpi_size_column);
        MPI_Comm_rank(comm2, &my_mpi_rank_column);
    }else{
        my_mpi_size_column=1;
        my_mpi_rank_column=0;
    }
    my_node=my_mpi_rank_row/my_cores_per_node_row;
    my_lrank_row=my_mpi_rank_row%my_cores_per_node_row;
    my_lrank_column=my_mpi_rank_column%my_cores_per_node_column;
    my_lrank_node=my_lrank_column*my_cores_per_node_row+my_lrank_row;
    my_mpi_size_global=my_mpi_size_row*my_cores_per_node_column;
    my_mpi_rank_global=my_mpi_rank_row*my_cores_per_node_column+my_mpi_rank_column%my_cores_per_node_column;
    new_counts_displs=(sdispls==NULL);
    if (new_counts_displs){
        sdispls=(int*)malloc(my_mpi_size_row*sizeof(int));
        recvcounts=(int*)malloc(my_mpi_size_row*sizeof(int));
        rdispls=(int*)malloc(my_mpi_size_row*sizeof(int));
        MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, comm);
        sdispls[0]=0;
        rdispls[0]=0;
        for (i=0; i<my_mpi_size_row-1; i++){
            sdispls[i+1]=sdispls[i]+sendcounts[i];
            rdispls[i+1]=rdispls[i]+recvcounts[i];
        }
    }
    setup_shared_memory(comm, my_cores_per_node_row, comm2, my_cores_per_node_column, my_mpi_size_row*my_cores_per_node_row*my_cores_per_node_column*2*sizeof(int), &lshmemid, (volatile char**)&lshmem);
    lshmem_sendcounts=lshmem;
    lshmem_recvcounts=lshmem+my_mpi_size_row*my_cores_per_node_row*my_cores_per_node_column;
    for (j=0; j<my_mpi_size_row/my_cores_per_node_row; j++){
        for (i=0; i<my_cores_per_node_row; i++){
            lshmem_sendcounts[(i+((my_mpi_size_row/my_cores_per_node_row+j-my_node)%(my_mpi_size_row/my_cores_per_node_row))*my_cores_per_node_row)*my_cores_per_node_row*my_cores_per_node_column+my_lrank_row*my_cores_per_node_column+my_lrank_column]=sendcounts[i+j*my_cores_per_node_row];
            lshmem_recvcounts[(my_lrank_row+((my_mpi_size_row/my_cores_per_node_row-j+my_node)%(my_mpi_size_row/my_cores_per_node_row))*my_cores_per_node_row)*my_cores_per_node_row*my_cores_per_node_column+i*my_cores_per_node_column+my_lrank_column]=recvcounts[i+j*my_cores_per_node_row];
        }
    }
    MPI_Barrier(comm);
    if (comm2!=MPI_COMM_NULL){
        MPI_Barrier(comm2);
        MPI_Barrier(comm);
    }
    my_size_shared_sendbuf=0;
    my_size_shared_recvbuf=0;
    for (i=0; i<my_mpi_size_row*my_cores_per_node_row*my_cores_per_node_column; i++){
        my_size_shared_sendbuf+=lshmem_sendcounts[i];
        my_size_shared_recvbuf+=lshmem_recvcounts[i];
    }
    my_size_shared_sendbuf*=type_size;
    my_size_shared_recvbuf*=type_size;
    i=my_size_shared_sendbuf+my_size_shared_recvbuf<=max_message_size_node;
    MPI_Allreduce(MPI_IN_PLACE, &i, 1, MPI_INT, MPI_MIN, comm);
    if (comm2!=MPI_COMM_NULL){
        MPI_Allreduce(MPI_IN_PLACE, &i, 1, MPI_INT, MPI_MIN, comm2);
    }
    if (!i) return(-1);
    if (!setup_shared_memory(comm, my_cores_per_node_row, comm2, my_cores_per_node_column, NUM_BARRIERS+max_message_size_node, &shmemid, &shmem)){
        MPI_Barrier(comm);
        if (comm2!=MPI_COMM_NULL){
            MPI_Barrier(comm2);
            MPI_Barrier(comm);
        }
        if ((my_mpi_rank_row%my_cores_per_node_row==0)&&(my_mpi_rank_column%my_cores_per_node_column==0)){
            memset((void *)shmem, 0, NUM_BARRIERS);
        }
        MPI_Barrier(comm);
        if (comm2!=MPI_COMM_NULL){
            MPI_Barrier(comm2);
            MPI_Barrier(comm);
        }
    }
    global_ranks=(int*)malloc(sizeof(int)*my_mpi_size_row*my_cores_per_node_column);
    setup_rank_translation(comm, my_cores_per_node_row, comm2, my_cores_per_node_column, global_ranks);
    my_shared_sendbuf=shmem+NUM_BARRIERS;
    my_shared_recvbuf=shmem+NUM_BARRIERS+my_size_shared_sendbuf;
    ip=0;
    if (my_mpi_size_row<=-my_cores_per_node_row){
        comm_code[handle][ip++]=(void *)(my_cores_per_node_row*my_cores_per_node_column);

        comm_code[handle][ip++]=NULL;

        add=0; m=0;
        for (i=0; i<my_lrank_row*my_cores_per_node_column+my_lrank_column; i++){
            add+=lshmem_sendcounts[m++];
        }
        for (i=0; i<my_cores_per_node_row; i++){
            if (i!=my_mpi_rank_row){
                comm_code[handle][ip++]=(void *)(my_shared_sendbuf+add*type_size);
                comm_code[handle][ip++]=(void *)(((char*)sendbuf)+sdispls[i]*type_size);
                comm_code[handle][ip++]=(void *)(sendcounts[i]*type_size);
            }
            for (j=0; j<my_cores_per_node_row*my_cores_per_node_column; j++){
                add+=lshmem_sendcounts[m++];
            }
        }
        comm_code[handle][ip++]=NULL;

        comm_code[handle][ip++]=NULL;

        comm_code[handle][ip++]=(void *)(((char *)recvbuf)+rdispls[my_mpi_rank_row]*type_size);
        comm_code[handle][ip++]=(void *)(((char *)sendbuf)+sdispls[my_mpi_rank_row]*type_size);
        comm_code[handle][ip++]=(void *)(sendcounts[my_mpi_rank_row]*type_size);
        for (i=0; i<my_cores_per_node_row; i++){
            add=0; m=0;
            for (j=0; j<i*my_cores_per_node_column+my_lrank_column; j++){
                add+=lshmem_recvcounts[m++];
            }
            for (k=0; k<my_lrank_row; k++){
                for (j=0; j<my_cores_per_node_row*my_cores_per_node_column; j++){
                    add+=lshmem_recvcounts[m++];
                }
            }
            if (i!=my_lrank_row){
                comm_code[handle][ip++]=(void *)(((char *)recvbuf)+rdispls[i]*type_size);
                comm_code[handle][ip++]=(void *)(my_shared_sendbuf+add*type_size);
                comm_code[handle][ip++]=(void *)(recvcounts[i]*type_size);
            }
        }
        comm_code[handle][ip++]=NULL;

        comm_code[handle][ip++]=NULL;

        comm_code[handle][ip++]=NULL;
    }else{
        comm_code[handle][ip++]=(void *)(my_cores_per_node_row*my_cores_per_node_column);

        for (port=my_lrank_node; port<my_mpi_size_row/my_cores_per_node_row-1; port+=my_cores_per_node_row*my_cores_per_node_column){
            add=0; m=0;
            for (i=0; i<my_cores_per_node_row*my_cores_per_node_row*my_cores_per_node_column*(port+1); i++){
                add+=lshmem_recvcounts[m++];
            }
            isize=0;
            for (i=0; i<my_cores_per_node_row*my_cores_per_node_row*my_cores_per_node_column; i++){
                isize+=lshmem_recvcounts[my_cores_per_node_row*my_cores_per_node_row*my_cores_per_node_column*(port+1)+i];
            }
            if (isize>0){
                comm_code[handle][ip++]=(void *)(((char*)my_shared_recvbuf)+add*type_size);
                comm_code[handle][ip++]=(void *)(global_ranks[(my_mpi_rank_global+my_mpi_size_global-(port+1)*my_cores_per_node_row*my_cores_per_node_column)%my_mpi_size_global]);
                comm_code[handle][ip++]=(void *)(isize*type_size);
            }
        }
        comm_code[handle][ip++]=NULL;

        add=0; m=0;
        for (i=0; i<my_lrank_row*my_cores_per_node_column+my_lrank_column; i++){
            add+=lshmem_sendcounts[m++];
        }
        for (k=0; k<my_mpi_size_row/my_cores_per_node_row; k++){
            for (i=0; i<my_cores_per_node_row; i++){
                j=(my_mpi_size_row/my_cores_per_node_row+k+my_node)%(my_mpi_size_row/my_cores_per_node_row);
                if (i+j*my_cores_per_node_row!=my_mpi_rank_row){
                    comm_code[handle][ip++]=(void *)(my_shared_sendbuf+add*type_size);
                    comm_code[handle][ip++]=(void *)(((char*)sendbuf)+sdispls[i+j*my_cores_per_node_row]*type_size);
                    comm_code[handle][ip++]=(void *)(sendcounts[i+j*my_cores_per_node_row]*type_size);
                }
                for (j=0; j<my_cores_per_node_row*my_cores_per_node_column; j++){
                    add+=lshmem_sendcounts[m++];
                }
            }
        }
        comm_code[handle][ip++]=NULL;

        for (port=my_lrank_node; port<my_mpi_size_row/my_cores_per_node_row-1; port+=my_cores_per_node_row*my_cores_per_node_column){
            add=0; m=0;
            for (i=0; i<my_cores_per_node_row*my_cores_per_node_row*my_cores_per_node_column*(port+1); i++){
                add+=lshmem_sendcounts[m++];
            }
            isize=0;
            for (i=0; i<my_cores_per_node_row*my_cores_per_node_row*my_cores_per_node_column; i++){
                isize+=lshmem_sendcounts[m++];
            }
            if (isize>0){
                comm_code[handle][ip++]=(void *)(((char*)my_shared_sendbuf)+add*type_size);
                comm_code[handle][ip++]=(void *)(global_ranks[(my_mpi_rank_global+(port+1)*my_cores_per_node_row*my_cores_per_node_column)%my_mpi_size_global]);
                comm_code[handle][ip++]=(void *)(isize*type_size);
            }
        }
        comm_code[handle][ip++]=NULL;

        comm_code[handle][ip++]=(void *)(((char *)recvbuf)+rdispls[my_mpi_rank_row]*type_size);
        comm_code[handle][ip++]=(void *)(((char *)sendbuf)+sdispls[my_mpi_rank_row]*type_size);
        comm_code[handle][ip++]=(void *)(sendcounts[my_mpi_rank_row]*type_size);
        for (i=0; i<my_cores_per_node_row; i++){
            if (i!=my_lrank_row){
                add=0; m=0;
                for (j=0; j<i*my_cores_per_node_column+my_lrank_column; j++){
                    add+=lshmem_recvcounts[m++];
                }
                for (k=0; k<my_lrank_row; k++){
                    for (j=0; j<my_cores_per_node_row*my_cores_per_node_column; j++){
                        add+=lshmem_recvcounts[m++];
                    }
                }
                for (j=0; j<my_cores_per_node_row*my_cores_per_node_row*my_cores_per_node_column*(i/my_cores_per_node_row); j++){
                    add+=lshmem_recvcounts[m++];
                }
                j=(my_mpi_size_row/my_cores_per_node_row-0+my_node)%(my_mpi_size_row/my_cores_per_node_row);
                comm_code[handle][ip++]=(void *)(((char *)recvbuf)+rdispls[i+j*my_cores_per_node_row]*type_size);
                comm_code[handle][ip++]=(void *)(my_shared_sendbuf+add*type_size);
                comm_code[handle][ip++]=(void *)(recvcounts[i+j*my_cores_per_node_row]*type_size);
            }
        }
        comm_code[handle][ip++]=NULL;

        for (k=1; k<my_mpi_size_row/my_cores_per_node_row; k++){
            for (l=0; l<my_cores_per_node_row; l++){
                add=0; m=0;
                for (i=0; i<l*my_cores_per_node_column+my_lrank_column; i++){
                    add+=lshmem_recvcounts[m++];
                }
                for (j=0; j<my_lrank_row; j++){
                    for (i=0; i<my_cores_per_node_row*my_cores_per_node_column; i++){
                        add+=lshmem_recvcounts[m++];
                    }
                }
                for (i=0; i<my_cores_per_node_row*my_cores_per_node_row*my_cores_per_node_column*k; i++){
                    add+=lshmem_recvcounts[m++];
                }
                j=(my_mpi_size_row/my_cores_per_node_row-k+my_node)%(my_mpi_size_row/my_cores_per_node_row);
                comm_code[handle][ip++]=(void *)(((char *)recvbuf)+rdispls[l+j*my_cores_per_node_row]*type_size);
                comm_code[handle][ip++]=(void *)(((char *)my_shared_recvbuf)+add*type_size);
                comm_code[handle][ip++]=(void *)(recvcounts[l+j*my_cores_per_node_row]*type_size);
            }
        }
        comm_code[handle][ip++]=NULL;
        comm_code[handle][ip++]=NULL;
    }
    if (new_counts_displs){
        free(rdispls);
        free(recvcounts);
        free(sdispls);
    }
    shmdt((const void*)lshmem);
    if (lshmemid!=-1){
        shmctl(lshmemid, IPC_RMID, NULL);
    }
    free(global_ranks);
    return(handle);
}

int MY_Alltoall_done(int handle){
    int i;
    free(comm_code[handle]);
    comm_code[handle]=NULL;
    for (i=0; i<MY_HANDLE_MAX; i++){
        if (comm_code[i]!=NULL){
            return(0);
        }
    }
    shmdt((const void*)shmem);
    if (shmemid!=-1){
        shmctl(shmemid, IPC_RMID, NULL);
    }
    shmem=NULL;
    shmemid=-1;
    return(0);
}

int barrier_count=0;

void node_barrier(int num_cores){
    __sync_fetch_and_add(shmem+barrier_count, 1);
    while (shmem[barrier_count]!=num_cores) {;}
    shmem[(barrier_count+NUM_BARRIERS-1)%NUM_BARRIERS]=0;
    barrier_count=(barrier_count+1)%NUM_BARRIERS;
}

int MY_Alltoall(int handle){
    MPI_Status my_status[MY_REQUEST_MAX];
    MPI_Request my_request[MY_REQUEST_MAX];
    void **lcomm_code = comm_code[handle];
    int num_cores;
    num_cores=*((int*)(lcomm_code));
    while (*lcomm_code!=NULL){
        lcomm_code++;
        int num_comm=0;
        while (*lcomm_code!=NULL){
            MPI_Irecv(*lcomm_code, *((int*)(lcomm_code+2)), MPI_CHAR, *((int*)(lcomm_code+1)), 0, MPI_COMM_WORLD, my_request+num_comm);
            lcomm_code+=3;
            num_comm++;
        }
        lcomm_code++;
        while (*lcomm_code!=NULL){
            memcpy(*lcomm_code, *(lcomm_code+1), *((int*)(lcomm_code+2)));
            lcomm_code+=3;
        }
        lcomm_code++;
        node_barrier(num_cores);
        while (*lcomm_code!=NULL){
            MPI_Isend(*lcomm_code, *((int*)(lcomm_code+2)), MPI_CHAR, *((int*)(lcomm_code+1)), 0, MPI_COMM_WORLD, my_request+num_comm);
            lcomm_code+=3;
            num_comm++;
        }
        lcomm_code++;
        while (*lcomm_code!=NULL){
            memcpy(*lcomm_code, *(lcomm_code+1), *((int*)(lcomm_code+2)));
            lcomm_code+=3;
        }
        lcomm_code++;
        MPI_Waitall(num_comm, my_request, my_status);
        node_barrier(num_cores);
        while (*lcomm_code!=NULL){
            memcpy(*lcomm_code, *(lcomm_code+1), *((int*)(lcomm_code+2)));
            lcomm_code+=3;
        }
        lcomm_code++;
    }
    return(0);
}
