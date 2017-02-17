/*! \file LATfield2_save_hdf5.h
 \brief LATfield2_save_hdf5.h contains the definition of the function used for hdf5 i/o.
 \author David Daverio
*/

#ifndef LATFIELD2_SAVE_HDF5
#define LATFIELD2_SAVE_HDF5

extern "C"{

#include <math.h>
#include <hdf5.h>
#include <stdlib.h>
#include <string.h>
#include "int2string.hpp"


int save_hdf5_externC(char **data,
                      long file_offset[2],
                      int *size,
                      int * sizeLocal,
                      int halo,
                      int lat_dim,
                      int comp,
                      int slice_offset,
                      int slice_thickness,
                      hid_t array_type,
                      int array_size,
                      string  filename_str,
                      string dataset_name_str)
{
  hid_t file_id, plist_id,filespace,memspace,dset_id,dtype_id,dtbase_id,root_id;
  hsize_t  Asize;

  char * filename;
  filename = (char*)malloc((filename_str.size()+1)*sizeof(char));
  for(int i = 0;i<filename_str.size();i++)filename[i]=filename_str[i];
  filename[filename_str.size()] = '\0';

  char  dataset_name[128];
  for(int i = 0;i<filename_str.size();i++)dataset_name[i]=dataset_name_str[i];
  dataset_name[dataset_name_str.size()] = '\0';

  string dname_str;
  dname_str = "/"+dataset_name_str;

  herr_t status;

  hsize_t * sizeGlobal;
  sizeGlobal = new hsize_t[lat_dim];
  hsize_t * localSize;
  localSize = new hsize_t[lat_dim];
  hsize_t * offset;
  offset = new hsize_t[lat_dim];
  hsize_t * offsetf;
  offsetf = new hsize_t[lat_dim];
  hsize_t * count;
  count = new hsize_t[lat_dim];

  for(int i=0;i<lat_dim;i++)
  {
    sizeGlobal[i]=size[lat_dim -1 - i];
    localSize[i]=sizeLocal[lat_dim -1 - i]+2*halo;
    count[i]=sizeLocal[lat_dim -1 -i];
    offset[i]=halo;
    offsetf[i]=0;

  }
  offsetf[0]=file_offset[0];
  offsetf[1]=file_offset[1];
  offset[lat_dim -1]=halo+slice_offset;
  count[lat_dim -1]=slice_thickness;


  if(array_size ==1)
  {
    dtype_id = H5Tcopy(array_type);
  }
  else
  {
    Asize = array_size;
    dtbase_id = H5Tcopy(array_type);
    dtype_id = H5Tarray_create(dtbase_id,1,&Asize);
  }


#ifdef H5_HAVE_PARALLEL //Parallel version, H5_HAVE_PARALLEL definition is needed by hdf5 to run in parallel too

MPI_Comm comm  = parallel.lat_world_comm();
MPI_Info info  = MPI_INFO_NULL;


plist_id = H5Pcreate(H5P_FILE_ACCESS);
H5Pset_fapl_mpio(plist_id, comm, info);

file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

H5Pclose(plist_id);
if(file_id<0) return -1;





for(long c=0;c<comp;c++)
{
  if(comp!=1)dname_str = "/"+dataset_name_str+"_"+int2string(c,999);

  filespace = H5Screate_simple(lat_dim,sizeGlobal,NULL);
  memspace = H5Screate_simple(lat_dim,localSize,NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);;
  dset_id = H5Dcreate1(file_id, dname_str.c_str(), dtype_id, filespace,plist_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  if(dset_id<0)
  {
    H5Sclose(memspace);
    H5Fclose(file_id);
    return -2;
  }

  filespace = H5Dget_space(dset_id);
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, count, NULL);
  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data[c]);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);


}
H5Fclose(file_id);

#else // serial version, without H5_HAVE_PARALLEL definition hdf5 will run in serial !

int mpi_size,mpi_rank,p;
mpi_size = parallel.size();
mpi_rank = parallel.rank();

//cout<<"rank: "<<mpi_rank<<" , calling save HDF5 extern c serial"<<endl;

if(mpi_rank==0)
{


  plist_id = H5Pcreate(H5P_FILE_ACCESS);

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

  H5Pclose(plist_id);
  if(file_id<0) return -1;

  for(long c=0;c<comp;c++)
  {
    if(comp!=1)dname_str = "/"+dataset_name_str+"_"+int2string(c,999);

    //dname_str = "/test";

    filespace = H5Screate_simple(lat_dim,sizeGlobal,NULL);
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, lat_dim, sizeGlobal);
    dset_id = H5Dcreate1(file_id, dname_str.c_str(), dtype_id, filespace, plist_id);
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    if(dset_id<0) return -2;
  }

  H5Fclose(file_id);
}

MPI_Barrier(parallel.lat_world_comm());


for(p=0;p < mpi_size;p++)
{

  MPI_Barrier(parallel.lat_world_comm());
  if(mpi_rank==p)
  {

    //cout<<"rank: "<<mpi_rank<<" , writting data"<<endl;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fopen(filename,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);
    if(file_id<0) return -1;

    for(long c=0;c<comp;c++)
    {
      if(comp!=1)dname_str = "/"+dataset_name_str+"_"+int2string(c,999);

      dset_id = H5Dopen(file_id, dname_str.c_str(), H5P_DEFAULT);
      if(dset_id<0) return -3;
      filespace = H5Dget_space(dset_id);
      dtype_id = H5Dget_type(dset_id);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      memspace = H5Screate_simple(lat_dim,localSize,NULL);

      status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL,count, NULL);
      status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL,count, NULL);
      status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data[c]);

      H5Pclose(plist_id);
      H5Sclose(filespace);
      H5Sclose(memspace);
      H5Dclose(dset_id);
    }
    H5Fclose(file_id);

  }
  MPI_Barrier(parallel.lat_world_comm());

}



#endif

free(filename);
return 1;
}


int load_hdf5_externC(char **data,
                      long file_offset[2],
                      int *size,
                      int *sizeLocal,
                      int halo,
                      int lat_dim,
                      int comp,
                      string  filename_str,
                      string dataset_name_str)
{
  hid_t file_id, plist_id,filespace,memspace,dset_id,dtype_id,dtbase_id,group_id,root_id;


  char * filename;
  filename = (char*)malloc((filename_str.size()+1)*sizeof(char));
  for(int i = 0;i<filename_str.size();i++)filename[i]=filename_str[i];
  filename[filename_str.size()] = '\0';

  string dname_str;
  dname_str = "/"+dataset_name_str;

  herr_t status;

  hsize_t * sizeGlobal;
  sizeGlobal = new hsize_t[lat_dim];
  hsize_t * localSize;
  localSize = new hsize_t[lat_dim];
  hsize_t * offset;
  offset = new hsize_t[lat_dim];
  hsize_t * offsetf;
  offsetf = new hsize_t[lat_dim];
  hsize_t * count;
  count = new hsize_t[lat_dim];

  hsize_t haloSize = 2*halo;

  for(int i=0;i<lat_dim;i++)
  {
    sizeGlobal[i]=size[lat_dim -1 - i];
    localSize[i]=sizeLocal[lat_dim -1 - i]+haloSize;
    count[i]=sizeLocal[lat_dim -1 -i];
    offset[i]=halo;
    offsetf[i]=0;
  }
  offsetf[0]=file_offset[0];
  offsetf[1]=file_offset[1];

#ifdef H5_HAVE_PARALLEL	//Parallel version, H5_HAVE_PARALLEL definition is needed by hdf5 to run in parallel too


  MPI_Comm comm  = parallel.lat_world_comm();
  MPI_Info info  = MPI_INFO_NULL;


  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);

  file_id = H5Fopen(filename,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);


  for(long c=0;c<comp;c++)
  {
    if(comp!=1)dname_str = "/"+dataset_name_str+"_"+int2string(c,999);

    dset_id = H5Dopen(file_id, dname_str.c_str(), H5P_DEFAULT);
  	filespace = H5Dget_space(dset_id);
  	dtype_id = H5Dget_type(dset_id);
  	plist_id = H5Pcreate(H5P_DATASET_XFER);

  	memspace = H5Screate_simple(lat_dim,localSize,NULL);
    //// verifier si l "offset" ne doit pas tenire compte du ....
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, count, NULL);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
  	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  	status = H5Dread(dset_id, dtype_id, memspace, filespace, plist_id, data[c]);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

  }
  H5Fclose(file_id);

#else // serial version, without H5_HAVE_PARALLEL definition hdf5 will run in serial !
  int mpi_size,mpi_rank,p;
  mpi_size = parallel.size();
  mpi_rank = parallel.rank();

  for(p=0;p < mpi_size;p++)
  {
    if(mpi_rank==p)
    {
      plist_id = H5Pcreate(H5P_FILE_ACCESS);

      file_id = H5Fopen(filename,H5F_ACC_RDWR,plist_id);
      H5Pclose(plist_id);

      for(long c=0;c<comp;c++)
      {
        if(comp!=1)dname_str = "/"+dataset_name_str+"_"+int2string(c,999);

        dset_id = H5Dopen(file_id, dname_str.c_str(), H5P_DEFAULT);
        filespace = H5Dget_space(dset_id);
        dtype_id = H5Dget_type(dset_id);
        plist_id = H5Pcreate(H5P_DATASET_XFER);

        memspace = H5Screate_simple(lat_dim,localSize,NULL);

        int status;

        status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL,count, NULL);
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL,count, NULL);
        status = H5Dread(dset_id, dtype_id, memspace, filespace, plist_id, data[c]);


        H5Pclose(plist_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Dclose(dset_id);

      }


      H5Fclose(file_id);

    }
    MPI_Barrier(MPI_COMM_WORLD);

  }

#endif

free(filename);
return 1;




}


}

#endif
