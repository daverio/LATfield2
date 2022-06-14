
/*! \file LATfield2_save_hdf5.h
 \brief LATfield2_save_hdf5.h contains the definition of the function used for hdf5 i/o.
 \author David Daverio
 */
extern "C"{

#include <math.h>
#include <hdf5.h>
#include <stdlib.h>
#include <string.h>


   int save_hdf5_externC(char *data,long file_offset[2],int *size,int * sizeLocal,int halo, int lat_dim,int comp,hid_t array_type,int array_size,string  filename_str, string dataset_name_str)
   {

	   hid_t file_id, plist_id,filespace,memspace,dset_id,dtype_id,dtbase_id,root_id;
	   hsize_t * components;

	   char * filename;
	   filename = (char*)malloc((filename_str.size()+1)*sizeof(char));
     for(int i = 0;i<filename_str.size();i++)filename[i]=filename_str[i];
     filename[filename_str.size()] = '\0';

	   char  dataset_name[128];
	   for(int i = 0;i<filename_str.size();i++)dataset_name[i]=dataset_name_str[i];
	   dataset_name[dataset_name_str.size()] = '\0';

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

	   ///////////////////////////////
	   // creat datatype
	   ///////////////////////////////
	   if(comp == 1 && array_size ==1)
	   {
		   dtype_id = H5Tcopy(array_type);
		   status = H5Tset_order(dtype_id, DATA_ORDER);
		   components = new hsize_t[1]; //to be sure is allocated when freed
	   }
	   if(comp == 1 && array_size !=1)
	   {
		   components = new hsize_t[1];
		   components[0] = array_size;
		   dtbase_id = H5Tcopy(array_type);
		   status = H5Tset_order(dtbase_id, DATA_ORDER);
		   dtype_id = H5Tarray_create(dtbase_id,1,components);
	   }
	   if(comp != 1 && array_size ==1)
	   {
		   components = new hsize_t[1];
		   components[0] = comp;
		   dtbase_id = H5Tcopy(array_type);
		   status = H5Tset_order(dtbase_id, DATA_ORDER);
		   dtype_id = H5Tarray_create(dtbase_id,1,components);
	   }
	   if(comp != 1 && array_size !=1)
	   {
		   components = new hsize_t[2];
		   components[0] = array_size;
		   components[1] = comp;
		   dtbase_id = H5Tcopy(array_type);
		   status = H5Tset_order(dtbase_id, DATA_ORDER);
		   dtype_id = H5Tarray_create(dtbase_id,2,components);
	   }
	   ///////////////////////////////
	   ///////////////////////////////

#ifdef H5_HAVE_PARALLEL //Parallel version, H5_HAVE_PARALLEL definition is needed by hdf5 to run in parallel too !


	   MPI_Comm comm  = parallel.lat_world_comm();
	   MPI_Info info  = MPI_INFO_NULL;


	   //creat the file
	   plist_id = H5Pcreate(H5P_FILE_ACCESS);
	   H5Pset_fapl_mpio(plist_id, comm, info);

	   file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	   H5Pclose(plist_id);


	   filespace = H5Screate_simple(lat_dim,sizeGlobal,NULL);
	   memspace = H5Screate_simple(lat_dim,localSize,NULL);

	   plist_id = H5Pcreate(H5P_DATASET_CREATE);
	   //H5Pset_chunk(plist_id, lat_dim,count);
	   dset_id = H5Dcreate1(file_id, dataset_name, dtype_id, filespace,plist_id);

	   H5Pclose(plist_id);
	   H5Sclose(filespace);

	   //save data

	   filespace = H5Dget_space(dset_id);
	   status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, count, NULL);
	   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

	   plist_id = H5Pcreate(H5P_DATASET_XFER);
	   H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	   status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data);


	   H5Dclose(dset_id);
	   H5Sclose(filespace);
	   H5Sclose(memspace);
	   H5Pclose(plist_id);
	   H5Fclose(file_id);
	   free(filename);

	   return 1;

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

		   filespace = H5Screate_simple(lat_dim,sizeGlobal,NULL);

		   plist_id = H5Pcreate(H5P_DATASET_CREATE);
		   H5Pset_chunk(plist_id, lat_dim, sizeGlobal);
		   dset_id = H5Dcreate1(file_id, dataset_name, dtype_id, filespace,plist_id);

		   H5Pclose(plist_id);
		   H5Sclose(filespace);
		   H5Dclose(dset_id);
		   H5Fclose(file_id);


	   }

	   MPI_Barrier(parallel.lat_world_comm());

	   for(p=0;p < mpi_size;p++)
	   {


		   if(mpi_rank==p)
		   {
         //cout<<"rank: "<<mpi_rank<<" , writting data"<<endl;
			   plist_id = H5Pcreate(H5P_FILE_ACCESS);

			   file_id = H5Fopen(filename,H5F_ACC_RDWR,plist_id);
			   H5Pclose(plist_id);
			   root_id = H5Gopen(file_id,"/",H5P_DEFAULT);
			   dset_id = H5Dopen(root_id, dataset_name, H5P_DEFAULT);
			   filespace = H5Dget_space(dset_id);
			   dtype_id = H5Dget_type(dset_id);
			   plist_id = H5Pcreate(H5P_DATASET_XFER);

			   memspace = H5Screate_simple(lat_dim,localSize,NULL);

			   int status;

			   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL,count, NULL);

			   status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL,count, NULL);

			   status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data);

			   H5Pclose(plist_id);
			   H5Sclose(filespace);
         H5Sclose(memspace);
			   H5Dclose(dset_id);
         H5Gclose(root_id);
			   H5Fclose(file_id);

		   }
		   MPI_Barrier(parallel.lat_world_comm());

	   }

	   free(filename);
	   return 1;
#endif

	   return -1;


 }


	int load_hdf5_externC(char *data,long file_offset[2],int *size,int * sizeLocal,int halo, int lat_dim,string  filename_str, string dataset_name_str)
	{



	    hid_t file_id, plist_id, plistxfer_id,filespace,memspace,dset_id,dtype_id,dtbase_id,group_id,root_id;


		char * filename;
		filename = (char*)malloc((filename_str.size()+1)*sizeof(char));
		for(int i = 0;i<filename_str.size();i++)filename[i]=filename_str[i];
		filename[filename_str.size()] = '\0';

		char  dataset_name[128];
		for(int i = 0;i<filename_str.size();i++)dataset_name[i]=dataset_name_str[i];
		dataset_name[dataset_name_str.size()] = '\0';

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

		file_id = H5Fopen(filename,H5F_ACC_RDONLY,plist_id);
		//H5Pclose(plist_id);

		root_id = H5Gopen(file_id,"/",H5P_DEFAULT);
		dset_id = H5Dopen(root_id, dataset_name, H5P_DEFAULT);
		filespace = H5Dget_space(dset_id);
		dtype_id = H5Dget_type(dset_id);
		//plist_id = H5Pcreate(H5P_DATASET_XFER);

		memspace = H5Screate_simple(lat_dim,localSize,NULL);



        //// verifier si l "offset" ne doit pas tenire compte du ....
		status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, count, NULL);
		status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

		plistxfer_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plistxfer_id, H5FD_MPIO_COLLECTIVE);
		status = H5Dread(dset_id, dtype_id, memspace, filespace, plistxfer_id, data);

		H5Dclose(dset_id);
		H5Gclose(root_id);
		H5Sclose(filespace);
		H5Sclose(memspace);
		H5Pclose(plistxfer_id);
		H5Fclose(file_id);
		H5Pclose(plist_id);


		return 1;



#else // serial version, without H5_HAVE_PARALLEL definition hdf5 will run in serial !
		int mpi_size,mpi_rank,p;
		//MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
		//MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		mpi_size = parallel.size();
		mpi_rank = parallel.rank();

		for(p=0;p < mpi_size;p++)
		{
			if(mpi_rank==p)
			{

				plist_id = H5Pcreate(H5P_FILE_ACCESS);

				file_id = H5Fopen(filename,H5F_ACC_RDONLY,plist_id);
        if(file_id<0)
        {
          cout<<"process "<<p << ", cant open file: "<<filename_str<<endl;
          parallel.abortForce();
        }
				H5Pclose(plist_id);

				root_id = H5Gopen(file_id,"/",H5P_DEFAULT);
				dset_id = H5Dopen(root_id, dataset_name, H5P_DEFAULT);
        if(dset_id<0)
        {
          cout<<"process "<<p << ", cant open dataset : "
              <<dataset_name_str<<" in file "<<filename_str<<endl;
          parallel.abortForce();
        }
				filespace = H5Dget_space(dset_id);
				dtype_id = H5Dget_type(dset_id);
				plist_id = H5Pcreate(H5P_DATASET_XFER);

				memspace = H5Screate_simple(lat_dim,localSize,NULL);

				int status;

				status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL,count, NULL);
				status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL,count, NULL);
				status = H5Dread(dset_id, dtype_id, memspace, filespace, plist_id, data);


        H5Pclose(plist_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Dclose(dset_id);
        H5Gclose(root_id);
        H5Fclose(file_id);

			}
			MPI_Barrier(MPI_COMM_WORLD);

		}

		return  1;

#endif
        return -1;
	}

}

template<class fieldType>
int save_hdf5(fieldType *data,hid_t type_id,int array_size,long file_offset[2],int *size,int * sizeLocal,int halo, int lat_dim,int comp,string  filename_str, string dataset_name_str)
{

	return save_hdf5_externC((char*)data, file_offset, size, sizeLocal, halo, lat_dim, comp, type_id, array_size, filename_str, dataset_name_str);
}
template<class fieldType>
int load_hdf5(fieldType *data,long file_offset[2],int *size,int * sizeLocal,int halo, int lat_dim,int comp,string  filename_str, string dataset_name_str)
{
    return load_hdf5_externC( (char*) data, file_offset, size, sizeLocal, halo, lat_dim,  filename_str, dataset_name_str);
}
