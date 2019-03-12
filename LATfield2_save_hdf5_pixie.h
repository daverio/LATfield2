#ifdef BIG_ENDIAN_ORDER
#define DATA_ORDER H5T_ORDER_BE
#else
#define  DATA_ORDER H5T_ORDER_LE
#endif

/*! \file LATfield2_save_hdf5_pixie.h
 \brief LATfield2_save_hdf5_pixie.h contains the definition of the function used for hdf5 i/o compatible to pixie reader.
 \author David Daverio
 */


extern "C"{


#include <math.h>
#include <hdf5.h>
#include <stdlib.h>
#include <string.h>
#include "int2string.hpp"
#include <sys/stat.h>

inline bool file_exists(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}


int save_hdf5_externC(char *data,long file_offset[2],int *size,int * sizeLocal,int halo, int lat_dim,int comp,hid_t array_type,int array_size,string  filename_str, string arg_dataset_name_str)
{


		bool fexist;


	   hid_t file_id, plist_id,filespace,memspace,dset_id,dtype_id,dtbase_id,plist_id_dataset,plist_id_write;
	   hsize_t  components;
		 string compName;

		 string dataset_name_str;


	   char * filename;
           filename = (char*)malloc((filename_str.size()+1)*sizeof(char));
           for(int i = 0;i<filename_str.size();i++)filename[i]=filename_str[i];
           filename[filename_str.size()] = '\0';

	   herr_t status;

	   hsize_t * sizeGlobal;
	   sizeGlobal = new hsize_t[lat_dim+1];
	   hsize_t * localSize;
	   localSize = new hsize_t[lat_dim+1];


       hsize_t * offset;
       offset = new hsize_t[lat_dim+1];
	   hsize_t * offsetf;
	   offsetf = new hsize_t[lat_dim+1];
	   hsize_t * count;
	   count = new hsize_t[lat_dim+1];


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

       sizeGlobal[lat_dim]=comp;
       localSize[lat_dim]=comp;
       count[lat_dim]=1;

       ///////////////////////////////
	   // creat datatype in big endian
	   ///////////////////////////////

       if(array_size ==1)
       {
       dtype_id = H5Tcopy(array_type);
		   status = H5Tset_order(dtype_id, DATA_ORDER);
       }
       else if(array_size > 1)
       {
		   components = array_size;
		   dtbase_id = H5Tcopy(array_type);
		   status = H5Tset_order(dtbase_id, DATA_ORDER);
		   dtype_id = H5Tarray_create(dtbase_id,1,&components);
       }
       ///////////////////////////////
	   ///////////////////////////////


#ifdef H5_HAVE_PARALLEL //Parallel version, H5_HAVE_PARALLEL definition is needed by hdf5 to run in parallel too !


		MPI_Comm comm  = parallel.lat_world_comm();
		MPI_Info info  = MPI_INFO_NULL;

		//verify if the file exists:
		if(parallel.rank()==0)fexist = file_exists(filename_str);
		parallel.broadcast(fexist,0);

		plist_id = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_id, comm, info);

		if(fexist) file_id = H5Fopen(filename_str.c_str(),H5F_ACC_RDWR,plist_id);
		else file_id = H5Fcreate(filename_str.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
		H5Pclose(plist_id);

		if(comp==1)
		{
			compName = "/"+arg_dataset_name_str;

			if(H5Lexists(file_id, compName.c_str(),H5P_DEFAULT)<=0)
			{
				filespace = H5Screate_simple(lat_dim,sizeGlobal,NULL);
				dset_id = H5Dcreate(file_id, compName.c_str(), dtype_id, filespace,
													  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				H5Sclose(filespace);
			}
			else
			{
				dset_id = H5Dopen(file_id, compName.c_str(), H5P_DEFAULT);
			}

			filespace = H5Dget_space(dset_id);
			memspace = H5Screate_simple(lat_dim,localSize,NULL);

			status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, count, NULL);
			status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

			plist_id = H5Pcreate(H5P_DATASET_XFER);
			H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
			status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data);

			H5Dclose(dset_id);
	    H5Sclose(filespace);
	    H5Sclose(memspace);
	    H5Pclose(plist_id);

		}
		else if(comp > 1)
		{
			for(int c = 0;c<comp;c++)
			{
				compName = "/"+arg_dataset_name_str+"_"+int2string(c,999);
				offset[lat_dim]=c;

				if(H5Lexists(file_id, compName.c_str(),H5P_DEFAULT)<=0)
				{
					filespace = H5Screate_simple(lat_dim,sizeGlobal,NULL);
					dset_id = H5Dcreate(file_id, compName.c_str(), dtype_id, filespace,
														  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
					H5Sclose(filespace);
				}
				else
				{
					dset_id = H5Dopen(file_id, compName.c_str(), H5P_DEFAULT);
				}

				filespace = H5Dget_space(dset_id);
				memspace = H5Screate_simple(lat_dim+1,localSize,NULL);

				status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, count, NULL);
				status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

				plist_id = H5Pcreate(H5P_DATASET_XFER);
				H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
				status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data);

				H5Dclose(dset_id);
		    H5Sclose(filespace);
		    H5Sclose(memspace);
		    H5Pclose(plist_id);
			}
		}

    H5Fclose(file_id);

	  return 1;

#else // serial version, without H5_HAVE_PARALLEL definition hdf5 will run in serial !

    int mpi_size,mpi_rank,p;
    MPI_Comm_size(parallel.lat_world_comm(), &mpi_size);
    MPI_Comm_rank(parallel.lat_world_comm(), &mpi_rank);

    //create the file

    if(mpi_rank==0)
    {
				fexist = file_exists(filename_str);

				plist_id = H5Pcreate(H5P_FILE_ACCESS);
				if(fexist) file_id = H5Fopen(filename,H5F_ACC_RDWR,plist_id);
				else file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
				H5Pclose(plist_id);

        filespace = H5Screate_simple(lat_dim,sizeGlobal,NULL);

        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, lat_dim, sizeGlobal);

        if(comp==1)
        {
						compName = "/"+arg_dataset_name_str;
						if(H5Lexists(file_id, compName.c_str(),H5P_DEFAULT)<=0)
						{
							dset_id = H5Dcreate1(file_id, compName.c_str(), dtype_id, filespace,plist_id);
							H5Dclose(dset_id);
						}
        }
        else if(comp > 1)
        {
            for(int c = 0;c<comp;c++)
            {
                compName = "/"+arg_dataset_name_str+"_"+int2string(c,999);
								if(H5Lexists(file_id, compName.c_str(),H5P_DEFAULT)<=0)
								{
									dset_id = H5Dcreate1(file_id, compName.c_str(), dtype_id, filespace,plist_id);
                	H5Dclose(dset_id);
								}
            }
        }
        H5Pclose(plist_id);
        H5Sclose(filespace);
        H5Fclose(file_id);
    }
    MPI_Barrier(parallel.lat_world_comm());
    for(p=0;p < mpi_size;p++)
    {
        MPI_Barrier(parallel.lat_world_comm());
        if(mpi_rank==p)
        {
            plist_id = H5Pcreate(H5P_FILE_ACCESS);

            file_id = H5Fopen(filename,H5F_ACC_RDWR,plist_id);
            H5Pclose(plist_id);

            if(comp==1)
            {
								compName = "/"+arg_dataset_name_str;
                dset_id = H5Dopen(file_id, compName.c_str(), H5P_DEFAULT);
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
                H5Dclose(dset_id);
            }
            else if(comp>1)
            {
                for(int c=0;c<comp;c++)
                {
										compName = "/"+arg_dataset_name_str+"_"+int2string(c,999);
                    dset_id = H5Dopen(file_id, compName.c_str(), H5P_DEFAULT);
                    filespace = H5Dget_space(dset_id);
                    dtype_id = H5Dget_type(dset_id);
                    plist_id = H5Pcreate(H5P_DATASET_XFER);

                    memspace = H5Screate_simple(lat_dim+1,localSize,NULL);

                    offset[lat_dim]=c;

                    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL,count, NULL);
                    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL,count, NULL);

                    status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data);

                    H5Pclose(plist_id);
                    H5Sclose(filespace);
                    H5Dclose(dset_id);
                }
            }

            H5Fclose(file_id);

        }
        MPI_Barrier(parallel.lat_world_comm());



    }
		MPI_Barrier(parallel.lat_world_comm());

           delete[] filename;


           delete[] sizeGlobal;
           delete[] localSize;
           delete[] offset;
           delete[] offsetf;
           delete[] count;

	   return 1;

#endif

   }



	int load_hdf5_externC(char *data,long file_offset[2],int *size,int * sizeLocal,int comp,int halo, int lat_dim,string  filename_str, string arg_dataset_name_str)
	{

        hid_t file_id, plist_id,filespace,memspace,dset_id,dtype_id,dtbase_id;


	char * filename;
	filename = (char*)malloc((filename_str.size()+1)*sizeof(char));
	for(int i = 0;i<filename_str.size();i++)filename[i]=filename_str[i];
	filename[filename_str.size()] = '\0';

	string compName;


        hsize_t * sizeGlobal;
        sizeGlobal = new hsize_t[lat_dim+1];
        hsize_t * localSize;
        localSize = new hsize_t[lat_dim+1];

        herr_t status;

        hsize_t * offset;
        offset = new hsize_t[lat_dim+1];
        hsize_t * offsetf;
        offsetf = new hsize_t[lat_dim+1];
        hsize_t * count;
        count = new hsize_t[lat_dim+1];

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

        sizeGlobal[lat_dim]=comp;
        localSize[lat_dim]=comp;
        count[lat_dim]=1;



#ifdef H5_HAVE_PARALLEL	//Parallel version, H5_HAVE_PARALLEL definition is needed by hdf5 to run in parallel too

      MPI_Comm comm  = parallel.lat_world_comm();
		  MPI_Info info  = MPI_INFO_NULL;

		plist_id = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_id, comm, info);

		file_id = H5Fopen(filename,H5F_ACC_RDWR,plist_id);
		H5Pclose(plist_id);


        if(comp==1)
        {
						compName = "/"+arg_dataset_name_str;
            dset_id = H5Dopen(file_id, compName.c_str(), H5P_DEFAULT);
            filespace = H5Dget_space(dset_id);
            dtype_id = H5Dget_type(dset_id);
            plist_id = H5Pcreate(H5P_DATASET_XFER);

            memspace = H5Screate_simple(lat_dim,localSize,NULL);

            status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, count, NULL);
            status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            status = H5Dread(dset_id, dtype_id, memspace, filespace, plist_id, data);

            H5Dclose(dset_id);
            H5Sclose(filespace);
            H5Sclose(memspace);
            H5Pclose(plist_id);
            H5Fclose(file_id);


        }
        else if(comp > 1)
        {
            //plist_id = H5Pcreate(H5P_DATASET_CREATE);
            for(int c = 0;c<comp;c++)
            {
								compName = "/"+arg_dataset_name_str+"_"+int2string(c,999);
                dset_id = H5Dopen(file_id, compName.c_str(), H5P_DEFAULT);
                filespace = H5Dget_space(dset_id);
                dtype_id = H5Dget_type(dset_id);
                plist_id = H5Pcreate(H5P_DATASET_XFER);

                memspace = H5Screate_simple(lat_dim+1,localSize,NULL);

                offset[lat_dim]=c;

                status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, count, NULL);
                status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
                status = H5Dread(dset_id, dtype_id, memspace, filespace, plist_id, data);

                H5Dclose(dset_id);
                H5Sclose(filespace);
                H5Sclose(memspace);
                H5Pclose(plist_id);



            }
            H5Fclose(file_id);

        }

		free(filename);
		return 1;



#else // serial version, without H5_HAVE_PARALLEL definition hdf5 will run in serial !


        int mpi_size,mpi_rank,p;
		MPI_Comm_size(parallel.lat_world_comm(), &mpi_size);
		MPI_Comm_rank(parallel.lat_world_comm(), &mpi_rank);


		for(p=0;p < mpi_size;p++)
		{
			if(mpi_rank==p)
			{
				plist_id = H5Pcreate(H5P_FILE_ACCESS);

				file_id = H5Fopen(filename,H5F_ACC_RDWR,plist_id);
				H5Pclose(plist_id);

                if(comp==1)
                {
										compName = "/"+arg_dataset_name_str;
                    dset_id = H5Dopen(file_id, compName.c_str(), H5P_DEFAULT);
                    filespace = H5Dget_space(dset_id);
                    dtype_id = H5Dget_type(dset_id);
                    plist_id = H5Pcreate(H5P_DATASET_XFER);

                    memspace = H5Screate_simple(lat_dim,localSize,NULL);


                    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL,count, NULL);

                    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL,count, NULL);

                    status = H5Dread(dset_id, dtype_id, memspace, filespace, plist_id, data);

                    H5Pclose(plist_id);
                    H5Sclose(filespace);
                    H5Sclose(memspace);
                    H5Dclose(dset_id);
                    H5Fclose(file_id);

                }
                else if(comp > 1)
                {
                    //plist_id = H5Pcreate(H5P_DATASET_CREATE);
                    for(int c = 0;c<comp;c++)
                    {
												compName = "/"+arg_dataset_name_str+"_"+int2string(c,999);
                        dset_id = H5Dopen(file_id, compName.c_str(), H5P_DEFAULT);
                        filespace = H5Dget_space(dset_id);


                        dtype_id = H5Dget_type(dset_id);
                        plist_id = H5Pcreate(H5P_DATASET_XFER);



                        memspace = H5Screate_simple(lat_dim+1,localSize,NULL);

                        offset[lat_dim]=c;


                        status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL,count, NULL);

                        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL,count, NULL);

                        status = H5Dread(dset_id, dtype_id, memspace, filespace, plist_id, data);

                        H5Pclose(plist_id);
                        H5Sclose(filespace);
                        H5Sclose(memspace);
                        H5Dclose(dset_id);


                    }
                    H5Fclose(file_id);

                }
                MPI_Barrier(parallel.lat_world_comm());

            }

        }



		return 1;



#endif

    delete[] filename;


    delete[] sizeGlobal;
    delete[] localSize;
    delete[] offset;
    delete[] offsetf;
    delete[] count;
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
    return load_hdf5_externC( (char*) data, file_offset, size, sizeLocal, comp, halo, lat_dim,  filename_str, dataset_name_str);
}
