#ifndef LATFIELD_PARTICLESIO_H
#define LATFIELD_PARTICLESIO_H


#ifndef RealC
#ifdef SINGLE
#define RealC float
#define MPI_RealC MPI_FLOAT
#else
#define RealC double
#define MPI_RealC MPI_DOUBLE
#endif
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct fileDsc
{
  int fileNumber;
  int numProcPerFile;
  int world_size;
  int grid_size[2];
  int latSize[3];
  RealC boxSize[3];
  RealC fileBoxSize;
  RealC fileBoxOffset;
  long numParts;
  RealC localBoxOffset[3];
  RealC localBoxSize[3];
};
#endif


template<typename part_info,typename parts_datatype>
void get_partInfo(string filename, part_info &partInfo, parts_datatype partdatatype)
{
#ifdef H5_HAVE_PARALLEL
hid_t plist_id,file_id,attr_id,root_id,dataset_id;


MPI_Info info  = MPI_INFO_NULL;


plist_id = H5Pcreate(H5P_FILE_ACCESS);
H5Pset_fapl_mpio(plist_id,parallel.lat_world_comm(),info);
file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);

if(file_id<0)cout<< "get_partInfo: cant open file: "<<filename<<endl;


if(parallel.rank()==0)
{
  dataset_id = H5Dopen(file_id, "/part_info", H5P_DEFAULT);
  H5Dread(dataset_id, partdatatype.part_info_memType, H5S_ALL, H5S_ALL, H5P_DEFAULT,&partInfo);
  H5Dclose(dataset_id);
}
H5Fclose(file_id);
H5Pclose(plist_id);


#else
    if(parallel.rank()==0)
    {
	hid_t plist_id,file_id,dataset_id;

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);

  if(file_id<0)cout<< "get_partInfo: cant open file: "<<filename<<endl;

	H5Pclose(plist_id);

	dataset_id = H5Dopen(file_id, "/part_info", H5P_DEFAULT);

	H5Dread(dataset_id, partdatatype.part_info_memType, H5S_ALL, H5S_ALL, H5P_DEFAULT,&partInfo);
	H5Dclose(dataset_id);
	H5Fclose(file_id);
    }

  #endif
  parallel.broadcast<part_info>(partInfo,0);
}

void get_fileDsc_global(string filename,fileDsc &fd)
{

  #ifdef H5_HAVE_PARALLEL

  hid_t plist_id,file_id,attr_id,root_id;


  MPI_Info info  = MPI_INFO_NULL;


  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,parallel.lat_world_comm(),info);
  file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);

  if(file_id<0)cout<< "get_fileDsc_global: cant open file: "<<filename<<endl;

  if(parallel.rank()==0)
  {
  	attr_id = H5Aopen_name(file_id, "fileNumber");
  	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.fileNumber));
  	H5Aclose(attr_id);

          attr_id = H5Aopen_name(file_id, "numProcPerFile");
  	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.numProcPerFile));
  	H5Aclose(attr_id);

  	attr_id = H5Aopen_name(file_id, "world_size");
  	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.world_size));
  	H5Aclose(attr_id);

  	attr_id = H5Aopen_name(file_id, "grid_size");
  	H5Aread(attr_id, H5T_NATIVE_INT, fd.grid_size);
  	H5Aclose(attr_id);

  	attr_id = H5Aopen_name(file_id, "boxSize");
  	H5Aread(attr_id, REAL_TYPE, fd.boxSize);
  	H5Aclose(attr_id);



  	attr_id = H5Aopen_name(file_id, "fileBoxSize");
  	H5Aread(attr_id, REAL_TYPE, &(fd.fileBoxSize));
  	H5Aclose(attr_id);

  	attr_id = H5Aopen_name(file_id, "fileBoxOffset");
  	H5Aread(attr_id, REAL_TYPE, &(fd.fileBoxOffset));
  	H5Aclose(attr_id);
  }
  H5Pclose(plist_id);
	H5Fclose(file_id);


  #else

    if(parallel.rank()==0)
    {
        hid_t plist_id,file_id,attr_id,root_id;

        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);




	attr_id = H5Aopen_name(file_id, "fileNumber");
	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.fileNumber));
	H5Aclose(attr_id);

        attr_id = H5Aopen_name(file_id, "numProcPerFile");
	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.numProcPerFile));
	H5Aclose(attr_id);

	attr_id = H5Aopen_name(file_id, "world_size");
	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.world_size));
	H5Aclose(attr_id);

	attr_id = H5Aopen_name(file_id, "grid_size");
	H5Aread(attr_id, H5T_NATIVE_INT, fd.grid_size);
	H5Aclose(attr_id);

	attr_id = H5Aopen_name(file_id, "boxSize");
	H5Aread(attr_id, REAL_TYPE, fd.boxSize);
	H5Aclose(attr_id);



	attr_id = H5Aopen_name(file_id, "fileBoxSize");
	H5Aread(attr_id, REAL_TYPE, &(fd.fileBoxSize));
	H5Aclose(attr_id);

	attr_id = H5Aopen_name(file_id, "fileBoxOffset");
	H5Aread(attr_id, REAL_TYPE, &(fd.fileBoxOffset));
	H5Aclose(attr_id);
	//cout<<"bosize loaded is "<< fd.boxSize[0] <<endl;
	H5Pclose(plist_id);
	H5Fclose(file_id);
    }

#endif
  parallel.broadcast<fileDsc>(fd,0);
  //cout<<"bosize loaded is "<< fd.boxSize[0] <<endl;


}
void get_fileDsc_local(string filename,long * numParts, RealC * localBoxOffset, RealC * localBoxSize, int numProcPerfile)
{
  #ifdef H5_HAVE_PARALLEL

  hid_t plist_id,file_id,attr_id,root_id,dataset_id;


  MPI_Info info  = MPI_INFO_NULL;


  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,parallel.lat_world_comm(),info);
  file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);
  if(file_id<0)cout<< "get_fileDsc_local: cant open file: "<<filename<<endl;

  if(parallel.rank()==0)
  {
    dataset_id = H5Dopen(file_id, "/numParts", H5P_DEFAULT);
  	H5Dread(dataset_id,H5T_NATIVE_LONG , H5S_ALL, H5S_ALL, H5P_DEFAULT,numParts);
  	H5Dclose(dataset_id);
  	dataset_id = H5Dopen(file_id, "/localBoxOffset", H5P_DEFAULT);
  	H5Dread(dataset_id,REAL_TYPE , H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxOffset);
  	H5Dclose(dataset_id);
  	dataset_id = H5Dopen(file_id, "/localBoxSize", H5P_DEFAULT);
  	H5Dread(dataset_id,REAL_TYPE , H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxSize);
  	H5Dclose(dataset_id);

  }
  H5Pclose(plist_id);
  H5Fclose(file_id);
  #else
    if(parallel.rank()==0)
    {
	hid_t plist_id,file_id,dataset_id;

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);

  if(file_id<0)cout<< "get_fileDsc_local: cant open file: "<<filename<<endl;

	dataset_id = H5Dopen(file_id, "/numParts", H5P_DEFAULT);
	H5Dread(dataset_id,H5T_NATIVE_LONG , H5S_ALL, H5S_ALL, H5P_DEFAULT,numParts);
	H5Dclose(dataset_id);
	dataset_id = H5Dopen(file_id, "/localBoxOffset", H5P_DEFAULT);
	H5Dread(dataset_id,REAL_TYPE , H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxOffset);
	H5Dclose(dataset_id);
	dataset_id = H5Dopen(file_id, "/localBoxSize", H5P_DEFAULT);
	H5Dread(dataset_id,REAL_TYPE , H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxSize);
	H5Dclose(dataset_id);
	H5Pclose(plist_id);
	H5Fclose(file_id);
    }

    #endif
    parallel.broadcast(numParts,numProcPerfile,0);
    parallel.broadcast(localBoxOffset,3*numProcPerfile,0);
    parallel.broadcast(localBoxSize,3*numProcPerfile,0);

}
template<typename parts,typename parts_datatype>
void get_part_sublist(string filename, long offset, long nparts, parts * parList, parts_datatype partdatatype )
{
  #ifdef H5_HAVE_PARALLEL
  hid_t plist_id,file_id,dataset_id,filespace_id,memspace_id;
  hsize_t offsetf,localNumParts;
  offsetf=offset;
  localNumParts=nparts;

  MPI_Info info  = MPI_INFO_NULL;


  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,parallel.lat_world_comm(),info);
  file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);
  if(file_id<0)cout<< "get_fileDsc_local: cant open file: "<<filename<<endl;

  dataset_id = H5Dopen(file_id, "/data", H5P_DEFAULT);
  filespace_id = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &offsetf, NULL,&localNumParts, NULL);

  memspace_id  = H5Screate_simple(1,&localNumParts,NULL);

  H5Dread(dataset_id, partdatatype.part_memType, memspace_id, filespace_id, H5P_DEFAULT,parList);

  H5Sclose(filespace_id);
  H5Dclose(dataset_id);
  H5Fclose(file_id);
  H5Pclose(plist_id);


  #else
    hid_t plist_id,file_id,dataset_id,filespace_id,memspace_id;
    hsize_t offsetf,localNumParts;
    offsetf=offset;
    localNumParts=nparts;

    plist_id = H5Pcreate(H5P_FILE_ACCESS);
	  file_id = H5Fopen(filename.c_str(),H5F_ACC_RDWR,plist_id);

    if(file_id<0)cout<< "get_part_sublist: cant open file: "<<filename<<endl;

	  H5Pclose(plist_id);

    dataset_id = H5Dopen(file_id, "/data", H5P_DEFAULT);
    filespace_id = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &offsetf, NULL,&localNumParts, NULL);

    memspace_id  = H5Screate_simple(1,&localNumParts,NULL);

    H5Dread(dataset_id, partdatatype.part_memType, memspace_id, filespace_id, H5P_DEFAULT,parList);

    H5Sclose(filespace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    #endif
}


template<typename parts,typename part_info,typename parts_datatype>
int save_hdf5_particles(string filename,
			parts * partlist,
			part_info partInfo,
			parts_datatype partdatatype,
			fileDsc fd,
			MPI_Comm comm){


  int mpi_size,mpi_rank;
  hsize_t one=1;
  hsize_t two=2;

  long *numParts;
    RealC * localBoxOffset;
    RealC * localBoxSize;

  hsize_t globalNumParts,localNumParts,offset,offsetf,numPartsSize[2];
  hid_t file_id,plist_id,filespace_id,memspace_id,dataset_id,dataspace_id,attribute_id,plist_id_file;

  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Info info  = MPI_INFO_NULL;

  if(mpi_size != fd.numProcPerFile)
    {
      COUT << "save particles internal error: error 1"<< endl;
      COUT << "please report this bug to developer@LATfield.org"<<endl;
      COUT << "Exiting executable..."<<endl;
      exit(1111);
    }

    numPartsSize[0]=mpi_size;
    numPartsSize[1]=3;

  numParts = new long[mpi_size];
    localBoxOffset = new RealC[mpi_size*3];
    localBoxSize  = new RealC[mpi_size*3];


  numParts[mpi_rank] = fd.numParts;
    for(int i=0;i<3;i++)
    {
        localBoxOffset[3*mpi_rank+i]=fd.localBoxOffset[i];
        localBoxSize[3*mpi_rank+i]=fd.localBoxSize[i];
    }


  for(int i=0;i<mpi_size;i++)
    {
        MPI_Bcast(&numParts[i],1,MPI_LONG,i,comm);
        MPI_Bcast(&localBoxOffset[i*3],3,MPI_RealC,i,comm);
        MPI_Bcast(&localBoxSize[i*3],3,MPI_RealC,i,comm);
        /*
        if(mpi_rank==i){
            cout<< localBoxSize[i*3] << " "<< localBoxSize[i*3+1] << " "<< localBoxSize[i*3+2] << endl;
            cout<< localBoxOffset[i*3] << " "<< localBoxOffset[i*3+1] << " "<< localBoxOffset[i*3+2] << endl;
        }*/
    }


  localNumParts = numParts[mpi_rank];

  globalNumParts=0;
  for(int i=0;i<mpi_size;i++)globalNumParts += numParts[i];

  offset=0;
  offsetf=0;
  for(int i=0;i<mpi_rank;i++)offsetf+=numParts[i];




#ifdef H5_HAVE_PARALLEL


  plist_id_file = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id_file,comm,info);

  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id_file);









  filespace_id = H5Screate_simple(1,&globalNumParts,NULL);
  memspace_id  = H5Screate_simple(1,&localNumParts,NULL);

  dataset_id = H5Dcreate(file_id, "data", partdatatype.part_fileType, filespace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, &offset, NULL,&localNumParts, NULL);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &offsetf, NULL,&localNumParts, NULL);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dataset_id,partdatatype.part_memType, memspace_id, filespace_id, plist_id,partlist);

  H5Pclose(plist_id);
  H5Sclose(memspace_id);
  H5Dclose(dataset_id);
  H5Sclose(filespace_id);


  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

  dataspace_id = H5Screate_simple(1,numPartsSize,NULL);
  dataset_id = H5Dcreate(file_id, "numParts", H5T_NATIVE_LONG, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(mpi_rank==0)H5Dwrite(dataset_id, H5T_NATIVE_LONG, dataspace_id, dataspace_id, plist_id,numParts);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate_simple(2,numPartsSize,NULL);
  dataset_id = H5Dcreate(file_id, "localBoxOffset", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(mpi_rank==0)H5Dwrite(dataset_id, REAL_TYPE, dataspace_id, dataspace_id, plist_id,localBoxOffset);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate_simple(2,numPartsSize,NULL);
  dataset_id = H5Dcreate(file_id, "localBoxSize", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(mpi_rank==0)H5Dwrite(dataset_id, REAL_TYPE, dataspace_id, dataspace_id, plist_id,localBoxSize);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);

  filespace_id = H5Screate_simple(1,&one,NULL);
  dataset_id = H5Dcreate(file_id, "part_info", partdatatype.part_info_fileType , filespace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if(mpi_rank==0)H5Dwrite(dataset_id, partdatatype.part_info_memType, filespace_id, filespace_id, plist_id,&partInfo);
  H5Dclose(dataset_id);
  H5Sclose(filespace_id);

  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate (file_id, "fileNumber", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id,H5T_NATIVE_INT , &(fd.fileNumber));
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);


        dataspace_id = H5Screate(H5S_SCALAR);
        attribute_id = H5Acreate (file_id, "numProcPerFile", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id,H5T_NATIVE_INT , &(fd.numProcPerFile));
        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);

        dataspace_id = H5Screate(H5S_SCALAR);
        attribute_id = H5Acreate (file_id, "world_size", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id,H5T_NATIVE_INT , &(fd.world_size));
        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(1,&two,NULL);
        attribute_id = H5Acreate (file_id, "grid_size", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id,H5T_NATIVE_INT , fd.grid_size);
        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(1,&numPartsSize[1],NULL);
        attribute_id = H5Acreate (file_id, "boxSize", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id,REAL_TYPE , fd.boxSize);
        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);

        dataspace_id = H5Screate(H5S_SCALAR);
        attribute_id = H5Acreate (file_id, "fileBoxSize", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id,REAL_TYPE , &(fd.fileBoxSize));
        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);

        dataspace_id = H5Screate(H5S_SCALAR);
        attribute_id = H5Acreate (file_id, "fileBoxOffset", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id,REAL_TYPE , &(fd.fileBoxOffset));
        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(1,&numPartsSize[1],NULL);
        attribute_id = H5Acreate (file_id, "latSize", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id,H5T_NATIVE_INT , fd.latSize);
        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);



  H5Pclose(plist_id_file);
  H5Fclose(file_id);



#else

  if(mpi_rank==0){

    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    filespace_id = H5Screate_simple(1,&globalNumParts,NULL);
    dataset_id = H5Dcreate(file_id, "data", partdatatype.part_fileType, filespace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dclose(dataset_id);
    H5Sclose(filespace_id);

    filespace_id = H5Screate_simple(1,&one,NULL);
    dataset_id = H5Dcreate(file_id, "part_info", partdatatype.part_info_fileType , filespace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(dataset_id, partdatatype.part_info_memType, H5S_ALL, H5S_ALL, H5P_DEFAULT,&partInfo);
    H5Dclose(dataset_id);
    H5Sclose(filespace_id);

      dataspace_id = H5Screate(H5S_SCALAR);
      attribute_id = H5Acreate (file_id, "fileNumber", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attribute_id,H5T_NATIVE_INT , &(fd.fileNumber));
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate(H5S_SCALAR);
      attribute_id = H5Acreate (file_id, "numProcPerFile", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attribute_id,H5T_NATIVE_INT , &(fd.numProcPerFile));
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate(H5S_SCALAR);
      attribute_id = H5Acreate (file_id, "world_size", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attribute_id,H5T_NATIVE_INT , &(fd.world_size));
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(1,&two,NULL);
      attribute_id = H5Acreate (file_id, "grid_size", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attribute_id,H5T_NATIVE_INT , fd.grid_size);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);


      dataspace_id = H5Screate_simple(1,&numPartsSize[1],NULL);
      attribute_id = H5Acreate (file_id, "boxSize", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attribute_id,REAL_TYPE , fd.boxSize);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate(H5S_SCALAR);
      attribute_id = H5Acreate (file_id, "fileBoxSize", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attribute_id,REAL_TYPE , &(fd.fileBoxSize));
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate(H5S_SCALAR);
      attribute_id = H5Acreate (file_id, "fileBoxOffset", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attribute_id,REAL_TYPE , &(fd.fileBoxOffset));
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);


      dataspace_id = H5Screate_simple(1,&numPartsSize[1],NULL);
      attribute_id = H5Acreate (file_id, "latSize", H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attribute_id,H5T_NATIVE_INT , fd.latSize);
      H5Aclose(attribute_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(1,numPartsSize,NULL);
      dataset_id = H5Dcreate(file_id, "numParts", H5T_NATIVE_LONG, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,numParts);
      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(2,numPartsSize,NULL);
      dataset_id = H5Dcreate(file_id, "localBoxOffset", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, REAL_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxOffset);
      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);

      dataspace_id = H5Screate_simple(2,numPartsSize,NULL);
      dataset_id = H5Dcreate(file_id, "localBoxSize", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, REAL_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxSize);
      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);

      //

    H5Fclose(file_id);
  }
  MPI_Barrier(comm);

  for(int p=0;p<mpi_size;p++){

      MPI_Barrier(comm);
      if(mpi_rank==p){

          plist_id = H5Pcreate(H5P_FILE_ACCESS);
          file_id = H5Fopen(filename.c_str(),H5F_ACC_RDWR,plist_id);
          H5Pclose(plist_id);

          dataset_id = H5Dopen(file_id, "/data", H5P_DEFAULT);
          filespace_id = H5Dget_space(dataset_id);

          memspace_id = H5Screate_simple(1,&localNumParts,NULL);

          H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, &offset, NULL,&localNumParts, NULL);
          H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &offsetf, NULL,&localNumParts, NULL);

          plist_id = H5Pcreate(H5P_DATASET_XFER);
          H5Dwrite(dataset_id,partdatatype.part_memType, memspace_id, filespace_id, plist_id,partlist);

          H5Pclose(plist_id);
          H5Sclose(memspace_id);
          H5Sclose(filespace_id);
          H5Dclose(dataset_id);

          H5Fclose(file_id);


      }
      MPI_Barrier(comm);
  }


#endif

    return 0;
}


#endif
