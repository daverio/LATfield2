#ifndef LATFIELD2_FIELD_DECL_HPP
#define LATFIELD2_FIELD_DECL_HPP

/*! \file LATfield2_Field.hpp
 \brief LATfield2_Field.hpp contain the class Field definition.
 \author David Daverio, Neil Bevis
 */


//Matrix-style component symmetry

inline int getCompSym(int i,int j)
{
  int component;
  if (i>j) component = i + j * 3 - (j * (1 + j)) / 2;
  else component = j + i * 3 - (i * (1 + i)) / 2;
  return component;
}

inline int gcs(int i,int j)
{
  int component;
  if (i>j) component = i + j * 3 - (j * (1 + j)) / 2;
  else component = j + i * 3 - (i * (1 + i)) / 2;
  return component;
}


extern  int symmetric;
extern  int unsymmetric;

//PROTO-TYPEs============================================
template<class FieldType>
void defaultFieldSave(fstream& file, FieldType* siteData, int components);

template<class FieldType>
void defaultFieldLoad(fstream& file, FieldType* siteData, int components);


/*! \class Field

 \brief The Field class represent a field on a given lattice.


 It stores the description of the field i.e. the datatype (given as template), the number of components, and the pointer to the field array in the memory.


 The datatype is under user control; a field of structur or class can be also used (to be able to use the update halo method, the = operator must be defined or overloaded). However, the I/O support only native datatype and 1d array of them.

 A field can be a single element, a vector of element or a 2d matrix of elements. In the case of a matrix it is possible to define the symmetry of the matrix:

    LATfield2d::unsymmetric  : no symmetry.

    LATfield2d::symmetric    : symmetric matrix (Tij = Tji)

    Antisymmetric field are not yet implemented.

 */
//FIELD CLASS DEFINITION==================================
template<class FieldType>
class Field
	{
	public:

        //! Constructor.
		Field();

        /*!
         Constructor of a "vector" field with initialization and allocation.
         \sa initialize(Lattice& lattice, int components = 1)
         \sa alloc()
         \param lattice    : lattice on which the field is defined
         \param components : number of components. The default is 1.
         */
		Field(Lattice& lattice, int components = 1);

        /*!
         Constructor of a "matrix" field with initialization and allocation.
         \sa initialize(Lattice& lattice, int rows, int cols, int symmetry=unsymmetric)
         \sa alloc()
         \param lattice    : lattice on which the field is defined
         \param matrixRows : matrix number of row .
         \param matrixCols : matrix number of colomn.
         \param symmetry   : symmetry of the matrix, default is unsymmetric. LATfield2d::symmetric can be passed to specify the symmetry, reducing memory usage.
         */
		Field(Lattice& lattice, int rows, int cols, int symmetry=unsymmetric);

        /*!
         Constructor of a vector of "matrix" field with initialization and allocation.
         \sa initialize(Lattice& lattice, int rows, int cols, int symmetry=unsymmetric)
         \sa alloc()
         \param lattice    : lattice on which the field is defined
         \param nMatrix    : number of matrix.
         \param matrixRows : matrix number of row .
         \param matrixCols : matrix number of colomn.
         \param symmetry   : symmetry of the matrix, default is unsymmetric. LATfield2d::symmetric can be passed to specify the symmetry, reducing memory usage.
         */
        Field(Lattice& lattice, int nMatrix, int rows, int cols, int symmetry);

        //!Destructor.
		~Field();

		//INITIALIZATION-TYPE FUNCTIONS
        /*!
         Initialization of a "vector" field. Without allocation.
         \param lattice    : lattice on which the field is defined
         \param components : number of components. Default is 1.
         */
		void initialize(Lattice& lattice, int components = 1);

        /*!
         Initialization of a "matrix" field. Without allocation.
         \param lattice    : lattice on which the field is defined
         \param matrixRows : matrix number of row .
         \param matrixCols : matrix number of colomn.
         \param symmetry   : symmetry of the matrix, default is unsymmetric. LATfield2d::symmetric  can be pass to specify the symmetry.
         */
		void initialize(Lattice& lattice, int rows, int cols, int symmetry=unsymmetric);

        /*!
         Initialization of a vector of "matrix" field. Without allocation.
         \param lattice    : lattice on which the field is defined
         \param nMatrix    : number of matrix.
         \param matrixRows : matrix number of row .
         \param matrixCols : matrix number of colomn.
         \param symmetry   : symmetry of the matrix, default is unsymmetric. LATfield2d::symmetric  can be pass to specify the symmetry.
         */
        void initialize(Lattice& lattice,int nMatrix ,int rows, int cols, int symmetry);

        /*!
         Memory allocation. Allocate the data_ array of this field. It allocated "components_*lattice_->sitesLocalGross()*sizeof(FieldType)" bytes. This method use malloc() to allocate the memory, in case the pointer is not allocated it will return a error message but not exiting the executable.
         */
		void alloc();

        /*!
         Memory allocation. Allocate the data_ array of this field. It allocated "size" bytes if "size" > "components_*lattice_->sitesLocalGross()*sizeof(FieldType)", if not it call this->alloc(). This method use malloc() to allocate the memory, in case the pointer is not allocated it will return a error message but not exiting the executable.
         */

        void alloc(long size);
        /*!
         Free the data_ array.
         */
		void dealloc();

    /* just to be able to check its status */
    bool IsAllocated() { return (status_ & allocated) > 0; }
    bool IsInitialized() { return (status_ & initialized) > 0; }

		//INDEXING


        /*!
         Returns the value of the field stored in data_[index]. User should used operator()(const Site& site) to refer and access to the value of the field.

         \param index: displacment on the data_ array.
         */
		FieldType& operator()(long index);

        /*!
         Returns the value of the field stored in data_[component + index*components_]. User should used operator()(const Site& site, int component) to refer and access to the value of the field.

         \param index    : number of site to skip.
         \param component: index of the desired component.
         */
		FieldType& operator()(long index, int component);

        /*!
         Returns the value of the field stored in data_[j*rows_ + i + index*components_]. In the symmetric case, it returns data_[abs(i-j) + min(i,j)*(rows_+0.5-0.5*min(i,j)) + index*components_]. User should used operator()(const Site& site, int i, int j) to refer and access to the value of the field.

         \param index : number of site to skip.
         \param i     : index of the row
         \param j     : index of the column
         */
		FieldType& operator()(long index, int i, int j);




        FieldType& operator()(long index, int k, int i, int j);


        /*!
         Returns the value of the field at the position pointed by the Site object (data_[site.index()]). Can be used only for field with one component!

         \param site: a site instance which points to the desired lattice site.

         \sa To have more description see the Site class documentation.
         */
		FieldType& operator()(const Site& site);

        /*!
         Returns the value of a (vector) field's component at the position pointed by the Site object (data_[component + site.index()*components_]).


         \param site: a site instance which points to the desired lattice site.
         \param component: index of the desired component.

         \sa To have more description see the Site class documentation.
         */
		FieldType& operator()(const Site& site, int component);

        /*!
         Returns the value of the (i,j) matrix component of the field at the position pointed by the Site object (data_[j*rows_ + i + site.index*components_]). In the symmetric case, it returns data_[abs(i-j) + min(i,j)*(rows_+0.5-0.5*min(i,j)) + site.index()*components_].

         \param site: a site instance which points to the desired lattice site.
         \param i     : index of the row
         \param j     : index of the column

         \sa To have more description see the Site class documentation.
         */
		FieldType& operator()(const Site& site, int i, int j);

        FieldType& operator()(const Site& site, int k, int i, int j);

#ifdef FFT3D

        /*!
         Equivalent to FieldType& operator()(const Site& site) for cKsite
         */
		FieldType& operator()(const cKSite& site);

        /*!
         Equivalent to FieldType& operator()(const Site& site, int component) for cKsite
         */
		FieldType& operator()(const cKSite& site, int component);

        /*!
         Equivalent to FieldType& operator()(const Site& site, int i, int j) for cKsite
         */
		FieldType& operator()(const cKSite& site, int i, int j);

        FieldType& operator()(const cKSite& site,int k, int i, int j);

        /*!
         Equivalent to FieldType& operator()(const Site& site) for rKsite
         */
		FieldType& operator()(const rKSite& site);

        /*!
         Equivalent to FieldType& operator()(const Site& site, int component) for rKsite
         */

		FieldType& operator()(const rKSite& site, int component);

        /*!
         Equivalent to FieldType& operator()(const Site& site, int i, int j) for rKsite
         */
		FieldType& operator()(const rKSite& site, int i, int j);

        FieldType& operator()(const rKSite& site,int k, int i, int j);
#endif

		//BOUNDARY UPDATE

        /*!
         Update the halo sites (ghost cells) of the field. This method us the operator = to assign values, thefor be sure that this operator is defined or correctly overloaded when using field of class or struct.
         */
		void updateHalo();

		//FILE I/O FUNCTIONS

        /*!
         Method to write a field in binary. This method use serial I/O so can be very slow. Should never be used during production, but can be usefull during development.

         \param filename: path to the file, from the executable folder.
         */
	    void write(const string filename);

        /*!
         Method to read a field in binary which have been writen by the void write(const string filename) method.

         \param filename: path to the file, from the executable folder.
         */
        void read(const string filename);

        /*!
         Method to write a field in Binary. This method uses serial I/O so can be very slow, but is faster than void save(const string filename). Should never be used! but it can be usefull on some architectures, where HDF5 is not installed and/or crashes the filesystem. There is no method to read back such a file. The file structure is dependent of the local geometry. This function dumps serially (in the paralle.lat_world_rank order) the data stored in each MPI process.

         \param filename: path to the file, from the executable folder.
         */
        void fastwrite(const string filename);

        /*!
         Method to write a field in ASCII. This method use serial I/O so can be very slow. Should never be used during production, but can be useful during development.

         \param filename      : path to the file, from the executable folder.
         \param FormatFunction: format used for the writting procedure.
         */
	    void save(const string filename,
		      void (*FormatFunction)(fstream&,FieldType*,int) = defaultFieldSave<FieldType>);

        /*!
         Method to read a field in ASCII which have been written by the void write(const string filename) method.

         \param filename      : path to the file, from the executable folder.
         \param FormatFunction: format used for the writting procedure.
         */
	    void load(const string filename,
		      void (*FormatFunction)(fstream&,FieldType*,int) = defaultFieldLoad<FieldType>);

        /*!
         Method to write a field in ASCII. This method use serial I/O so can be very slow, but is faster than void write(const string filename). Should never be used! but it can be usefull on some architectures, where HDF5 is not installed and/or crashes the filesystem. There is no method to read back such a file. The file structur is dependent of the local geometry. This function dumps serially (in the paralle.lat_world_rank order) the data stored in each MPI process.
         */
        void fastsave(const string filename,
                      void (*FormatFunction)(fstream&,FieldType*,int) = defaultFieldSave<FieldType>);
        /*!
         Method to write a field with HDF5. To be able to use this method the flag HDF5 need to be set at compilation (-DHDF5). This method use serial HDF5 by default. For parallel HDF5 the flag -DH5_HAVE_PARALLEL must be used at compilation.

         This methods will write 1 dataset (named "/field") which contain all components of the field. If one want to use a dataset per components (named "/comp0" to "/compN") the flag -DH5_HAVE_PIXIE need to be set at compilation.

         \param filename : path to the file, from the executable folder.
         */

	    void saveHDF5(string filename, string dataset_name);
	    void saveHDF5(string filename){this->saveHDF5(filename, "data");}

	    /*!
         Method to load a field with HDF5. To be able to use this method the flag HDF5 need to be set at compilation (-DHDF5). This method use serial HDF5 by default. For parallel HDF5 the flag -DH5_HAVE_PARALLEL must be set at compilation.

         This methods will expect 1 dataset named field which contain all component of the field. If one want to use a dataset per components (named comp0 to compN) the flag -DH5_HAVE_PIXIE need to be set at compilation.

         \param filename : path to the file, from the executable folder.
         */


        void loadHDF5(string filename,string dataset_name);
	    	void loadHDF5(string filename){this->loadHDF5(filename, "data");}


        /*!
         A way to save coarse grained version of the fields. To be able to use this method the flag HDF5 need to be set at compilation (-DHDF5). Work only for 3D lattice!!!

         This methods will write 1 dataset (named "/field") which contain all component of the field. If one want to use a dataset per components (named "/comp0" to "/compN") the flag -DH5_HAVE_PIXIE need to be set at compilation.

         \param filename : path to the file, from the executable folder.
         \param ration   : ration of the coarse graining. Must be an integer divider of the size of each dimension of this->lattice()
         */

	    void saveHDF5_coarseGrain3D(string filename,string dataset_name ,int ratio);
	    void saveHDF5_coarseGrain3D(string filename,int ratio){this->saveHDF5_coarseGrain3D(filename,"data",ratio);}


        /*!
         Save a slice perpendicular to the first coordinate, at xcoord. To be able to use this method the flag HDF5 need to be set at compilation (-DHDF5).

          This methods will write 1 dataset (named "/field") which contain all component of the field. If one want to use a dataset per components (named "/comp0" to "/compN") the flag -DH5_HAVE_PIXIE need to be set at compilation.

         \param filename : path to the file, from the executable folder.
         \param xcoord   : coordinate of the slice on the first dimension of the lattice.
         \param thickness: thickness of the slice, the default is 1.
         */

	    void saveSliceHDF5(string filename,string dataset_name , int xcoord, int thickness = 1);
	    void saveSliceHDF5(string filename, int xcoord, int thickness = 1){this->saveSliceHDF5(filename,"data",xcoord,thickness);}





#ifdef EXTERNAL_IO
        void saveHDF5_server_open(string filename_base,int offset=-1,int thickness=-1);
        void saveHDF5_server_write(int number_of_message = 4, string filename_base = "defaultfilename",int offset=-1,int thickness=-1);
#endif


	    //MISCELLANEOUS
	    /*!
         Returns a pointer to the lattice on which the field is defined.
	    */
		Lattice& lattice();
        /*!
         Returns the number of components of the field at each sites.
         */
		int   components();
        /*!
         Returns the number of rows of the component matrix at each sites.
         */
		int   rows();
        /*!
         Returns the number of columns of the component matrix at each sites.
         */
		int   cols();
    /*!
     Returns the number of component matrix at each sites.
     */
    int   nMatrix();
        /*!
         returns the symmetry of the component matrix at each sites.
         */
		int   symmetry();

        /*!
         Returns the pointer to the data_ array of the field.
         */
		FieldType*& data();
#ifdef HDF5
    hid_t type_id();
    int array_size();
#endif



	private:
		//PRIVATE FUNCTIONS
		void updateHaloComms();

#ifdef HDF5
        void get_h5type();
#endif
    public:
        FieldType* data_;
	protected:
		//MEMBER DATA
		Lattice*   lattice_;

		int        components_;
		int        rows_;
		int        cols_;
    int        nMatrix_;
    int        matrixSize_;
		int        symmetry_;
		unsigned int sizeof_fieldType_;

		int        status_;
		static int initialized;
		static int allocated;

    unsigned long long data_memSize_;
#ifdef HDF5
        hid_t type_id_;
        int array_size_;
#endif
#ifdef EXTERNAL_IO
        ioserver_file io_file_;
        int iof_offset_;
        int iof_thickness_;
#endif

	};

template <class FieldType>
int Field<FieldType>::initialized = 1;        //Status flag for initialized
template <class FieldType>
int Field<FieldType>::allocated = 2;          //Status flag for allocated memory

//CONSTRUCTORS=================

template <class FieldType>
Field<FieldType>::Field() : data_memSize_(0), data_(NULL) {
    status_=0;
#ifdef HDF5
    this->get_h5type();
#endif
}

template <class FieldType>
Field<FieldType>::Field(Lattice& lattice, int components) : data_memSize_(0), data_(NULL)
{
	status_=0;
	this->initialize(lattice, components);
	this->alloc();
#ifdef HDF5
    this->get_h5type();
#endif

}

template <class FieldType>
Field<FieldType>::Field(Lattice& lattice, int rows, int cols, int symmetry) : data_memSize_(0), data_(NULL)
{
	status_=0;
	this->initialize(lattice, rows, cols, symmetry);
	this->alloc();
#ifdef HDF5
    this->get_h5type();
#endif

}

template <class FieldType>
Field<FieldType>::Field(Lattice& lattice, int nMatrix, int rows, int cols, int symmetry) : data_memSize_(0), data_(NULL)
{
    status_=0;
    this->initialize(lattice,nMatrix, rows, cols, symmetry);
    this->alloc();
#ifdef HDF5
    this->get_h5type();
#endif

}


#ifdef HDF5
template <class FieldType>
void Field<FieldType>::get_h5type()
{
    string type_name;
	type_name = typeid(FieldType).name();
	char  str_array_size_[20];
	//string str_array_size_;
	array_size_ = 1;
	int nt;

	nt = 0;
	//get the array size
	if(type_name[0]=='A')
	{

		while(type_name[nt+1]!='_')
		{
			str_array_size_[nt]=type_name[nt+1];
			nt++;
		}
		str_array_size_[nt]='\0';

		nt = nt + 2;
		if(type_name[nt]=='A')
		{
			COUT<<"LATField2d::Field::save_hdf5  :  Cannot recognize type of field: "<< type_name <<endl;
			COUT<<"LATField2d::Field::save_hdf5  :  ---------------------------------------------------------------------------------" <<endl;
			COUT<<"LATField2d::Field::save_hdf5  :  List of accepted type : short, unsigned short, int, unsigned int, long," <<endl;
			COUT<<"LATField2d::Field::save_hdf5  :  unsigned long, long long, unsigned long long, float, double, long double, Imag ." <<endl;
			COUT<<"LATField2d::Field::save_hdf5  :  1 dimensional array of those type are also accepted" <<endl;
			COUT<<"LATField2d::Field::save_hdf5  :  ---------------------------------------------------------------------------------" <<endl;
			return;

		}

		array_size_=atoi(str_array_size_);
	}
	//get the type

	if(type_name[nt]=='s')
	{
		//COUT << " type : short ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_SHORT;
	}
	else if(type_name[nt]=='t')
	{
		//COUT << " type : unsigned short ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_USHORT;
	}

	else if(type_name[nt]=='i')
	{
		//COUT << " type : int ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_INT;
	}
	else if(type_name[nt]=='j')
	{
		//COUT << " type : unsigned int ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_UINT;
	}

	else if(type_name[nt]=='l')
	{
		//COUT << " type : long ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_LONG;
	}
	else if(type_name[nt]=='m')
	{
		//COUT << " type : unsigned long ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_ULONG;
	}
	else if(type_name[nt]=='x')
	{
		//COUT << " type : long long ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_LLONG;
	}
	else if(type_name[nt]=='y')
	{
		//COUT << " type : unsigned long long ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_ULLONG;
	}
	else if(type_name[nt]=='f')
	{
		//COUT << " type : float ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_FLOAT;
	}
	else if(type_name[nt]=='d')
	{
		//COUT << " type : double ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_DOUBLE;
	}

	else if(type_name[nt]=='e')
	{
		//COUT << " type : long double ; size : "<<array_size_<<endl;
		type_id_ = H5T_NATIVE_LDOUBLE;
	}
	else if (type_name[nt]=='N'&&type_name[nt+12]=='I'&&type_name[nt+13]=='m'&&type_name[nt+14]=='a'&&type_name[nt+15]=='g')
	{
        type_id_ = H5Tcreate (H5T_COMPOUND, sizeof (Imag));
#ifdef SINGLE
        H5Tinsert (type_id_, "real", 0,H5T_NATIVE_FLOAT);
        H5Tinsert (type_id_, "imaginary", sizeof(Real),H5T_NATIVE_FLOAT);
#else
	    H5Tinsert (type_id_, "real", 0,H5T_NATIVE_DOUBLE);
        H5Tinsert (type_id_, "imaginary", sizeof(Real),H5T_NATIVE_DOUBLE);
#endif
    }
    else if (type_name[nt+12]=='p' && type_name[nt+13]=='a' && type_name[nt+14]=='r' && type_name[nt+15]=='t' &&
             type_name[nt+16]=='L' && type_name[nt+17]=='i' && type_name[nt+18]=='s' && type_name[nt+19]=='t')
    {
        //COUT<< "initializing particle field... unable to call Field.saveHDF5"<<endl;
    }
	else
	{
		COUT<<"LATField2d::Field::save_hdf5  :  Cannot recognize type of field: " << type_name <<endl;
		COUT<<"LATField2d::Field::save_hdf5  :  ---------------------------------------------------------------------------------" <<endl;
		COUT<<"LATField2d::Field::save_hdf5  :  List of accepted type : short, unsigned short, int, unsigned int, long," <<endl;
		COUT<<"LATField2d::Field::save_hdf5  :  unsigned long, long long, unsigned long long, float, double, long double, Imag ." <<endl;
		COUT<<"LATField2d::Field::save_hdf5  :  1 dimensional array of those type are also accepted" <<endl;
		COUT<<"LATField2d::Field::save_hdf5  :  ---------------------------------------------------------------------------------" <<endl;
        /*for(int i = 0;i<9 ;i++)
        {
            COUT<<type_name[nt+i]<<endl;
        }*/
	}


}
#endif

//DESTRUCTOR=======================

template <class FieldType>
Field<FieldType>::~Field()
{
  if(status_ & allocated) this->dealloc();

}

template <class FieldType>
void Field<FieldType>::initialize(Lattice& lattice, int components)
{
	if(status_ & allocated) { this->dealloc(); }

	sizeof_fieldType_ = sizeof(FieldType);
	status_= initialized;
	lattice_=&lattice;
	components_=components;
	rows_=components_;
	cols_=1;
    nMatrix_ = 1;
	symmetry_=unsymmetric;
    matrixSize_ = components;


}

template <class FieldType>
void Field<FieldType>::initialize(Lattice& lattice, int rows, int cols, int symmetry)
{
	int components;

	if(status_ & allocated) { this->dealloc(); }

	sizeof_fieldType_ = sizeof(FieldType);
	status_= initialized;
	lattice_=&lattice;
	rows_=rows;
	cols_=cols;
    nMatrix_ =1;
	symmetry_=symmetry;
	if(symmetry_==symmetric) { components = ( rows_ * (rows_+1) ) / 2; }
	else { components = rows*cols; }
	components_=components;
    matrixSize_ = components;


}
template <class FieldType>
void Field<FieldType>::initialize(Lattice& lattice,int nMatrix, int rows, int cols, int symmetry)
{
    int components;

    if(status_ & allocated) { this->dealloc(); }

    sizeof_fieldType_ = sizeof(FieldType);
    status_= initialized;
    lattice_=&lattice;
    rows_=rows;
    cols_=cols;
    nMatrix_ = nMatrix;
    symmetry_=symmetry;
    if(symmetry_==symmetric) { components = ( rows_ * (rows_+1) ) / 2; }
    else { components = rows*cols; }
    matrixSize_ = components;
    components_=components * nMatrix;


}

template <class FieldType>
void Field<FieldType>::alloc()
{
	if((status_ & allocated) == 0)
    {

    data_memSize_ = components_*lattice_->sitesLocalGross();
		data_= new FieldType[data_memSize_];
		status_ = status_ | allocated;

        if(data_==NULL)
        {
            cout<<"LATField2d::Field::alloc(long size)  :process "<< parallel.rank() <<" cannot allocate the field data array."<<endl;
            data_memSize_ = 0;
            allocated = 0;
        }

    }
}

template <class FieldType>
void Field<FieldType>::alloc(long size)
{
    long alloc_number;

    if(size < lattice_->sitesLocalGross()) alloc_number =  lattice_->sitesLocalGross();
    else alloc_number= size;

	if((status_ & allocated) == 0)
    {

    data_memSize_ = components_* size;
		data_= new FieldType[data_memSize_];
		status_ = status_ | allocated;

        if(data_==NULL)
        {
            cout<<"LATField2d::Field::alloc(long size)  :process "<< parallel.rank() <<" cannot allocate the field data array."<<endl;
            data_memSize_ = 0;
            allocated = 0;
        }
    }
    else if(size > lattice_->sitesLocalGross())
    {
        this->dealloc();
        this->alloc(size);
    }


}



template <class FieldType>
void Field<FieldType>::dealloc()
{
	if((status_ & allocated) > 0)
    {
	//std::cerr << "Going to deallocate, size: " << data_memSize_ << "\n";// Commented out to descrease outputs....
      delete[] data_;
      //std::cerr << "Finished deallocate, size: " << data_memSize_ << "\n";// Commented out to descrease outputs....
    data_memSize_ = 0;
		status_= status_ ^ allocated;
    }
}


//FIELD INDEXING===============

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(long index)
{
	return data_[index];
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(long index, int component)
{
	return data_[index*components_ + component];
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(long index, int i, int j)
{
	int component;
	if(symmetry_==symmetric)
    {
		if (i>j) component = i + j * rows_ - (j * (1 + j)) / 2;
		else component = j + i * rows_ - (i * (1 + i)) / 2;
    }
	else { component = j*rows_ + i; }
	return data_[index*components_ + component];
}


template <class FieldType>
inline FieldType& Field<FieldType>::operator()(long index, int k, int i, int j)
{
    int component;
    if(symmetry_==symmetric)
    {
        if (i>j) component = i + j * rows_ - (j * (1 + j)) / 2;
        else component = j + i * rows_ - (i * (1 + i)) / 2;
    }
    else { component = j*rows_ + i; }
    component += matrixSize_ * k;
    return data_[index*components_ + component];
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const Site& site)
{
	return this->operator()(site.index());
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const Site& site, int component)
{
	return this->operator()(site.index(),component);
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const Site& site, int i, int j)
{
	return this->operator()(site.index(),i,j);
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const Site& site,int k, int i, int j)
{
    return this->operator()(site.index(),k,i,j);
}

#ifdef FFT3D

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const cKSite& site)
{
	return this->operator()(site.index());
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const cKSite& site, int component)
{
	return this->operator()(site.index(),component);
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const cKSite& site, int i, int j)
{
	return this->operator()(site.index(),i,j);
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const cKSite& site, int k, int i, int j)
{
    return this->operator()(site.index(),k,i,j);
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const rKSite& site)
{
	return this->operator()(site.index());
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const rKSite& site, int component)
{
	return this->operator()(site.index(),component);
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const rKSite& site, int i, int j)
{
	return this->operator()(site.index(),i,j);
}

template <class FieldType>
inline FieldType& Field<FieldType>::operator()(const rKSite& site, int k, int i, int j)
{
    return this->operator()(site.index(),k,i,j);
}

#endif

//FIELD BOUNDARY UPDATE=========
////////////////////////////////

template <class FieldType>
void Field<FieldType>::updateHalo()
{

	Site site(*lattice_);
	int copyfrom;

	int  dim=lattice_->dim();
	int* jump=new int[dim];
	int* size=new int[dim];

	for(int i=0; i<dim; i++)
	{
		jump[i]=lattice_->jump(i);
		size[i]=lattice_->sizeLocal(i);
	}

	for(site.haloFirst(); site.haloTest(); site.haloNext())
	{
		//Work out where to copy from
		copyfrom=site.index();
		for(int i=0; i<dim; i++)
		{
			if( site.coordLocal(i)<0 )
			{
				copyfrom += jump[i] * size[i];
			}
			else if( site.coordLocal(i) >= size[i] )
			{
				copyfrom -= jump[i] * size[i];
			}
		}
		//Copy data
		for(int i=0; i<components_; i++)
		{
			data_[site.index()*components_+i] = data_[copyfrom*components_+i];
			//memcpy(data_[site.index()*components_+i],data_[copyfrom*components_+i],sizeof_fieldType_);

		}
	}

	delete[] jump;
	delete[] size;

	if( parallel.size()>1 ) { updateHaloComms(); }
}

template <class FieldType>
void Field<FieldType>::updateHaloComms()
{

	int buffer_size0, buffer_size1,temp;
	int i,j;

	//Size of buffer : max size between 2 scatered dimension;
	buffer_size0 = lattice_->halo() * components_*lattice_->jump(lattice_->dim()-1);
	buffer_size1 = lattice_->halo() * components_*lattice_->jump(lattice_->dim()-2) *lattice_->sizeLocal(lattice_->dim()-1);

	if(buffer_size0>buffer_size1)temp=buffer_size0;
	else temp=buffer_size1;

	FieldType* buffer_send = new FieldType[ temp ];
	FieldType* buffer_rec = new FieldType[ temp ];

	FieldType* pointer_send_up;
	FieldType* pointer_send_down;
	FieldType* pointer_rec_up;
	FieldType* pointer_rec_down;


	pointer_send_up = data_ + ((lattice_->halo()+1)*lattice_->jump(lattice_->dim()-1) - 2*lattice_->halo()*lattice_->jump(lattice_->dim()-2))*components_;
	pointer_rec_up = data_ + ((lattice_->halo()+1)*lattice_->jump(lattice_->dim()-1) - lattice_->halo()*lattice_->jump(lattice_->dim()-2))*components_;

	pointer_send_down = data_ + buffer_size0 + lattice_->jump(lattice_->dim()-2)*lattice_->halo()*components_   ;
	pointer_rec_down = data_ + buffer_size0;

	if(parallel.grid_rank()[1]%2==0)
	{


		for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
		{
			for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
			{
				buffer_send[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)]= pointer_send_up[i+j*lattice_->jump(lattice_->dim()-1)*components_];

			}
		}

		if(parallel.grid_rank()[1]!=parallel.grid_size()[1]-1)
		{
			parallel.send_dim1( buffer_send, buffer_size1, parallel.grid_rank()[1]+1);
			parallel.receive_dim1( buffer_rec, buffer_size1, parallel.grid_rank()[1]+1);
		}

		for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
		{
			for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
			{
				pointer_rec_up[i+j*lattice_->jump(lattice_->dim()-1)*components_] = buffer_rec[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)];
			}
		}


		for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
		{
			for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
			{
				buffer_send[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)]= pointer_send_down[i+j*lattice_->jump(lattice_->dim()-1)*components_];
			}
		}


		if(parallel.grid_rank()[1] != 0)
		{
			parallel.send_dim1( buffer_send, buffer_size1, parallel.grid_rank()[1]-1);
			parallel.receive_dim1( buffer_rec, buffer_size1,  parallel.grid_rank()[1]-1);
		}
		else if(parallel.grid_size()[1]%2==0)
		{
			parallel.send_dim1( buffer_send, buffer_size1, parallel.grid_size()[1]-1);
			parallel.receive_dim1( buffer_rec, buffer_size1,  parallel.grid_size()[1]-1);
		}

		for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
		{
			for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
			{
				pointer_rec_down[i+j*lattice_->jump(lattice_->dim()-1)*components_] = buffer_rec[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)];
			}
		}


	}
	else
	{
        for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
		{
			for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
			{
				buffer_send[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)]= pointer_send_down[i+j*lattice_->jump(lattice_->dim()-1)*components_];
			}
		}

		parallel.receive_dim1( buffer_rec, buffer_size1, parallel.grid_rank()[1]-1);
		parallel.send_dim1( buffer_send, buffer_size1, parallel.grid_rank()[1]-1);


		for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
		{
			for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
			{
				pointer_rec_down[i+j*lattice_->jump(lattice_->dim()-1)*components_] = buffer_rec[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)];
			}
		}


		for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
		{
			for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
			{
				buffer_send[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)]= pointer_send_up[i+j*lattice_->jump(lattice_->dim()-1)*components_];
			}
		}

		if(parallel.grid_rank()[1]!=parallel.grid_size()[1]-1)
		{
			parallel.receive_dim1( buffer_rec, buffer_size1, parallel.grid_rank()[1]+1);
			parallel.send_dim1( buffer_send, buffer_size1, parallel.grid_rank()[1]+1);
		}
		else
		{
			parallel.receive_dim1( buffer_rec, buffer_size1,0);
			parallel.send_dim1( buffer_send, buffer_size1, 0);
		}

		for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
		{
			for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
			{
				pointer_rec_up[i+j*lattice_->jump(lattice_->dim()-1)*components_] = buffer_rec[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)];
			}
		}


	}


	if(parallel.grid_size()[1]%2!=0)
	{
		if(parallel.grid_rank()[1]==0)
		{
			for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
			{
				for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
				{
					buffer_send[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)]= pointer_send_down[i+j*lattice_->jump(lattice_->dim()-1)*components_];

				}
			}

			parallel.send_dim1( buffer_send, buffer_size1, parallel.grid_size()[1]-1);
			parallel.receive_dim1( buffer_rec, buffer_size1,  parallel.grid_size()[1]-1);

			for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
			{
				for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
				{
					pointer_rec_down[i+j*lattice_->jump(lattice_->dim()-1)*components_] = buffer_rec[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)];
				}
			}

		}
		if(parallel.grid_rank()[1]==parallel.grid_size()[1]-1)
		{
			for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
			{
				for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
				{
					buffer_send[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)]= pointer_send_up[i+j*lattice_->jump(lattice_->dim()-1)*components_];
				}
			}
			parallel.receive_dim1( buffer_rec, buffer_size1,0);
			parallel.send_dim1( buffer_send, buffer_size1, 0);

			for(j=0;j<lattice_->sizeLocal(lattice_->dim()-1);j++)
			{
				for(i=0;i<buffer_size1/lattice_->sizeLocal(lattice_->dim()-1);i++)
				{
					pointer_rec_up[i+j*lattice_->jump(lattice_->dim()-1)*components_] = buffer_rec[i+j*buffer_size1/lattice_->sizeLocal(lattice_->dim()-1)];

				}
			}
		}
	}



	pointer_send_up = data_ + lattice_->sitesLocalGross() * components_ - 2*buffer_size0;
	pointer_rec_up = data_ + lattice_->sitesLocalGross() * components_ - buffer_size0;

	pointer_send_down = data_ + buffer_size0;
	pointer_rec_down = data_;


	if(parallel.grid_rank()[0]%2==0)
	{

		if(parallel.grid_rank()[0]!=parallel.grid_size()[0]-1)
		{
			parallel.send_dim0( pointer_send_up, buffer_size0, parallel.grid_rank()[0]+1);
			parallel.receive_dim0( pointer_rec_up, buffer_size0, parallel.grid_rank()[0]+1);
		}




		if(parallel.grid_rank()[0] != 0)
		{
			parallel.send_dim0( pointer_send_down, buffer_size0, parallel.grid_rank()[0]-1);
			parallel.receive_dim0( pointer_rec_down, buffer_size0,  parallel.grid_rank()[0]-1);
		}
		else if(parallel.grid_size()[0]%2==0)
		{
			parallel.send_dim0( pointer_send_down, buffer_size0, parallel.grid_size()[0]-1);
			parallel.receive_dim0( pointer_rec_down, buffer_size0,  parallel.grid_size()[0]-1);
		}


	}
	else
	{

		parallel.receive_dim0( pointer_rec_down, buffer_size0, parallel.grid_rank()[0]-1);
		parallel.send_dim0( pointer_send_down, buffer_size0, parallel.grid_rank()[0]-1);


		if(parallel.grid_rank()[0]!=parallel.grid_size()[0]-1)
		{
			parallel.receive_dim0( pointer_rec_up, buffer_size0, parallel.grid_rank()[0]+1);
			parallel.send_dim0( pointer_send_up, buffer_size0, parallel.grid_rank()[0]+1);
		}
		else
		{
			parallel.receive_dim0( pointer_rec_up, buffer_size0,0);
			parallel.send_dim0( pointer_send_up, buffer_size0, 0);
		}



	}

	if(parallel.grid_size()[0]%2!=0)
	{
		if(parallel.grid_rank()[0]==0)
		{
			parallel.send_dim0( pointer_send_down, buffer_size0, parallel.grid_size()[0]-1);
			parallel.receive_dim0( pointer_rec_down, buffer_size0,  parallel.grid_size()[0]-1);
		}
		if(parallel.grid_rank()[0]==parallel.grid_size()[0]-1)
		{
			parallel.receive_dim0( pointer_rec_up, buffer_size0,0);
			parallel.send_dim0( pointer_send_up, buffer_size0, 0);
		}
	}




	delete[] buffer_send;
	delete[] buffer_rec;


}



	////////////////////////////////

	//FILE BINARY I/O FUNCTIONS=============
template <class FieldType>
void Field<FieldType>::write(const string filename)
{
	fstream file;
	Site    x(*lattice_);
	int     p,j,i,k;
	int start,stop;
	int haloreverse = 0;
	for(i=0;i<lattice_->dim()-2;i++) haloreverse+=lattice_->jump(i);


	//world_rank 0 creat the file , or if exist trunc it with nothing
	if(parallel.rank()==0)
	{
		file.open(filename.c_str(), fstream::out | fstream::trunc);
		file.close();
	}
	MPI_Barrier(parallel.lat_world_comm());

	for(p=0;p<parallel.grid_size()[0];p++)
	{
		if(parallel.grid_rank()[0]==p)
		{
			for(k=0;k<lattice_->sizeLocal(lattice_->dim()-1);k++)
			{
				for(j=0;j<parallel.grid_size()[1];j++)
				{
					if(parallel.grid_rank()[1]==j)
					{
						file.open(filename.c_str(), fstream::out | fstream::app);

						if(!file.is_open())
						{
							cout<<"Latfield::Field::save - Could not open file for writing"<<endl;
							cout<<"Latfield::Field::save - File: "<<filename<<endl;
							parallel.abortRequest();
						}

						start = lattice_->siteFirst()  + lattice_->jump(lattice_->dim()-1) * k; //oky
						stop = start + (lattice_->jump(lattice_->dim()-2))*lattice_->sizeLocal(lattice_->dim()-2) -2 * haloreverse * lattice_->halo();

						for(x.setIndex(start);x.index()<stop;x.next())
						{
							file.write( (char*)(data_+x.index()*components_),
									   components_*sizeof(FieldType)*lattice_->sizeLocal(0) );
							x.indexAdvance(lattice_->sizeLocal(0)-1);
						}
						file.close();


					}
					MPI_Barrier(parallel.dim1_comm()[p]);
				}
			}
		}
		MPI_Barrier(parallel.lat_world_comm());
	}
	/* cout file size
	if(parallel.rank()==0)
	{
		file.open(filename.c_str(), fstream::out | fstream::app);
		file.seekg(0,ios::end);
		cout<<"file size : "<<file.tellg()<<endl;
		file.close();
	}*/
}

template <class FieldType>
void Field<FieldType>::fastwrite(const string filename)
{
	fstream file;
	Site    x(*lattice_);
	int     i;
	int     p;
	const string fn = filename + "fast";
	//Loop over processes
	for( p=0; p<parallel.size(); p++ )
	{
		if( parallel.rank()==p )
		{
			if(p==0)
			{
				//Truncate file if first process
				file.open(fn.c_str(), fstream::out | fstream::binary | fstream::trunc);
			}
			else
			{
				//Append to file if another process
				file.open(fn.c_str(), fstream::out | fstream::binary | fstream::app);
			}
			if(!file.is_open())
			{
				cerr<<"Latfield::Field::write - Could not open file for writing"<<endl;
				cerr<<"Latfield::Field::write - File: "<<filename<<endl;
				parallel.abortRequest();
			}

			for( x.first(); x.test(); x.next() )
			{
				file.write( (char*)(data_+x.index()*components_),
						   components_*sizeof(FieldType)*lattice_->sizeLocal(0) );
				x.indexAdvance(lattice_->sizeLocal(0)-1);
			}
			file.close();
		}
		parallel.barrier();
	}

}

template <class FieldType>
void Field<FieldType>::read(const string filename)
{
	fstream file;
	Site    x(*lattice_);
	int p,i;
	int sel = components_*sizeof(FieldType);
	int start,stop,skip;
	int haloreverse = 0;
	for(i=0;i<lattice_->dim()-2;i++) haloreverse+=lattice_->jump(i);

	for( p=0; p<parallel.size(); p++ )
	{
		if( parallel.rank()==p )
		{
			file.open(filename.c_str(), fstream::in | fstream::binary);

			if(!file.is_open())
			{
				cerr<<"Latfield::Field::read - Could not open file for reading"<<endl;
				cerr<<"Latfield::Field::read - File: "<<filename<<endl;
				parallel.abortRequest();
			}

			file.seekg( (lattice_->sitesSkip2d())*sel, fstream::beg );

			start = lattice_->siteFirst();
			stop = start + (lattice_->jump(lattice_->dim()-2))*lattice_->sizeLocal(lattice_->dim()-2) - 2*haloreverse * lattice_->halo();
			skip = 1;
			for(i=0;i<lattice_->dim()-1;i++)skip*=lattice_->size()[i];

			for( x.setIndex(start); x.index()<stop ; x.next() )
			{
				file.read((char*)(data_+x.index()*components_),components_*sizeof(FieldType)*lattice_->sizeLocal(0) );
				x.indexAdvance(lattice_->sizeLocal(0)-1);
			}

			start += lattice_->jump(lattice_->dim()-1);
			stop = start + (lattice_->jump(lattice_->dim()-2))*lattice_->sizeLocal(lattice_->dim()-2) - 2*haloreverse * lattice_->halo();

			for(i=1;i<lattice_->sizeLocal(lattice_->dim()-1);i++)
			{

				file.seekg( sel*(lattice_->sitesSkip2d()+i*skip), fstream::beg );

				for( x.setIndex(start); x.index()<stop ; x.next() )
				{
					file.read((char*)(data_+x.index()*components_),components_*sizeof(FieldType)*lattice_->sizeLocal(0) );
					x.indexAdvance(lattice_->sizeLocal(0)-1);
				}

				start += lattice_->jump(lattice_->dim()-1);
				stop = start + (lattice_->jump(lattice_->dim()-2))*lattice_->sizeLocal(lattice_->dim()-2) - 2*haloreverse * lattice_->halo();
			}
			file.close();
		}
		parallel.barrier();
	}
	this->updateHalo();

}

//FILE ASCII I/O FUNCTIONS=======================
template <class FieldType>
void Field<FieldType>::save(const string filename, void (*FormatFunction)(fstream&,FieldType*, int))
{
	fstream file;
	Site    x(*lattice_);
	int     p,i,j,k;

	int haloreverse = 0;
	for(i=0;i<lattice_->dim()-2;i++) haloreverse+=lattice_->jump(i);

	if(parallel.rank()==0)
	{
		file.open(filename.c_str(), fstream::out | fstream::trunc);
		file.close();
	}
	MPI_Barrier(parallel.lat_world_comm());

	for(p=0;p<parallel.grid_size()[0];p++)
	{
		if(parallel.grid_rank()[0]==p)
		{
			for(k=0;k<lattice_->sizeLocal(lattice_->dim()-1);k++)
			{
				for(j=0;j<parallel.grid_size()[1];j++)
				{

					if(parallel.grid_rank()[1]==j && parallel.grid_rank()[0]==p)
					{

						file.open(filename.c_str(), fstream::out | fstream::app);

						if(!file.is_open())
						{
							cout<<"Latfield::Field::save - Could not open file for writing"<<endl;
							cout<<"Latfield::Field::save - File: "<<filename<<endl;
							parallel.abortRequest();
						}


						int start,stop;

						start = lattice_->siteFirst()  + lattice_->jump(lattice_->dim()-1) * k; //oky
						stop = start + (lattice_->jump(lattice_->dim()-2))*lattice_->sizeLocal(lattice_->dim()-2) - 2*haloreverse * lattice_->halo();


						for(x.setIndex(start);x.index()<stop;x.next())
						{
							(*FormatFunction)(file, data_+x.index()*components_, components_);
						}
						file.close();
					}
					MPI_Barrier(parallel.dim1_comm()[p]);
				}

			}


		}
		MPI_Barrier(parallel.lat_world_comm());
	}

}

template <class FieldType>
void Field<FieldType>::fastsave(const string filename,
							void (*FormatFunction)(fstream&,FieldType*, int))
{
	const string arch = filename + "arch";
	const string fn = filename + "fast";


	if(!lattice_->is_arch_saved())
	{
		lattice_->save_arch(arch);
	}

	{
		fstream file;
		Site    x(*lattice_);
		int     p;

		//Loop over processes write data
		for( p=0; p<parallel.size(); p++ )
		{
			if( parallel.rank()==p )
			{
				if( parallel.rank()==0)
				{
					//Truncate file if first process
					file.open(fn.c_str(), fstream::out | fstream::trunc);
				}
				else
				{
					//Append to file if another process
					file.open(fn.c_str(), fstream::out | fstream::app);
				}
				if(!file.is_open())
				{
					cerr<<"Latfield::Field::save - Could not open file for writing"<<endl;
					cerr<<"Latfield::Field::save - File: "<<fn<<endl;
					parallel.abortRequest();
				}

				for( x.first(); x.test(); x.next() )
				{
					(*FormatFunction)(file, data_+x.index()*components_, components_);
				}
				file.close();
			}
			MPI_Barrier(parallel.lat_world_comm());
		}
	}


}


template <class FieldType>
void Field<FieldType>::load(const string filename,
							void (*FormatFunction)(fstream&,FieldType*, int) )
{
	fstream file;
	Site    x(*lattice_);
	int     p,i;
	int     sel;   //site entry length

	int haloreverse = 0;
	for(i=0;i<lattice_->dim()-2;i++) haloreverse+=lattice_->jump(i);



	//Loop over processes
	for( p=0; p<parallel.size(); p++ )
	{
		if( parallel.rank()==p )
		{
            cout<<p<<endl;
			//Open file
			file.open(filename.c_str(), fstream::in);
			//Check for open
			if(!file.is_open())
			{
				cerr<<"Latfield::Field::load - Could not open file for reading"<<endl;
				cerr<<"Latfield::Field::load - File: "<<filename<<endl;
				parallel.abortRequest();
			}

			//Get site entry length
			x.first();
			(*FormatFunction)(file, data_+x.index()*components_, components_);
			sel=int(file.tellg())+1;

			//Move to start of local data
			file.seekg( sel*lattice_->sitesSkip2d(), fstream::beg );

			int start,stop,skip;

			start = lattice_->siteFirst();
			stop = start + (lattice_->jump(lattice_->dim()-2))*lattice_->sizeLocal(lattice_->dim()-2) - 2* haloreverse * lattice_->halo();
			skip = 1;
			for(i=0;i<lattice_->dim()-1;i++)skip*=lattice_->size()[i];

			for( x.setIndex(start); x.index()<stop ; x.next() )
			{
				(*FormatFunction)(file, data_+x.index()*components_, components_);
			}


			start += lattice_->jump(lattice_->dim()-1);
			stop = start + (lattice_->jump(lattice_->dim()-2))*lattice_->sizeLocal(lattice_->dim()-2) - 2 * haloreverse * lattice_->halo();

			for(i=1;i<lattice_->sizeLocal(lattice_->dim()-1);i++)
			{

				file.seekg( sel*(lattice_->sitesSkip2d()+i*skip), fstream::beg );

				for( x.setIndex(start); x.index()<stop ; x.next() )
				{
					(*FormatFunction)(file, data_+x.index()*components_, components_);
				}

				start += lattice_->jump(lattice_->dim()-1);
				stop = start + (lattice_->jump(lattice_->dim()-2))*lattice_->sizeLocal(lattice_->dim()-2) - 2 * haloreverse * lattice_->halo();



			}



			file.close();
		}
		parallel.barrier();
	}
	this->updateHalo();


}

template <class FieldType>
void  Field<FieldType>::saveHDF5(string filename, string dataset_name)
{
#ifdef HDF5

	save_hdf5(data_,type_id_,array_size_,lattice_->coordSkip(),lattice_->size(),lattice_->sizeLocal(),lattice_->halo(),lattice_->dim(),components_,filename,dataset_name);
#else
    COUT<<"LATfield2d must be compiled with HDF5 (flag HDF5 turn on!)"<<endl;
    COUT<<"to be able to use hdf5 data format!!!)"<<endl;
    COUT<<"saving file in binary: "<<filename<<"BIN"<<endl;
    this->write(filename+"BIN");
#endif
}
template <class FieldType>
void  Field<FieldType>::loadHDF5(string filename, string dataset_name)
{
#ifdef HDF5
	load_hdf5(data_,lattice_->coordSkip(),lattice_->size(),lattice_->sizeLocal(),lattice_->halo(),lattice_->dim(),components_,filename,dataset_name);
#else
    COUT<<"LATfield2d must be compiled with HDF5 (flag HDF5 turn on!)"<<endl;
    COUT<<"to be able to use hdf5 data format!!!)"<<endl;
    COUT<<"aborting...."<<endl;
#endif
}

template <class FieldType>

void  Field<FieldType>::saveSliceHDF5(string filename, string dataset_name,int xcoord, int thickness)
{
//write a slice which is in the y-z (done to insure same amount of data from each proc)

#ifdef HDF5

    Lattice slat;
    Field<FieldType> sfield;

    Site X(*lattice_);
    Site sX;

    int dim = this->lattice_->dim();

    int sSize[dim];
		int sSize_thin[dim];
    int r[dim];

    if(thickness>1)
    {

        sSize[0]=thickness;
        for(int i=1;i<dim;i++)sSize[i]=this->lattice_->size(i);
        slat.initialize(dim,sSize,0);
        sfield.initialize(slat,rows_,cols_,symmetry_);
        sfield.alloc();

        sX.initialize(slat);

        for(sX.first();sX.test();sX.next())
        {
            for(int l=0;l<dim;l++)r[l]=sX.coord(l);
            r[0]+=xcoord;

            X.setCoord(r);

            for(int i =0; i<components_;i++)sfield(sX,i)=this->operator()(X,i);
        }
    }
    else
    {
			for(int i=0;i<dim-1;i++)sSize_thin[i]=this->lattice_->size(i+1);
			slat.initialize(dim-1,sSize_thin,0);
			sfield.initialize(slat,rows_,cols_,symmetry_);
			sfield.alloc();

			sX.initialize(slat);

			for(sX.first();sX.test();sX.next())
			{
					for(int l=1;l<dim;l++)r[l]=sX.coord(l-1);
					r[0]=xcoord;
					X.setCoord(r);
					for(int i =0; i<components_;i++)sfield(sX,i)=this->operator()(X,i);
			}
    }

    sfield.saveHDF5(filename,dataset_name);
    sfield.dealloc();


#else
    COUT<<"LATfield2d must be compiled with HDF5 (flag HDF5 turn on!)"<<endl;
    COUT<<"to be able to use hdf5 data format!!!)"<<endl;
    COUT<<"aborting...."<<endl;
#endif

}

template <class FieldType>
void  Field<FieldType>::saveHDF5_coarseGrain3D(string filename, string dataset_name, int ratio)
{
#ifdef HDF5
    Lattice slat;
    Field<FieldType> sfield;

    int dim =lattice_->dim();
    long localsize[dim];

    int sSize[dim];
    int slocalsize[dim];

    long blocksize = array_size_*components_;
    long halo = lattice_->halo();

    int number_cg = ratio*ratio*ratio;
    long index_cg[number_cg];

    long index;
    long sindex;




    for(int i=0; i<dim;i++)
    {
        if(lattice_->sizeLocal(i)%ratio != 0 )
        {
            cout<<"process "<<parallel.rank()<<" have wrong ratio aborting coarse grain write"<<endl;
            return;
        }
        sSize[i]=(lattice_->size(i))/ratio;
        localsize[i]=lattice_->sizeLocal(i);
        slocalsize[i]=lattice_->sizeLocal(i)/ratio;
    }

    for(int k=0;k<ratio;k++)
    {
        for(int j=0;j<ratio;j++)
        {
            for(int i=0;i<ratio;i++)
            {
                index_cg[i+ratio*(j+ratio*k)]=((long)i+(localsize[0]+2l*halo)*((long)j + (localsize[1]+2l*halo)*(long)k))*blocksize;
            }
        }
    }
    slat.initialize(dim,sSize,0);
    sfield.initialize(slat,rows_,cols_,symmetry_);
    sfield.alloc();


    index= halo + (localsize[0]+2l*halo)*(halo+(localsize[1]+2l*halo)*halo);
    sindex=0;

    for(int k=0;k<slocalsize[2];k++)
    {
        for(int j=0;j<slocalsize[1];j++)
        {
            for(int i=0;i<slocalsize[0];i++)
            {
                for(int i_block=0;i_block<blocksize;i_block++)sfield.data()[sindex*blocksize+i_block] = data_[index*blocksize + i_block];
								for(int s=1;s<number_cg;s++)
                {
                    for(int i_block=0;i_block<blocksize;i_block++)sfield.data()[sindex*blocksize+i_block] += data_[index*blocksize + index_cg[s] + i_block];
                }
                for(int i_block=0;i_block<blocksize;i_block++)sfield.data()[sindex*blocksize+i_block]/=number_cg;

                index+=ratio;
                sindex+=1;
            }
            index+=(localsize[0]+(2l*halo))*ratio-localsize[0];
        }
        index += ( (localsize[1]+(2*halo) )*ratio-localsize[1] ) * (localsize[0]+(2*halo) );
    }

    sfield.saveHDF5(filename,dataset_name);

    sfield.dealloc();

#else
    COUT<<"LATfield2d must be compiled with HDF5 (flag HDF5 turn on!!)"<<endl;
    COUT<<"to be able to use hdf5 data format!!!)"<<endl;
    COUT<<"aborting.... "<<endl;
#endif
}


#ifdef EXTERNAL_IO
template <class FieldType>
void Field<FieldType>:: saveHDF5_server_open(string filename_base,int offset,int thickness)
{
    Lattice lat;
    int size[lattice_->dim()];
    for(int i=0;i<lattice_->dim();i++)size[i]=lattice_->size(i);
    iof_offset_ = offset;
    if(iof_offset_!=-1 && thickness == -1)   iof_thickness_ = 1;
    else iof_thickness_ = thickness;
    if(lattice_->dim()==2)iof_thickness_=-1;
    if(iof_thickness_!=-1)size[0]=iof_thickness_;
    lat.initialize(lattice_->dim(),size,lattice_->halo());
    io_file_ = ioserver.openFile(filename_base.c_str(),STRUCTURED_H5_FILE,type_id_,type_id_,&lat,components_,array_size_);
}
template <class FieldType>
void Field<FieldType>:: saveHDF5_server_write(int number_of_message, string filename_base, int offset,int thickness)
{

    if(!(io_file_.is_open))this->saveHDF5_server_open(filename_base,offset,thickness);



    if(!(io_file_.is_open)){COUT<<"Field<FieldType>::saveHDF5_server_write cannot open file: "<<filename_base<<endl;return;}


    int nom;
    int npt_nom;
    hsize_t size[lattice_->dim()];
    hsize_t offsetM[lattice_->dim()];
    int coord[lattice_->dim()];
    Site x(*lattice_);
    FieldType * buffer;
    long bufSize;


    if(lattice_->sizeLocal(lattice_->dim()-1) % number_of_message   == 0)
    {
        nom = number_of_message;
    }
    else
    {
        //cout<<"need to compute nom"<<endl;

        int nom_min,nom_max,nom_temp;



        bool flag=true;
        nom_temp = number_of_message;
        while(flag)
        {
            nom_temp -= 1;
            if(lattice_->sizeLocal(lattice_->dim()-1) % nom_temp == 0)
            {
                nom_min=nom_temp;
                flag=false;
            }
        }

        flag=true;
        nom_temp = number_of_message;
        while(flag)
        {
            if(nom_temp > lattice_->sizeLocal(lattice_->dim()-1))
            {
                nom_max=-1;
                flag=false;
            }
            nom_temp += 1;
            if(lattice_->sizeLocal(lattice_->dim()-1) % nom_temp == 0)
            {
                nom_max=nom_temp;
                flag=false;
            }
        }
        //cout<< lattice_->sizeLocal(lattice_->dim()-1) <<" , "<<  nom_min<<" , "<<nom_max<<endl;


        if(abs(number_of_message - nom_max) <=  abs(number_of_message - nom_min))nom = nom_max;
        else nom =nom_min;

        //cout<< lattice_->sizeLocal(lattice_->dim()-1) <<" , "<<  nom_min<<" , "<<nom_max<<" , "<< nom<<endl
    }
    npt_nom = lattice_->sizeLocal(lattice_->dim()-1) / nom;


    for(int i=0;i<lattice_->dim();i++)coord[i]=0;

    coord[lattice_->dim()-1]=lattice_->coordSkip()[0];
    coord[lattice_->dim()-2]=lattice_->coordSkip()[1];

    for(int i=0;i<lattice_->dim();i++)
    {
        offsetM[i]=coord[i];
        size[i] = lattice_->sizeLocal(i);
    }
    size[lattice_->dim()-1]=npt_nom;

    if(iof_thickness_ != -1)
    {
        size[0] = iof_thickness_;
        coord[0] = iof_offset_;
    }

    if(iof_thickness_ != -1 || lattice_->halo()!=0)
    {
        bufSize = 1;
        for(int i=0;i<lattice_->dim();i++)bufSize *= size[i];
        buffer = new FieldType[bufSize*components_];
    }

    for(int i=0;i<nom;i++)
    {
        if(!x.setCoord(coord))cout<<"internal error: ioserver send field out of local space..."<<endl;

        if(lattice_->halo()==0 && iof_thickness_ == -1)
        {
            //copy directly from data_
            ioserver.sendData(io_file_,(char*)(data_ + x.index()*components_),size,offsetM);

        }
        else
        {
            //have to make a buffer to send;

            if(iof_thickness_ == -1)
            {
                for(int i=0;i<bufSize;x.next(),i++)
                {
                    for(int c=0;c<components_;c++)buffer[i*components_ + c] = data_[x.index()*components_ + c];
                }
            }
            else
            {
                for(int i=0;i<bufSize;x.nextInSlice(iof_offset_,iof_thickness_),i++)
                {
                    for(int c=0;c<components_;c++)buffer[i*components_ + c] = data_[x.index()*components_ + c];
                }
            }
            ioserver.sendData(io_file_,(char*)buffer,size,offsetM);

        }

        coord[lattice_->dim()-1] += npt_nom;
        offsetM[lattice_->dim()-1] += npt_nom;
    }


    ioserver.closeFile(io_file_);
}
#endif



//MISCELLANEOUS=================

template <class FieldType>
Lattice& Field<FieldType>::lattice() { return *lattice_; }

template <class FieldType>
int Field<FieldType>::components() { return components_; }

template <class FieldType>
int Field<FieldType>::rows() { return rows_; }

template <class FieldType>
int Field<FieldType>::cols() { return cols_; }

template <class FieldType>
int Field<FieldType>::nMatrix() { return nMatrix_; }

template <class FieldType>
int Field<FieldType>::symmetry() { return symmetry_; }

template <class FieldType>
FieldType*& Field<FieldType>::data() { return data_; }

#ifdef HDF5
template <class FieldType>
hid_t Field<FieldType>::type_id(){return type_id_;}
template <class FieldType>
int Field<FieldType>::array_size(){return array_size_;}
#endif

//DEFAULT ASCII FIELD I/O FUNCTIONS======================
template<class FieldType>
void defaultFieldSave(fstream& file, FieldType* siteData, int components)
{
	file.precision(7);
	file.width(file.precision()+7);
	file.setf(fstream::scientific | fstream::showpos);

	file<<siteData[0];
	for(int i=1; i<components; i++)
    {
		file<<" "<<siteData[i];
    }
	file<<endl;
}

template<class FieldType>
void defaultFieldLoad(fstream& file, FieldType* siteData, int components)
{
	for(int i=0; i<components; i++) { file>>siteData[i]; }
}


#endif
