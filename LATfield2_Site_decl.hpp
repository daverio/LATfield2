#ifndef LATFIELD2_SITE_DECL_HPP
#define LATFIELD2_SITE_DECL_HPP



/*! \class Site
 \brief A class for referencing values of an instance of the field class at a given point on a lattice.

 A class which simplify the map of the field data array index. This class allow to get coordinate on the lattice, loop over each site of the lattice and perform displacment on the lattice.

 The site class encapsulates the mapping between the coordinate on the lattice and the index of the Field::data_ array which store the value of an instance of the Field class. It also contain method to loop over each site of the lattice, and to perform spacial displacement on the lattice.
 */
class Site
{
public:
		//CONSTRUCTORS=================

  //! Constructor.
		Site();

  /*!
   Constructor with initialization.

   \param lattice : the lattice on which the Site is defined.

   \sa initialize(Lattice& lattice)

   */
		Site(Lattice& lattice);

  /*!
   Constructor with initialization.

   \param lattice : the lattice on which the Site is defined.
   \param index   : set the current index of the field.

   \sa initialize(Lattice& lattice, long index)

   */
		Site(Lattice& lattice, long index);

		//INITIALIZATION=================
		/*!
     Initialization.
     \param lattice : the lattice on which the Site is defined.
     */
		void initialize(Lattice& lattice);

  /*!
   Constructor with initialization.
   \param lattice : the lattice on which the Site is defined.
   \param index   : set the current index of the field.
   */
		void initialize(Lattice& lattice, long index);

		//LOOPING OPERATIONS==============
  /*!
   Method to set the Site to the first site which is not within the halo. This method is used for loopping over the all lattice sites:

   for(site.first();site.test();site.next());

   \sa test()
   \sa next()
   */
		void first();

  /*!
   Method to test if the Site have a smaller or equal index than the last index not within the halo. This method is used for loopping over the all lattice sites:

   for(site.first();site.test();site.next());

   \sa first()
   \sa next()
   */
		bool test();

  /*!
   Method to jump to the next index which is not in the halo. This method is used for loopping over the all lattice sites:

   for(site.first();site.test();site.next());

   \sa first()
   \sa test()
   */
		void next();

  void nextInSlice(int offset,int thickness);

		//HALO OPERATIONS==============

  /*!
   Method to set the Site to the first site which is within the halo. This method is used for loopping over the all halo sites:

   for(site.haloFirst();site.haloTest();site.haloNext());

   \sa haloTest()
   \sa haloNext()
   */
		void haloFirst();

  /*!
   Method to test if the Site have a smaller or equal index than the last index within the halo. This method is used for loopping over the all halo sites:

   for(site.haloFirst();site.haloTest();site.haloNext());

   \sa haloFirst()
   \sa haloNext()
   */
		bool haloTest();

  /*!
   Method to jump to the next index which is in the halo. This method is used for loopping over the all halo sites:

   for(site.haloFirst();site.haloTest();site.haloNext());

   \sa haloFirst()
   \sa haloTest()
   */
		void haloNext();

		//NEIGHBOURING SITE OPERATORS==

  /*!
   Overloaded operator +
   The + operator is used to make a displacement of +1 site the the asked direction.
   \param direction : direction of the displacement
   */
		Site operator+(int direction);
  /*!
   Overloaded operator -
   The - operator is used to make a displacement of -1 site the the asked direction.
   \param direction : direction of the displacement
   */
		Site operator-(int direction);

		Site move(int direction);
    Site move(int direction, int step);
		Site move(int * steps);
    Site move3d(int sx, int sy, int sz);

		//SITE INDEX ADVANCE===========
  /*!
   Method which add "number" to the current index.
   */
		void indexAdvance(long number);

		//MISCELLANEOUS================
  /*!
   \return this method return the current index pointed by the site.
   */
		long index() const;

  /*!
   Method to set the current index of the site.
   \param new_index: the site index is set to new_index.
   */
		void setIndex(long new_index);

  /*!
   Method which return the site coordinate of a give dimension
   \param direction : label of the coordinate.
   \return site coordinate of the "direction" dimension
   */
		int coord(int direction);
  /*!
   Method which return the local site coordinate of a give dimension
   \param direction : label of the coordinate.
   \return site local coordinate of the "direction" dimension
   */
		int coordLocal(int direction);
  /*!
   Method to set the site to a given coordinate.
   \param r : array which contain the coordinate. The array size must be equal to the number of dimension of the lattice
   \return True: if the coordinate is local.
   False: if the local part of the lattice does not have this coordinate.
   */
		bool setCoord(int* r);
  /*!
   Method to set the site to a given coordinate for 3d lattices.
   \param x : coordinate of the 0 dimension.
   \param y : coordinate of the 1 dimension.
   \param z : coordinate of the 2 dimension.
   \return True: if the coordinate is local.
   False: if the local part of the lattice does not have this coordinate.
   */
		bool setCoord(int x, int y, int z);

  bool setCoordLocal(int *r);

  /*!
   \return Returns the pointer to the lattice on which the site is defined.
   */
		Lattice& lattice();

protected:
		Lattice* lattice_;
		long index_;
};



#ifdef FFT3D

/*! \class cKSite
 \brief A child of Site, built to work with the Fourier space Lattice for complex to complex transforms.

 A class which simplify the map of the field data array index. This class allow to get coordinate on the lattice, loop over each site of the lattice and perform displacment on the lattice.

 WARNING: this site class must be used only on lattice initialized using initializeComplexFFT() method of the Lattice class.

 This class have same binding that the Site class, so one can refer to the Site class for the documentation.

 */
class cKSite:public Site{
public:

  cKSite(){;}
  cKSite(Lattice& lattice){initialize(lattice);}
  cKSite(Lattice& lattice, long index){initialize(lattice, index);}

  void initialize(Lattice& lattice);
  void initialize(Lattice& lattice, long index);


  cKSite operator+(int asked_direction);
  cKSite operator-(int asked_direction);


  int coordLocal(int asked_direction);
  int coord(int asked_direction) ;
  int latCoord(int direction);
  int latCoordLocal(int direction);

  bool setCoord(int* r_asked);
  bool setCoord(int x, int y, int z);

private:
  int directions_[3];
};


/*! \class rKSite
 \brief A child of Site, built to work with the Fourier space Lattice for real to complex transforms.

 A class which simplifies the map of the field data array index. This class allow to get coordinate on the Lattice, loop over each site of the Lattice and access neighboring lattices sites

 WARNING: the rKSite class must be used only on Lattice initialized using initializeRealFFT() method of the Lattice class.

 This class has same binding as the Site class, please refer to the Site class for the documentation.

 */
class rKSite:public Site{
public:

  rKSite(){;}
  rKSite(Lattice& lattice){initialize(lattice);}
  rKSite(Lattice& lattice, long index){initialize(lattice, index);}

  void initialize(Lattice& lattice);
  void initialize(Lattice& lattice, long index);


  rKSite operator+(int asked_direction);
  rKSite operator-(int asked_direction);


  int coordLocal(int asked_direction);
  int coord(int asked_direction) ;
  int latCoord(int direction);
  int latCoordLocal(int direction);

  bool setCoord(int* r_asked);
  bool setCoord(int x, int y, int z);

private:
  int directions_[3];
};

#endif


#endif
