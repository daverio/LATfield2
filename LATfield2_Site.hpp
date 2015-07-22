#ifndef LATFIELD2_SITE_HPP
#define LATFIELD2_SITE_HPP
/*! \file LATfield2_Site.hpp
 \brief Site class definition
 
 LATfield2_Site.hpp contain the class Site definition.
 
 */ 


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
         Mehtod to set the current index of the site.
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
        /*!
         \return Returns the pointer to the lattice on which the site is defined.
         */
		Lattice& lattice();
		
	protected:
		Lattice* lattice_;
		long index_;
	};

//CONSTRUCTORS===================

Site::Site() {;}
Site::Site(Lattice& lattice) { initialize(lattice); }
Site::Site(Lattice& lattice, long index) { initialize(lattice, index); }

//INITIALIZATION=================

void Site::initialize(Lattice& lattice) { lattice_=&lattice; index_=0l;}
void Site::initialize(Lattice& lattice, long index) { lattice_ = &lattice; index_ = index; }

//NEIGHBOURING SITE OPERATORS==

Site Site::operator+(int direction)
{
	return Site( *lattice_, index_ + lattice_->jump(direction) );
}

Site Site::operator-(int direction)
{
	return Site( *lattice_, index_ - lattice_->jump(direction) );
}

//LOOPING OPERATIONS====================

void Site::first() { index_=lattice_->siteFirst(); }

bool Site::test() { return index_ <= lattice_->siteLast(); }

void Site::next()
{
	index_++;
	//If coordLocal(0) != sizeLocal(0) then next site reached
	if( coordLocal(0) != lattice_->sizeLocal(0) ) { return; }
	else
	{
		index_ -= lattice_->sizeLocal(0);
		for(int i=1; i<lattice_->dim(); i++)
		{
			index_ += lattice_->jump(i);
			//If coordLocal(i) !=sizeLocal(0) then next site reached
			if( coordLocal(i) != lattice_->sizeLocal(i) ) { return; }
			index_ -= lattice_->sizeLocal(i) * lattice_->jump(i);
		}
		index_ = lattice_->siteLast() + 1;
	}
}

//HALO OPERATIONS====================

void Site::haloFirst() { index_ =  0; }
bool Site::haloTest() { return index_ < lattice_->sitesLocalGross(); }

void Site::haloNext()
{
	index_++;
	
	//Can only leave boundary by reaching coord(0)=0, so otherwise done
	if(coordLocal(0)==0)
	{
		bool is_in_halo = 0;
		for(int i=1; i<lattice_->dim(); i++)
		{
			if( coordLocal(i)<0 || coordLocal(i)>=lattice_->sizeLocal(i) ) { is_in_halo = 1; break; }
		}
		if(!is_in_halo) index_ += lattice_->sizeLocal(0);
	}
}

//INDEX ADVANCE======================

void Site::indexAdvance(long number) { index_ += number; }

//MISCELLANEOUS======================

long Site::index() const { return index_; }

void Site::setIndex(long new_index) {index_ = new_index;}

int Site::coord(int direction) ////////sensible a quelle dim est scatter (seul modif a faire ici)
{
	if(direction<lattice_->dim()-2) { return coordLocal(direction); }
	else if (direction==lattice_->dim()-2) {return coordLocal(direction)+lattice_->coordSkip()[1]; }
	else {return coordLocal(direction)+lattice_->coordSkip()[0]; }
}

int Site::coordLocal(int direction)
{
	if(direction==lattice_->dim()-1)
	{
		return index_/lattice_->jump(direction) - lattice_->halo();
	}
	else if(direction==0)
	{
		return index_%lattice_->jump(1) - lattice_->halo();
	}
	else
	{
		return (index_%lattice_->jump(direction+1)) / lattice_->jump(direction) - lattice_->halo();
	}
}

bool Site::setCoord(int* r)
{
	
	this->first();
	//Check site is local
	if(r[lattice_->dim()-1]<this->coord(lattice_->dim()-1) || r[lattice_->dim()-1]>=this->coord(lattice_->dim()-1)+lattice_->sizeLocal(lattice_->dim()-1)
	   || r[lattice_->dim()-2]<this->coord(lattice_->dim()-2) || r[lattice_->dim()-2]>=this->coord(lattice_->dim()-2)+lattice_->sizeLocal(lattice_->dim()-2) )
	{
		return false;
		//COUT<<"LATfield::Site::setCoord(int*) - Site desired non-local!"<<endl;
	} 
	else
	{
		
		int jump=0;
		for(int i=0; i<lattice_->dim(); i++)
		{
			jump+=(r[i]-coord(i))*lattice_->jump(i);
		} 
		
		this->indexAdvance(jump);
		return true;
	}
}    

bool Site::setCoord(int x, int y=0, int z=0)
{
	int* r = new int[3];
	r[0]=x;
	r[1]=y;
	r[2]=z;
	return this->setCoord(r);
	delete[] r;
}    
Lattice& Site::lattice() { return *lattice_ ; }


#ifdef FFT3D

/*! \class cKSite  
 \brief A child of Site, built to work with the Fourier space lattices for complex to complex transforms.
 
 A class which simplify the map of the field data array index. This class allow to get coordinate on the lattice, loop over each site of the lattice and perform displacment on the lattice. 
 
 WARNING: this site class must be used only on lattices initialized using initializeComplexFFT() method of the Lattice class.
 
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

void cKSite::initialize(Lattice& lattice) { lattice_=&lattice; directions_[0]=1; directions_[1]=2; directions_[2]=0;}
void cKSite::initialize(Lattice& lattice, long index) { lattice_ = &lattice; index_ = index; directions_[0]=1; directions_[1]=2; directions_[2]=0; }

cKSite cKSite::operator+(int asked_direction)
{
    return cKSite( *lattice_, index_ + lattice_->jump(directions_[asked_direction]));
}

cKSite cKSite::operator-(int asked_direction)
{
    return cKSite( *lattice_, index_ - lattice_->jump(directions_[asked_direction]));
}
int cKSite::coord(int asked_direction) ////////sensible a quelle dim est scatter (seul modif a faire ici)
{
	int direction= directions_[asked_direction] ;
    
	if(direction<lattice_->dim()-2) { return latCoordLocal(direction); }
	else if (direction==lattice_->dim()-2) {return latCoordLocal(direction)+lattice_->coordSkip()[1]; }
	else {return latCoordLocal(direction)+lattice_->coordSkip()[0]; }
}

int cKSite::latCoord(int direction) ////////sensible a quelle dim est scatter (seul modif a faire ici)
{
	
	if(direction<lattice_->dim()-2) { return coordLocal(direction); }
	else if (direction==lattice_->dim()-2) {return coordLocal(direction)+lattice_->coordSkip()[1]; }
	else {return coordLocal(direction)+lattice_->coordSkip()[0]; }
}

int cKSite::coordLocal(int asked_direction)
{
	int direction= directions_[asked_direction] ;
	
	if(direction==lattice_->dim()-1)
	{
		return index_/lattice_->jump(direction) - lattice_->halo();
	}
	else if(direction==0)
	{
		return index_%lattice_->jump(1) - lattice_->halo();
	}
	else
	{
		return (index_%lattice_->jump(direction+1)) / lattice_->jump(direction) - lattice_->halo();
	}
}
int cKSite::latCoordLocal(int direction)
{
    
	if(direction==lattice_->dim()-1)
	{
		return index_/lattice_->jump(direction) - lattice_->halo();
	}
	else if(direction==0)
	{
		return index_%lattice_->jump(1) - lattice_->halo();
	}
	else
	{
		return (index_%lattice_->jump(direction+1)) / lattice_->jump(direction) - lattice_->halo();
	}
}

bool cKSite::setCoord(int* r_asked)
{
    int r[3];
	r[0]=r_asked[2];
	r[1]=r_asked[0];
	r[2]=r_asked[1];
	this->first();
	//Check site is local
	if(r[lattice_->dim()-1]<this->latCoord(lattice_->dim()-1) || r[lattice_->dim()-1]>=this->latCoord(lattice_->dim()-1)+lattice_->sizeLocal(lattice_->dim()-1)
	   || r[lattice_->dim()-2]<this->latCoord(lattice_->dim()-2) || r[lattice_->dim()-2]>=this->latCoord(lattice_->dim()-2)+lattice_->sizeLocal(lattice_->dim()-2) )
	{
		return false;
		//COUT<<"LATfield::Site::setCoord(int*) - Site desired non-local!"<<endl;
	} 
	else
	{
		
		int jump=0;
		for(int i=0; i<lattice_->dim(); i++)
		{
			jump+=(r[i]-latCoord(i))*lattice_->jump(i);
		} 
		
		this->indexAdvance(jump);
		return true;
	}
}
bool cKSite::setCoord(int x, int y=0, int z=0)
{
	int* r = new int[3];
	r[0]=x;
	r[1]=y;
	r[2]=z;
	return this->setCoord(r);
	delete[] r;
}   

/*! \class rKSite  
 \brief A child of Site, built to work with the Fourier space lattices for real to complex transforms.
 
 A class which simplifies the map of the field data array index. This class allow to get coordinate on the lattice, loop over each site of the lattice and access neighboring lattices sites
 
 WARNING: the rKSite class must be used only on lattices initialized using initializeRealFFT() method of the Lattice class.
 
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

void rKSite::initialize(Lattice& lattice) { lattice_=&lattice; directions_[0]=0; directions_[1]=2; directions_[2]=1;}
void rKSite::initialize(Lattice& lattice, long index) { lattice_ = &lattice; index_ = index; directions_[0]=0; directions_[1]=2; directions_[2]=1; }

rKSite rKSite::operator+(int asked_direction)
{
    return rKSite( *lattice_, index_ + lattice_->jump(directions_[asked_direction]));
}

rKSite rKSite::operator-(int asked_direction)
{
    return rKSite( *lattice_, index_ - lattice_->jump(directions_[asked_direction]));
}
int rKSite::coord(int asked_direction) ////////sensible a quelle dim est scatter (seul modif a faire ici)
{
	int direction= directions_[asked_direction] ;
    
	if(direction<lattice_->dim()-2) { return latCoordLocal(direction); }
	else if (direction==lattice_->dim()-2) {return latCoordLocal(direction)+lattice_->coordSkip()[1]; }
	else {return latCoordLocal(direction)+lattice_->coordSkip()[0]; }
}

int rKSite::latCoord(int direction) ////////sensible a quelle dim est scatter (seul modif a faire ici)
{
	
	if(direction<lattice_->dim()-2) { return coordLocal(direction); }
	else if (direction==lattice_->dim()-2) {return coordLocal(direction)+lattice_->coordSkip()[1]; }
	else {return coordLocal(direction)+lattice_->coordSkip()[0]; }
}

int rKSite::coordLocal(int asked_direction)
{
	int direction= directions_[asked_direction] ;
	
	if(direction==lattice_->dim()-1)
	{
		return index_/lattice_->jump(direction) - lattice_->halo();
	}
	else if(direction==0)
	{
		return index_%lattice_->jump(1) - lattice_->halo();
	}
	else
	{
		return (index_%lattice_->jump(direction+1)) / lattice_->jump(direction) - lattice_->halo();
	}
}
int rKSite::latCoordLocal(int direction)
{
    
	if(direction==lattice_->dim()-1)
	{
		return index_/lattice_->jump(direction) - lattice_->halo();
	}
	else if(direction==0)
	{
		return index_%lattice_->jump(1) - lattice_->halo();
	}
	else
	{
		return (index_%lattice_->jump(direction+1)) / lattice_->jump(direction) - lattice_->halo();
	}
}

bool rKSite::setCoord(int* r_asked)
{
    int r[3];
	r[0]=r_asked[0];
	r[1]=r_asked[2];
	r[2]=r_asked[1];
	this->first();
	//Check site is local
	if(r[lattice_->dim()-1]<this->latCoord(lattice_->dim()-1) || r[lattice_->dim()-1]>=this->latCoord(lattice_->dim()-1)+lattice_->sizeLocal(lattice_->dim()-1)
	   || r[lattice_->dim()-2]<this->latCoord(lattice_->dim()-2) || r[lattice_->dim()-2]>=this->latCoord(lattice_->dim()-2)+lattice_->sizeLocal(lattice_->dim()-2) )
	{
		return false;
		//COUT<<"LATfield::Site::setCoord(int*) - Site desired non-local!"<<endl;
	} 
	else
	{
		
		int jump=0;
		for(int i=0; i<lattice_->dim(); i++)
		{
			jump+=(r[i]-latCoord(i))*lattice_->jump(i);
		} 
		
		this->indexAdvance(jump);
		return true;
	}
}  
bool rKSite::setCoord(int x, int y=0, int z=0)
{
	int* r = new int[3];
	r[0]=x;
	r[1]=y;
	r[2]=z;
	return this->setCoord(r);
	delete[] r;
} 


#endif


#endif


