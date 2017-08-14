#ifndef LATFIELD2_SITE_HPP
#define LATFIELD2_SITE_HPP
/*! \file LATfield2_Site.hpp
 \brief LATfield2_Site.hpp contains the Site, rKSite, and cKSite definition.
 \author David Daveio, Neil Bevis, with modifications by Wessel Valkenburg

 */

#include "LATfield2_Site_decl.hpp"

//CONSTRUCTORS===================

Site::Site() {index_=0; lattice_ = NULL;}
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

Site Site::move(int direction, int step)
{
	return Site(*lattice_, index_ + step * lattice_->jump(direction));
}
Site Site::move(int direction)
{
	return Site( *lattice_, index_ + lattice_->jump(direction) );
}
Site Site::move(int * steps)
{
	double index = index_;
	for(int i=0;i<lattice_->dim();i++)index += steps[i]*lattice_->jump(i);
	return Site(*lattice_,index);
}
Site Site::move3d(int sx, int sy, int sz)
{
	return Site(*lattice_,index_ + (sx*lattice_->jump(0)) + (sy*lattice_->jump(2)) + (sz*lattice_->jump(2)) );
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

void Site::nextInSlice(int offset,int thickness)
{
    index_++;
    //If coordLocal(0) != sizeLocal(0) then next site reached
    if( coordLocal(0) != offset+thickness ) { return; }
    else
    {
        index_ -= thickness;
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
	int r[3];
	r[0]=x;
	r[1]=y;
	r[2]=z;
	return this->setCoord(r);
}
bool Site::setCoordLocal(int *r)
{
    this->first();
    long jump=0;
    for(int i=0; i<lattice_->dim(); i++)
    {
        jump+=r[i]*lattice_->jump(i);
    }

    this->indexAdvance(jump);
    return true;
}
Lattice& Site::lattice() { return *lattice_ ; }


ostream& operator<<(ostream& os,  Site& x)
{
    os << " (";
    int i=0;
    for(;i < (x.lattice().dim()-1); )
    {
	os << x.coord(i) << " , ";
	i++;
    }
    os << x.coord(i) << ") ";
    return os;
}




#ifdef FFT3D

/* ckSite implmentation */

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
	int r[3];
	r[0]=x;
	r[1]=y;
	r[2]=z;
	return this->setCoord(r);
}

/* rkSite implmentation */
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
	int r[3];
	r[0]=x;
	r[1]=y;
	r[2]=z;
	return this->setCoord(r);
}


#endif


#endif
