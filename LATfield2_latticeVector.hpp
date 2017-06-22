#ifndef LATFIELD2_LATTICEVECTOR_HPP
#define LATFIELD2_LATTICEVECTOR_HPP


/*! \file LATfield2_Site.hpp
 \brief LATfield2_Site.hpp contains the Site, rKSite, and cKSite definition.
 \author David Daveio, Neil Bevis, with modifications by Wessel Valkenburg

 */

#include "LATfield2_latticeVector_decl.hpp"

latticeVector::latticeVector()
{
	lattice_ = NULL;
}
latticeVector::latticeVector(Lattice& lattice)
{
	lattice_ = &lattice;
}
latticeVector::jump()
{
	return lattice_->jump(direction_);
}

inline latticeVector& latticeVector::operator=(const int dir)
{
	direction_ = dir;
	return *this;
}
inline latticeVector& latticeVector::operator++()
{
	direction_++;
	return *this;
}
inline latticeVector& latticeVector::operator--()
{
	direction_--;
	return *this;
}


inline bool operator< (const latticeVector& lhs, const int rhs)
{
	return lhs.direction_ < rhs;
}

inline bool operator> (const latticeVector& lhs, const int rhs)
{
	return lhs.direction_ > rhs;
}

inline long operator+ (const long& i, const latticeVector&  v)
{
	return i + v.jump();
}

inline long operator- (const long& i, const latticeVector&  v)
{
	return i - v.jump();
}

inline long operator+ (const Site& x, const latticeVector&  v)
{
	return x.index() + v.jump();
}

inline long operator- (const Site x, const latticeVector&  v)
{
	return x.index() - v.jump();
}

#endif
