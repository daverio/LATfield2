#ifndef LATFIELD2_LATTICEVECTOR_DECL_HPP
#define LATFIELD2_LATTICEVECTOR_DECL_HPP


#include "LATfield2_Lattice_decl.hpp"

class latticeVector
{
public:

		latticeVector();
		latticeVector(Lattice& lattice);
		latticeVector& operator=(const int dir);
		latticeVector& operator++();
		latticeVector& operator--();
		bool operator<(const int lim);
		bool operator>(const int lim);
		long jump();

		int direction_;

protected:
		Lattice* lattice_;

};




#endif
