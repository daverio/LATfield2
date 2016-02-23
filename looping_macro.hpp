#ifndef LOOPING_MACRO_HPP
#define LOOPING_MACRO_HPP

/*! \file looping_macro.hpp
 \brief looping over the lattice site macro
 \author David Daverio,Neil Bevis

 */

#define forallboundary_start(latttt,xxxx) \
Site xxxx(latttt); \
for(xxxx.haloFirst();xxxx.haloTest();xxxx.haloNext()) \
{ \
bool inBoundary = false; \
for(int i=0; i<lat.dim(); i++)if(x.coord(i)>=lat.size(i) || x.coord(i)<0)inBoundary=true; \
if(inBoundary==true) \
{



#define forallboundary_stop }};


#endif
