#ifndef LOOPING_MACRO_HPP
#define LOOPING_MACRO_HPP

/*! \file looping_macro.hpp
 \brief looping over the lattice site macro
 \author David Daverio

 */

#define forallboundary_start(latttt,xxxx) \
Site xxxx(latttt); \
for(xxxx.haloFirst();xxxx.haloTest();xxxx.haloNext()) \
{ \
bool inBoundary = false; \
for(int i=0; i<latttt.dim(); i++)if(xxxx.coord(i)>=latttt.size(i) || xxxx.coord(i)<0)inBoundary=true; \
if(inBoundary==true) \
{



#define forallboundary_stop }};


#endif
