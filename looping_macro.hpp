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
for(int iii=0; iii<latttt.dim(); iii++)if(xxxx.coord(iii)>=latttt.size(iii) || xxxx.coord(iii)<0)inBoundary=true; \
if(inBoundary==true) \
{



#define forallboundary_stop }};


#endif
