#ifndef INT2STRING_HPP
#define INT2STRING_HPP

/* \file int2string.hpp
   \brief a small function to convert integer to string
 
 */

//function: int2string
//
//Generates a string object containing an integer padded with zeros.

#include <string>
using std::string;

string int2string(int number, int max = 999, bool zeropad = true)
{
  string output;
  char c;
  int i;

  //Get number of digits needed for max
  int digits=1;
  for(i=10; (i-1)<max; i*=10) digits++;

  //Create string of length digits (padded with zeros)
  for(i=0; i<digits; i++)
    {
      c = '0' + number%10;
      if(c!='0' || zeropad)
	{
	  output = c + output;
	}
      number /= 10;
    }

  return output;
}

#endif
