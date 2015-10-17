//   file: render.cpp
//
//   Program to render a 3D scene using Monte Carlo ray tracing algorithm.
//   
//   Programmer:   Matthew Acosta acosta.62@osu.edu
//
//   Revision history:
//	03/21/2015 Original
//	   
//
//
//   Notes:
//	* Method implements Monte Carlo ray tracing to render a scene
//	* Code base is based on <http://www.kevinbeason.com/
//	  smallpt/#moreinfo> 
//********************************************************************

//***************************Include Files****************************

#include <cstdio>
#include <cstdlib>

#include "raytrace_func.h"


//***************************Main Program*****************************

int main()
{
   int w = 1024, h = 768;   //image size
   int samps = 4;   //number of samples to take   

   render(w, h, samps);

return 0;
}

