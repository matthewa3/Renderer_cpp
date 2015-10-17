//   file: raytrace_func.cpp
//
//   Definitions for functions and Vec3 class operations
//
//   Programmer:   Matthew Acosta   acosta.62@osu.edu
//
//   Revision history:
//	03/21/2015   Original version created
//
//  Git Test
//
//   Notes: 
//	* Much of definitions have been redefined for clarity and ease
//	  of use
//	* To see if a ray hits the sphere we start with the equation of
//	  a sphere in vector form (p - c) * (p - c) - r ^ 2 = 0, where 
//	  p is a point along the ray, c is the center of a circle, 
//	  and r is the radius of this circle. With the equation for
//	  a point in parametric terms, p(t) = o + t * d (o is origin
//	  and d is direction), we plug this into sphere equation. We
//	  get the quadratic equation: 
//	  (d * d) * t ^ 2 + 2 * d * (o - c) * t + (o - c) ^ 2 - r ^ 2,
//	  which leads to the following solution:
//	  t = (- b +- sqrt(b ^ 2 - 4 * a * c))/(2 * a). 
//	  In the above a = (d * d), b = 2 * d * (o - c), 
//	  c = (o - c) ^ 2 - r ^ 2.
//	  If the discriminant is negative, the ray misses the sphere.
//	  If both solutions are negative then the sphere is behind the 
//	  ray.
//	* Gamma Correction: <https://en.wikipedia.org/wiki/
//	  		    Gamma_correction>
//	* Source code uses erand48 as to not screw up parallelization
//	  http://www-01.ibm.com/support/knowledgecenter/
//	       SSLTBW_2.1.0/com.ibm.zos.v2r1.bpxbd00/rernd4.htm
//	* PPM: https://en.wikipedia.org/wiki/Netpbm_format
//
//***********************Include Files**********************************

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <iostream>

#include "raytrace_func.h"	  //include header for class/functions

using namespace std;

//************************ Definitions *********************************

/******************** Beginning of Vec3 routines **********************/

Vec3 Vec3::operator+ (const Vec3& v)const 
{
   return Vec3(x + v.x, y + v.y, z + v.z);
}

Vec3 Vec3::operator- (const Vec3& v)const 
{
   return Vec3(x - v.x, y - v.y, z - v.z);
}

Vec3 Vec3::operator* (double v)const
{
   return Vec3(x * v, y * v, z * v);
}

Vec3 Vec3::mult(const Vec3& v)const
{
   return Vec3(x * v.x, y * v.y, z * v.z);
}

double Vec3::dot (const Vec3& v)const
{
   return x * v.x + y * v.y + z * v.z;
}

Vec3 Vec3::cross(const Vec3& v)const
{
   return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
} 

/* "this" is used to refer to Vec3*/
Vec3& Vec3::normalize ()
{
   return *this = *this * (1/sqrt(x * x + y * y + z * z));  
}
  

/*********************** End of Vec3 routines *************************/



/*********************** Begin Sphere routines ************************/

//See notes to look at how the routine below was determined
double Sphere::intersect(const Ray &r) const 
{
   //pos is the center of circle and op connects origin of ray to center
   //of circle
   Vec3 op = pos - r.o;
   double t, eps = 1e-4;   //eps is our tolerance for hitting 
   double b = op.dot(r.d);   //1/2 b from quad. eq.
   double det = b * b - op.dot(op) + rad * rad;   //(b^2-4ac):a=1 by norm
   
   if (det < 0) return 0;
   else det = sqrt(det);
   //return smallest positive t since this is the point the camera can see
   return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);  
}

/*********************** End Sphere routines **************************/



/******************** Begin Scene Description *************************/

//scene: radius, position, emission, color, material

/************************* Cornell Box *****************************/
   //consists of a reflective sphere and transparent sphere in a room with
   //light above
   Sphere spheres[] = 
   {
   //Left Wall
   Sphere(1e5, Vec3( 1e5+1,40.8,81.6), Vec3(),Vec3(.75,.25,.25),DIFF),
   //Right Wall 
   Sphere(1e5, Vec3(-1e5+99,40.8,81.6),Vec3(),Vec3(.25,.25,.75),DIFF),
   //Back Wall 
   Sphere(1e5, Vec3(50,40.8, 1e5),     Vec3(),Vec3(.75,.75,.75),DIFF),
   //Front Wall
   Sphere(1e5, Vec3(50,40.8,-1e5+170), Vec3(),Vec3(),           DIFF),
   //Floor
   Sphere(1e5, Vec3(50, 1e5, 81.6),    Vec3(),Vec3(.75,.75,.75),DIFF),
   //Ceiling
   Sphere(1e5, Vec3(50,-1e5+81.6,81.6),Vec3(),Vec3(.75,.75,.75),DIFF),
   //Mirror Ball
   Sphere(16.5,Vec3(27,16.5,47),       Vec3(),Vec3(1,1,1)*.999, SPEC),
   //Glass Ball
   Sphere(16.5,Vec3(73,16.5,78),       Vec3(),Vec3(1,1,1)*.999, REFR),
   //Light 
   Sphere(1.5, Vec3(50,81.6-16.5,81.6),Vec3(4,4,4)*100,  Vec3(), DIFF)
   }; 


int numSpheres = sizeof(spheres)/sizeof(Sphere); 
/******************** End Scene Description ***************************/



/******************** Begin Intersect Scene Function ******************/

//check to see if our ray intersects with the spheres in our scene,
//if ray does hit, then keep closest intersection since again this is
//the point our camera sees.

inline bool intersect(const Ray &r, double &t, int &id)
{
   double n = numSpheres;   //num of spheres
   double d;   //distance to potential intersection
   double inf = t = 1e20;   //cutoff for intersection points

   for (int i = int (n); i--;)
   {
      if((d = spheres[i].intersect(r)) && d < t)
      {
         t = d;    //distance to intersection
         id = i;   //id of intersected object
      }
   }
    
   return t < inf;
}

/****************** End Intersect Scene Function **********************/



/************************ Begin Radiance Function**********************/

Vec3 radiance(const Ray &r, int depth, unsigned short *Xi, int E = 1)
{
   /*intersection*/
   double t;	// distance to intersection
   int id = 0;	// id of intersected object
   if (!intersect(r, t, id)) return Vec3();	//if missed, return black
   const Sphere &obj = spheres[id];	//hit object
   if (depth > 10) return Vec3();   //too deep

   /*surface properties*/
   Vec3 x = r.o + r.d * t;	//ray intersection point
   Vec3 n = (x - obj.pos).normalize();	//sphere normal
   //determine if inside or outside glass object
   //dot product of normal and ray will determine this
   Vec3 nl = n.dot(r.d) < 0 ? n : n *-1;   //properly oriented surface norm
   Vec3 f = obj.col;  // object color

   /*Russian Roulette for reflectivity*/
   //stops recusion based on surface reflectivity
   //uses max component of RGB surface color and only activates after depth 5

   double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; //max refl
   if (++depth > 5 || !p)
   {
     if (erand48(Xi) < p) f = f * (1/p); 
     else return obj.emi*E;
   }

   /*Diffuse Reflection*/
   if (obj.refl == DIFF)
   {
      double r1 = 2 * M_PI * erand48(Xi);   //random angle 
      double r2 = erand48(Xi), r2s = sqrt(r2);	//random distance from center
      //use normal to create orthonormal coordinates w, u, v
      Vec3 w = nl;   // w = normal
      Vec3 u = ((fabs(w.x) > 0.1 ? Vec3(0,1) : Vec3 (1)).cross(w)).normalize();
      Vec3 v = w.cross(u);
      //d is random reflection ray
      Vec3 d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalize(); 
      
      /* Loop over any lights*/
      Vec3 e;
      for (int i = 0; i < numSpheres; i++)
      {
         const Sphere &s = spheres[i];
	 if (s.emi.x <= 0 && s.emi.y <= 0 && s.emi.z <= 0) continue; //skip non lights
         
	 // Create random direction towards sphere using method from realistic ray 
	 // tracing text
	 //create coordinate system for sampling by solid angle
	 Vec3 sw = s.pos-x; 
	 Vec3 su = ((fabs(sw.x) > 0.1 ? Vec3 (0,1) : Vec3(1)).cross(sw)).normalize();
	 Vec3 sv = sw.cross(su);

	 //determine max angle of sphere
	 double cos_a_max = sqrt(1 - s.rad * s.rad / (x - s.pos).dot(x - s.pos));

	 // calculate sample direction based on random numbers according to 
	 // equation from text Realistic Ray Tracing
	 double eps1 = erand48(Xi), eps2 = erand48(Xi);
	 double cos_a = 1 - eps1 + eps1 * cos_a_max;
	 double sin_a = sqrt(1 - cos_a * cos_a);
	 double phi = 2 * M_PI * eps2;
	 Vec3 l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
	 l.normalize();

	 /* Create shadow ray*/
	 
	 // check for occlusion with shadow ray or see if shadow ray is hidden
	 // occlusion is when something you want to see is hidden 
	 if (intersect(Ray(x,l), t, id) && id == i)
	 {
	    double omega = 2 * M_PI * (1 - cos_a_max);
	    // calculate lighting and add to current value
	    e = e + f.mult(s.emi * l.dot(nl) * omega) * M_1_PI;   //M_1_Pi = 1 / pi
	 }
      }
      // make recursive call with random ray direction computed
      // 0 parameter below turns off emmision terms for next recursion
      return obj.emi * E + e + f.mult(radiance(Ray(x,d),depth,Xi,0));
   }  
   /* Ideal Specular Reflection*/ 
   else if (obj.refl == SPEC)
   {
      return obj.emi + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)),depth,Xi));
	// Angle of incidence = Angle of reflection^^^^^^^^^^^^^^^^^^^
   }
   /* Dielectric Surface(Glass)*/
   // Glass is reflective and refractive so compute reflected ray
   Ray reflRay(x, r.d - n * 2 * n.dot(r.d));	//Ideal reflection
   bool into = n.dot(nl) > 0;	// Check to see if ray is outside to inside
   // Index of Refraction for air is 1 and for glass is 1.5
   // nnt is either 1.5 or 1/1.5
   double nc = 1, nt = 1.5, nnt = into ? nc/nt:nt/nc, ddn = r.d.dot(nl), cos2t;
   
   // If total internal reflection, Reflect
   // Occurs if light ray tries to leave at to shallow an angle
   if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)	//Total int refl
   {   
      return obj.emi + f.mult(radiance(reflRay,depth,Xi));
   }

   // Else, choose to reflect or refract
   
   // calculate refracted ray
   Vec3 tdir = (r.d * nnt - n * ((into ? 1 : -1)  
	       	* (ddn * nnt + sqrt(cos2t)))).normalize();

   double a = nt - nc;
   double b = nt + nc; 
   double R0 = a * a / (b * b);	 // reflectance at norm incidence based on IOR
   double c = 1 - (into ? -ddn:tdir.dot(n));   // 1-cos(theta)
   double Re = R0 + (1 - R0) * c * c * c * c * c;   //fresnel reflectance
   double Tr = 1 - Re;
   double P = 0.25 + 0.5 * Re;   //Prob of reflecting
   double RP = Re/P;
   double TP = Tr / (1 - P);

   // Russian roulette so that 1 recursive call randomly if depth > 2
   // and 2 calls if depth <= 2
   return obj.emi + f.mult(depth > 2 ? (erand48(Xi) < P ?
          radiance(reflRay, depth, Xi) * RP : radiance(Ray(x, tdir), depth, Xi) * TP) :
	  radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);

}

/********************** End Radiance Function *************************/



/***************** Begin Color Correction Functions *******************/

//functions needed to convert colors from radiance to something 
//displayable in the rbg range(0-255).
//utilizes a gamma correction of 2.2 to correct brightness(see notes) 

inline double clamp(double x) 
{
   return x < 0 ? 0 : x > 1 ? 1 : x;   //get values between 0 and 1
}

//convert doubles to integers for ppm file
inline int toInt(double x)
{
   return int(pow(clamp(x), 1. / 2.2) * 255 + 0.5);
}

/**************** End Color Correction Functions **********************/

/******************** Begin Render Function ***************************/
void render( int w_, int h_, int samps_)
{
   /*set up image*/

   int w = w_, h = h_;   //image size
   int samps = samps_ > 0 ? samps_ : 1;   //# of samples(default of 1)
   Vec3 r;   //used for colors of samples
   Vec3 *c = new Vec3[w * h];   //array stores image
   
   /*set up camera*/
 
   //Look from and gaze direction(camera pos, dir)
   /* not 100% sure how these numbers were determined*/
   Ray cam(Vec3(50, 52, 295.6), Vec3(0, -0.042612, -1).normalize());
 
   //horizontal(x) camera direction with field of view 
   //FoV angle of 30 deg(0.5135 rad). 
   //uses implicit 0 for y and z
   Vec3 camx = Vec3(w * 0.5135 / h);
  
   //vertical(y) camera direction
   //uses cross product to get vector perp. to both cx and gaze dir.
   Vec3 camy = camx.cross(cam.d).normalize() * 0.5135;

   //OpenMP for parallelizing code
   //each loop iteration in own thread 
   #pragma omp parallel for schedule(dynamic, 1) private(r)

   //Loop over all image pixels
   for (int y = 0; y < h; y++)   //loop over image rows
   {
     // print out the progress of image
     fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps,100.*y/(h-1));        unsigned short y3 = y * y * y;
     unsigned short Xi[3] = {0, 0, y3};   //seed for random num gen
     for (unsigned short x = 0; x < w; x++)   //loop over image columns
     {
       //For each pixel use 2x2 subpixels and average thier colors

       //2x2 subpixel rows
       for (int sy = 0; sy < 2; sy++)
       {
          int i = (h - y - 1) * w + x;  //array index for pixel

          //2x2 subixel columns
          for (int sx = 0; sx < 2; sx++)
          {
             r = Vec3();   //radiance

	     for (int s = 0; s < samps; s++)
             {
                //Create a tent filter to determine location
		//of sample within pixel
		double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1: 1 - sqrt(2 - r1);
		double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1: 1 - sqrt(2 - r2);

		//compute ray direction using cam.d, camx, and camy
		Vec3 d = camx * ( ( (sx + 0.5 + dx) / 2 + x)/w - 0.5) +
			 camy * ( ( (sy + 0.5 + dy) / 2 + y)/h - 0.5) + cam.d;
		
		//Use radiance function to calculate radiance
		r = r + radiance(Ray(cam.o + d * 140, d.normalize()), 0, Xi) * (1./samps);
             }
	     	    // Camera rays are pushed ^^^^^^ forward to start in interior
	     // gamma correction for subpixel color for pixel color 
   	     c[i] = c[i] + Vec3(clamp(r.x), clamp(r.y), clamp(r.z))*0.25;
          }
       }
     }  

   }
   // Write to a PPM file(see notes for wiki page on PPM)
   FILE *f = fopen("image.ppm", "w");   // write image to PPM file
   fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
   for (int i = 0; i < w * h; i++)
   {
      fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
   }    
}

/*********************** End Render Function **************************/




