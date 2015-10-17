//   file: raytrace_func.h
//
//   Header file for render.cpp
//
//   Programmer:   Matthew Acosta   acosta.62@osu.edu
//
//   Revision history:
//	03/21/2015   Original Version Created
//
//   Notes:
//	* Must define my own vector class and operations since vector 
//	  math is not as straight forward as just using predefined 
//	  operators.
//	* Vec3 redefined from V1 to make more trasparent how Vec3 
//	  is defined
//	* Some vector methods have also been changed for clarity

class Vec3
{
public:
   //used double since double will be unnessary precision for coordinates
   double x, y, z;	
  
   Vec3 (double x_= 0, double y_= 0, double z_= 0) : x(x_), y(y_), z(z_) {}   
  
   // declare new vector operations
   Vec3 operator + (const Vec3&)const;
   Vec3 operator - (const Vec3&)const;
   Vec3 operator * (double)const;   //mult vector and scalar
   Vec3 mult(const Vec3&)const;   //mult two vectors
   double dot(const Vec3&)const;
   Vec3 cross(const Vec3&)const;
   /* & after Vec3 for "this" pointer */
   Vec3& normalize();
};

class Ray
{
public:
   //specify the ray as a parametic line with origin (o) and 
   //direction (d)
   //a point along the ray can then be specified as p(t) = o + t * d,
   //where t is our parameter
   Vec3 o, d;
   Ray(Vec3 o_, Vec3 d_) : o(o_), d(d_) {}
};

//enum of material types to be used in radiance function
//DIFFuse(most objects), SPECular(mirror like), REFRactive(glass like)
enum Refl_t { DIFF, SPEC, REFR};

class Sphere
{ 
public:
   double rad;   //radius
   Vec3 pos, emi, col;   //position, emision, color
   Refl_t refl;   //reflection type DIFFuse, SPECular, REFRactive
   
   Sphere(double rad_, Vec3 pos_, Vec3 emi_, Vec3 col_, Refl_t refl_) :
      rad(rad_), pos(pos_), emi(emi_), col(col_), refl(refl_) {}
   //function determines distance to object or 0 if miss
   double intersect(const Ray &) const;
};


//main render function
extern void render(int w, int h, int samps);


