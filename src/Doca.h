#ifndef Doca_H
#define Doca_H

#include "geometry/Ray.h"

/** @class Doca

  @brief does some calculations about the DOCA between two Rays 
  
  A utility for calculating the distance of closest approach between two
  (assumed straight) tracks and the point on each line where this occurs.
    
  The method is taken from "Distance between Lines and Segments with their
  Closest Point of Approach" found at 
  http://www.geometryalgorithms.com/Archive/algorithm_0106/algorithm_0106.htm
            
 @author Tracy Usher 03/05/02

 $Header$
          
*/
class Doca
{
public:
    Doca(const Ray& ray1, const Ray& ray2);
    ~Doca() {}
    
    /// return doca between two Rays 
    double docaRay1Ray2()   {return doca;}
    /// return distance of doca from origin of 1st Ray
    double arcLenRay1()     {return s;}
    /// return distance of doca from origin of 2nd Ray
    double arcLenRay2()     {return t;}
    // I think the next 2 methods are broken...
    //    anyway, no-one calls them.
    //double docaPoint1Ray2();
    //double docaRay1Point2();
    /// point of doca on 1st Ray
    Point  docaPointRay1();
    /// point of doca on 2nd Ray
    Point  docaPointRay2();
    
private:
    /// origin of 1st Ray
    Point  P;
    /// direction of 1st Ray
    Vector u;
    /// origin of 2nd Ray
    Point  Q;
    /// direction of 2nd Ray
    Vector v;
    /// DOCA between the two Rays
    double doca;
    /// distance of DOCA point from origin 1
    double s;
    /// distance of DOCA point from origin 2
    double t;
};

#endif