/** @file Doca.cxx
@brief calculates the DOCA between two straight lines, and related quantities
@author Tracy Usher, Leon Rochester

$Header$
*/

#include "Doca.h"

//
// Constructor is where all the real work gets done...
//


Doca::Doca(const Ray& ray1, const Ray& ray2)
{
    //Copy the input rays to the class
    P = ray1.position();
    u = ray1.direction().unit();
    Q = ray2.position();
    v = ray2.direction().unit();

    //Determine vector from start point track 1 to start of track 2
    Vector w     = P - Q;

    //Projections of of tracks along vector between start points
    double d     = u.dot(w);
    double e     = v.dot(w);

    //Dot product between two tracks to check if parallel
    double b     = u.dot(v);
    double denom = 1. - b*b;

    //Lines are not parallel
    if (fabs(b) < 1.)
    {
        s    = (b*e - d  ) / denom;
        t    = (e   - b*d) / denom;
        w    = w + s * u - t * v;
        doca = w.magnitude();
    }
    //Lines are parallel
    else
    {
        s    = 0;
        t    = d / b;
        w    = w - t * v;
        doca = w.magnitude();
    }
}


Point Doca::docaPointRay1()
{
    Point w = P + s * u;

    return w;
}

Point Doca::docaPointRay2()
{
    Point w = Q + t * v;

    return w;
}

/*
double Doca::docaPoint1Ray2()
{
    Vector w = P - Q;
    Vector c = v.cross(w);

    return c.magnitude();
}
double Doca::docaRay1Point2()
{
    Vector w = Q - P;
    Vector c = u.cross(w.unit());

    return c.magnitude();
}
*/