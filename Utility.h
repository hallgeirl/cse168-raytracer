#ifndef __UTILITY__H_
#define __UTILITY__H_

#include <cstdlib>
#include <cmath>
#include "Vector3.h"

double getTime();
void getEigenVector(const float (&A)[3][3], float (&outV)[3], float lambda);

//Returns random number between 0 and 1
inline float frand()
{
    return (float)rand() / (float)RAND_MAX;
}

inline float sigmoid(float x)
{
    return 1/(1+exp(-x));
}

//Get two tangents to the normal vector that is supplied. Store the results in t1 and t2.
inline void getTangents(const Vector3& normal, Vector3& t1, Vector3& t2)
{
    //determine one surface tangent from cross of v and one of the unit vectors 
    t1 = cross(Vector3(0,0,1), normal);
    if (t1.length2() < 1e-6) t1 = cross(Vector3(0, 1, 0), normal);
    t2 = cross(t1, normal);
}

//Aligns a hemisphere to the vector v, and returns the vector formed by the spherical coordinates theta and phi.
inline Vector3 alignHemisphereToVector(const Vector3& v, float theta, float phi) 
{
    //convert spherical coords to cartesian coords
    float u1 = sin(phi) * cos(theta);
    float u2 = sin(phi) * sin(theta);
    float u3 = cos(phi);

    //determine one surface tangent from cross of v and one of the unit vectors 
    Vector3 t1 = cross(Vector3(0,0,1), v);
    if (t1.length2() < 1e-6) t1 = cross(Vector3(0, 1, 0), v);

    //Aligned direction determined by combination of cartesian coordiantes along tangents
    Vector3 aligned_d(u1*t1 + u2*cross(t1, v) + u3*v);

    aligned_d.normalize();
    return aligned_d;
}

#endif
