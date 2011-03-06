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

//Returns a random direction in the hemisphere that is oriented in the direction specified
inline Vector3 sampleHemisphereDirection(const Vector3& hemisphereOrientation)
{
    //bias to the surface normal
    float x, y, z;
    do
    { 
        x = 2*frand() - 1;
        y = 2*frand() - 1;
        z = 2*frand() - 1;
    } while (x*x + y*y + z*z > 1.0f && dot(Vector3(x, y, z), hemisphereOrientation) < 0);

    return Vector3(x,y,z);
}

inline VectorR2 sampleDisc(float radius)
{
	float x_rand, y_rand;
	Vector3 new_eye;
	do {
		x_rand = (2*frand() - 1) * radius;
		y_rand = (2*frand() - 1) * radius;
	} while (x_rand*x_rand + y_rand*y_rand > radius*radius);

    VectorR2 v;
    v.x = x_rand; v.y = y_rand;

    return v;
}

#endif
