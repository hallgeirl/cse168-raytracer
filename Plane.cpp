#include <iostream>
#include "Plane.h"
#include "Ray.h"

using namespace std;

Plane::Plane()
{
}


Plane::~Plane()
{

}


void
Plane::renderGL()
{
/*    glBegin(GL_QUADS);
        glVertex3f(-
        , v0.y, v0.z);
        glVertex3f(v1.x, v1.y, v1.z);
        glVertex3f(v2.x, v2.y, v2.z);
    glEnd();*/
}


bool
Plane::intersect(HitInfo& result, const Ray& r,float tMin, float tMax)
{
    float ndotd = dot(normal, r.d);
    if (fabs(ndotd) < 1e-6) return false;
    
    float t = dot(normal, (origin-r.o)) / ndotd;
    if (t < tMin || t > tMax)
        return false;

    result.P = r.o + t*r.d;
    result.t = t;
    result.N = normal;
    result.material = m_material;

    return true;
}
