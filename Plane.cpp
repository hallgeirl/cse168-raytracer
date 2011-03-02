#include <iostream>
#include "Plane.h"
#include "Ray.h"

using namespace std;

Plane::Plane()
{
    m_normal = Vector3(0,1,0);
    m_origin = Vector3(0,0,0);
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
    float ndotd = dot(m_normal, r.d);
    if (fabs(ndotd) < 1e-6) return false;
    
    float t = dot(m_normal, (m_origin-r.o)) / ndotd;
    if (t < tMin || t > tMax)
        return false;

    result.P = r.o + t*r.d;
    result.t = t;
    result.N = m_normal;
    result.material = m_material;

    return true;
}

tex_coord2d_t Plane::toUVCoordinates(const Vector3 & xyz) const
{
    float inX = xyz.x, inY = xyz.z;
    // Wrap coordinates to [0,1]
/*    inX -= (int)inX; 
    inY -= (int)inY;
    if (inX < 0) inX++;
    if (inY < 0) inY++;*/
    
    return tex_coord2d_t(inX, inY);
}
