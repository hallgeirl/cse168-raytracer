#include "Plane.h"
#include "Ray.h"


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
    float t = dot(normal, (origin-r.o)) / (dot(normal, r.d));
    if (t < 0)
        return false;

    result.P = r.o + t*r.d;
    result.t = t;
    result.N = normal;
    result.material = m_material;

    return true;
}
