#ifndef CSE168_Plane_H_INCLUDED
#define CSE168_Plane_H_INCLUDED

#include "Object.h"

/*
    The Plane class stores a pointer to a mesh and an index into its
    Plane array. The mesh stores all data needed by this Plane.
*/
class Plane : public Object
{
public:
    Plane();
    virtual ~Plane();

    void setNormal(Vector3 normal) { this->normal = normal; }
    void setOrigin(Vector3 origin) { this->origin = origin; }

    virtual void renderGL();
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX);

protected:
    Vector3 normal, origin;
};

#endif // CSE168_Plane_H_INCLUDED
