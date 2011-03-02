#ifndef CSE168_Plane_H_INCLUDED
#define CSE168_Plane_H_INCLUDED

#include <limits>
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

    //Object boundaries used with bounding box creation
    virtual Vector3 coordsMin() const { return -Vector3(infinity); }
    virtual Vector3 coordsMax() const { return Vector3(infinity); }
    virtual Vector3 center() const { return this->m_origin; } 

    void setNormal(Vector3 normal) { this->m_normal = normal; }
    void setOrigin(Vector3 origin) { this->m_origin = origin; }
	virtual tex_coord2d_t toUVCoordinates(const Vector3 & xyz) const;

    virtual void renderGL();
 
    virtual bool isBounded() const { return false; }
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX);

protected:
    Vector3 m_normal, m_origin;
};

#endif // CSE168_Plane_H_INCLUDED
