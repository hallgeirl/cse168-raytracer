#ifndef CSE168_SPHERE_H_INCLUDED
#define CSE168_SPHERE_H_INCLUDED

#include "Vector3.h"
#include "Object.h"

class Sphere : public Object
{
public:
    Sphere();
    virtual ~Sphere();

    void setCenter(const Vector3& v)    {m_center = v;}
    void setRadius(const float f)       {m_radius = f;}

    float radius() const                {return m_radius;}
    
    //Object boundaries used with bounding box creation
    virtual Vector3 coordsMin() const { return m_center - m_radius; }
    virtual Vector3 coordsMax() const { return m_center + m_radius; }
    virtual Vector3 center() const { return m_center; }
	
	virtual float getArea(const Vector3& lightPos) {return PI*m_radius*m_radius;}
	virtual Vector3 samplePosition(const Vector3& lightPos) const;
    virtual void renderGL();

    //For spherical texture mapping
    virtual tex_coord2d_t toUVCoordinates(const Vector3 & xyz) const;

    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX);

protected:
    Vector3 m_center;
    float m_radius;
};

#endif // CSE168_SPHERE_H_INCLUDED
