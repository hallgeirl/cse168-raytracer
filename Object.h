#ifndef CSE168_OBJECT_H_INCLUDED
#define CSE168_OBJECT_H_INCLUDED

#include <vector>
#include "Miro.h"
#include "Material.h"



class Object
{
public:
    Object() {}
    virtual ~Object() {}

    void setMaterial(const Material* m) {m_material = m;}
	const Material* getMaterial() const {return m_material;}

    virtual void renderGL() {}
    virtual void preCalc() {}

    //Object boundaries and center used with bounding box creation.
    virtual Vector3 coordsMin() const = 0;
    virtual Vector3 coordsMax() const = 0;
    virtual Vector3 center() const = 0;

    //Returns the area and a position sample as seen from some point in space (e.g. a light source)
	virtual float getArea(const Vector3& lightPos) {return 0;}
	virtual Vector3 samplePosition(const Vector3& asSeenFrom) const { return samplePosition(); }
	
	virtual Vector3 samplePosition() const {return center();}
    
    //Unbounded objects like planes should override this and return false.
    virtual bool isBounded() const { return true; }

	//Returns the (u,v) coordinates corresponding to a (x,y,z) coordinate.
	virtual tex_coord2d_t toUVCoordinates(const Vector3 & xyz) const { return tex_coord2d_t(xyz.x, xyz.z); }


    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX) = 0;

protected:
    const Material* m_material;
};

typedef std::vector<Object*> Objects;

#endif // CSE168_OBJECT_H_INCLUDED
