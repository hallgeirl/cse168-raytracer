#ifndef CSE168_OBJECT_H_INCLUDED
#define CSE168_OBJECT_H_INCLUDED

#include <vector>
#include "Miro.h"
#include "Material.h"

typedef struct tex_coord_s
{
	tex_coord_s() : u(0), v(0) {}
	tex_coord_s(float _u, float _v) : u(_u), v(_v) {}
	float u, v;
} tex_coord_t;

class Object
{
public:
    Object() {}
    virtual ~Object() {}

    void setMaterial(const Material* m) {m_material = m;}

    virtual void renderGL() {}
    virtual void preCalc() {}

	//Returns the (u,v) coordinates corresponding to a (x,y,z) coordinate.
	virtual tex_coord_t toTextureCoordinates(const Vector3 & xyz) const { return tex_coord_t(xyz.x, xyz.z); }

    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX) = 0;

protected:
    const Material* m_material;
};

typedef std::vector<Object*> Objects;

#endif // CSE168_OBJECT_H_INCLUDED
