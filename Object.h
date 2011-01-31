#ifndef CSE168_OBJECT_H_INCLUDED
#define CSE168_OBJECT_H_INCLUDED

#include <vector>
#include "Miro.h"
#include "Material.h"

typedef struct tex_coord2d_s
{
	tex_coord2d_s() : u(0), v(0) {}
	tex_coord2d_s(float _u, float _v) : u(_u), v(_v) {}
	float u, v;
} tex_coord2d_t;

typedef struct tex_coord3d_s
{
	tex_coord3d_s() : u(0), v(0), w(0) {}
	tex_coord3d_s(float _u, float _v, float _w) : u(_u), v(_v), w(_w) {}
	float u, v, w;
} tex_coord3d_t;

class Object
{
public:
    Object() {}
    virtual ~Object() {}

    void setMaterial(const Material* m) {m_material = m;}

    virtual void renderGL() {}
    virtual void preCalc() {}

	//Returns the (u,v) coordinates corresponding to a (x,y,z) coordinate.
	virtual tex_coord2d_t toUVCoordinates(const Vector3 & xyz) const { return tex_coord2d_t(xyz.x, xyz.z); }

    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX) = 0;

protected:
    const Material* m_material;
};

typedef std::vector<Object*> Objects;

#endif // CSE168_OBJECT_H_INCLUDED
