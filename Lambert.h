#ifndef CSE168_LAMBERT_H_INCLUDED
#define CSE168_LAMBERT_H_INCLUDED

#include <cmath>
#include "Object.h"
#include "Material.h"

class Lambert : public Material
{
public:
    Lambert(const Vector3 & kd = Vector3(1),
            const Vector3 & ka = Vector3(0));
    virtual ~Lambert();

	//For lambert we just return the color we have set in the constructor (or with setKd).
	//For other materials, we might return different values based on the texture coordinates given.
    virtual Vector3 diffuse2D(const tex_coord2d_t & texcoords) const {return m_kd;}
    virtual Vector3 diffuse3D(const tex_coord3d_t & texcoords) const { return m_kd; }
    virtual Vector3 ka(const tex_coord2d_t & position) const {return m_ka;}

    void setKd(const Vector3 & kd) {m_kd = kd;}
    void setKa(const Vector3 & ka) {m_ka = ka;}

    virtual void preCalc() {}

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;
private:
    Vector3 m_kd;
    Vector3 m_ka;
};

#endif // CSE168_LAMBERT_H_INCLUDED
