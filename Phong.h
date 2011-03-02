#ifndef CSE168_PHONG_H_INCLUDED
#define CSE168_PHONG_H_INCLUDED

#include "Object.h"
#include "Material.h"

class Phong : public Material
{
public:
	Phong(const Vector3 & diffuseColor = Vector3(1),       
			const Vector3 & specularColor = Vector3(0),    //Reflectivity for each color. 1 is fully reflective, 0 is fully non-reflective.
			const Vector3 & transparentColor = Vector3(0), //Transparency for each color. 1 is fully transparent (refracting according to refractIndex), 0 is fully opaque.
			const float shinyness = 1.f,
			const float refractIndex = 1);
	virtual ~Phong();

	virtual bool isDiffuse() const;
	
	//For other materials, we might return different values based on the texture coordinates given.
    virtual Vector3 diffuse2D(const tex_coord2d_t & texcoords) const { return m_diffuse; }
    virtual Vector3 diffuse3D(const tex_coord3d_t & texcoords) const { return m_diffuse; }
	virtual Vector3 getDiffuse() const {return m_diffuse;}

	void setDiffuse(const Vector3 & kd) {m_diffuse = kd;}

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;

protected:
	Vector3 m_diffuse;
};

#endif		// CSE168_PHONG_H_INCLUDED
