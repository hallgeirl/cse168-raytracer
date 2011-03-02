#ifndef CSE168_MATERIAL_H_INCLUDED
#define CSE168_MATERIAL_H_INCLUDED

#include "Miro.h"
#include "Vector3.h"

//Determines how the shader should query the material for its color.
//UV is a regular 2D lookup using u,v coordinates (e.g. 2D textures),
//UVW is a 3D lookup using 3 coordinates (e.g. 3D textures)
//For UV, the object's texture mapping function is called to get the UV coordinates.
enum LookupCoordinates
{
	UV = 0,
	UVW = 1
};

class Material
{
public:
    Material();
    virtual ~Material();

	//Default to 2D lookup
	virtual LookupCoordinates GetLookupCoordinates() const { return UV; }

    virtual float   bumpHeight2D(const tex_coord2d_t & texture_coords) const { return 0; }
    virtual float   bumpHeight3D(const tex_coord3d_t & texture_coords) const { return 0; }

	bool         isReflective() const { return (m_specular.x > 0.f || m_specular.y > 0.f || m_specular.z > 0.f); }
	bool         isRefractive() const { return (m_transmission.x > 0.f || m_transmission.y > 0.f || m_transmission.z > 0.f);}
	virtual bool isDiffuse() const    { return true; }

	Vector3         getReflection() const       {return m_specular;}
	Vector3         getRefraction() const       {return m_transmission;}
	virtual Vector3 getDiffuse() const          {return Vector3(1.f);}
	float           getRefractionIndex() const  {return m_refractIndex;}

	void setReflection(const Vector3 reflect) { m_specular = reflect; }
	void setRefraction(const Vector3 refract, const float refractIndex) { m_refractIndex = refractIndex; m_transmission = refract; }
	void setRefractionIndex(const float refractIndex) {m_refractIndex = refractIndex;}

	void  setShininess(const float shininess) {m_shininess = shininess;}
    float getShininess() const { return m_shininess; }	

    virtual void preCalc() {}

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit, const Scene& scene) const;

protected:
	Vector3 m_specular;     //Reflection
	Vector3 m_transmission; //Refraction
	float m_refractIndex;
    float m_gloss;
	float m_shininess;
};

#endif // CSE168_MATERIAL_H_INCLUDED
