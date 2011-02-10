#ifndef CSE168_PHONG_H_INCLUDED
#define CSE168_PHONG_H_INCLUDED

#include "Object.h"
#include "Material.h"

class Phong : public Material
{
public:
	Phong(const Vector3 & diffuseColor = Vector3(1),       
			const Vector3 & ambientColor = Vector3(0),     //Magnitude of ambient lighting for this material. ambient*diffuse=ambient contribution.
			const Vector3 & specularColor = Vector3(0),    //Reflectivity for each color. 1 is fully reflective, 0 is fully non-reflective.
			const Vector3 & transparentColor = Vector3(0), //Transparency for each color. 1 is fully transparent (refracting according to refractIndex), 0 is fully opaque.
			const float shinyness = 1.f,
			const float refractIndex = 1);
	virtual ~Phong();

	
	virtual bool IsReflective() const;
	virtual bool IsRefractive() const;

	//For other materials, we might return different values based on the texture coordinates given.
    virtual Vector3 diffuse2D(const tex_coord2d_t & texcoords) const {return m_kd;}
    virtual Vector3 diffuse3D(const tex_coord3d_t & texcoords) const { return m_kd; }
    virtual Vector3 ka(const tex_coord2d_t & position) const {return m_ka;}
	const float a() const {return m_a;}
	virtual Vector3 GetReflection() const {return m_ks;}
	virtual Vector3 GetRefraction() const {return m_kt;}
	virtual Vector3 GetAbsorption() const {return Vector3(1.f) - m_ks - m_kt;}
	virtual float GetRefractionIndex() const {return m_refractIndex;}

	void setKd(const Vector3 & kd) {m_kd = kd;}
    void setKa(const Vector3 & ka) {m_ka = ka;}
	void setKs(const Vector3 & ks) {m_ks = ks;}
	void setKt(const Vector3 & kt) {m_kt = kt;}
	void setA(const float a) {m_a = a;}
	void SetRefractionIndex(const float refractIndex) {m_refractIndex = refractIndex;}

    virtual void preCalc() {}

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;

protected:
    Vector3 m_kd;
    Vector3 m_ka;
	Vector3 m_ks;
	Vector3 m_kt;
	Vector3 m_absorb;
	float m_a;
	float m_refractIndex;
};

#endif		// CSE168_PHONG_H_INCLUDED
