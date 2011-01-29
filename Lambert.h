#ifndef CSE168_LAMBERT_H_INCLUDED
#define CSE168_LAMBERT_H_INCLUDED

#include "Material.h"

class Lambert : public Material
{
public:
    Lambert(const Vector3 & kd = Vector3(1),
            const Vector3 & ka = Vector3(0),
			const float m_reflect = 0,
			const float m_refract = 0,
			const float m_refractIndex = 1);
    virtual ~Lambert();

    virtual const Vector3 & kd(const Vector3 & position) const {return m_kd;} //For lambert we just return the color we have set in the constructor (or with setKd).
    virtual const Vector3 & ka(const Vector3 & position) const {return m_ka;}
	virtual float GetReflection() const {return m_reflect;}
	virtual float GetRefraction() const {return m_refract;}
	virtual float GetRefractionIndex() const {return m_refractIndex;}

    void setKd(const Vector3 & kd) {m_kd = kd;}
    void setKa(const Vector3 & ka) {m_ka = ka;}
	virtual void SetReflection(const float reflect) {m_reflect = reflect;};
	virtual void SetRefraction(const float refract, const float refractIndex) 
		{m_refract = refract; m_refractIndex = refractIndex;}

    virtual void preCalc() {}
    
    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;
private:
    Vector3 m_kd;
    Vector3 m_ka;

protected:
	float m_reflect;
	float m_refract;
	float m_refractIndex;
};

#endif // CSE168_LAMBERT_H_INCLUDED
