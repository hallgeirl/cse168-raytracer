#ifndef CSE168_PHONG_H_INCLUDED
#define CSE168_PHONG_H_INCLUDED

#include "Material.h"

class Phong : public Material
{
public:
	Phong(const Vector3 & kd = Vector3(1), 
			const Vector3 & ka = Vector3(0),
			const Vector3 & ks = Vector3(1),
			const float a = 0.f);
	virtual ~Phong();

	const Vector3 & kd() const {return m_kd;}
	const Vector3 & ka() const {return m_ka;}
	const Vector3 & ks() const {return m_ks;}
	const float a() const {return m_a;}

	void setKd(const Vector3 & kd) {m_kd = kd;}
	void setKa(const Vector3 & ka) {m_ka = ka;}
	void setKs(const Vector3 & ks) {m_ks = ks;}
	void setA(const float a) {m_a = a;}

    virtual void preCalc() {}
    
    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;

private:
	Vector3 m_kd;
	Vector3 m_ka;
	Vector3 m_ks;
	float m_a;
};

#endif		// CSE168_PHONG_H_INCLUDED
