#ifndef CSE168_PHONG_H_INCLUDED
#define CSE168_PHONG_H_INCLUDED

#include "Lambert.h"

class Phong : public Lambert
{
public:
	Phong(const Vector3 & kd = Vector3(1), 
			const Vector3 & ka = Vector3(0),
			const Vector3 & ks = Vector3(1),
			const float a = 0.f,
			const float m_reflect = 0,
			const float m_refract = 0,
			const float m_refractIndex = 1);
	virtual ~Phong();

	const Vector3 & ks() const {return m_ks;}
	const float a() const {return m_a;}	

	void setKs(const Vector3 & ks) {m_ks = ks;}
	void setA(const float a) {m_a = a;}

    virtual void preCalc() {}
    
    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;

protected:
	Vector3 m_ks;
	float m_a;
};

#endif		// CSE168_PHONG_H_INCLUDED
