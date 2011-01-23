#ifndef CSE168_MATERIAL_H_INCLUDED
#define CSE168_MATERIAL_H_INCLUDED

#include "Miro.h"
#include "Vector3.h"

class Material
{
public:
    Material();
    virtual ~Material();

	virtual float GetReflection() const {return 0;}
	virtual float GetRefraction() const {return 0;}
	virtual float GetRefractionIndex() const {return 1.0f;}

	virtual void SetReflection(const float reflect){};
	virtual void SetRefraction(const float refract, const float refractIndex){};

    virtual void preCalc() {}
    
    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;
};

#endif // CSE168_MATERIAL_H_INCLUDED
