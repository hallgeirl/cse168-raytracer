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
