#include "Material.h"

Material::Material()
{
    setReflection(Vector3(0));
    setRefraction(Vector3(0), 1);
    setShininess(infinity); //Default to a perfectly reflective surface if the specular component is non-zero
}

Material::~Material()
{
}

Vector3
Material::shade(const Ray&, const HitInfo&, const Scene&) const
{
    return Vector3(1.0f, 1.0f, 1.0f);
}
