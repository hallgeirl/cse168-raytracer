#include "Material.h"

Material::Material()
{
    m_emittance = Vector3(0);
}

Material::~Material()
{
}

Vector3
Material::shade(const Ray&, const HitInfo&, const Scene&) const
{
    return Vector3(1.0f, 1.0f, 1.0f);
}
