#include <iostream>
#include "Phong.h"
#include "Ray.h"
#include "Scene.h"

using namespace std;

Phong::Phong(const Vector3 &kd, const Vector3 &ks, const Vector3 &kt,
				const float shinyness, const float refractIndex)
	: m_kd(kd), m_ks(ks), m_kt(kt), m_a(shinyness), m_refractIndex(refractIndex)
{
	//Keep the energy equation balanced (that is, don't refract+reflect+absorb more than 100% of the ray)

	//refraction is between 0.f and 1.f - reflection 
	m_kt.x = std::max(std::min(m_kt.x, 1.0f-m_ks.x), 0.f);
	m_kt.y = std::max(std::min(m_kt.y, 1.0f-m_ks.y), 0.f);
	m_kt.z = std::max(std::min(m_kt.z, 1.0f-m_ks.z), 0.f);

	//absorption is between 0.f and 1.f - reflection - refraction
	m_diffuse.x = std::max(1.0f-m_ks.x-m_kt.x, 0.f);
	m_diffuse.y = std::max(1.0f-m_ks.y-m_kt.y, 0.f);
	m_diffuse.z = std::max(1.0f-m_ks.z-m_kt.z, 0.f);
}

Phong::~Phong()
{

}

bool Phong::IsDiffuse() const
{
	return (m_diffuse.x > 0.f || m_diffuse.y > 0.f || m_diffuse.z > 0.f);
}

bool Phong::IsReflective() const
{
	return (m_ks.x > 0.f || m_ks.y > 0.f || m_ks.z > 0.f);
}

bool Phong::IsRefractive() const
{
	return (m_kt.x > 0.f || m_kt.y > 0.f || m_kt.z > 0.f);
}

Vector3
Phong::shade(const Ray &ray, const HitInfo &hit, const Scene &scene) const
{
	Vector3 L = Vector3(0.0f, 0.0f, 0.0f);

	Vector3 e = -ray.d;

		//Look up the diffuse color
	Vector3 diffuseColor;
	if (GetLookupCoordinates() == UV)
		diffuseColor = diffuse2D(hit.object->toUVCoordinates(hit.P));
	else
		diffuseColor = diffuse3D(tex_coord3d_t(hit.P.x, hit.P.y, hit.P.z));


	const Lights *lightlist = scene.lights();

	// loop over all of the lights
    Lights::const_iterator lightIter;
    for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
    {
        PointLight* pLight = *lightIter;

        Vector3 l = pLight->position() - hit.P;

		Vector3 r = - l + 2 * dot(l, hit.N) * hit.N;
		r.normalize();

        // the inverse-squared falloff
        float falloff = l.length2();

        // normalize the light direction
        l /= sqrt(falloff);

		// No light contribution if Ray hits an object 
#ifndef DISABLE_SHADOWS
		Ray Shadow(hit.P+(l*epsilon), l);
		HitInfo hitInfo;
		if (scene.trace(hitInfo, Shadow, 0.f, sqrt(falloff)))
		{
			continue;
		}
#endif

        // get the diffuse component
        float nDotL = dot(hit.N, l);

		// get the specular component
        float eDotr = dot(e, r);

		Vector3 result = pLight->color();

		//removed m_ks from specular highlight 
		//specular highlight should be dependent on shinyness rather than reflective component
		L += result * (std::max(0.0f, nDotL/falloff * pLight->wattage() / (4 * PI)) * diffuseColor * m_diffuse + (pow(std::max(0.0f, eDotr/falloff * pLight->wattage() / (4 * PI)), m_a)));
    }

    // add the ambient component
    L += getEmittance();

    return L;
}
