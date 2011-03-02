#include <iostream>
#include "Phong.h"
#include "Ray.h"
#include "Scene.h"

#ifdef STATS
#include "Stats.h"
#endif

using namespace std;

Phong::Phong(const Vector3 &kd, const Vector3 &ks, const Vector3 &kt,
				const float shininess, const float refractIndex)
{
	//Keep the energy equation balanced (that is, don't refract+reflect+absorb more than 100% of the ray)
    m_diffuse = kd;
    m_specular = ks; 
    m_transmission = kt;
    m_shininess = shininess;
    m_refractIndex = refractIndex;

	//refraction is between 0.f and 1.f - reflection 
	m_transmission.x = std::max(std::min(m_transmission.x, 1.0f-m_specular.x), 0.f);
	m_transmission.y = std::max(std::min(m_transmission.y, 1.0f-m_specular.y), 0.f);
	m_transmission.z = std::max(std::min(m_transmission.z, 1.0f-m_specular.z), 0.f);

	//absorption is between 0.f and 1.f - reflection - refraction
	m_diffuse.x = std::max(1.0f-m_specular.x-m_transmission.x, 0.f);
	m_diffuse.y = std::max(1.0f-m_specular.y-m_transmission.y, 0.f);
	m_diffuse.z = std::max(1.0f-m_specular.z-m_transmission.z, 0.f);
}

Phong::~Phong()
{

}

bool Phong::isDiffuse() const
{
	return (m_diffuse.x > 0.f || m_diffuse.y > 0.f || m_diffuse.z > 0.f);
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
#ifdef STATS 
		Stats::Shadow_Rays++;
#endif
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

		//removed m_specular from specular highlight 
		//specular highlight should be dependent on shinyness rather than reflective component
		L += result * (std::max(0.0f, nDotL/falloff * pLight->wattage() / (4 * PI)) * diffuseColor * m_diffuse);
        //Old specular highlights component. /*+ (pow(std::max(0.0f, eDotr/falloff * pLight->wattage() / (4 * PI)), m_a))*/
    }

    return L;
}
