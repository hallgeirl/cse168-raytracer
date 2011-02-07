#include <iostream>
#include "Phong.h"
#include "Ray.h"
#include "Scene.h"

using namespace std;

Phong::Phong(const Vector3 &kd, const Vector3 &ka, const Vector3 &ks, const Vector3 &kt,
				const float shinyness, const float refractIndex)
	: m_kd(kd), m_ka(ka), m_ks(ks), m_kt(kt), m_a(shinyness), m_refractIndex(refractIndex)
{
	//Keep the energy equation balanced (that is, don't refract+reflect+absorb more than 100% of the ray)
	m_ks.x = std::min(m_ks.x, 1.0f-m_kt.x);
	m_ks.y = std::min(m_ks.y, 1.0f-m_kt.y);
	m_ks.z = std::min(m_ks.z, 1.0f-m_kt.z);
}

Phong::~Phong()
{

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
		Ray Shadow(hit.P+(l*epsilon), l);
		HitInfo hitInfo;
		if (scene.trace(hitInfo, Shadow))
		{
			if (hitInfo.t < sqrt(falloff))
				continue;
		}

        // get the diffuse component
        float nDotL = dot(hit.N, l);

		// get the specular component
        float eDotr = dot(e, r);

		Vector3 result = pLight->color();

		//Look up the diffuse color
		Vector3 diffuse;
		if (GetLookupCoordinates() == UV)
			diffuse = diffuse2D(hit.object->toUVCoordinates(hit.P));
		else
			diffuse = diffuse3D(tex_coord3d_t(hit.P.x, hit.P.y, hit.P.z));

		if (m_ks.x > 0 && eDotr > 0.9)
			int count = 0;

		L += result * (std::max(0.0f, nDotL/falloff * pLight->wattage() / (4 * PI)) * diffuse + (pow(std::max(0.0f, eDotr/falloff * pLight->wattage() / (4 * PI)), m_a))*m_ks);
    }

    // add the ambient component
    L += ka(hit.object->toUVCoordinates(hit.P));

    return L;
}
