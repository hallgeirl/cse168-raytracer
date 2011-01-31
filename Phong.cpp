#include <iostream>
#include "Phong.h"
#include "Ray.h"
#include "Scene.h"

using namespace std;

Phong::Phong(const Vector3 &kd, const Vector3 &ka, const Vector3 &ks,
				const float shinyness, const float reflect, const float refract, const float refractIndex)
	:Lambert(kd, ka, reflect, refract, refractIndex), m_ks(ks), m_a(shinyness)
{

}

Phong::~Phong()
{

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

		//should specular component use material specular color?
		L += result * (std::max(0.0f, nDotL/falloff * pLight->wattage() / (4 * PI)) * diffuse + (pow(std::max(0.0f, eDotr/falloff * pLight->wattage() / (4 * PI)), m_a)) * m_ks);
    }

    // add the ambient component
    L += ka(hit.object->toUVCoordinates(hit.P));

    return L;
}
