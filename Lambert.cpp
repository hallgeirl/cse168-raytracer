#include "Lambert.h"
#include "Ray.h"
#include "Scene.h"

Lambert::Lambert(const Vector3 & kd) :
    m_kd(kd)
{
}

Lambert::~Lambert()
{
}

Vector3
Lambert::shade(const Ray& ray, const HitInfo& hit, const Scene& scene) const
{
    Vector3 L = Vector3(0.0f, 0.0f, 0.0f);

    const Vector3 viewDir = -ray.d; // d is a unit vector

    const Lights *lightlist = scene.lights();

    // loop over all of the lights
    Lights::const_iterator lightIter;
    for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
    {
        PointLight* pLight = *lightIter;

        Vector3 l = pLight->position() - hit.P;

        // the inverse-squared falloff
        float falloff = l.length2();

        // normalize the light direction
        l /= sqrt(falloff);

		Ray Shadow(hit.P+(l*epsilon), l);
		HitInfo hitInfo;
		if (scene.trace(hitInfo, Shadow, 0.f, sqrt(falloff)))
		{
			continue;
		}

        // get the diffuse component
        float nDotL = dot(hit.N, l);
        Vector3 result = pLight->color();
        result *= m_kd;

        L += std::max(0.0f, nDotL/falloff * pLight->wattage() / PI) * result;
    }

    return L;
}
