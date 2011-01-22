#include "Phong.h"
#include "Ray.h"
#include "Scene.h"

Phong::Phong(const Vector3 &kd, const Vector3 &ka, const Vector3 &ks, const float a)
	:m_kd(kd), m_ka(ka), m_ks(ks), m_a(0.f)
{

}

Phong::~Phong()
{

}

Vector3 
Phong::shade(const Ray &ray, const HitInfo &hit, const Scene &scene) const
{
	/*Vector3 L = Vector3(0.0f, 0.0f, 0.0f);

	Vector3 viewDir = -ray.d;

	const Lights *lightlist = scene.lights();

	// loop over all of the lights
    Lights::const_iterator lightIter;
    for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
    {
        PointLight* pLight = *lightIter;
    
        Vector3 surfaceToLight = pLight->position() - hit.P;
        
        // the inverse-squared falloff
        float falloff = surfaceToLight.length2();
        
        // normalize the light direction
        surfaceToLight /= sqrt(falloff);

		Vector3 lightToSurface = -surfaceToLight;

		Vector3 reflectRay = lightToSurface - 2 * cross(lightToSurface, hit.N) * hit.N;
	}*/

	return Vector3();
}