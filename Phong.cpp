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
	Vector3 specular = Vector3(0.0f, 0.0f, 0.0f);
    Vector3 reflectResult;
	Vector3 refractResult;
    HitInfo hitInfo;
	Vector3 e = -ray.d;

	const Lights *lightlist = scene.lights();
    float refraction = std::min(GetRefraction(), 1.0f);
    float reflection = std::min(GetReflection(), 1.0f);
    
    //Keep the energy equation balanced (that is, don't refract+reflect+absorb more than 100% of the ray)
    reflection = std::min(reflection, 1.0f-refraction);
    float diffuse = 1-reflection-refraction;
   	Vector3 diffuseColor;
    
    if (diffuse > epsilon)
    {
        //Look up the diffuse color
	    if (GetLookupCoordinates() == UV)
    		diffuseColor = diffuse2D(hit.object->toUVCoordinates(hit.P));
	    else
		    diffuseColor = diffuse3D(tex_coord3d_t(hit.P.x, hit.P.y, hit.P.z));
    }
    //if reflective material, send trace with ReflectRay
    if (reflection > 0.0f)
    {
        Ray reflectRay = ray.Reflect(hit);

        if (!scene.traceScene(reflectRay, reflectResult, hit.depth))
        {
            reflection = 0;
            reflectResult.set(0);
        }
    }


    //if refractive material, send trace with RefractRay
    if (refraction > 0.0f)
    {
        Ray	refractRay = ray.Refract(hit);
        if (!scene.traceScene(refractRay, refractResult, hit.depth))
        {
            refraction = 0;
            refractResult.set(0);
        }
    }

    
	
    // loop over all of the lights
    Lights::const_iterator lightIter;
    for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
    {
        PointLight* pLight = *lightIter;
        Vector3 l = pLight->position() - hit.P;

        //Cast shadow ray
        {
            Vector3 d = l;
            d.normalize();

            Ray r(hit.P + d*epsilon, d);
            HitInfo h;
            if (scene.trace(h, r, 0.0f, l.length()))
                continue;
        }


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

		//should specular component use material specular color?
		L += result * (std::max(0.0f, nDotL/falloff * pLight->wattage() / (4 * PI)) * diffuseColor);
        specular += result * (pow(std::max(0.0f, eDotr/falloff * pLight->wattage() / (4 * PI)), m_a)) * m_ks;
    }

    // add the ambient component
    L += diffuse * 0.1;

    return refraction * refractResult + reflection * reflectResult + diffuse * L + specular;
}
