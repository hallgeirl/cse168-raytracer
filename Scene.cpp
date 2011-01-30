#include <cmath>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"

Scene * g_scene = 0;

void
Scene::openGL(Camera *cam)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    cam->drawGL();

    // draw objects
    for (size_t i = 0; i < m_objects.size(); ++i)
        m_objects[i]->renderGL();

    glutSwapBuffers();
}

void
Scene::preCalc()
{
    Objects::iterator it;
    for (it = m_objects.begin(); it != m_objects.end(); it++)
    {
        Object* pObject = *it;
        pObject->preCalc();
    }
    Lights::iterator lit;
    for (lit = m_lights.begin(); lit != m_lights.end(); lit++)
    {
        PointLight* pLight = *lit;
        pLight->preCalc();
    }

    m_bvh.build(&m_objects);
}

void
Scene::raytraceImage(Camera *cam, Image *img)
{
    Ray ray;
    Vector3 shadeResult;
	int depth = 0;

    // loop over all pixels in the image
    #ifdef OPENMP
    #pragma omp parallel for private(ray, shadeResult, depth)
    #endif
    for (int i = 0; i < img->height(); ++i)
    {
        for (int j = 0; j < img->width(); ++j)
        {
            ray = cam->eyeRay(j, i, img->width(), img->height());
            if (traceScene(ray, shadeResult, depth))
            {
				img->setPixel(j, i, shadeResult);
            }
        }
        #ifndef NO_GFX //If not rendering graphics to screen, don't draw scan lines (it will segfault in multithreading mode)
        img->drawScanline(i);
        #endif
        if (i % 10 == 0)
        {
            printf("Rendering Progress: %.3f%%\r", i/float(img->height())*100.0f);
            fflush(stdout);
        }
    }

    printf("Rendering Progress: 100.000%%\n");
    debug("done Raytracing!\n");
}

bool
Scene::trace(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    return m_bvh.intersect(minHit, ray, tMin, tMax);
}

bool
Scene::traceScene(const Ray& ray, Vector3& shadeResult, int depth)
{
    HitInfo hitInfo;
	Vector3 reflectResult;
	Vector3 refractResult;

    if (depth < TRACE_DEPTH)
    {
		if (trace(hitInfo, ray))
		{
			//shadeResult = hitInfo.material->shade(ray, hitInfo, *this);
			++depth;

			//if reflective material, send trace with ReflectRay
			float reflection = fmin(hitInfo.material->GetReflection(), 1);
			if (reflection > 0.0f)
			{
				Ray reflectRay = ray.Reflect(hitInfo);
				//fudge factor for now
				reflectRay.o += reflectRay.d * 0.001;
				if (!traceScene(reflectRay, reflectResult, depth))
				{
					reflection = 0;
					reflectResult.set(0);
				}
			}

			float refraction = hitInfo.material->GetRefraction();

			//if refractive material, send trace with RefractRay
			if (refraction > 0.0f)
			{
				Ray	refractRay = ray.Refract(hitInfo);
				refractRay.o += refractRay.d * 0.0005;
				if (!traceScene(refractRay, refractResult, depth))
				{
					refraction = 0;
					refractResult.set(0);
				}
			}
			//Keep the energy equation balanced (that is, don't refract+reflect+absorb more than 100% of the ray)
			reflection = fmin(reflection, 1-refraction);
			float absorb = 1-reflection-refraction;
			shadeResult = refraction * refractResult + reflection * reflectResult + absorb * hitInfo.material->shade(ray, hitInfo, *this);
		}
		else
		{
			//Environment mapping here
			if (m_environment != 0)
			{
				tex_coord_t coords;
				//Calculate texture coordinates for where the ray hits the "sphere"
				coords.u = (atan2(ray.d.x, ray.d.z)) / (2.0f * PI) + 0.5;
				coords.v = (asin(ray.d.y)) / PI + 0.5;
				//And just look up the shading value in the texture.
				shadeResult = m_environment->lookup(coords);
			}
			else
			{
				shadeResult.x = 0.5; shadeResult.y = 0.5; shadeResult.z = 0.5;
			}

		}
		return true;
	}
	return false;
}
