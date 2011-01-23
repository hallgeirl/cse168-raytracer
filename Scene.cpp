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
    #pragma omp parallel for private(ray, hitInfo, shadeResult)
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
        img->drawScanline(i);
        printf("Rendering Progress: %.3f%%\r", i/float(img->height())*100.0f);
        fflush(stdout);
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

    if (depth < TRACE_DEPTH && trace(hitInfo, ray))
    {
        shadeResult = hitInfo.material->shade(ray, hitInfo, *this);
		++depth;

		//if reflective material, send trace with ReflectRay
		float reflection = hitInfo.material->GetReflection(); 
		if (reflection > 0.0f)
		{
			Ray reflectRay = ray.Reflect(hitInfo.P, hitInfo.N);
			//fudge factor for now
			reflectRay.o += reflectRay.d * 0.0005;
			if (traceScene(reflectRay, reflectResult, depth))
			{
				shadeResult = reflection * reflectResult + (1 - reflection) * shadeResult;
			}
		}

		//if refractive material, send trace with RefractRay
		if (hitInfo.material->GetRefraction() > 0.0f)
		{
			Ray	refractRay;
			if (traceScene(refractRay, refractResult, depth))
			{

			}
		}

		return true;
	}
	return false;
}