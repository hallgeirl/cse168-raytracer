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
    #pragma omp parallel for private(ray, shadeResult)
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
			/*else
			{
				img->setPixel(j, i, getEnvironmentMap(ray));
			}*/
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
    bool result = m_bvh.intersect(minHit, ray, tMin, tMax);
    if (result)
    {
        //Bump mapping
        float delta = 0.00001;
        if (minHit.material->GetLookupCoordinates() == UV)
        {
            //Take a few samples to calculate the derivative
            tex_coord2d_t center = minHit.object->toUVCoordinates(minHit.P);
            float u = center.u, v = center.v;
            float u1 = minHit.material->bumpHeight2D(tex_coord2d_t(u-delta, v)), 
                  u2 = minHit.material->bumpHeight2D(tex_coord2d_t(u+delta, v)),
                  v1 = minHit.material->bumpHeight2D(tex_coord2d_t(u, v-delta)),
                  v2 = minHit.material->bumpHeight2D(tex_coord2d_t(u, v+delta));
                  
            //Approximate derivatives using central finite differences
            float dx = (u2-u1)/(2*delta),
                  dy = (v2-v1)/(2*delta);
            
            minHit.N += dx*(cross(minHit.N, Vector3(0,0,1)))-dy*(cross(minHit.N, Vector3(1,0,0)));
            minHit.N.normalize();
            
        }   
        //Todo: implement for 3D
        //bumpHeight = minHit.material->bumpHeight3D(tex_coord3d_t(minHit.P.x, minHit.P.y, minHit.P.z));
            
        
    }
    return result;
}

bool
Scene::traceScene(const Ray& ray, Vector3& shadeResult, int depth)
{
    HitInfo hitInfo;
	shadeResult = Vector3(0.f);

    if (depth < TRACE_DEPTH)
    {
		if (trace(hitInfo, ray))
		{
			//shadeResult = hitInfo.material->shade(ray, hitInfo, *this);
			++depth;
			
			//min checks should be done in material constructor, not every trace
		
			//if reflective material, send trace with ReflectRay
			if (hitInfo.material->IsReflective())
			{
				Vector3 reflectResult;
				Ray reflectRay = ray.Reflect(hitInfo);

				//Nudge Ray along normal to avoid Acne
				reflectRay.o += reflectRay.d * epsilon;
				if (traceScene(reflectRay, reflectResult, depth))
				{
					shadeResult += hitInfo.material->GetReflection()* reflectResult;
				}
			}

			//if refractive material, send trace with RefractRay
			if (hitInfo.material->IsRefractive())
			{
				Vector3 refractResult;
				Ray	refractRay = ray.Refract(hitInfo);

				//Nudge Ray along normal to avoid Acne
				refractRay.o += refractRay.d * epsilon;
				if (traceScene(refractRay, refractResult, depth))
				{
					shadeResult += hitInfo.material->GetRefraction()* refractResult;
				}
			}

			shadeResult += hitInfo.material->shade(ray, hitInfo, *this);
		}
		else
		{
			shadeResult = getEnvironmentMap(ray);
		}
		return true;
	}
	return false;
}

Vector3
Scene::getEnvironmentMap(const Ray & ray)
{
	Vector3 envResult;
	//Environment mapping here
	if (m_environment != 0)
	{
		tex_coord2d_t coords;
		//Calculate texture coordinates for where the ray hits the "sphere"
		coords.u = (atan2(ray.d.x, ray.d.z)) / (2.0f * PI) + 0.5;
		coords.v = (asin(ray.d.y)) / PI + 0.5;
		//And just look up the shading value in the texture.
		envResult = m_environment->lookup2D(coords);
	}
	else
	{
		envResult.x = 0.5; envResult.y = 0.5; envResult.z = 0.5;
	}
	return envResult;
}
