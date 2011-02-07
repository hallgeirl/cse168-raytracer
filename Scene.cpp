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
    #pragma omp parallel for private(ray, shadeResult) schedule(dynamic, 10)
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
    bool result = m_bvh.intersect(minHit, ray, tMin, tMax);
    if (result)
    {
        //Bump mapping
        float delta = 0.0001;
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
	Vector3 reflectResult;
	Vector3 refractResult;

    if (depth < TRACE_DEPTH)
    {
		if (trace(hitInfo, ray))
		{
			//shadeResult = hitInfo.material->shade(ray, hitInfo, *this);
			++depth;

			//if reflective material, send trace with ReflectRay
			float reflection = std::min(hitInfo.material->GetReflection(), 1.0f);
			if (reflection > 0.0f)
			{
				Ray reflectRay = ray.Reflect(hitInfo);
				//fudge factor for now
				//reflectRay.o += reflectRay.d * 0.001;
				if (!traceScene(reflectRay, reflectResult, depth))
				{
					reflection = 0;
					reflectResult.set(0);
				}
			}

			float refraction = std::min(hitInfo.material->GetRefraction(), 1.0f);

			//if refractive material, send trace with RefractRay
			if (refraction > 0.0f)
			{
				Ray	refractRay = ray.Refract(hitInfo);
				//refractRay.o += refractRay.d * 0.0005;
				if (!traceScene(refractRay, refractResult, depth))
				{
					refraction = 0;
					refractResult.set(0);
				}
			}
			
			//Keep the energy equation balanced (that is, don't refract+reflect+absorb more than 100% of the ray)
			reflection = std::min(reflection, 1.0f-refraction);
			float diffuse = 1-reflection-refraction;
			shadeResult = refraction * refractResult + reflection * reflectResult + diffuse * hitInfo.material->shade(ray, hitInfo, *this);
		}
		else
		{
    		//Environment mapping here
			if (m_environment != 0)
			{
				tex_coord2d_t coords;
				//Calculate texture coordinates for where the ray hits the "sphere"
				coords.u = (atan2(ray.d.x, ray.d.z)) / (2.0f * PI) + 0.5;
				coords.v = (asin(ray.d.y)) / PI + 0.5;
				//And just look up the shading value in the texture.
				shadeResult = m_environment->lookup2D(coords);
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
