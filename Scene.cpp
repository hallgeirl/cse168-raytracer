#include "Utility.h"
#include <cmath>
#include <iostream>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"

#ifdef OPENMP
#include <omp.h>
#endif 

using namespace std;

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
    debug("Precalcing objects\n");
    double t1 = -getTime();
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
    t1 += getTime();
    printf("Time spent preprocessing objects and lights: %lf\n", t1);
    debug("Building BVH...\n");
    t1 = -getTime();
    m_bvh.build(&m_objects);
    t1 += getTime();
    debug("Done building BVH.\n");
    printf("Time spent building BVH: %lf\n", t1);
}


inline float tonemapValue(float value, float maxIntensity)
{
    return std::min(pow(value / maxIntensity, 0.85f)*1.5f, 1.0f);

}

void
Scene::raytraceImage(Camera *cam, Image *img)
{
	int depth = TRACE_DEPTH;
    float minIntensity = infinity, maxIntensity = -infinity;
    printf("Rendering Progress: %.3f%%\r", 0.0f);
    fflush(stdout);

    //For tone mapping. The Image class stores the pixels internally as 1 byte integers. We want to store the actual values first.
    int width = img->width(), height = img->height();
    Vector3 *tempImage = new Vector3[height*width];

    double t1 = -getTime();

    // loop over all pixels in the image
    #ifdef OPENMP
    #pragma omp parallel for schedule(dynamic, 2)
    #endif
    for (int i = 0; i < height; ++i)
    {
        float localMaxIntensity = -infinity,
              localMinIntensity = infinity;

        for (int j = 0; j < width; ++j)
        {
            Ray ray;
            Vector3 tempShadeResult;
            Vector3 shadeResult(0.f);

			#ifdef PATH_TRACING
			for (int k = 0; k < TRACE_SAMPLES; ++k)
			{
                ray = cam->eyeRay(j, i, width, height, true);
				if (traceScene(ray, tempShadeResult, depth))
				{
					shadeResult += tempShadeResult;
				}
			}
			shadeResult /= TRACE_SAMPLES; 
            tempImage[i*width+j] = shadeResult;
			#else
            ray = cam->eyeRay(j, i, width, height, false);
			if (traceScene(ray, shadeResult, depth))
			{
				tempImage[i*width+j] = shadeResult;
			}
			#endif // PATH_TRACING
            for (int k = 0; k < 3; k++)
            {
                if (shadeResult[k] > localMaxIntensity)
                    localMaxIntensity = shadeResult[k];
                if (shadeResult[k] < localMinIntensity)
                    localMinIntensity = shadeResult[k];
            }
            #ifdef OPENMP
            #pragma omp critical
            {
                if (localMinIntensity < minIntensity) minIntensity = localMinIntensity;
                if (localMaxIntensity > maxIntensity) maxIntensity = localMaxIntensity;
            }
            #else
            minIntensity = localMinIntensity
            #endif
        }
        #ifdef OPENMP
        #pragma omp master
        #endif
        {
            printf("Rendering Progress: %.3f%%\r", i/float(img->height())*100.0f);
            fflush(stdout);
        }
    }
    debug("Performing tone mapping...");

    #ifdef OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            Vector3 finalColor = tempImage[i*width+j];

            #pragma unroll(3)
            for (int k = 0; k < 3; k++)
                finalColor[k] = tonemapValue(finalColor[k], maxIntensity);

            img->setPixel(j, i, finalColor);
        }
        #ifndef NO_GFX //If not rendering graphics to screen, don't draw scan lines (it will segfault in multithreading mode)
        img->drawScanline(i);
        #endif
    }
    t1 += getTime();

    printf("Rendering Progress: 100.000%%\n");
    debug("Done raytracing!\n");
    printf("Time spent raytracing image: %lf seconds.\n", t1);
}

bool
Scene::trace(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    bool result = m_bvh.intersect(minHit, ray, tMin, tMax);
    
    //Trace the unbounded objects (like planes)
    for (int i = 0; i < m_unboundedObjects.size(); i++)
    {
        HitInfo tempMinHit;
        bool ubresult = m_unboundedObjects[i]->intersect(tempMinHit, ray, tMin, tMax);
        if (ubresult && (!result || tempMinHit.t < minHit.t))
        {
            result = true;
            minHit = tempMinHit;
            minHit.object = m_unboundedObjects[i];
        }
    }

    if (result)
    {
        //Bump mapping
        float delta = 0.0001;
//        printf("object %p, material from %p\n", minHit.object, minHit.material);

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
            
            //Find two tangents
            float n[3] = { minHit.N.x, minHit.N.y, minHit.N.z };

            int m = 0;
            if (n[1] > n[0]) m = 1;
            if (n[2] > n[m]) m = 2;
            Vector3 randomVec(m == 2 ? -n[2] : 0, m == 0 ? -n[0] : 0, m == 1 ? -n[1] : 0);
           
            Vector3 t1 = cross(minHit.N, randomVec);
            minHit.N += dx*(cross(minHit.N, t1))-dy*(cross(minHit.N, cross(minHit.N, t1)));
            minHit.N.normalize();
        }   
        //Todo: implement for 3D
        //bumpHeight = minHit.material->bumpHeight3D(tex_coord3d_t(minHit.P.x, minHit.P.y, minHit.P.z));
            
        
    }
    return result;
}

bool Scene::traceScene(const Ray& ray, Vector3& shadeResult, int depth)
{
    HitInfo hitInfo;
	shadeResult = Vector3(0.f);
    bool hit = false;
    
    if (depth >= 0)
    {
		if (trace(hitInfo, ray))
		{
            hit = true;

			--depth;
			
			shadeResult = hitInfo.material->shade(ray, hitInfo, *this);
			
            #ifdef PATH_TRACING
			//if diffuse material, send trace with RandomRay generate by Monte Carlo
			if (hitInfo.material->IsDiffuse())
			{
				Vector3 diffuseResult;
				Ray diffuseRay = ray.Random(hitInfo);

				if (traceScene(diffuseRay, diffuseResult, depth))
				{
					shadeResult += (hitInfo.material->GetDiffuse() * diffuseResult);
				}
			}
			#endif

			
			//if reflective material, send trace with ReflectRay
			if (hitInfo.material->IsReflective())
			{
				Vector3 reflectResult;
				Ray reflectRay = ray.Reflect(hitInfo);

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

				if (traceScene(refractRay, refractResult, depth))
				{
					shadeResult += hitInfo.material->GetRefraction()* refractResult;
				}
			}
		}
		else
		{
            shadeResult = getEnvironmentMap(ray);
            hit = true;
		}
	}
    
    return hit;
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
		envResult = m_bgColor; 
	}
	return envResult;
}
