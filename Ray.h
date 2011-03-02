#ifndef CSE168_RAY_H_INCLUDED
#define CSE168_RAY_H_INCLUDED

#include <stdlib.h>
#include "Vector3.h"
#include "Material.h"
#include "Utility.h"

#ifdef STATS
#include "Stats.h"
#endif

#include "SSE.h"
    
    //! Contains information about a ray hit with a surface.
    /*!
        HitInfos are used by object intersection routines. They are useful in
        order to return more than just the hit distance.
    */
    class HitInfo
    {
    public:
        float t;                            //!< The hit distance
        Vector3 P;                          //!< The hit point
        Vector3 N;                          //!< Shading normal vector
        const Material* material;           //!< Material of the intersected object
        const Object  * object;             //!< Material of the intersected object
    
        //! Default constructor.
        explicit HitInfo(float t = 0.0f,
                         const Vector3& P = Vector3(),
                         const Vector3& N = Vector3(0.0f, 1.0f, 0.0f)) :
            t(t), P(P), N(N), material (0)
        {
            // empty
        }
    };
    
    class Ray
    {
    
    public:
        bool diffuse;
        Vector3 o,      //!< Origin of ray
                d;      //!< Direction of ray
    
    
        #ifdef __SSE4_1__
        __m128 d_SSE, o_SSE, d_SSE_rcp;
    //    SSEVectorTuple3 o_SSE, d_SSE;
    
        void setupSSE()
        {
            static const __m128 _zero = _mm_setzero_ps();
            float _o[4] = {o.x, o.y, o.z, o.x};
            float _d[4] = {d.x, d.y, d.z, d.x};

            d_SSE = _mm_sub_ps(_zero, _mm_loadu_ps(_d));
            d_SSE_rcp = _mm_rcp_ps(_mm_loadu_ps(_d));
            o_SSE = _mm_loadu_ps(_o);


        /*__m128 _d_sse = _mm_loadu_ps(_d);
        __m128 _o_sse = _mm_loadu_ps(_o);
        #pragma unroll(3)
        for (int i = 0; i < 3; i++)
        {
            o_SSE.v[i] = _mm_shuffle_ps(_o_sse, _o_sse, _MM_SHUFFLE(i,i,i,i));
            d_SSE.v[i] = _mm_sub_ps(_zero, _mm_shuffle_ps(_d_sse, _d_sse, _MM_SHUFFLE(i,i,i,i)));
        }*/
    }
    #endif

    Ray() : o(), d(Vector3(0.0f,0.0f,1.0f))
    {
        diffuse = false;
#ifdef STATS
		Stats::Rays++;
#endif

        #ifdef __SSE4_1__
        setupSSE();
        #endif
    }

    Ray(const Vector3& o, const Vector3& d) : o(o), d(d)
    {
        diffuse = false;
#ifdef STATS
		Stats::Rays++;
#endif
        #ifdef __SSE4_1__
        setupSSE();
        #endif
    }

	Ray Random(const HitInfo & hitInfo) const
	{
		//bias to the surface normal
	    float theta = asin(sqrt((float) rand() / (float)RAND_MAX));
		float phi = 2.0f * PI * ((float) rand() / (float)RAND_MAX);

		//convert spherical coords to cartesian coords
		float u1 = sin(theta) * cos(phi);
		float u2 = sin(theta) * sin(phi);
		float u3 = cos(theta);

		//pick random vector
		float n[3] = {hitInfo.N.x, hitInfo.N.y, hitInfo.N.z};
		int maxDimIndex = 0;
		if (n[1] > n[0]) 
			maxDimIndex = 1;
		if (n[2] > n[maxDimIndex])
			maxDimIndex = 2;
		
		//determine one surface tangent from cross of normal and random vec
		Vector3 randomVec(maxDimIndex == 2 ? -n[2] : 0, maxDimIndex == 0 ? -n[0] : 0, maxDimIndex == 1 ? -n[1] : 0);
		Vector3 t1 = cross(randomVec, hitInfo.N);

		//Random direction determined by combination of random values by surface axes
		Vector3 random_d(u1*t1 + u2*cross(t1, hitInfo.N) + u3*hitInfo.N);

		random_d.normalize();
		Ray random(hitInfo.P + random_d * epsilon, random_d);
        random.diffuse = true;
#ifdef STATS
		Stats::Secondary_Rays++;
#endif
		return random;
	}
	
    //Shoots a reflection ray. If path tracing is enabled, shoot a random ray according to the glossyness of the material.
	Ray Reflect(const HitInfo & hitInfo) const
	{
#ifdef PATH_TRACING
        //Generate randomized reflection ray based on glossyness
        //bias to the perfectly reflected ray
	    float phi = acos(pow(frand(), 1/(1+hitInfo.material->getShininess())));
		float theta = 2.0f * PI * frand();

        //Direction of perfect reflection
	    Vector3 d_reflect = d - 2 * dot(hitInfo.N, d) * hitInfo.N;

		//convert spherical coords to cartesian coords
		float u1 = sin(phi) * cos(theta);
		float u2 = sin(phi) * sin(theta);
		float u3 = cos(phi); 

		//determine one surface tangent from cross of normal and one of the unit vectors 
		Vector3 t1 = cross(Vector3(0,0,1), d_reflect);
        if (t1.length2() < 1e-6) t1 = cross(Vector3(0, 1, 0), d_reflect);

		//Random direction determined by combination of random values by surface axes
		Vector3 random_d(u1*t1 + u2*cross(t1, d_reflect) + u3*d_reflect);

		random_d.normalize();
		Ray reflect(hitInfo.P + random_d * epsilon, random_d);
#else
	    Vector3 d_r = d - 2 * dot(hitInfo.N, d) * hitInfo.N;
		Ray reflect(hitInfo.P + d_r * epsilon, d_r);
#endif
        

#ifdef STATS
		Stats::Secondary_Rays++;
#endif

		return reflect;
	}

	Ray Refract(const HitInfo & hitInfo) const
	{
		float n1, n2;
		Vector3 n;

		// if ray enters object, else ray exits object
		if ( dot(d, hitInfo.N) < 0)
		{
			n1 = 1.0f;
			n2 = hitInfo.material->getRefractionIndex();
			n = hitInfo.N;
		}
		else 
		{
			n1 = hitInfo.material->getRefractionIndex();
			n2 = 1.0f;
			n = -hitInfo.N;
		}

		float energy = 1 - (pow(n1, 2) * (1 - pow(dot(d, n), 2)) / pow(n2, 2));

		// Total internal reflection: all of the energy is reflected
		if (energy < 0)
        {
            return Reflect(hitInfo);
        }

        Vector3 d_r = n1 * (d - n * dot(d, n)) / n2 - n * sqrt(energy);
		Ray refract(hitInfo.P + d_r * epsilon, d_r);

#ifdef STATS
		Stats::Secondary_Rays++;
#endif

		return refract;
	}

};


#endif // CSE168_RAY_H_INCLUDED
