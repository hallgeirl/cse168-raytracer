#ifndef CSE168_RAY_H_INCLUDED
#define CSE168_RAY_H_INCLUDED

#include "Vector3.h"
#include "Material.h"

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
    Vector3 o,      //!< Origin of ray
            d;      //!< Direction of ray


    #ifdef __SSE4_1__
    SSEVectorTuple3 o_SSE, d_SSE;

    void setupSSE()
    {
        static const __m128 _zero = _mm_setzero_ps();
        float _o[4] = {o.x, o.y, o.z, 0};
        float _d[4] = {d.x, d.y, d.z, 0};
        __m128 _d_sse = _mm_loadu_ps(_d);
        __m128 _o_sse = _mm_loadu_ps(_o);
        #pragma unroll(3)
        for (int i = 0; i < 3; i++)
        {
            o_SSE.v[i] = _mm_shuffle_ps(_o_sse, _o_sse, _MM_SHUFFLE(i,i,i,i));
            d_SSE.v[i] = _mm_sub_ps(_zero, _mm_shuffle_ps(_d_sse, _d_sse, _MM_SHUFFLE(i,i,i,i)));
        }
    }
    #endif

    Ray() : o(), d(Vector3(0.0f,0.0f,1.0f))
    {
        #ifdef __SSE4_1__
        setupSSE();
        #endif
    }

    Ray(const Vector3& o, const Vector3& d) : o(o), d(d)
    {
        #ifdef __SSE4_1__
        setupSSE();
        #endif
    }

	Ray Reflect(const HitInfo & hitInfo) const
	{
	    Vector3 d_r = d - 2 * dot(hitInfo.N, d) * hitInfo.N;
		Ray reflect(hitInfo.P + d_r * epsilon, d_r);
        
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
			n2 = hitInfo.material->GetRefractionIndex();
			n = hitInfo.N;
		}
		else 
		{
			n1 = hitInfo.material->GetRefractionIndex();
			n2 = 1.0f;
			n = -hitInfo.N;
		}

		float energy = 1 - (pow(n1, 2) * (1 - pow(dot(d, n), 2)) / pow(n2, 2));

		// Total internal reflection: all of the energy is reflected
		if (energy < 0)
			return Reflect(hitInfo);

        Vector3 d_r = n1 * (d - n * dot(d, n)) / n2 - n * sqrt(energy);
		Ray refract(hitInfo.P + d_r * epsilon, d_r);

		return refract;
	}

};


#endif // CSE168_RAY_H_INCLUDED
