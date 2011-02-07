#include "Triangle.h"
#include "TriangleMesh.h"
#include "Ray.h"

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

using namespace std;

Triangle::Triangle(TriangleMesh * m, unsigned int i) :
    m_mesh(m), m_index(i)
{

}


Triangle::~Triangle()
{

}


void
Triangle::renderGL()
{
    TriangleMesh::TupleI3 ti3 = m_mesh->vIndices()[m_index];
    const Vector3 & v0 = m_mesh->vertices()[ti3.x]; //vertex a of triangle
    const Vector3 & v1 = m_mesh->vertices()[ti3.y]; //vertex b of triangle
    const Vector3 & v2 = m_mesh->vertices()[ti3.z]; //vertex c of triangle

    glBegin(GL_TRIANGLES);
        glVertex3f(v0.x, v0.y, v0.z);
        glVertex3f(v1.x, v1.y, v1.z);
        glVertex3f(v2.x, v2.y, v2.z);
    glEnd();
}

#ifdef __SSE4_1__
inline __m128 
_mm_cross_ps( __m128 a , __m128 b ) {
	__m128 ea , eb;
	//account for reverse ordering

	// set to a[1][2][0][3] , b[2][0][1][3]
	ea = _mm_shuffle_ps( a, a, _MM_SHUFFLE(2,1,3,0) ); //3,0,2,1
	eb = _mm_shuffle_ps( b, b, _MM_SHUFFLE(1,3,2,0) ); //3,1,0,2
	// multiply
	__m128 xa = _mm_mul_ps( ea , eb );
	// set to a[2][0][1][3] , b[1][2][0][3]
	a = _mm_shuffle_ps( a, a, _MM_SHUFFLE(1,3,2,0) ); //3,1,0,2
	b = _mm_shuffle_ps( b, b, _MM_SHUFFLE(2,1,3,0) );	//3,0,2,1
	// multiply
	__m128 xb = _mm_mul_ps( a , b );
	// subtract
	return _mm_sub_ps( xa , xb );
}

inline __m128 
_mm_dot_ps(__m128 a, __m128 b)
{
	__m128 mul0 = _mm_mul_ps(a, b);
	__m128 swp0 = _mm_shuffle_ps(mul0, mul0, _MM_SHUFFLE(2,1,3,0)); //2, 3, 0, 1
	__m128 add0 = _mm_add_ps(mul0, swp0);
	__m128 swp1 = _mm_shuffle_ps(mul0, mul0, _MM_SHUFFLE(1,3,2,0)); //0, 1, 2, 3
	__m128 add1 = _mm_add_ps(add0, swp1);
	return add1;
}

inline __m128 
_mm_negate_ps( __m128 a )
{
  return _mm_sub_ps(_mm_setzero_ps(),a);
}

#endif

bool
Triangle::intersect(HitInfo& result, const Ray& r,float tMin, float tMax)
{
    TriangleMesh::TupleI3 ti3 = m_mesh->vIndices()[m_index];
    TriangleMesh::TupleI3 ni3 = m_mesh->nIndices()[m_index];

#ifdef __SSE4_1__
    const __m128 _A = m_mesh->SSEvertices()[ti3.x];
    const __m128 _B = m_mesh->SSEvertices()[ti3.y];
    const __m128 _C = m_mesh->SSEvertices()[ti3.z];
    const __m128 _nA = m_mesh->SSEnormals()[ni3.x];
    const __m128 _nB = m_mesh->SSEnormals()[ni3.y];
    const __m128 _nC = m_mesh->SSEnormals()[ni3.z];

	//Load the ray origin and direction
	__m128 _rd = _mm_sub_ps(_mm_setzero_ps(), r.d_SSE),
	       _ro = r.o_SSE;
	       
	__m128 _BmA = _mm_sub_ps(_B, _A);
	__m128 _CmA = _mm_sub_ps(_C, _A);
	__m128 _normal = _mm_cross_ps(_BmA, _CmA);
	__m128 _ddotn = _mm_dp_ps(_rd, _normal, 0xEE);

    //Calculate t, beta and gamma
    __m128 _beta_gamma_t = _mm_div_ps(_mm_add_ps(_mm_dp_ps(_mm_sub_ps(_ro, _A), _normal, 0xE2), //(O-A) * N = t -> r0
                                                 _mm_add_ps(_mm_dp_ps(_rd, _mm_cross_ps(_mm_sub_ps(_ro, _A), _CmA), 0xE4), //d * (O-A)x(C-A) = beta -> r1
                                                            _mm_dp_ps(_rd, _mm_cross_ps(_BmA, _mm_sub_ps(_ro, _A)), 0xE8))), // = gamma -> r2
                                      _ddotn);
    
    //Compare to see if we have a hit 
    float betagammat[4];
    _mm_store_ps(betagammat, _beta_gamma_t);
    
    if (betagammat[2] < -epsilon || betagammat[3] < -epsilon || betagammat[2]+betagammat[3] > 1+epsilon || betagammat[1] < tMin || betagammat[1] > tMax) 
    {
		return false;
	}


    __m128 _beta  = _mm_shuffle_ps(_beta_gamma_t, _beta_gamma_t, _MM_SHUFFLE(2, 2, 2, 2)),
           _gamma = _mm_shuffle_ps(_beta_gamma_t, _beta_gamma_t, _MM_SHUFFLE(3, 3, 3, 3));
	__m128 _alpha = _mm_sub_ps(_mm_sub_ps(_mm_set1_ps(1.0f), _beta), _gamma);
           
    //_mm_test_all_zeros(_beta, _gamma);
    __m128 _P = _mm_add_ps(_A, _mm_add_ps(_mm_mul_ps(_beta, _BmA), _mm_mul_ps(_gamma, _CmA)));
	__m128 _N = _mm_add_ps(_mm_mul_ps(_alpha, _nA), _mm_add_ps(_mm_mul_ps(_beta, _nB), _mm_mul_ps(_gamma, _nC)));
	
    float P[4], N[4];
    
    _mm_store_ps(P, _P);
    _mm_store_ps(N, _N);

	result.P.x = P[3];
	result.P.y = P[2];
	result.P.z = P[1];

	result.N.x = N[3];
	result.N.y = N[2];
	result.N.z = N[1];
	
	result.t = betagammat[1];

#else

    const Vector3 & A = m_mesh->vertices()[ti3.x];
    const Vector3 & B = m_mesh->vertices()[ti3.y];
    const Vector3 & C = m_mesh->vertices()[ti3.z];
    const Vector3 & nA = m_mesh->normals()[ni3.x];
    const Vector3 & nB = m_mesh->normals()[ni3.y];
    const Vector3 & nC = m_mesh->normals()[ni3.z];

    Vector3 BmA = B-A, CmA = C-A;
    Vector3 normal = cross(BmA, CmA);
    float ddotn = (dot(-r.d, normal));

    float t = dot(r.o-A, normal) / ddotn;
    float beta = dot(-r.d, cross(r.o-A, CmA)) / ddotn;
    float gamma = dot(-r.d, cross(BmA, r.o-A)) / ddotn;

    if (beta < -epsilon || gamma < -epsilon || beta+gamma > 1+epsilon || t < tMin || t > tMax) return false;

    result.P = A + beta*BmA + gamma*CmA;
    result.t = t;
    result.N = (1-beta-gamma)*nA + beta*nB + gamma*nC;

#endif

    result.material = m_material;

    return true;
}
