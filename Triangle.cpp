#include "Triangle.h"
#include "TriangleMesh.h"
#include "Ray.h"

#ifdef __SSE4_1__
#include <emmintrin.h>
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
    const Vector3 & A = m_mesh->vertices()[ti3.x];
    const Vector3 & B = m_mesh->vertices()[ti3.y];
    const Vector3 & C = m_mesh->vertices()[ti3.z];
    const Vector3 & nA = m_mesh->normals()[ni3.x];
    const Vector3 & nB = m_mesh->normals()[ni3.y];
    const Vector3 & nC = m_mesh->normals()[ni3.z];

	HitInfo result2;

#ifdef __SSE4_1__	
	//this actually loads in reverse order...
	__m128 _rd = _mm_set_ps(-r.d.x, -r.d.y, -r.d.z, 0.0f);
	__m128 _ro = _mm_set_ps(r.o.x, r.o.y, r.o.z, 0.0f);

	__m128 _A = _mm_set_ps(A.x, A.y, A.z, 0.0f);
	__m128 _B = _mm_set_ps(B.x, B.y, B.z, 0.0f);
	__m128 _C = _mm_set_ps(C.x, C.y, C.z, 0.0f);

	__m128 _nA = _mm_set_ps(nA.x, nA.y, nA.z, 0.0f);
	__m128 _nB = _mm_set_ps(nB.x, nB.y, nB.z, 0.0f);
	__m128 _nC = _mm_set_ps(nC.x, nC.y, nC.z, 0.0f);

	__m128 _BmA = _mm_sub_ps(_B, _A);
	__m128 _CmA = _mm_sub_ps(_C, _A);
	__m128 _normal = _mm_cross_ps(_BmA, _CmA);
	__m128 _ddotn = _mm_dot_ps(_rd, _normal);

	__m128 _t = _mm_div_ps(_mm_dot_ps(_mm_sub_ps(_ro, _A), _normal), _ddotn);
	__m128 _beta = _mm_div_ps(_mm_dot_ps(_rd, _mm_cross_ps(_mm_sub_ps(_ro, _A), _CmA)), _ddotn);
	__m128 _gamma = _mm_div_ps(_mm_dot_ps(_rd, _mm_cross_ps(_BmA, _mm_sub_ps(_ro, _A))), _ddotn);

	//can clean up with a compare function
    if (_beta.m128_f32[3] < -epsilon || _gamma.m128_f32[3] < -epsilon || _beta.m128_f32[3]+_gamma.m128_f32[3] > 1+epsilon || _t.m128_f32[3] < tMin || _t.m128_f32[3] > tMax) 
		return false;

	__m128 _P = _mm_add_ps(_A, _mm_add_ps(_mm_mul_ps(_beta, _BmA), _mm_mul_ps(_gamma, _CmA)));
	__m128 _alpha = _mm_sub_ps(_mm_sub_ps(_mm_set1_ps(1.0f), _beta), _gamma);
	__m128 _N = _mm_add_ps(_mm_mul_ps(_alpha, _nA), _mm_add_ps(_mm_mul_ps(_beta, _nB), _mm_mul_ps(_gamma, _nC)));

	result.P.x = _P.m128_f32[3];
	result.P.y = _P.m128_f32[2];
	result.P.z = _P.m128_f32[1];

	result.N.x = _N.m128_f32[3];
	result.N.y = _N.m128_f32[2];
	result.N.z = _N.m128_f32[1];

	result.t = _t.m128_f32[3];

#else

    Vector3 BmA = B-A, CmA = C-A;
    Vector3 normal = cross(BmA, CmA);
    float ddotn = (dot(-r.d, normal));

    float t = dot(r.o-A, normal) / ddotn;
    float beta = dot(-r.d, cross(r.o-A, CmA)) / ddotn;
    float gamma = dot(-r.d, cross(BmA, r.o-A)) / ddotn;

    if (beta < -epsilon || gamma < -epsilon || beta+gamma > 1+epsilon || t < tMin || t > tMax) return false;

    result2.P = A + beta*BmA + gamma*CmA;
    result2.t = t;
    result2.N = (1-beta-gamma)*nA + beta*nB + gamma*nC;

#endif

    result.material = m_material;

    return true;
}
