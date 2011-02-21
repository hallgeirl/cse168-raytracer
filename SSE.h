#ifndef __SSE_H__
#define __SSE_H__

#ifdef __SSE4_1__

#include <iostream>
#include <smmintrin.h>

struct SSEVectorTuple3
{
    __m128 v[3];
};

//Do the cross product of 4 pairs of vectors
inline SSEVectorTuple3 SSEmultiCross(const SSEVectorTuple3 &a, const SSEVectorTuple3 &b)
{
    SSEVectorTuple3 out;
    //Do the cross product
    out.v[0] = _mm_sub_ps(_mm_mul_ps(a.v[1], b.v[2]), _mm_mul_ps(a.v[2], b.v[1]));
    out.v[1] = _mm_sub_ps(_mm_mul_ps(a.v[2], b.v[0]), _mm_mul_ps(a.v[0], b.v[2]));
    out.v[2] = _mm_sub_ps(_mm_mul_ps(a.v[0], b.v[1]), _mm_mul_ps(a.v[1], b.v[0]));
    return out;
}


inline void SSEprintVec(const __m128 &v)
{
    float tmp[4];
    _mm_storeu_ps(tmp, v);
    std::cout << "__m128(v3=" << tmp[0] << ",\tv2=" << tmp[1] << ",\tv1=" << tmp[2] << ",\tv0=" << tmp[3] << ")" << std::endl;
}

inline void SSEprintVecTuple(const SSEVectorTuple3 &v)
{
    std::cout << "SSEVectorTuple3(" << std::endl;
    for (int i = 0; i < 3; i++)
        SSEprintVec(v.v[i]);
    std::cout << ")" << std::endl;
}

inline __m128 SSEmultiDot13(const __m128 &a, const SSEVectorTuple3 &b)
{
    return _mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 0, 0, 0)), b.v[0]) ,_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 1, 1, 1)), b.v[1]), _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(2, 2, 2, 2)), b.v[2])));
}

inline __m128 SSEmultiDot33(const SSEVectorTuple3 &a, const SSEVectorTuple3 &b)
{
    return _mm_add_ps(_mm_mul_ps(a.v[0], b.v[0]) ,_mm_add_ps(_mm_mul_ps(a.v[1], b.v[1]), _mm_mul_ps(a.v[2], b.v[2])));
}

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

#endif

#endif
