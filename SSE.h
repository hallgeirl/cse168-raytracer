#ifndef __SSE_H__

#ifdef __SSE4_1__

#include <smmintrin.h>

typedef struct SSEVectorTuple3
{
    __m128 v[3];
};

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
