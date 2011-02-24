#ifndef __UTILITY__H_
#define __UTILITY__H_

#include <cstdlib>

double getTime();
void getEigenVector(const float (&A)[3][3], float (&outV)[3], float lambda);

//Returns random number between 0 and 1
inline float frand()
{
    return (float)rand() / (float)RAND_MAX;
}

#endif
