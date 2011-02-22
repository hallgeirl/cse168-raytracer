#ifndef __UTILITY__H_
#define __UTILITY__H_

//typedef float[3][3] Matrix3x3;

double getTime();
void getEigenVector(const float (&A)[3][3], float (&outV)[3], float lambda);

#endif
