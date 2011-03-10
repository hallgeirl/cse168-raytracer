#ifndef WIN32
#include <sys/time.h>
#else
#include <windows.h>
#endif
#include <stdio.h>
#include <cmath>
#include <iostream>

#include "includes.h"

using namespace std;

void addModel(const char* filename, Material *mat, Scene* scene, Vector3 position, float rotY, Vector3 scale)
{
	TriangleMesh * mesh = new TriangleMesh();
    Matrix4x4 m_rot(Vector4(cos(rotY),0,-sin(rotY),0), Vector4(0,1,0,0), Vector4(sin(rotY),0,cos(rotY),0), Vector4(0, 0, 0, 1));
    Matrix4x4 m_trans(Vector4(1,0,0,0), Vector4(0,1,0,0), Vector4(0,0,1,0), Vector4(position.x, position.y, position.z, 1));
    Matrix4x4 m_scale(Vector4(scale.x,0,0,0), Vector4(0,scale.y,0,0), Vector4(0,0,scale.z,0), Vector4(0, 0, 0, 1));
	mesh->load(filename, m_trans*m_rot*m_scale);
    for (int i = 0; i < mesh->numTris(); i++)
    {
		Triangle *tri = new Triangle();
		tri->setMesh(mesh);
		tri->setIndex(i);
		tri->setMaterial(mat);
		g_scene->addObject(tri);
    }
}

//Returns the time in seconds
double getTime()
{
    #ifndef WIN32
    struct timeval tv;
    int ret = gettimeofday(&tv, 0);
    double t = (double)tv.tv_sec + 1e-6 * (double)tv.tv_usec;
    if (ret == -1)
    {
        printf("Something bad happened when calling gettimeofday().\n");
        return -1;
    }
    return t;

    #else
	return GetTickCount()/1000.f;
	#endif
}

void printMat(float A[3][3] )
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

//Returns one eigen vector for a 3x3 matrix A given the eigen value lambda
//Easily expandable to arbitrary dimensions but we probably won't need to.
void getEigenVector(const float (&A)[3][3], float (&outV)[3], float lambda)
{
    int n = 3;
    float E[2][3];

    //Reduce A to echelon form. Note that the last row is always 0 if lambda is actually an eigenvalue of A, so don't bother calculate it.
    //Copy the array (so we leave the original one intact)
    for (int i = 0; i < n-1; i++)
    {
        for (int j = 0; j < n; j++)
        {
            E[i][j] = A[i][j];
        }
        E[i][i] -= lambda;
    }

    //Do the row reductions (for all except the last row)
    for (int j = 0; j < n-2; j++)
    {
        for (int i = j+1; i < n; i++)
        {
            float m = E[i][j]/E[j][j];
            for (int k = 0; k < n; k++)
            {
                E[i][k] = E[i][k] - m*E[j][k];
            }
        }
    }


    //Now, get it to reduced echelon form by going backwards.
    for (int j = n-2; j >= 1; j--)
    {
        for (int i = j-1; i >= 0; i--)
        {
            float m = E[i][j]/E[j][j];
            for (int k = j; k < n; k++)
            {
                E[i][k] = E[i][k] - m*E[j][k];
            }
        }
    }

    //Get the eigenvectors from E
    outV[n-1] = 1;
    float length = 1;
    for (int i = 0; i < n-1; i++)
    {
        outV[i] = E[i][n-1] / E[i][i];
        length += (outV[i]*outV[i]);
    }
    
    //Normalize it
    for (int i = 0; i < 3; i++)
    {
        outV[i] /= sqrt(length);
    }
}
