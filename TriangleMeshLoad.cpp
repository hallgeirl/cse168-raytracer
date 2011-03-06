#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "TriangleMesh.h"
#include "Console.h"
#include <vector>

#ifndef __GNUC__
// disable useless warnings
#pragma warning(disable:4996)
#endif

using namespace std;

void
TriangleMesh::createSingleTriangle()
{
    m_normals = new Vector3[3];
    m_vertices = new Vector3[3];
    m_texCoords = new VectorR2[3];

    m_texCoords[0].x = 0.0f;
    m_texCoords[0].y = 0.0f;
    m_texCoords[1].x = 1.0f;
    m_texCoords[1].y = 0.0f;
    m_texCoords[2].x = 0.0f;
    m_texCoords[2].y = 1.0f;

    m_normalIndices = new TupleI3[1];
    m_vertexIndices = new TupleI3[1];
    m_texCoordIndices = new TupleI3[1];

    m_vertexIndices[0].v[0] = 0;
    m_vertexIndices[0].v[1] = 1;
    m_vertexIndices[0].v[2] = 2;

    m_normalIndices[0].v[0] = 0;
    m_normalIndices[0].v[1] = 1;
    m_normalIndices[0].v[2] = 2;

    m_texCoordIndices[0].v[0] = 0;
    m_texCoordIndices[0].v[1] = 1;
    m_texCoordIndices[0].v[2] = 2;

    #ifdef __SSE4_1__
#ifdef WIN32
    m_SSEvertices = (__m128*)_aligned_malloc(3 * sizeof(__m128), 16);
    m_SSEnormals = (__m128*)_aligned_malloc(3 * sizeof(__m128), 16);
#else
    m_SSEvertices = new __m128[3];
    m_SSEnormals = new __m128[3];
#endif
	#endif    

    m_numTris = 1;
}

//************************************************************************
// You probably don't want to modify the following functions
// They are for loading .obj files
//************************************************************************

bool
TriangleMesh::load(const char* file, const Matrix4x4& ctm)
{
    FILE *fp = fopen(file, "rb");
    if (!fp)
    {
        error("Cannot open \"%s\" for reading\n",file);
        return false;
    }
    debug("Loading \"%s\"...\n", file);

    loadObj(fp, ctm);
    debug("Loaded \"%s\" with %d triangles\n",file,m_numTris);
    fclose(fp);

    return true;
}

void
getIndices(char *word, int *vindex, int *tindex, int *nindex)
{
    char *null = (char*)" ";
    char *ptr;
    char *tp;
    char *np;

    // by default, the texture and normal pointers are set to the null string
    tp = null;
    np = null;

    // replace slashes with null characters and cause tp and np to point
    // to character immediately following the first or second slash
    for (ptr = word; *ptr != '\0'; ptr++)
    {
        if (*ptr == '/')
        {
            if (tp == null)
                tp = ptr + 1;
            else
                np = ptr + 1;

            *ptr = '\0';
        }
    }

    *vindex = atoi (word);
    *tindex = atoi (tp);
    *nindex = atoi (np);
}


void
TriangleMesh::loadObj(FILE* fp, const Matrix4x4& ctm)
{
    int nv=0, nt=0, nn=0, nf=0;
    char line[81];
    while (fgets(line, 80, fp) != 0)
    {
        if (line[0] == 'v')
        {
            if (line[1] == 'n')
                nn++;
            else if (line[1] == 't')
                nt++;
            else
                nv++;
        }
        else if (line[0] == 'f')
        {
            nf++;
        }
    }
    fseek(fp, 0, 0);


    m_normals = new Vector3[std::max(nv,nf*3)];
    m_vertices = new Vector3[nv];
    m_numVertices = nv;    
    #ifdef __SSE4_1__

#ifdef WIN32
    m_SSEvertices = (__m128*)_aligned_malloc(nv * sizeof(__m128), 16);
    m_SSEnormals = (__m128*)_aligned_malloc(std::max(nv,nf*3) * sizeof(__m128), 16);
#else
    m_SSEnormals = new __m128[std::max(nv,nf*3)];
    m_SSEvertices = new __m128[nv];
#endif

    #endif   

    if (nt)
    {   // got texture coordinates
        m_texCoords = new VectorR2[nt];
        m_texCoordIndices = new TupleI3[nf];
    }
	m_numTextCoords = nt;
    
    //For normal averaging.
    //Index i contains a list of all normals that are neighbors to vertex i.
    vector<int>* neighboringNormals = new vector<int>[nv];
    
    // For indicating wether or not the normal must be "fixed" (averaged).
    // For instance, if no normal is specified in the .obj file, it must be averaged.
    bool * fixNormal = new bool[std::max(nv,nf*3)];
   
    m_normalIndices = new TupleI3[nf]; // always make normals
    m_vertexIndices = new TupleI3[nf]; // always have vertices

    m_numTris = 0;
    int nvertices = 0;
    int nnormals = 0;
    int ntextures = 0;

    Matrix4x4 nctm = ctm;
    nctm.invert();
    nctm.transpose();
    
    while (fgets(line, 80, fp) != 0)
    {
        if (line[0] == 'v')
        {
            if (line[1] == 'n')
            {
                float x, y, z;
                sscanf(&line[2], "%f %f %f\n", &x, &y, &z);
                Vector3 n(x, y, z);
                m_normals[nnormals] = nctm*n;
                m_normals[nnormals].normalize();
                fixNormal[nnormals] = false;
                
                #ifdef __SSE4_1__
                m_SSEnormals[nnormals]  = _mm_set_ps(m_normals[nnormals].x, m_normals[nnormals].y, m_normals[nnormals].z, 0.0f);
                #endif   
                nnormals++;
            }
            else if (line[1] == 't')
            {
                float x, y;
                sscanf(&line[2], "%f %f\n", &x, &y);
                m_texCoords[ntextures].x = x;
                m_texCoords[ntextures].y = y;
                ntextures++;
            }
            else
            {
                float x, y, z;
                sscanf(&line[1], "%f %f %f\n", &x, &y, &z);
                Vector3 v(x, y, z);
                m_vertices[nvertices] = ctm*v;
                #ifdef __SSE4_1__
                m_SSEvertices[nvertices]  = _mm_set_ps(m_vertices[nvertices].x, m_vertices[nvertices].y, m_vertices[nvertices].z, 0.0f);
                #endif   
                nvertices++;
            }
        }
        else if (line[0] == 'f')
        {
            char s1[32], s2[32], s3[32];
            int v, t, n;
            sscanf(&line[1], "%s %s %s\n", s1, s2, s3);

            getIndices(s1, &v, &t, &n);
            m_vertexIndices[m_numTris].v[0] = v-1;
            if (n)
            {
                m_normalIndices[m_numTris].v[0] = n-1;
                neighboringNormals[v-1].push_back(n-1);
            }
            if (t)
                m_texCoordIndices[m_numTris].v[0] = t-1;
            getIndices(s2, &v, &t, &n);
            m_vertexIndices[m_numTris].v[1] = v-1;
            if (n)
            {
                m_normalIndices[m_numTris].v[1] = n-1;
                neighboringNormals[v-1].push_back(n-1);
            }
            if (t)
                m_texCoordIndices[m_numTris].v[1] = t-1;
            getIndices(s3, &v, &t, &n);
            m_vertexIndices[m_numTris].v[2] = v-1;
            if (n)
            {
                m_normalIndices[m_numTris].v[2] = n-1;
                neighboringNormals[v-1].push_back(n-1);
            }
            if (t)
                m_texCoordIndices[m_numTris].v[2] = t-1;

            if (!n)
            {   // if no normal was supplied
                Vector3 e1 = m_vertices[m_vertexIndices[m_numTris].v[1]] -
                             m_vertices[m_vertexIndices[m_numTris].v[0]];
                Vector3 e2 = m_vertices[m_vertexIndices[m_numTris].v[2]] -
                             m_vertices[m_vertexIndices[m_numTris].v[0]];

                for (int i = 0; i < 3; i++)
                {
                    m_normals[nnormals] = cross(e1, e2);
                    m_normals[nnormals].normalize();
                    fixNormal[nnormals] = true;

                    #ifdef __SSE4_1__
                    m_SSEnormals[nnormals]  = _mm_set_ps(m_normals[nnormals].x, m_normals[nnormals].y, m_normals[nnormals].z, 0.0f);
                    #endif   
                    
                    nnormals++;
                }
                
                m_normalIndices[m_numTris].v[0] = nnormals-3;
                m_normalIndices[m_numTris].v[1] = nnormals-2;
                m_normalIndices[m_numTris].v[2] = nnormals-1;
                nn++;
                
                //Add the normals as neighbors to the vertices
                neighboringNormals[m_vertexIndices[m_numTris].v[0]].push_back(m_normalIndices[m_numTris].v[0]);
                neighboringNormals[m_vertexIndices[m_numTris].v[1]].push_back(m_normalIndices[m_numTris].v[1]);
                neighboringNormals[m_vertexIndices[m_numTris].v[2]].push_back(m_normalIndices[m_numTris].v[2]);
            }

            m_numTris++;
        } //  else ignore line
    }

    //Average normals that needs to be averaged 
    for (int i = 0; i < nvertices; i++)
    {
        if (neighboringNormals[i].size() == 0) continue;
        
        Vector3 avg;
        for (int j = 0; j < neighboringNormals[i].size(); j++)
        {
            avg += m_normals[neighboringNormals[i][j]];
        }
        avg /= neighboringNormals[i].size();
        avg.normalize();

        for (int j = 0; j < neighboringNormals[i].size(); j++)
        {
            if (!fixNormal[neighboringNormals[i][j]]) continue;
            m_normals[neighboringNormals[i][j]] = avg;
            #ifdef __SSE4_1__
            m_SSEnormals[neighboringNormals[i][j]]  = _mm_set_ps(m_normals[neighboringNormals[i][j]].x, m_normals[neighboringNormals[i][j]].y, m_normals[neighboringNormals[i][j]].z, 0.0f);
            #endif
        }
    }
    
    delete [] neighboringNormals;
}

