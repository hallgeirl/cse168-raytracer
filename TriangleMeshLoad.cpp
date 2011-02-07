#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "TriangleMesh.h"
#include "Console.h"

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

    m_vertexIndices[0].x = 0;
    m_vertexIndices[0].y = 1;
    m_vertexIndices[0].z = 2;

    m_normalIndices[0].x = 0;
    m_normalIndices[0].y = 1;
    m_normalIndices[0].z = 2;

    m_texCoordIndices[0].x = 0;
    m_texCoordIndices[0].y = 1;
    m_texCoordIndices[0].z = 2;

    #ifdef __SSE4_1__
    m_SSEvertices = new __m128[3];
    m_SSEnormals = new __m128[3];
    
    for (int i = 0; i < 3; i++)
    {
        m_SSEvertices[i] = _mm_set_ps(m_vertices[i].x, m_vertices[i].y, m_vertices[i].z, 0.0f);
        m_SSEnormals[i]  = _mm_set_ps(m_normals[i].x, m_normals[i].y, m_normals[i].z, 0.0f);
    }
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


    m_normals = new Vector3[std::max(nv,nf)];
    m_vertices = new Vector3[nv];
    
    #ifdef __SSE4_1__
    m_SSEnormals = new __m128[std::max(nv,nf)];
    m_SSEvertices = new __m128[nv];
    #endif   

    if (nt)
    {   // got texture coordinates
        m_texCoords = new VectorR2[nt];
        m_texCoordIndices = new TupleI3[nf];
    }
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
            m_vertexIndices[m_numTris].x = v-1;
            if (n)
                m_normalIndices[m_numTris].x = n-1;
            if (t)
                m_texCoordIndices[m_numTris].x = t-1;
            getIndices(s2, &v, &t, &n);
            m_vertexIndices[m_numTris].y = v-1;
            if (n)
                m_normalIndices[m_numTris].y = n-1;
            if (t)
                m_texCoordIndices[m_numTris].y = t-1;
            getIndices(s3, &v, &t, &n);
            m_vertexIndices[m_numTris].z = v-1;
            if (n)
                m_normalIndices[m_numTris].z = n-1;
            if (t)
                m_texCoordIndices[m_numTris].z = t-1;

            if (!n)
            {   // if no normal was supplied
                Vector3 e1 = m_vertices[m_vertexIndices[m_numTris].y] -
                             m_vertices[m_vertexIndices[m_numTris].x];
                Vector3 e2 = m_vertices[m_vertexIndices[m_numTris].z] -
                             m_vertices[m_vertexIndices[m_numTris].x];

                m_normals[nn] = cross(e1, e2);
                m_normals[nn].normalize();
                
                #ifdef __SSE4_1__
                m_SSEnormals[nn]  = _mm_set_ps(m_normals[nn].x, m_normals[nn].y, m_normals[nn].z, 0.0f);
                #endif   
                
                m_normalIndices[nn].x = nn;
                m_normalIndices[nn].y = nn;
                m_normalIndices[nn].z = nn;
                nn++;
            }

            m_numTris++;
        } //  else ignore line
    }
}

