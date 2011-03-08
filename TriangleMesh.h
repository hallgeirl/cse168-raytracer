#ifndef CSE168_TRIANGLE_MESH_H_INCLUDED
#define CSE168_TRIANGLE_MESH_H_INCLUDED

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif
#include "Matrix4x4.h"

class TriangleMesh
{
public:
    TriangleMesh();
    ~TriangleMesh();

    // load from an OBJ file
    bool load(const char* file, const Matrix4x4& ctm = Matrix4x4());

    // for single triangles
    void createSingleTriangle();
    inline void setV1(const Vector3& v) { setVertex(0, v); }
    inline void setV2(const Vector3& v) { setVertex(1, v); }
    inline void setV3(const Vector3& v) { setVertex(2, v); }
    inline void setN1(const Vector3& n) { setNormal(0, n); }
    inline void setN2(const Vector3& n) { setNormal(1, n); }
    inline void setN3(const Vector3& n) { setNormal(2, n); }

    struct TupleI3
    {
        unsigned int v[3];
    };


    Vector3* vertices()     {return m_vertices;}
    Vector3* normals()      {return m_normals;}
    VectorR2* texCoords()   {return m_texCoords;}
    
    #ifdef __SSE4_1__
    __m128* SSEvertices()     {return m_SSEvertices;}
    __m128* SSEnormals()      {return m_SSEnormals;}
    #endif
    
    TupleI3* vIndices()     {return m_vertexIndices;}
    TupleI3* nIndices()     {return m_normalIndices;}
	TupleI3* tIndices()		{return m_texCoordIndices;}
    int numTris()           {return m_numTris;}
    int numTextCoords()		{return m_numTextCoords;}

protected:
    void setVertex(int index, const Vector3 &v);
    void setNormal(int index, const Vector3 &n);

    void loadObj(FILE* fp, const Matrix4x4& ctm);

    Vector3* m_normals;
    Vector3* m_vertices;
    #ifdef __SSE4_1__
    //We use these with the SSE intersection code since _mm_set_ps is terribly slow to do for each test.
    __m128* m_SSEnormals;
    __m128* m_SSEvertices;
    #endif
    VectorR2* m_texCoords;
    
    TupleI3* m_normalIndices;
    TupleI3* m_vertexIndices;
    TupleI3* m_texCoordIndices;
    unsigned int m_numVertices;
    unsigned int m_numTris;
	unsigned int m_numTextCoords;
};


#endif // CSE168_TRIANGLE_MESH_H_INCLUDED
