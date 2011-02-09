#include "TriangleMesh.h"
#include "Triangle.h"
#include "Scene.h"


TriangleMesh::TriangleMesh() :
    m_normals(0),
    m_vertices(0),
    m_texCoords(0),
    m_normalIndices(0),
    m_vertexIndices(0),
    m_texCoordIndices(0)
{

}

TriangleMesh::~TriangleMesh()
{
    delete [] m_normals;
    delete [] m_vertices;
    delete [] m_texCoords;
    delete [] m_normalIndices;
    delete [] m_vertexIndices;
    delete [] m_texCoordIndices;
}

void TriangleMesh::setVertex(int index, const Vector3 &v)
{
    m_vertices[index] = v;

    #ifdef __SSE4_1__
    m_SSEvertices[index] = _mm_set_ps(v.x, v.y, v.z, 0.0f);
    #endif
}

void TriangleMesh::setNormal(int index, const Vector3 &n)
{
    m_normals[index] = n;

    #ifdef __SSE4_1__
    m_SSEnormals[index] = _mm_set_ps(n.x, n.y, n.z, 0.0f);
    #endif
}
