#include "Triangle.h"
#include "TriangleMesh.h"
#include "Ray.h"


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

    Vector3 BmA = B-A, CmA = C-A;
    Vector3 normal = cross(BmA, CmA);
    float ddotn = (dot(-r.d, normal));

    float t = dot(r.o-A, normal) / ddotn;
    float beta = dot(-r.d, cross(r.o-A, CmA)) / ddotn;
    float gamma = dot(-r.d, cross(BmA, r.o-A)) / ddotn;

    if (beta < 0 || gamma < 0 || beta+gamma > 1 || t < tMin || t > tMax) return false;

    result.P = A + beta*BmA + gamma*CmA;
    result.t = t;
    result.N = (1-beta-gamma)*nA + beta*nB + gamma*nC;
    result.material = m_material;

    return true;
}
