#include "SSE.h"
#include "Triangle.h"
#include "TriangleMesh.h"
#include "Ray.h"
#include "Vector3.h"
#include "Utility.h"

#ifdef STATS
#include "Stats.h"
#endif

using namespace std;

namespace
{
	float Det(float a, float b, float c, float d)
	{
		return a * d - b * c;
	}
}

Triangle::Triangle(TriangleMesh * m, unsigned int i)
{
    setMesh(m);
    setIndex(i);
}


Triangle::~Triangle()
{

}

void Triangle::preCalc()
{
    //Find the min and max boundaries
    updateMinMax();
}


Vector3 Triangle::center() const
{
    TriangleMesh::TupleI3 ti3 = m_mesh->vIndices()[m_index];
    Vector3 verts[3] = {m_mesh->vertices()[ti3.v[0]], m_mesh->vertices()[ti3.v[1]], m_mesh->vertices()[ti3.v[2]]};
    Vector3 BmA = verts[1]-verts[0], CmA = verts[2]-verts[0];
 
    return verts[0] + BmA/3 + CmA/3;
}

float Triangle::getArea(const Vector3& lightPos)
{
	//need to cache center
	Vector3 l_dir = (lightPos - center()).normalize();
	Vector3 ut, vt;
	getTangents(l_dir, ut, vt);

    TriangleMesh::TupleI3 ti3 = m_mesh->vIndices()[m_index];
    Vector3 verts[3] = {m_mesh->vertices()[ti3.v[0]], m_mesh->vertices()[ti3.v[1]], m_mesh->vertices()[ti3.v[2]]};
	float minBase = infinity;
	float minHeight = infinity;
	float maxBase = -infinity;
	float maxHeight = -infinity;

	for (int i = 0; i < 3; ++i)
	{
		float projU = dot(ut, verts[i]);
		float projV = dot(vt, verts[i]);
		
		if (projU < minBase)
			minBase = projU;
		if (projU > maxBase)
			maxBase = projU;

		if (projV < minHeight)
			minHeight = projV;
		if (projV > maxHeight)
			maxHeight = projV;
	}

	return (0.5f * (maxHeight-minHeight) * (maxBase-minBase));
}

Vector3 Triangle::samplePosition() const
{
	float u1 = frand();
	float beta = 1.0f - sqrt(u1);
	float gamma = sqrt(frand()) * u1;	

	TriangleMesh::TupleI3 ti3 = m_mesh->vIndices()[m_index];
    Vector3 verts[3] = {m_mesh->vertices()[ti3.v[0]], m_mesh->vertices()[ti3.v[1]], m_mesh->vertices()[ti3.v[2]]};
    Vector3 BmA = (verts[1]-verts[0]), CmA = (verts[2]-verts[0]);

	return (verts[0] + beta * BmA + gamma * CmA);
}


void Triangle::updateMinMax()
{
    if (m_mesh != 0 && m_index < m_mesh->numTris() && m_index >= 0)
    {
        TriangleMesh::TupleI3 ti3 = m_mesh->vIndices()[m_index];
        Vector3 verts[3] = {m_mesh->vertices()[ti3.v[0]], m_mesh->vertices()[ti3.v[1]], m_mesh->vertices()[ti3.v[2]]};

        m_cachedMin = verts[0];
        m_cachedMax = verts[0];


        for (int i = 1; i < 3; i++)
        {
            if (verts[i].x < m_cachedMin.x) m_cachedMin.x = verts[i].x;
            if (verts[i].y < m_cachedMin.y) m_cachedMin.y = verts[i].y;
            if (verts[i].z < m_cachedMin.z) m_cachedMin.z = verts[i].z;
            if (verts[i].x > m_cachedMax.x) m_cachedMax.x = verts[i].x;
            if (verts[i].y > m_cachedMax.y) m_cachedMax.y = verts[i].y;
            if (verts[i].z > m_cachedMax.z) m_cachedMax.z = verts[i].z;
        }
    }
}


void
Triangle::renderGL()
{
    TriangleMesh::TupleI3 ti3 = m_mesh->vIndices()[m_index];
    const Vector3 & v0 = m_mesh->vertices()[ti3.v[0]]; //vertex a of triangle
    const Vector3 & v1 = m_mesh->vertices()[ti3.v[1]]; //vertex b of triangle
    const Vector3 & v2 = m_mesh->vertices()[ti3.v[2]]; //vertex c of triangle

    glBegin(GL_TRIANGLES);
        glVertex3f(v0.x, v0.y, v0.z);
        glVertex3f(v1.x, v1.y, v1.z);
        glVertex3f(v2.x, v2.y, v2.z);
    glEnd();
}

bool
Triangle::intersect(HitInfo& result, const Ray& r,float tMin, float tMax)
{
#ifndef __SSE4_1__
    TriangleMesh::TupleI3 ti3 = m_mesh->vIndices()[m_index];
    TriangleMesh::TupleI3 ni3 = m_mesh->nIndices()[m_index];

    const Vector3 & A = m_mesh->vertices()[ti3.v[0]];
    const Vector3 & B = m_mesh->vertices()[ti3.v[1]];
    const Vector3 & C = m_mesh->vertices()[ti3.v[2]];
    const Vector3 & nA = m_mesh->normals()[ni3.v[0]];
    const Vector3 & nB = m_mesh->normals()[ni3.v[1]];
    const Vector3 & nC = m_mesh->normals()[ni3.v[2]];

    Vector3 BmA = B-A, CmA = C-A;
    Vector3 normal = cross(BmA, CmA);
    float ddotn = (dot(-r.d, normal));

    float t = dot(r.o-A, normal) / ddotn;
    float beta = dot(-r.d, cross(r.o-A, CmA)) / ddotn;
    float gamma = dot(-r.d, cross(BmA, r.o-A)) / ddotn;

    if (beta < -epsilon || gamma < -epsilon || beta+gamma > 1+epsilon || t < tMin || t > tMax) return false;

    result.P = A + beta*BmA + gamma*CmA;
    result.t = t;
    result.N = (1-beta-gamma)*nA + beta*nB + gamma*nC;

#endif

    result.material = m_material;

    return true;
}

//We planned to use the triangle texturing for the leaf and stem, but could not find a suitable texture and did not have time to make one.  
tex_coord2d_t Triangle::toUVCoordinates(const Vector3 & xyz) const
{
	if (m_mesh->numTextCoords() == 0)
		return tex_coord2d_t();

    TriangleMesh::TupleI3 vi3 = m_mesh->vIndices()[m_index];
    TriangleMesh::TupleI3 ti3 = m_mesh->tIndices()[m_index];

    const Vector3 & vA = m_mesh->vertices()[vi3.v[0]];
    const Vector3 & vB = m_mesh->vertices()[vi3.v[1]];
    const Vector3 & vC = m_mesh->vertices()[vi3.v[2]];
    const VectorR2 & tA = m_mesh->texCoords()[ti3.v[0]];
    const VectorR2 & tB = m_mesh->texCoords()[ti3.v[1]];
    const VectorR2 & tC = m_mesh->texCoords()[ti3.v[2]];

	// if no texture coords exist
	//return tex_coord2d_t(xyz.x, xyz.z);

	//discard largest normal component
    Vector3 BmA = vB-vA, CmA = vC-vA;
    Vector3 normal = cross(BmA, CmA);
	int i = 0;
	int j = 1;

	if (normal.x > normal.z)
		i = 2;
	else if (normal.y > normal.z)
		j = 2;

	// convert to float arrays to alleviate calculations
	float p[3] = {xyz.x-vA.x, xyz.y-vA.y, xyz.z-vA.z};
	float B[3] = {vB.x-vA.x, vB.y-vA.y, vB.z-vA.z};
	float C[3] = {vC.x-vA.x, vC.y-vA.y, vC.z-vA.z};

	//Cramer's rule to determine barycentric coords
	float detPC = Det(p[i], C[i], p[j], C[j]);
	float detBP = Det(B[i], p[i], B[j], p[j]);
	float detBC = Det(B[i], C[i], B[j], C[j]);

	float beta = max(detPC / detBC, 0.f);
	float gamma = max(detBP / detBC, 0.f);
	//this shouldn't happen...but just in case
	float alpha = max(1 - (beta + gamma), 0.f);

	//interpolate vertices texture coords
	tex_coord2d_t UV;
	UV.u = alpha * tA.x + beta * tB.x + gamma * tC.x;
	UV.v = alpha * tA.y + beta * tB.y + gamma * tC.y;

	return UV;
}
