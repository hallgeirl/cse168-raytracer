#include <algorithm>
#include "SSE.h"
#include "BVH.h"
#include "Ray.h"
#include "Console.h"

using namespace std;

void getCornerPoints(Vector3 (&outCorners)[2], Objects * objs)
{
    //Find the bounds of the root node
    outCorners[0].set(infinity, infinity, infinity),
    outCorners[1].set(infinity, -infinity, -infinity);

    for (size_t i = 0; i < objs->size(); ++i)
    {
        //Ignore infinitely spanning objects like planes
        if (!(*objs)[i]->isBounded()) continue;

        Vector3 objMax = (*objs)[i]->coordsMax(),
                objMin = (*objs)[i]->coordsMin();
    
        for (int j = 0; j < 3; j++)
        {
            if (outCorners[1][j] < objMax[j])
                outCorners[1][j] = objMax[j];
            if (outCorners[0][j] > objMin[j])
                outCorners[0][j] = objMin[j];
        }
    }
}

float getArea(const Vector3 (&corners)[2])
{
    float area = corners[1][0]-corners[0][0];
    for (int i = 1; i < 3; i++)
    {
        area *= (corners[1][i]-corners[0][i]);
    }
    return area;
}

void
BVH::build(Objects * objs, int depth)
{
    // construct the bounding volume hierarchy
    m_objects = objs;
    return;
    //Find the bounds of this node 
    getCornerPoints(m_corners, objs); 

    //Check if we're done
    if (objs->size() < 4 || depth >= MAX_TREE_DEPTH)
    {
        m_objects = objs;
        m_isLeaf = true;
    }
    else
    {
        //Split the node somehow
        //and recurse into children
        float bestCost = infinity, bestPosition = infinity;
        int bestDim = 0;
        const int maxSearchDepth = 16; //Max search depth for binary search.

        //Do a binary search for the best splitting position for each dimension. Pick the dimension with the best split.
        for (int dim = 0; dim < 3; dim++)
        {
            float current = m_corners[1][dim] - m_corners[0][dim], beg = m_corners[0][dim], end = m_corners[1][dim];
            bool done = false;
            Objects left, right;

            //Put objects into the "left" and "right" node respectively
            for (int i = 0; i < objs->size(); i++)
            {
                if ((*objs)[i]->center()[dim] < current)
                    left.push_back((*objs)[i]);
                else
                    right.push_back((*objs)[i]);
            }

            for (int searchDepth = 0; searchDepth < maxSearchDepth; searchDepth++)
            {
                int n1 = left.size(), n2 = right.size();
                Vector3 cornersLeft[2], cornersRight[2];
                getCornerPoints(cornersLeft, &left);
                getCornerPoints(cornersRight, &right);

                float costLeft = n1 * getArea(cornersLeft), costRight = n2 * getArea(cornersRight);
                if (costLeft+costRight < bestCost)
                {
                    bestCost = costLeft+costRight;
                    bestDim = dim;
                    bestPosition = current;
                }
                
                //If the left node has a higher cost, we want to reduce that area. If not, we want to reduce the right area.
                if (costLeft > costRight)
                {
                    end = current;
                    current = (beg+end)/2;
                    
                    //Move the objects from left to right
                    for (int i = left.size() - 1; i >= 0; i--)
                    {
                        if (left[i]->center()[dim] > current) 
                        {
                            right.push_back(left[i]);
                            left.erase(left.begin()+i);
                        }
                    }
                }
                else
                {
                    beg = current;
                    current = (beg+end)/2;

                    //Move the objects from right to left
                    for (int i = right.size() - 1; i >= 0; i--)
                    {
                        if (right[i]->center()[dim] < current) 
                        {
                            left.push_back(left[i]);
                            right.erase(right.begin()+i);
                        }
                    }
                }
            }
        }    
        //Add child nodes
        m_children = new vector<BVH*>;
        Objects* left, * right;
        left = new Objects; right = new Objects;

        //Split the object array according to the best splitting plane we found
        for (int i = 0; i < objs->size(); i++)
        {
            if ((*objs)[i]->center()[bestDim] < bestPosition)
                left->push_back((*objs)[i]);
            else
                right->push_back((*objs)[i]);
        }
        
        for (int i = 0; i < 2; i++)
        {
            Objects* current = (i == 0 ? left : right);
            m_children->push_back(new BVH);
            (*m_children)[i]->build(current, depth+1);

            // If the new node is an internal one, free the object list since it wasn't used.
            if (!(*m_children)[i]->m_isLeaf) delete current;
        }
    }
}

#ifdef __SSE4_1__

//Do the cross product of 4 pairs of vectors
inline SSEVectorTuple3 SSEmultiCross(const SSEVectorTuple3 &a, const SSEVectorTuple3 &b)
{
    SSEVectorTuple3 out;
    //Do the cross product
    out.v[0] = _mm_sub_ps(_mm_mul_ps(a.v[1], b.v[2]), _mm_mul_ps(a.v[2], b.v[1]));
    out.v[1] = _mm_sub_ps(_mm_mul_ps(a.v[2], b.v[0]), _mm_mul_ps(a.v[0], b.v[2]));
    out.v[2] = _mm_sub_ps(_mm_mul_ps(a.v[0], b.v[1]), _mm_mul_ps(a.v[1], b.v[0]));
}

inline __m128 SSEmultiDot(const SSEVectorTuple3 &a, const SSEVectorTuple3 &b)
{
    return _mm_add_ps(_mm_mul_ps(a.v[0], b.v[0]) ,_mm_add_ps(_mm_mul_ps(a.v[1], b.v[1]), _mm_mul_ps(a.v[2], b.v[2])));
}
/*
inline bool SSEintersectTriangles(Triangle *triangles[], int nTriangles, HitInfo& hitInfo, const Ray &ray, float tMin, float tMax)
{
    float verts[3][12];   //3 dimensions, up to 4*3 vertices. For each dimension, store (A1x, A2x, A3x, A4x, B1x, B2x, B3x...)
    float normals[3][12]; 
    float outT[4], outGamma[4], outBeta[4];

    #pragma unroll(4)
    for (int t = 0; t < nTriangles; j++)
    {
        const TriangleMesh* m = triangles[i]->getMesh();
        TriangleMesh::TupleI3 vInd = m->vIndices()[triangles[t]->getIndex()],
                              nInd = m->vIndices()[triangles[t]->getIndex()]; 

        //for each dimension
        #pragma unroll(3)
        for (int i = 0; i < 3; i++)
        {
            //for each point (A, B, C)
            #pragma unroll(3)
            for (int k = 0; k < 3; k++)
            {
                verts[i][k*4+t] = m->vertices()[vInd.v[k]][i];
                normals[i][k*4+t] = m->normals()[nInd.v[k]][i];
            }
        }
    }

    //Each __m128 contains the coordinates for 4 triangles and there is 3 dimensions.
    SSEVectorTuple3 A, B, C, nA, nB, nC, BmA, CmA, normal, rayD, rayO, RomA;
    __m128 ddotn, t, beta, gamma; 

    //For each dimension, load the A, B and C
    #pragma unroll(3)
    for (int i = 0; i < 3; i++)
    {
        A.v[i] = _mm_loadu(&verts[i][0]);
        B.v[i] = _mm_loadu(&verts[i][4]);
        C.v[i] = _mm_loadu(&verts[i][8]);
        nA.v[i] = _mm_loadu(&normals[i][0]);
        nB.v[i] = _mm_loadu(&normals[i][4]);
        nC.v[i] = _mm_loadu(&normals[i][8]);
        BmA.v[i] = _mm_sub_ps(B[i], A[i]);
        CmA.v[i] = _mm_sub_ps(C[i], A[i]);
        rayO.v[i] = _mm_shuffle_ps(r.o_SSE, r.o_SSE, _MM_SHUFFLE(i,i,i,i));
        rayD.v[i] = _mm_sub_ps(_mm_setzero_ps(), _mm_shuffle_ps(r.d_SSE, r.d_SSE, _MM_SHUFFLE(i,i,i,i)));
        RomA.v[i] = _mm_sub_ps(rayO.v[i], A.v[i]);
    }

    //Do the cross product
    normal = SSEmultiCross(BmA, CmA);

    ddotn = SSEmultiDot(rayD, normal); 
    
    //Calculate t, beta and gamma
    t = _mm_mul_ps(SSEmultiDot(rayO, normal), _mm_rcp_ps(ddotn)); 
    beta = _mm_mul_ps(SSEmultiDot(rayD, SSEmultiCross(RomA, CmA)), _mm_rcp_ps(ddotn)); 
    gamma = _mm_mul_ps(SSEmultiDot(rayD, SSEmultiCross(BmA, RomA)), _mm_rcp_ps(ddotn)); 

    _mm_storeu_ps(outT, t);
    _mm_storeu_ps(outBeta, beta);
    _mm_storeu_ps(outGamma, gamma);

    //Find the lowest t > tMin
    int best = -1;
    for (int i = 0; i < nTriangles; i++)
    {
        if (outT[i] > tMin && (best == -1 || outT[i] < outT[best]))
            best = i;
    }

    //Compare to see if we have a hit 
	float betagammat[4];
    _mm_storeu_ps(betagammat, _beta_gamma_t);
    
    if (betagammat[2] < -epsilon || betagammat[3] < -epsilon || betagammat[2]+betagammat[3] > 1+epsilon || betagammat[1] < tMin || betagammat[1] > tMax) 
    {
		return false;
	}


    __m128 _beta  = _mm_shuffle_ps(_beta_gamma_t, _beta_gamma_t, _MM_SHUFFLE(2, 2, 2, 2)),
           _gamma = _mm_shuffle_ps(_beta_gamma_t, _beta_gamma_t, _MM_SHUFFLE(3, 3, 3, 3));
	__m128 _alpha = _mm_sub_ps(_mm_sub_ps(_mm_set1_ps(1.0f), _beta), _gamma);
           
    //_mm_test_all_zeros(_beta, _gamma);
    __m128 _P = _mm_add_ps(_A, _mm_add_ps(_mm_mul_ps(_beta, _BmA), _mm_mul_ps(_gamma, _CmA)));
	__m128 _N = _mm_add_ps(_mm_mul_ps(_alpha, _nA), _mm_add_ps(_mm_mul_ps(_beta, _nB), _mm_mul_ps(_gamma, _nC)));
	
    float P[4], N[4];
    
    _mm_storeu_ps(P, _P);
    _mm_storeu_ps(N, _N);

	result.P.x = P[3];
	result.P.y = P[2];
	result.P.z = P[1];

	result.N.x = N[3];
	result.N.y = N[2];
	result.N.z = N[1];
	
	result.t = betagammat[1];
    
}*/
#endif

bool
BVH::intersect(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    // Here you would need to traverse the BVH to perform ray-intersection
    // acceleration. For now we just intersect every object.

    bool hit = false;
    HitInfo tempMinHit;
    minHit.t = MIRO_TMAX;

    for (size_t i = 0; i < m_objects->size(); ++i)
    {
        if ((*m_objects)[i]->intersect(tempMinHit, ray, tMin, tMax))
        {
            hit = true;
            if (tempMinHit.t < minHit.t)
            {
            	minHit = tempMinHit;

            	//Update object reference
            	minHit.object = (*m_objects)[i];
            }
        }
    }

    return hit;
}
