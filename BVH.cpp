#include <algorithm>
#include "Triangle.h"
#include "SSE.h"
#include "BVH.h"
#include "Ray.h"
#include "Console.h"

using namespace std;

void getCornerPoints(Vector3 (&outCorners)[2], Objects * objs)
{
    //Find the bounds of the root node
    outCorners[0].set(infinity, infinity, infinity),
    outCorners[1].set(-infinity, -infinity, -infinity);

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
    //Find the bounds of this node 
    getCornerPoints(m_corners, objs); 

    //Check if we're done
    if (objs->size() <= OBJECTS_PER_LEAF || depth >= MAX_TREE_DEPTH)
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
        m_isLeaf = false;

        //Do a binary search for the best splitting position for each dimension. Pick the dimension with the best split.
        for (int dim = 0; dim < 3; dim++)
        {
            float current = (m_corners[1][dim] + m_corners[0][dim])/2.0f, beg = m_corners[0][dim], end = m_corners[1][dim];
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
                            left.push_back(right[i]);
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


int SSEintersectTriangles(Triangle *triangles[], int nTriangles, HitInfo& hitInfo, const Ray &ray, float tMin, float tMax)
{
    //Define some constants that are used throughout.
    static const __m128 _one = _mm_set1_ps(1.0f);
    static const __m128 _zero = _mm_setzero_ps();
    static const __m128 _minus_epsilon = _mm_set1_ps(-epsilon);
    static const __m128 _one_plus_epsilon = _mm_set1_ps(1.0f+epsilon);

    //3 dimensions, up to 4*3 floats (3 components per vertex). For each dimension, store (A1x, A2x, A3x, A4x, B1x, B2x, B3x...)
    float verts[3][12];       
    float normals[3][12]; 
    float outT[4], outP[3][4], outN[3][4];

    #pragma unroll(4)
    for (int t = 0; t < nTriangles; t++)
    {
        TriangleMesh* m = triangles[t]->getMesh();
        TriangleMesh::TupleI3 vInd = m->vIndices()[triangles[t]->getIndex()],
                              nInd = m->nIndices()[triangles[t]->getIndex()];


        //for each vertex (A, B, C)
        #pragma unroll(3)
        for (int k = 0; k < 3; k++)
        {
            //for each dimension
            #pragma unroll(3)
            for (int i = 0; i < 3; i++)
            {
                {
                    verts[i][k*4+t] = m->vertices()[vInd.v[k]][i];
                    normals[i][k*4+t] = m->normals()[nInd.v[k]][i];
                }
            }
        }
    }

    //Each __m128 contains the coordinates for 4 triangles and there is 3 dimensions.
    SSEVectorTuple3 A, B, C, nA, nB, nC, BmA, CmA, normal, rayD, rayO, RomA;
    __m128 ddotn, t, beta, gamma, alpha, P[4], N[4]; 

    //For each dimension, load the A, B and C
    #pragma unroll(3)
    for (int i = 0; i < 3; i++)
    {
        A.v[i] = _mm_loadu_ps(&verts[i][0]);
        B.v[i] = _mm_loadu_ps(&verts[i][4]);
        C.v[i] = _mm_loadu_ps(&verts[i][8]);
        nA.v[i] = _mm_loadu_ps(&normals[i][0]);
        nB.v[i] = _mm_loadu_ps(&normals[i][4]);
        nC.v[i] = _mm_loadu_ps(&normals[i][8]);
        BmA.v[i] = _mm_sub_ps(B.v[i], A.v[i]);
        CmA.v[i] = _mm_sub_ps(C.v[i], A.v[i]);
        rayO.v[i] = _mm_shuffle_ps(ray.o_SSE, ray.o_SSE, _MM_SHUFFLE(3-i,3-i,3-i,3-i));
        rayD.v[i] = _mm_sub_ps(_zero, _mm_shuffle_ps(ray.d_SSE, ray.d_SSE, _MM_SHUFFLE(3-i,3-i,3-i,3-i)));
        RomA.v[i] = _mm_sub_ps(rayO.v[i], A.v[i]);
    }

    //Do the cross product
    normal = SSEmultiCross(BmA, CmA);
    ddotn = SSEmultiDot(rayD, normal);

    //Calculate t, beta and gamma
    t = _mm_mul_ps(SSEmultiDot(RomA, normal), _mm_rcp_ps(ddotn)); 
    beta = _mm_mul_ps(SSEmultiDot(rayD, SSEmultiCross(RomA, CmA)), _mm_rcp_ps(ddotn)); 
    gamma = _mm_mul_ps(SSEmultiDot(rayD, SSEmultiCross(BmA, RomA)), _mm_rcp_ps(ddotn)); 

    //Test t, beta and gamma
    int mask = _mm_movemask_ps(_mm_and_ps(_mm_cmpgt_ps(beta, _minus_epsilon),
                          _mm_and_ps(_mm_cmpgt_ps(gamma, _minus_epsilon),
                                     _mm_and_ps(_mm_cmplt_ps(_mm_add_ps(gamma, beta), _one_plus_epsilon),
                                                _mm_and_ps(_mm_cmpgt_ps(t, _mm_set1_ps(tMin)), _mm_cmplt_ps(t, _mm_set1_ps(tMax)))))));

    if (mask == 0) return -1;
  
    
    _mm_storeu_ps(outT, t);

    //Find the lowest t > tMin
    int best = -1;
    for (int i = 0; i < nTriangles; i++)
    {
        if ((mask & (1 << i)) == 0)
            continue;

        if (best == -1 || outT[i] < outT[best]) best = i;
    }

    alpha = _mm_sub_ps(_mm_sub_ps(_one, beta), gamma);
    
    #pragma unroll(3)
    for (int i = 0; i < 3; i++)
    {
        _mm_storeu_ps(outP[i], _mm_add_ps(A.v[i], _mm_add_ps(_mm_mul_ps(beta, BmA.v[i]), _mm_mul_ps(gamma, CmA.v[i]))));
	    _mm_storeu_ps(outN[i], _mm_add_ps(_mm_mul_ps(alpha, nA.v[i]), _mm_add_ps(_mm_mul_ps(beta, nB.v[i]), _mm_mul_ps(gamma, nC.v[i]))));
    }


    hitInfo.t = outT[best];

    hitInfo.P = Vector3(outP[0][best], outP[1][best], outP[2][best]);
    hitInfo.N = Vector3(outN[0][best], outN[1][best], outN[2][best]);

    return best;
}

inline bool intersectTriangleList(Triangle* triangles[4], int nTriangles, HitInfo& minHit, const Ray& ray, float tMin, float tMax)
{
    HitInfo tempMinHit;

    int best = SSEintersectTriangles(triangles, nTriangles, tempMinHit, ray, tMin, tMax);
    if (best != -1)
    {       
        if (tempMinHit.t < minHit.t)
        {
            minHit = tempMinHit;

            //Just call intersect to get the material
            triangles[best]->intersect(minHit, ray, tMin, tMax);
             
            //Update object reference
            minHit.object = triangles[best];
        }
        return true;
    }
    return false;
}
#endif

bool
BVH::intersect(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    // Traverse the BVH to perform ray-intersection acceleration.
    bool hit = false;
    HitInfo tempMinHit;
    minHit.t = MIRO_TMAX;

	if (m_isLeaf)
	{
        //For SSE, we put up to 4 triangles in a list before testing them.
        #ifdef __SSE4_1__
        Triangle* triangles[4];
        int nTriangles = 0;
        #endif
		
        for (size_t i = 0; i < m_objects->size(); ++i)
		{
            //Build list for SSE
            #ifdef __SSE4_1__
            Triangle* tri = (Triangle*)(*m_objects)[i];
            if (tri != 0)
            {
                //If this is indeed a triangle, add it to the array
                triangles[nTriangles] = tri;
                nTriangles++;
                //If the array is full, test them
                if (nTriangles == 4)
                {
                    if (intersectTriangleList(triangles, nTriangles, minHit, ray, tMin, tMax))
                        hit = true;
                    nTriangles = 0;
                }
            }
            else
            {
            #endif
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
            #ifdef __SSE4_1__
            }
            #endif
        }
        #ifdef __SSE4_1__
        //Take care of any leftover triangles
        if (nTriangles > 0)
        {
            if (intersectTriangleList(triangles, nTriangles, minHit, ray, tMin, tMax))
                hit = true;
            nTriangles = 0;
        }
        #endif
        
		return hit;
	}

	// intersect with node bounding box
	Component t[3];
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			t[i].Bounds[j] = (m_corners[j][i] - ray.o[i]) / ray.d[i];
		}
		if (t[i].Bounds[0] > t[i].Bounds[1])
			std::swap(t[i].Bounds[0], t[i].Bounds[1]);
	}

	// sort with worst case 3 comparisons and 9 memory records
	// previous version had 2-3 comparisons and 4 memory records, but
	// significantly more branching
	if (t[0].Bounds[0] > t[1].Bounds[0])
		std::swap(t[0], t[1]);
	if (t[1].Bounds[0] > t[2].Bounds[0])
		std::swap(t[1], t[2]);
	if (t[0].Bounds[0] > t[1].Bounds[0])
		std::swap(t[0], t[1]);

	// return false if there is no overlap
	if (t[0].Bounds[1] < t[1].Bounds[0] && t[1].Bounds[1] < t[2].Bounds[0])
		return false;

	// References do not seem to conflict
	if ((*m_children)[0]->intersect(tempMinHit, ray, tMin, tMax))
	{
		minHit = tempMinHit;
		hit = true;
	}
	if ((*m_children)[1]->intersect(tempMinHit, ray, tMin, tMax))
	{
		if (tempMinHit.t < minHit.t)
		{	
			minHit = tempMinHit;
			hit = true;
		}
	}
	return hit;

}
