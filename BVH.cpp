#include <algorithm>
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


bool
BVH::intersect(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    // Here you would need to traverse the BVH to perform ray-intersection
    // acceleration. For now we just intersect every object.

    bool hit = false;
    HitInfo tempMinHit;
    minHit.t = MIRO_TMAX;

	if (m_isLeaf)
	{
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
	}
	// intersect with node bounding box
	else
	{
		Component t[3];
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				t[i].Bounds[j] = (m_corners[j][i] - ray.o[i]) / ray.d[0];
			}
			if (t[i].Bounds[0] > t[i].Bounds[1])
				std::swap(t[i].Bounds[0], t[i].Bounds[0]);
		}

		Component t_sorted[3];
		if (t[0].Bounds[0] < t[1].Bounds[0] && t[0].Bounds[0] < t[2].Bounds[0])
		{
			t_sorted[0] = t[0];
			if (t[1].Bounds[0] < t[2].Bounds[0])
			{
				t_sorted[1] = t[1];
				t_sorted[2] = t[2];
			}
			else
			{
				t_sorted[1] = t[2];
				t_sorted[2] = t[1];
			}
		}
		else if (t[1].Bounds[0] < t[2].Bounds[0])
		{
			t_sorted[0] = t[1];
			if (t[0].Bounds[0] < t[2].Bounds[0])
			{
				t_sorted[1] = t[0];
				t_sorted[2] = t[2];
			}
			else
			{
				t_sorted[1] = t[2];
				t_sorted[2] = t[0];
			}
		}
		else
		{
			t_sorted[0] = t[2];
			if (t[0].Bounds[0] < t[1].Bounds[0])
			{
				t_sorted[1] = t[0];
				t_sorted[2] = t[1];
			}
			else
			{
				t_sorted[1] = t[1];
				t_sorted[2] = t[0];
			}
		}

		if (t_sorted[0].Bounds[1] < t_sorted[1].Bounds[0] && t_sorted[1].Bounds[1] < t_sorted[2].Bounds[0])
			return hit;


	}
    return hit;
}
