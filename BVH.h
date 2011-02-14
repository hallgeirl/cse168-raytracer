#ifndef CSE168_BVH_H_INCLUDED
#define CSE168_BVH_H_INCLUDED

#include <vector>
#include <limits>
#include "SSE.h"
#include "Miro.h"
#include "Object.h"

struct Component
{
	float Bounds[2];
};

#ifdef __SSE4_1__
struct SSETriangleCache
{
    SSETriangleCache() { nTriangles = 0; }
    Triangle *triangles[4];
    int nTriangles;
    SSEVectorTuple3 A, B, C, nA, nB, nC, BmA, CmA, normal;
};
#endif


//Represents a node in the bounding volume hierarchy
class BVH
{
public:
    void build(Objects * objs, int depth = 0);

    bool intersect(HitInfo& result, const Ray& ray,
                   float tMin = 0.0f, float tMax = MIRO_TMAX) const;
protected:
    union
    {
        std::vector<BVH*> * m_children; //Child nodes of this BVH, which are also BVHs. Only applicable for inner nodes.
        Objects * m_objects;          //Objects contained in the BVH. Only applicable for child nodes. 
    };

    Vector3 m_corners[2]; //The min and max corner of the box.    
    bool m_isLeaf;
    static const int MAX_TREE_DEPTH = 32;
    #ifdef __SSE4_1__
    std::vector<SSETriangleCache> * m_triangleCache;
    static const int OBJECTS_PER_LEAF = 12; //For SSE, it is beneficial to have more objects in each leaf.
    #else
    static const int OBJECTS_PER_LEAF = 4;
    #endif
};

#endif // CSE168_BVH_H_INCLUDED
