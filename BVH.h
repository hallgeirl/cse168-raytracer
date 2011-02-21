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
    SSEVectorTuple3 A, BmA, CmA, normal, nA, nB, nC;
};
#endif

typedef float Corner[4];


//Represents a node in the bounding volume hierarchy
class BVH
{
public:
    BVH() { m_corners[0][0] = infinity; } 
    void build(Objects * objs, int depth = 0);

    bool intersect(HitInfo& result, const Ray& ray,
                   float tMin = 0.0f, float tMax = MIRO_TMAX) const;
    bool intersectChildren(HitInfo& result, const Ray& ray,
                   float tMin, float tMax) const;
protected:
    union
    {
        std::vector<BVH*> * m_children; //Child nodes of this BVH, which are also BVHs. Only applicable for inner nodes.
        Objects * m_objects;          //Objects contained in the BVH. Only applicable for child nodes. 
    };

    union
    {
        Corner m_corners[2]; //The min and max corner of the box. 
#ifdef __SSE4_1__
        __m128 m_corners_SSE[2];
#endif
    };

    bool m_isLeaf;
    static const int MAX_TREE_DEPTH = 32;
    #ifdef __SSE4_1__
    std::vector<SSETriangleCache> * m_triangleCache;
    //For SSE, it is beneficial to have more objects in each leaf. There's a sweet spot between having too many leaf objects, and having enough leaf objects so that we don't have too many half-empty vectors.
    static const int OBJECTS_PER_LEAF = 8;     
    #else
    static const int OBJECTS_PER_LEAF = 4;
    #endif
};

#endif // CSE168_BVH_H_INCLUDED
