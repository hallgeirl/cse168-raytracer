#ifndef CSE168_BVH_H_INCLUDED
#define CSE168_BVH_H_INCLUDED

#include <vector>
#include <limits>
#include "Miro.h"
#include "Object.h"

struct Component
{
	float Bounds[2];
};

//Represents a node in the bounding volume hierarchy
class BVH
{
public:
/*   BVH() { m_corners[0] = Vector3(infinity); m_corners[1] = Vector3(infinity);}
    BVH(Vector3 corners[2]) { m_corners[0] = corners[0]; m_corners[1] = corners[1]; }*/
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
    static const int OBJECTS_PER_LEAF = 12; //For SSE, it is beneficial to have more objects in each leaf.
    #else
    static const int OBJECTS_PER_LEAF = 4;
    #endif
};

#endif // CSE168_BVH_H_INCLUDED
