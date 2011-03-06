#ifndef __STATS_H__
#define __STATS_H__

class Stats
{
public:

	static int BVH_Nodes;
	static int BVH_LeafNodes;
	static int Rays;
	static int Primary_Rays;
	static int Secondary_Rays;
	static int Shadow_Rays;
	static int Photon_Bounces;
	static int Ray_Box_Intersect;
	static int Ray_Tri_Intersect;

	void PrintStats();
};

#endif