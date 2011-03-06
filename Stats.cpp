#include <iostream>
#include <cstdio>
#include "Stats.h"

int Stats::BVH_Nodes = 0;
int Stats::BVH_LeafNodes = 0;
int Stats::Rays = 0;
int Stats::Primary_Rays = 0;
int Stats::Secondary_Rays = 0;
int Stats::Shadow_Rays = 0;
int Stats::Photon_Bounces = 0;
int Stats::Ray_Box_Intersect = 0;
int Stats::Ray_Tri_Intersect = 0;

void Stats::PrintStats()
{
	printf("\n~~Ray Tracer Stats~~\n\n");
	printf("BVH Nodes: %d\n", BVH_Nodes);
	printf("BVH Leaf Nodes: %d\n", BVH_LeafNodes);
	printf("Total Rays: %d\n", Primary_Rays + Secondary_Rays + Shadow_Rays);
	printf("Primary Rays: %d\n", Primary_Rays);
	printf("Secondary Rays: %d\n", Secondary_Rays);
	printf("Shadow Rays: %d\n", Shadow_Rays);
	printf("Average Photon Bounces: %f\n", (float)Photon_Bounces/100000.0f);
	printf("Ray-Box Intersections: %d\n", Ray_Box_Intersect);
	printf("Ray-Triangle Intersections: %d\n", Ray_Tri_Intersect);
}
