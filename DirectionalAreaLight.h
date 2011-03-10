#ifndef _FOCUSEDAREALIGHT_H_
#define _FOCUSEDAREALIGHT_H_
#include "SquareLight.h"
#include <iostream>

//Circular light with one direction
class DirectionalAreaLight : public SquareLight
{
public:
    DirectionalAreaLight(float radius = 1) { m_radius = radius; }
   
    virtual float getRadius() { return m_radius; } 

	virtual float GetLightRatio(float objArea, const Vector3& objCenter) const 
	{
		return (objArea / (PI * m_radius * m_radius));
	}

	virtual Vector3 samplePhotonOrigin(int sampleNumber = 0, int totalSamples = 1) const  
    {
        VectorR2 discSample = sampleDisc(m_radius);
        return m_position + (discSample.x*m_tangent1 + discSample.y*m_tangent2);
    }
   
    virtual Vector3 getLightDirection(const Vector3 &origin, const Vector3 &directionOf) const
    {
        //The light direction is always the negative normal
        return -m_normal;
    }
   
    virtual Vector3 samplePhotonDirection() const
    {
        return m_normal;
    }
    
protected:
    float m_radius;
};


#endif
