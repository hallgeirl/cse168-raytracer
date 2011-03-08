#ifndef _SQUARELIGHT_H_
#define _SQUARELIGHT_H_
#include "PointLight.h"
#include <iostream>

class SquareLight : public PointLight
{
public:
    SquareLight() { m_normal = Vector3(0,1,0); }

    void setNormal(Vector3 n) 
    { 
        m_normal = n;
    }
    
    Vector3 getNormal()
    {
        return m_normal;
    }

    void setDimensions(float width, float height) { m_dimensions[0] = width; m_dimensions[1] = height; }
    
    virtual Vector3 samplePhotonOrigin(int sampleNumber = 0, int totalSamples = 1) const  
    {
        //Take samples within a subdivided rectangle. For simplicity we assume that the light is square so we have nxn cells.
        //First find the cell dimensions
        float sideLength = sqrt((float)totalSamples);
        float du = m_dimensions[0] / sideLength;
        float dv = m_dimensions[1] / sideLength;

        //Sample index in grid
        int sx = sampleNumber % int(sideLength);
        int sy = sampleNumber / int(sideLength);

        float u = (du*frand()) + sx * du - m_dimensions[0]/2.0f;
        float v = (dv*frand()) + sy * dv - m_dimensions[1]/2.0f;

        return m_position + u*m_tangent1 + v*m_tangent2;
    }

    virtual Vector3 samplePhotonDirection() const
    {
        //bias to the light normal
        float phi = asin(sqrt(frand()));
        float theta = 2.0f * PI * (frand());

        return alignHemisphereToVector(m_normal, theta, phi);
    }

    virtual void preCalc()
    {
        getTangents(m_normal, m_tangent1, m_tangent2);
    }
    
protected:
    Vector3 m_normal; Vector3 m_tangent1, m_tangent2;
    float m_dimensions[2];
};


#endif
