#ifndef CSE168_POINTLIGHT_H_INCLUDED
#define CSE168_POINTLIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"

class PointLight
{
public:
    void setPosition(const Vector3& v)  {m_position = v;}
    void setColor(const Vector3& v)     {m_color = v;}
    void setWattage(float f)            {m_wattage = f;}
    
    float wattage() const               {return m_wattage;}
    const Vector3 & color() const       {return m_color;}
    const Vector3& position() const     {return m_position;}

    void preCalc() {} // use this if you need to

protected:
    Vector3 m_position;
    Vector3 m_color;
    float m_wattage;
};

typedef std::vector<PointLight*> Lights;

#endif // CSE168_POINTLIGHT_H_INCLUDED
