#include <iostream>
#include <cmath>
#include "Sphere.h"
#include "Ray.h"
#include "Console.h"

using namespace std;

Sphere::Sphere()
{
}

Sphere::~Sphere()
{
}

void
Sphere::renderGL()
{
    glColor3f(1, 1, 1);
    glPushMatrix();
    glTranslatef(m_center.x, m_center.y, m_center.z);
    glutWireSphere(m_radius, 20, 20);
    glPopMatrix();
}

bool
Sphere::intersect(HitInfo& result, const Ray& ray,
                  float tMin, float tMax)
{
    const Vector3 toO = ray.o - m_center;

    const float a = ray.d.length2();
    const float b = dot(2*ray.d, toO);
    const float c = toO.length2() - m_radius*m_radius;

    const float discrim = b*b-4.0f*a*c;

    if (discrim < 0)
        return false;   // quadratic equation would yield imaginary numbers

    const float sqrt_discrim = sqrt(discrim);

    // solve the quadratic equation
    const float t[2] = {(-b-sqrt_discrim)/(2.0f*a), (-b+sqrt_discrim)/(2.0f*a)};

    // since we know that discrim >= 0, t[0] < t{1]
    // return the t closest to us that is within range
    if ((t[0] > tMin) && (t[0] < tMax))
    {
        result.t = t[0];
    }
    else if((t[1] > tMin) && (t[1] < tMax))
    {
        result.t = t[1];
    }
    else
    {
        // neither of the solutions are in the required range
        return false;
    }

    result.P = ray.o + result.t*ray.d;
    result.N = (result.P-m_center);
    result.N.normalize();
    result.material = this->m_material;

    return true;
}

Vector3 Sphere::samplePosition(const Vector3& lightPos) const
{
    Vector3 tang1, tang2, l = lightPos - m_center;
    l.normalize();
    getTangents(l, tang1, tang2);
    
    VectorR2 discSample = sampleDisc(m_radius);

    return m_center + tang1 * discSample.x + tang2 * discSample.y;
}

// Sphere mapping
tex_coord2d_t Sphere::toUVCoordinates(const Vector3 & xyz) const
{
	//Get vector from xyz to center
	Vector3 dir = xyz - m_center;
	dir.normalize();

	tex_coord2d_t coords;
 	coords.u = (atan2(dir.x, dir.z)) / (2.0f * PI) + 0.5;

	coords.v = (std::max(-1.0f, std::min(1.0f, std::asin(dir.y)))) / PI + 0.5;

	return coords;
}
