#ifndef CSE168_SCENE_H_INCLUDED
#define CSE168_SCENE_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include "PointLight.h"
#include "BVH.h"
#include "Texture.h"
#include "PhotonMap.h"

class Camera;
class Image;

class Scene
{
public:
	Scene() : m_photonMap(PhotonsPerLightSource) { m_environment = 0; m_bgColor = Vector3(0.0f); }
    void addObject(Object* pObj)        
    { 
        if (pObj->isBounded()) m_objects.push_back(pObj);
        else m_unboundedObjects.push_back(pObj);
    }
    const Objects* objects() const      {return &m_objects;}

    void addLight(PointLight* pObj)     {m_lights.push_back(pObj);}
    const Lights* lights() const        {return &m_lights;}

    void generatePhotonMap();

    void preCalc();
    void openGL(Camera *cam);

    void raytraceImage(Camera *cam, Image *img);
    bool trace(HitInfo& minHit, const Ray& ray,
               float tMin = 0.0f, float tMax = MIRO_TMAX) const;
	bool traceScene(const Ray& ray, Vector3& shadeResult, int depth);

    void tracePhotons();
    void tracePhoton(const Vector3& position, const Vector3& direction, const Vector3& power, int depth);

	void setEnvironment(Texture* environment) { m_environment = environment; }
	Vector3 getEnvironmentMap(const Ray & ray);

    void setBgColor(Vector3 color) { m_bgColor = color; }

protected:
    Objects m_objects;
    Objects m_unboundedObjects;
    Photon_map m_photonMap;
    BVH m_bvh;
    Lights m_lights;
    Texture * m_environment; //Environment map
    Vector3 m_bgColor;       //Background color (for when environment map is not available)

    static const int PhotonsPerLightSource =
    #ifdef PHOTONS_PER_LIGHT
    PHOTONS_PER_LIGHTSOURCE;
    #else
    100000;
    #endif
};

extern Scene * g_scene;

#endif // CSE168_SCENE_H_INCLUDED
