#ifndef CSE168_CAMERA_H_INCLUDED
#define CSE168_CAMERA_H_INCLUDED

#include "Vector3.h"
#include "Miro.h"
#include "Ray.h"

class Camera
{
public:
    Camera();
    virtual ~Camera();

    enum
    {
        RENDER_OPENGL   = 0,
        RENDER_RAYTRACE = 1
    };

    void click(Scene* pScene, Image* pImage);

    inline bool isOpenGL() const {return m_renderer == RENDER_OPENGL;}
    inline void setRenderer(int i) {m_renderer = i;}

    inline void setEye(float x, float y, float z);
    inline void setEye(const Vector3& eye);
    inline void setUp(float x, float y, float z);
    inline void setUp(const Vector3& up);
    inline void setViewDir(float x, float y, float z);
    inline void setViewDir(const Vector3& vd);
    inline void setLookAt(float x, float y, float z);
    inline void setLookAt(const Vector3& look);
    inline void setBGColor(float x, float y, float z);
    inline void setBGColor(const Vector3& color);
    inline void setFOV(float fov) {m_fov = fov; m_initialized = false;}

    inline float fov() const                {return m_fov;}
    inline const Vector3 & viewDir() const  {return m_viewDir;}
    inline const Vector3 & lookAt() const   {return m_lookAt;}
    inline const Vector3 & up() const       {return m_up;}
    inline const Vector3 & eye() const      {return m_eye;}
    inline const Vector3 & bgColor() const  {return m_bgColor;}

    Ray eyeRay(int x, int y, int imageWidth, int imageHeight, bool randomize);
    
    void drawGL();

private:

    void calcLookAt();
    bool m_initialized;
    Vector3 m_bgColor;
    int m_renderer;

    // main screen params
    Vector3 m_eye;
    Vector3 m_up;
    Vector3 m_viewDir;
    Vector3 m_lookAt;
    float m_fov;
};

extern Camera * g_camera;

//--------------------------------------------------------

inline void Camera::setEye(float x, float y, float z)
{
    m_eye.set(x, y, z);
    m_initialized = false;
}

inline void Camera::setEye(const Vector3& eye)
{
    setEye(eye.x, eye.y, eye.z);
    m_initialized = false;
}

inline void Camera::setUp(float x, float y, float z)
{
    m_up.set(x, y, z);
    m_up.normalize();
    m_initialized = false;
}

inline void Camera::setUp(const Vector3& up)
{
    setUp(up.x, up.y, up.z);
    m_initialized = false;
}

inline void Camera::setViewDir(float x, float y, float z)
{
    m_viewDir.set(x, y, z);
    m_viewDir.normalize();
    m_initialized = false;
}

inline void Camera::setViewDir(const Vector3& vd)
{
    setViewDir(vd.x, vd.y, vd.z);
    m_initialized = false;
}

inline void Camera::setLookAt(float x, float y, float z)
{
    Vector3 dir = Vector3(x, y, z) - m_eye;
    setViewDir(dir);
    m_initialized = false;
}

inline void Camera::setLookAt(const Vector3& vd)
{
    setLookAt(vd.x, vd.y, vd.z);
    m_initialized = false;
}

inline void Camera::setBGColor(float x, float y, float z)
{
    m_bgColor.set(x, y, z);
}

inline void Camera::setBGColor(const Vector3& vd)
{
    setBGColor(vd.x, vd.y, vd.z);
}

#endif // CSE168_CAMERA_H_INCLUDED
