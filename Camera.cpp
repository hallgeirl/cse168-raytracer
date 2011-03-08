#include <stdio.h>
#include <stdlib.h>
#include "Utility.h"
#include "Miro.h"
#include "Camera.h"
#include "Image.h"
#include "Scene.h"
#include "Console.h"
#include "OpenGL.h"

Camera * g_camera = 0;

//static bool firstRayTrace = true;

const float HalfDegToRad = DegToRad/2.0f;

Camera::Camera() :
    m_bgColor(0,0,0),
    m_renderer(RENDER_OPENGL),
    m_eye(0,0,0),
    m_up(0,1,0),
    m_viewDir(0,0,-1),
    m_lookAt(FLT_MAX, FLT_MAX, FLT_MAX),
    m_fov((45.)*(PI/180.))
{
    m_initialized = false;
    calcLookAt();
}


Camera::~Camera()
{

}


void
Camera::click(Scene* pScene, Image* pImage)
{
    calcLookAt();
    static bool firstRayTrace = false;

    if (m_renderer == RENDER_OPENGL)
    {
        glDrawBuffer(GL_BACK);
        pScene->openGL(this);
        firstRayTrace = true;
    }
    else if (m_renderer == RENDER_RAYTRACE)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glDrawBuffer(GL_FRONT);
        #ifndef NO_GFX
        if (firstRayTrace)
        {
        #endif
            pImage->clear(bgColor());
            pScene->raytraceImage(this, g_image);
            firstRayTrace = false;
        #ifndef NO_GFX
        }
        #endif

        g_image->draw();
    }
}


void
Camera::calcLookAt()
{
    // this is true when a "lookat" is not used in the config file
    if (m_lookAt.x != FLT_MAX)
    {
        setLookAt(m_lookAt);
        m_lookAt.set(FLT_MAX, FLT_MAX, FLT_MAX);
    }
}


void
Camera::drawGL()
{
    // set up the screen with our camera parameters
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fov(), g_image->width()/(float)g_image->height(),
                   0.01, 10000);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    Vector3 vCenter = eye() + viewDir();
    gluLookAt(eye().x, eye().y, eye().z,
              vCenter.x, vCenter.y, vCenter.z,
              up().x, up().y, up().z);
}


Ray
Camera::eyeRay(int x, int y, int imageWidth, int imageHeight, bool randomize)
{
    static bool initialized = false;
    static Vector3 wDir, uDir, vDir;
    static float top, left, bottom, right, aspectRatio;
    if (!initialized)
    {
		initialized = true;

        // first compute the camera coordinate system
        // wDir = e - (e+m_viewDir) = -m_vView
        wDir = Vector3(-m_viewDir).normalize();
        uDir = cross(m_up, wDir).normalize();
        vDir = cross(wDir, uDir);

        // next find the corners of the image plane in camera space
        aspectRatio = (float)imageWidth/(float)imageHeight;
        top     = tan(m_fov*HalfDegToRad);
        right   = aspectRatio*top;
        bottom  = -top;
        left    = -right;
    }

    float dx = 0.5, dy = 0.5;

    if (randomize)
    {
        dx = frand();
        dy = frand();
    }

	#ifdef DOF

	//randomize eye location around circle of confusion
    VectorR2 discSample = sampleDisc(DOF_APERTURE);

	Vector3 new_eye = m_eye + (discSample.x*uDir + discSample.y*vDir);

	Vector3 new_viewDir = m_eye + m_viewDir * DOF_FOCUS_PLANE - new_eye;

	//need to recalulate view plane
    Vector3 localwDir = Vector3(-new_viewDir).normalize();
    //uDir = cross(m_up, wDir).normalize();
    //vDir = cross(wDir, uDir);

	#else
    Vector3 new_eye = m_eye;
    Vector3 localwDir = wDir;
    // need stratified sampling
    // transform x and y into camera space
    // -----------------------------------
	#endif

    const float imPlaneUPos = left   + (right - left)*(((float)x+dx)/(float)imageWidth);
    const float imPlaneVPos = bottom + (top - bottom)*(((float)y+dy)/(float)imageHeight);

    return Ray(new_eye, (imPlaneUPos*uDir + imPlaneVPos*vDir - localwDir).normalize());
}
