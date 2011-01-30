#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include "Phong.h"
#include <string>
#include <FreeImage.h>
class Texture
{
public:
    virtual Vector3 lookup(const Vector3 & position) = 0; // Look up the color value for a specified position
};

class ProceduralTexture3D : public Texture
{

};

//Texture loaded from a file
class LoadedTexture : public Texture
{
public:
    LoadedTexture(std::string filename);
	~LoadedTexture();
    virtual Vector3 lookup(const Vector3 & position);
protected:
	FIBITMAP* m_bitmap;
};

//Shading model that also
class TexturedPhong : public Phong
{
public:
    TexturedPhong(Texture * texture,
			            const Vector3 & ka = Vector3(0),
			            const Vector3 & ks = Vector3(1),
			            const float shinyness = 0.f,
			            const float reflect = 0,
			            const float refract = 0,
			            const float refractIndex = 1);

    virtual Vector3 kd(const Vector3 & texture_coords) const;
protected:
    Texture * m_texture;
};

#endif
