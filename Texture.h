#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include "Phong.h"
#include <string>
#include <iostream>
#include <FreeImage.h>
#include <noise/noisegen.h>

class Texture
{
public:
	virtual LookupCoordinates GetLookupCoordinates() const = 0;
    virtual Vector3 lookup2D(const tex_coord2d_t & coords) { return Vector3(0,0,0); } // Look up the color value for a specified position
    virtual Vector3 lookup3D(const tex_coord3d_t & coords) { return Vector3(0,0,0); } // For 3D textures
};

class Texture2D : public Texture
{
public:
	virtual LookupCoordinates GetLookupCoordinates() const {return UV; }
};

class Texture3D : public Texture
{
public:
	virtual LookupCoordinates GetLookupCoordinates() const {return UVW; }
};

class TestTexture3D : public Texture3D
{
public:
	virtual Vector3 lookup3D(const tex_coord3d_t & coords)
	{
		float x = fabs(coords.u) * 1.0, y = fabs(coords.v) * 1.0, z = fabs(coords.w) * 1.0;
		//std::cout << x << " " << y << " " << z << std::endl;
		//std::cout << Vector3((float)(x%256) / 255.0, (float)(y%256) / 255.0, (float)(z%256) / 255.0) << std::endl;
		//return Vector3(0,0,0);
		float out = 0;
		float ampl = 1;
		float freq = 1;
		for (int i = 0; i < 6; i++)
		{
			out += ampl * (noise::GradientCoherentNoise3D (x*freq, y*freq, z*freq)/2.0 + 0.5);
			ampl *= 0.7;
			freq *= 2;
		}


		//std::cout << out << std::endl;
		return Vector3(out);
		//return Vector3((float)(x%256) / 255.0, (float)(y%256) / 255.0, (float)(z%256) / 255.0);
	}
};

//Texture loaded from a file
class LoadedTexture : public Texture2D
{
public:
    LoadedTexture(std::string filename);
	~LoadedTexture();

    virtual Vector3 lookup2D(const tex_coord2d_t & coords);
protected:
	FIBITMAP* m_bitmap;
};

//Shading model that also does textures
class TexturedPhong : public Phong
{
public:
    TexturedPhong(Texture * texture,
			            const Vector3 & ka = Vector3(0),
			            const Vector3 & ks = Vector3(1),
			            const float shinyness = 1.0f,
			            const float reflect = 0,
			            const float refract = 0,
			            const float refractIndex = 1);
	virtual LookupCoordinates GetLookupCoordinates() const { return m_texture->GetLookupCoordinates(); }

    virtual Vector3 diffuse2D(const tex_coord2d_t & texture_coords) const;
    virtual Vector3 diffuse3D(const tex_coord3d_t & texture_coords) const;

protected:
    Texture * m_texture;
};

#endif
