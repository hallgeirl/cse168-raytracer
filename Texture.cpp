#include "Texture.h"
#include <FreeImage.h>

using namespace std;


LoadedTexture::LoadedTexture(std::string filename)
{
	FIBITMAP* tmp = FreeImage_Load(FIF_HDR, filename.c_str());

	m_bitmap = FreeImage_TmoDrago03(tmp); //Do tonemapping to get pixel values in a sensible range
}
LoadedTexture::~LoadedTexture()
{
	free(m_bitmap->data);
	free(m_bitmap);
}

//Aux function to retrieve pixel values and convert them to vectors.
Vector3 getPixel(FIBITMAP* img, int x, int y)
{
	RGBQUAD color;
	FreeImage_GetPixelColor(img, x, y, &color);
	return Vector3((float)color.rgbRed/255.0f,(float)color.rgbGreen/255.0f,(float)color.rgbBlue/255.0f);
}

Vector3 LoadedTexture::lookup2D(const tex_coord2d_t & texture_coords)
{
	float u = texture_coords.u, v = texture_coords.v;

	//Image dimensions
	int w = FreeImage_GetWidth(m_bitmap), h = FreeImage_GetHeight(m_bitmap);

	//Do bilinear filtering
	float px_real = (float)w*u, py_real = (float)h*v; //Point in image we really want

	int x1 = (int)px_real, x2 = x1+1;
	x1 %= w; x2 %= w;
	float x1_error = px_real - (float)x1;

	int y1 = (int)py_real, y2 = y1+1;
	y1 %= h; y2 %= h;
	float y1_error = py_real - (float)y1;

	//Get pixel values
	Vector3 f = (getPixel(m_bitmap, x1, y1) * (1-x1_error) + getPixel(m_bitmap, x2, y1) * x1_error) * (1 - y1_error) + (getPixel(m_bitmap, x1, y2) * (1-x1_error) + getPixel(m_bitmap, x2, y2) * x1_error) * y1_error;

	return f;
	//int tx = (int)((float)w*texture_coords.u) % w, ty = (int)((float)h*texture_coords.v) % h;

	//TODO: Bilinear filtering.



	//return Vector3((float)color.rgbRed/255.0f,(float)color.rgbGreen/255.0f,(float)color.rgbBlue/255.0f);
}

TexturedPhong::TexturedPhong(Texture * texture, const Vector3 & ka, const Vector3 & ks, const float shinyness, const float reflect, const float refract, const float refractIndex)
 : Phong(Vector3(0,0,0), ka, ks, shinyness, reflect, refract, refractIndex), m_texture(texture)
{

}

Vector3 TexturedPhong::diffuse2D(const tex_coord2d_t & texture_coords) const
{
    return m_texture->lookup2D(texture_coords);
}

Vector3 TexturedPhong::diffuse3D(const tex_coord3d_t & texture_coords) const
{
    return m_texture->lookup3D(texture_coords);
}
