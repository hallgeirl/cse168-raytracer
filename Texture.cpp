#include "Texture.h"
#include <FreeImage.h>
using namespace std;


LoadedTexture::LoadedTexture(std::string filename)
{
	FIBITMAP* tmp = FreeImage_Load(FIF_HDR, filename.c_str());

	m_bitmap = FreeImage_TmoDrago03(tmp);
}
LoadedTexture::~LoadedTexture()
{
	free(m_bitmap->data);
	free(m_bitmap);
}

Vector3 LoadedTexture::lookup(const Vector3 & texture_coords)
{
	int w = FreeImage_GetWidth(m_bitmap), h = FreeImage_GetHeight(m_bitmap);
	int tx = abs((int)(texture_coords.x*100)) % w, ty = (int)texture_coords.y, tz = abs((int)(texture_coords.z*100)) % h;
	//cout << FreeImage_GetImageType(m_bitmap) << endl;
	//FIRGBF color = ((FIRGBF*)(m_bitmap->data))[tz*w+tx];
	RGBQUAD color;
	FreeImage_GetPixelColor(m_bitmap, tx, tz, &color);
	//cout << color.red << " " << color.green << " " << color.blue << endl;
	//cout << tx << "\t" << tz << endl;
	//return Vector3(color.red, color.green, color.blue);

	//return Vector3(1,1,1);
	return Vector3((float)color.rgbRed/256.0f,(float)color.rgbGreen/256.0f,(float)color.rgbBlue/256.0f);
}

TexturedPhong::TexturedPhong(Texture * texture, const Vector3 & ka, const Vector3 & ks, const float shinyness, const float reflect, const float refract, const float refractIndex)
 : Phong(Vector3(0,0,0), ka, ks, shinyness, reflect, refract, refractIndex), m_texture(texture)
{

}

Vector3 TexturedPhong::kd(const Vector3 & texture_coords) const
{
    return m_texture->lookup(texture_coords);
}
