#include <IL/il.h>
#include "Texture.h"

using namespace std;


LoadedTexture::LoadedTexture(std::string filename)
{
    
}

Vector3 LoadedTexture::lookup(const Vector3 & texture_coords)
{
    return Vector3(1, 0, 1);
}

TexturedPhong::TexturedPhong(Texture * texture, const Vector3 & ka, const Vector3 & ks, const float shinyness, const float reflect, const float refract, const float refractIndex)
 : Phong(Vector3(0,0,0), ka, ks, shinyness, reflect, refract, refractIndex), m_texture(texture)
{
    
}

const Vector3 & TexturedPhong::kd(const Vector3 & texture_coords) const
{
    return m_texture->lookup(texture_coords);
}
