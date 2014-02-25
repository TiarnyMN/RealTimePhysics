#include "texture.h"

Texture::Texture(GLenum TextureTarget, const std::string& FileName)
{
    m_textureTarget = TextureTarget;
    m_fileName      = FileName;
}

bool Texture::Load()
{	
	m_textureObj = ilutGLLoadImage((ILstring)m_fileName.c_str());
	iluFlipImage();

	if(m_textureObj == 0)
	{
		cout << "Error loading texture " << m_fileName.c_str() << "\n";
		return false;
	}

    return true;
}

void Texture::Bind(GLenum TextureUnit)
{
    glActiveTexture(TextureUnit);
    glBindTexture(m_textureTarget, m_textureObj);
}