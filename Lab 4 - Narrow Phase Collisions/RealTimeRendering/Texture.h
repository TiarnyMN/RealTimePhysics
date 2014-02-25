#pragma once

#ifdef WIN32
#include <Windows.h>
#endif

#include "shared.h"

#include "DevIL\include\IL\ilu.h"
#include "DevIL\include\IL\il.h"
#include "DevIL\include\IL\ilut.h"

class Texture
{
public:
    Texture(GLenum TextureTarget, const std::string& FileName);

    bool Load();
    void Bind(GLenum TextureUnit);

private:
    std::string m_fileName;
    GLenum m_textureTarget;
    GLuint m_textureObj;
};

