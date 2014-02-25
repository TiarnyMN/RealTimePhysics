#pragma once

#include <stdio.h>
#ifdef WIN32
#define _USE_MATH_DEFINES 
#include <cmath>
#else
#include <math.h>
#endif

#define ToRadian(x) (float)(((x) * M_PI / 180.0f))
#define ToDegree(x) (float)(((x) * 180.0f / M_PI))

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <glm.hpp>
#include <gtx/transform.hpp>
#include <gtc/type_ptr.hpp>
#include <iostream>
#include <vector>

#include "Vertex.h"
#include "Texture.h"

using namespace std;

// The Wonder Function
static glm::vec3 operator*(glm::mat4 mat, glm::vec3 param)
{
	return glm::vec3(mat * glm::vec4(param, 1.0f));
}
