#pragma once
#define ONE_RAD_DEG (2.0f * 3.14159265358979323846264338327950288f) / 360.0f

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <glm.hpp>
#include <gtx/transform.hpp>
#include <gtc/type_ptr.hpp>
#include <iostream>
#include <vector>

#include "particle.h"
#include "camera.h"
#include "shader.h"
#include "force.h"
#include "plane.h"

using namespace std;