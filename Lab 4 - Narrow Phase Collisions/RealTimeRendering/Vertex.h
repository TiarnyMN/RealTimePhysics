#pragma once
#include "shared.h"

struct Vertex
{
	glm::vec3 position;
	glm::vec2 texCoords;
	glm::vec3 normal;

	Vertex(){}

	Vertex(glm::vec3 pos, glm::vec2 uv, glm::vec3 norm)
	{
		position = pos;
		texCoords = uv;
		normal = norm;
	}
};

