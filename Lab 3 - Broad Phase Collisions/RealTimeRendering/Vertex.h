#pragma once
#include "shared.h"

struct Vertex
{
	glm::vec3 position;
	glm::vec2 texCoords;

	Vertex(){}

	Vertex(glm::vec3 pos, glm::vec2 uv)
	{
		position = pos;
		texCoords = uv;
	}
};

