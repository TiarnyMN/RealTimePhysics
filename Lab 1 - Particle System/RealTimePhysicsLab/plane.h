#pragma once
#include "shared.h"

class Plane
{
public:
	//glm::vec3 GetPos(void);
	//glm::vec3 GetNorm(void);

	glm::vec3 position, normal;

	void DrawPlane(void);
	void GenerateVertices(glm::vec3 pos, glm::vec3 norm);

private:
	glm::vec3 vOne, vTwo, vThree, vFour;
};

