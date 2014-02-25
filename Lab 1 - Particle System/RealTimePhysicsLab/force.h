#pragma once
#include "shared.h"

class Force
{
public:
	Force(void);
	Force(glm::vec3 pos, float stren);

	glm::vec3 position;
	//glm::vec3 impulse;
	float strength;
	
	bool forceEnabled;
	bool forcePushing;
};

