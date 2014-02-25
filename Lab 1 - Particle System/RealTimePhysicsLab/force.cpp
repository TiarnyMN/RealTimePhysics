#include "force.h"

Force::Force(void)
{
	position = glm::vec3(0.0f, -50.0f, 0.0f);
	//impulse = glm::vec3(0.0f, 50.0f, 0.0f);	//Is "impulse" used?
	strength = 500.0f;
	
	forceEnabled = true;
	forcePushing = false;
}

Force::Force(glm::vec3 pos, float stren)
{
	position = pos;
	//impulse = imp;
	strength = stren;
	
	forceEnabled = true;
	forcePushing = false;
}

