#pragma once
#include "shared.h"

class Particle
{
public:
	Particle(void);
	glm::vec3 position, velocity, force, color, forcePressure, acceleration, oldVelocity;
	float mass, density, pressure, forceViscosity;

public:
	void SetPosition(GLfloat x, GLfloat y, GLfloat z);
	glm::vec3 GetPosition(void);

	void SetVelocity(GLfloat x, GLfloat y, GLfloat z);
	glm::vec3 GetVelocity(void);

	void Update(GLfloat elapsedTime);
};

