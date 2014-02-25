#pragma once
#include "shared.h"

class Particle
{
public:
	Particle(void);
	glm::vec3 position, velocity;
	float mass, age;

public:
	void SetPosition(GLfloat x, GLfloat y, GLfloat z);
	glm::vec3 GetPosition(void);

	void SetVelocity(GLfloat x, GLfloat y, GLfloat z);
	glm::vec3 GetVelocity(void);

	void Update(GLfloat elapsedTime);
};

