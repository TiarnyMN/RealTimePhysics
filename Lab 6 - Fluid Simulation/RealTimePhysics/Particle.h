#pragma once
#include "shared.h"

using namespace std;

class Particle
{
public:
	Particle(float m);
	glm::vec3 position, velocity, force, color, forcePressure, forceViscosity, forceSurface, forceInternal, forceExternal, acceleration, oldVelocity;
	float mass, inverseMass, density, pressure;

	vector<Particle*> neighbours;

public:
	void SetPosition(GLfloat x, GLfloat y, GLfloat z);
	glm::vec3 GetPosition(void);

	void SetVelocity(GLfloat x, GLfloat y, GLfloat z);
	glm::vec3 GetVelocity(void);

	void Update(GLfloat elapsedTime);
};

