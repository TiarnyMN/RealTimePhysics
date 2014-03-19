#include "Particle.h"

Particle::Particle(void)
{
	mass = 1.0f;
	velocity = glm::vec3();
}

void Particle::SetPosition(GLfloat x, GLfloat y, GLfloat z)
{
	position = glm::vec3(x, y, z);
}

glm::vec3 Particle::GetPosition(void)
{
	return position;
}

void Particle::SetVelocity(GLfloat x, GLfloat y, GLfloat z)
{
	velocity = glm::vec3(x, y, z);
}

glm::vec3 Particle::GetVelocity(void)
{
	return velocity;
}

void Particle::Update(GLfloat elapsedTime)
{
	position += velocity * elapsedTime;
	velocity += velocity * elapsedTime;
}
