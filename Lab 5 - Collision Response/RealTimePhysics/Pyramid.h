#pragma once
#include "Polyhedron.h"

class Pyramid: public Polyhedron
{
public:
	Pyramid(void);
	Pyramid(glm::vec3 position, glm::vec3 initMomentum, glm::vec3 initAngularMomentum, glm::vec3 initRotation, float m, bool gravityEnabled);

	virtual void Update(float t, float dt);
	virtual void Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix, double interpolatedAlpha);

private:
	virtual void CreateVertexBuffer(void);
	virtual void CreateIndexBuffer(void);
};
