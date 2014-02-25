#pragma once
#include "Polyhedron.h"

class BasicCube: public Polyhedron
{
public:
	BasicCube(void);
	BasicCube(glm::vec3 position, glm::vec3 initMomentum, glm::vec3 initAngularMomentum, float m, bool gravityEnabled);

	virtual void Update(float t, float dt);
	virtual void Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix, double interpolatedAlpha);

private:
	virtual void CreateVertexBuffer(void);
	virtual void CreateIndexBuffer(void);
};

