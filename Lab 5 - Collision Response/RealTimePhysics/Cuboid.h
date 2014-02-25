#pragma once
#include "polyhedron.h"

class Cuboid : public Polyhedron
{
public:
	Cuboid(void);
	Cuboid(glm::vec3 position, glm::vec3 initMomentum, glm::vec3 initAngularMomentum, const float mass, const float width, const float height, const float depth, bool gravityEnabled);

	virtual void Update(float t, float dt);
	virtual void Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix, double interpolatedAlpha);

private:
	virtual void CreateVertexBuffer(const float width, const float height, const float depth);
	virtual void CreateVertexBuffer(void);
	virtual void CreateIndexBuffer(void);
};

