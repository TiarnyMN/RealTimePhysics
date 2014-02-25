#pragma once
#include "polyhedron.h"
class Plane : public Polyhedron
{
public:
	Plane(void);
	Plane(glm::vec3 initialPosition, const float size);

	virtual void Update(float t, float dt);
	virtual void Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix, double interpolatedAlpha);

private:
	virtual void CreateVertexBuffer(void);
	virtual void CreateVertexBuffer(const float size);
	virtual void CreateIndexBuffer(void);
};

