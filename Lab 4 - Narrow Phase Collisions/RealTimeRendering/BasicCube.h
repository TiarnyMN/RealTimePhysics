#pragma once
#include "shared.h"
#include "RigidBody.h"
#include "AxisAlignedBoundingBox.h"
#include "glm/gtx/euler_angles.hpp"

#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 720

class BasicCube
{
public:
	BasicCube(void);
	BasicCube(glm::vec3 position, glm::vec3 initMomentum, glm::vec3 initAngularMomentum);

	void Update(float t, float dt);
	//void Render(glm::mat4 &transformationMat);
	void Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix, double interpolatedAlpha);
	void RenderAABB(void);
	void TransformVertices(glm::mat4 transformationMat);

	BodyState* GetCurrentState(void);
	AxisAlignedBoundingBox* GetBoundingBox(void);

	vector<Vertex> GetTransformedVertices(void);
	vector<glm::vec3> GetTransformedCollisionVertices(void);

private:
	GLuint VBO;
	GLuint IBO;

	RigidBody* rigidBody;
	AxisAlignedBoundingBox* AABB;

	glm::vec3 scale;

	Texture* cubeTexture;

	vector<Vertex> verts;
	vector<Vertex> transformedVerts;
	
	vector<glm::vec3> collisionVerts;
	vector<glm::vec3> transformedCollisionVerts;

	void CreateVertexBuffer(void);
	void CreateIndexBuffer(void);
};

