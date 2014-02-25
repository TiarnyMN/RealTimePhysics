#pragma once
#include "shared.h"
#include "RigidBody.h"
#include "AxisAlignedBoundingBox.h"
#include "glm/gtx/euler_angles.hpp"

class Polyhedron
{
public:
	//Polyhedron(void);

	virtual void Update(float t, float dt) = 0;
	virtual void Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix, double interpolatedAlpha) = 0;
	virtual void RenderNormals(void);
	
	void TransformVertices(glm::mat4 transformationMat);
	void RenderAABB(void);

	BodyState* GetCurrentState(void);
	AxisAlignedBoundingBox* GetBoundingBox(void);

	vector<Vertex> GetTransformedVertices(void);
	vector<glm::vec3> GetTransformedCollisionVertices(void);
	vector<int> GetFaceIndices(void);

	void SetGravity(const bool enabled);
	void SetDamping(const bool enabled);

protected:
	GLuint VBO;
	GLuint IBO;

	RigidBody* rigidBody;
	AxisAlignedBoundingBox* AABB;

	glm::vec3 scale;

	Texture* polyTexture;

	vector<Vertex> verts;
	vector<Vertex> transformedVerts;
	vector<int> faceIndices;
	vector<int> indices;
	
	vector<glm::vec3> collisionVerts;
	vector<glm::vec3> transformedCollisionVerts;

	int renderVertsCount;
	int colliderVertsCount;
	int faceIndicesCount;

	virtual void CreateVertexBuffer(void) = 0;
	virtual void CreateIndexBuffer(void) = 0;
};

