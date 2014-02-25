#pragma once
#include "shared.h"
#include "glm/gtx/euler_angles.hpp"

#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 720

class BasicCube
{
public:
	BasicCube(void);
	
	void Update(float elapsedTimeStep);
	void Render(glm::mat4 &transformationMat);
	void Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix);
	void TransformVertices(glm::mat4 transformationMat);
	Vertex GetTransformedVertex(GLuint vertexID);

private:
	GLuint VBO;
	GLuint IBO;

	glm::vec3 position;
	glm::vec3 rotation;
	glm::vec3 scale;

	Texture* cubeTexture;

	vector<Vertex> verts;
	vector<Vertex> transformedVerts;

	void CreateVertexBuffer(void);
	void CreateIndexBuffer(void);
};

