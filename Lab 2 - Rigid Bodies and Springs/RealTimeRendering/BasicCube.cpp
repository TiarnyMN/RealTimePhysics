#include "BasicCube.h"

BasicCube::BasicCube(void)
{
	position = glm::vec3(0.0f, 0.0f, 0.0f);
	rotation = glm::vec3(0.0f, 0.0f, 0.0f);
	scale = glm::vec3(1.0f, 1.0f, 1.0f);

	CreateVertexBuffer();
	CreateIndexBuffer();

	cubeTexture = new Texture(GL_TEXTURE_2D, "dogecoin.png");

	if (!cubeTexture->Load()) 
	{
		assert(cubeTexture == NULL);
    }
}

Vertex BasicCube::GetTransformedVertex(GLuint vertexID)
{
	if(vertexID < transformedVerts.size())
		return transformedVerts[vertexID];

	return Vertex();
}

void BasicCube::TransformVertices(glm::mat4 transformationMat)
{
	transformedVerts.clear();

	for(int i = 0; i < 24; i++)
		transformedVerts.push_back(Vertex(glm::vec3(transformationMat * glm::vec4(verts[i].position, 1.0f)), verts[i].texCoords));
}

void BasicCube::CreateVertexBuffer(void)
{
	Vertex cubeVerts[] = 
	{
		//Front
		Vertex(glm::vec3(1.0f, 1.0f, -1.0f), glm::vec2(0.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec2(1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, -1.0f), glm::vec2(0.0f, 1.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, -1.0f), glm::vec2(1.0f, 1.0f)),

		//Back
		Vertex(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec2(1.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, 1.0f), glm::vec2(0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, 1.0f), glm::vec2(1.0f, 1.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, 1.0f), glm::vec2(0.0f, 1.0f)),

		//Right
		Vertex(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec2(0.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, 1.0f), glm::vec2(1.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, -1.0f), glm::vec2(0.0f, 1.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, 1.0f), glm::vec2(1.0f, 1.0f)),

		//Left
		Vertex(glm::vec3(1.0f, 1.0f, -1.0f), glm::vec2(1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec2(0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, -1.0f), glm::vec2(1.0f, 1.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, 1.0f), glm::vec2(0.0f, 1.0f)),

		//Top
		Vertex(glm::vec3(-1.0f, 1.0f, 1.0f), glm::vec2(1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec2(0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, 1.0f, -1.0f), glm::vec2(0.0f, 1.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec2(1.0f, 1.0f)),

		//Bottom
		Vertex(glm::vec3(-1.0f, -1.0f, 1.0f), glm::vec2(1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, 1.0f), glm::vec2(0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, -1.0f), glm::vec2(0.0f, 1.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, -1.0f), glm::vec2(1.0f, 1.0f))
	};

	

	for(int i = 0; i < 24; i++)
	{
		//cubeVerts[i].texCoords.x = 1.0f - cubeVerts[i].texCoords.x;
		cubeVerts[i].texCoords.y = 1.0f - cubeVerts[i].texCoords.y;
		verts.push_back(cubeVerts[i]);
	}

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVerts), cubeVerts, GL_STATIC_DRAW);
}

void BasicCube::CreateIndexBuffer(void)
{
	unsigned int indices[] = 
	{
		// Front
		0, 1, 2,
		2, 1, 3,
		// Back
		5, 4, 7,
		7, 4, 6,
		// Right
		8, 9, 10,
		10, 9, 11,
		// Left
		13, 12, 15,
		15, 12, 14,
		// Top
		17, 16, 18,
		18, 16, 19,
		// Bottom
		22, 23, 21,
		21, 23, 20,
	};

	glGenBuffers(1, &IBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
}

void BasicCube::Update(float elapsedTimeStep)
{
	//rotation.y += 1.0f * elapsedTimeStep;
	//rotation.x += 1.0f * elapsedTimeStep;
	//rotation.z += 1.0f * elapsedTimeStep;
}

void BasicCube::Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix) //viewMatrix, glm::mat4 &projMatrix)
{
	GLint shaderID;
	glGetIntegerv(GL_CURRENT_PROGRAM, &shaderID);

	glm::mat4 scaleMat = glm::scale(scale);
	glm::mat4 rotationMat = glm::eulerAngleYXZ(rotation.y, rotation.x, rotation.z);
	glm::mat4 translationMat = glm::translate(position);

	glm::mat4 transformationMat = projMatrix * viewMatrix * translationMat * rotationMat * scaleMat;

    glUniformMatrix4fv(glGetUniformLocation(shaderID, "gWVP"), 1, GL_FALSE, &transformationMat[0][0]);

	glEnableVertexAttribArray(0);	//Why is this...?
	glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), 0);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)12);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);

	cubeTexture->Bind(GL_TEXTURE0);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}

void BasicCube::Render(glm::mat4 &transformationMat) //viewMatrix, glm::mat4 &projMatrix)
{
	GLint shaderID;
	glGetIntegerv(GL_CURRENT_PROGRAM, &shaderID);

	//glm::mat4 scaleMat = glm::scale(scale);
	//glm::mat4 rotationMat = glm::eulerAngleYXZ(rotation.y, rotation.x, rotation.z);
	//glm::mat4 translationMat = glm::translate(position);

	//glm::mat4 transformationMat = projMatrix * viewMatrix * translationMat * rotationMat * scaleMat;

    glUniformMatrix4fv(glGetUniformLocation(shaderID, "gWVP"), 1, GL_FALSE, &transformationMat[0][0]);

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), 0);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)12);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);

	cubeTexture->Bind(GL_TEXTURE0);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}