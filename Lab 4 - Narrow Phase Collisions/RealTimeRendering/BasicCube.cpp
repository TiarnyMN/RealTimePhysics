#include "BasicCube.h"

BasicCube::BasicCube(void)
{
	CreateVertexBuffer();
	CreateIndexBuffer();

	cubeTexture = new Texture(GL_TEXTURE_2D, "dogecoin.png");
	if (!cubeTexture->Load()) 
		assert(cubeTexture == NULL);

	rigidBody = new RigidBody();
	AABB = new AxisAlignedBoundingBox();

	rigidBody->SetGravity(false);
	rigidBody->SetDamping(false);
}

BasicCube::BasicCube(glm::vec3 initPosition, glm::vec3 initMomentum, glm::vec3 initAngularMomentum)
{
	CreateVertexBuffer();
	CreateIndexBuffer();

	cubeTexture = new Texture(GL_TEXTURE_2D, "dogecoin.png");
	if (!cubeTexture->Load()) 
		assert(cubeTexture == NULL);

	float height, depth, width, size, mass;
	height = depth = width = size = 2.0f;
	mass = 1.0f;

	glm::mat3 approxCubeTensor = glm::mat3(
		(1.0f /12.0f) * mass * (height * height + depth * depth), 0.0f, 0.0f,
		0.0f, (1.0f / 12.0f) * mass * (width * width + depth * depth), 0.0f,
		0.0f, 0.0f, (1.0f / 12.0f) * mass * (width * width + height * height));

	rigidBody = new RigidBody(mass, initPosition, initMomentum, glm::quat(), initAngularMomentum, approxCubeTensor);
	AABB = new AxisAlignedBoundingBox();

	rigidBody->SetGravity(false);
	rigidBody->SetDamping(false);
}

BodyState* BasicCube::GetCurrentState(void)
{
	return rigidBody->GetCurrentState();
}

AxisAlignedBoundingBox* BasicCube::GetBoundingBox(void)
{
	return AABB;
}

//Vertex BasicCube::GetTransformedVertex(GLuint vertexID)
//{
//	if(vertexID < transformedVerts.size())
//		return transformedVerts[vertexID];
//
//	return Vertex();
//}

vector<Vertex> BasicCube::GetTransformedVertices(void)
{
	return transformedVerts;
}

vector<glm::vec3> BasicCube::GetTransformedCollisionVertices(void)
{
	return transformedCollisionVerts;
}

void BasicCube::TransformVertices(glm::mat4 transformationMat)
{
	transformedVerts.clear();

	for(int i = 0; i < 24; i++)
		transformedVerts.push_back(Vertex(glm::vec3(transformationMat * glm::vec4(verts[i].position, 1.0f)), verts[i].texCoords, glm::normalize(glm::vec3(transformationMat * glm::vec4(verts[i].normal, 0.0f)))));

	transformedCollisionVerts.clear();

	for(int i = 0; i < 8; i++)
		transformedCollisionVerts.push_back(glm::vec3(transformationMat * glm::vec4(collisionVerts[i], 1.0f)));
}

void BasicCube::CreateVertexBuffer(void)
{
	Vertex cubeVerts[] = 
	{
		//Front
		Vertex(glm::vec3(1.0f, 1.0f, -1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec2(1.0f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, -1.0f), glm::vec2(1.0f, 1.0f), glm::vec3(0.0f, 0.0f, -1.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, -1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, 0.0f, -1.0f)),

		//Back
		Vertex(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec2(1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, 1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, 1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, 0.0f, 1.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, 1.0f), glm::vec2(1.0f, 1.0f), glm::vec3(0.0f, 0.0f, 1.0f)),

		//Right
		Vertex(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(-1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, 1.0f), glm::vec2(1.0f, 0.0f), glm::vec3(-1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, 1.0f), glm::vec2(1.0f, 1.0f), glm::vec3(-1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, -1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(-1.0f, 0.0f, 0.0f)),	

		//Left
		Vertex(glm::vec3(1.0f, 1.0f, -1.0f), glm::vec2(1.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, 1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, -1.0f), glm::vec2(1.0f, 1.0f), glm::vec3(1.0f, 0.0f, 0.0f)),		

		//Top
		Vertex(glm::vec3(-1.0f, 1.0f, 1.0f), glm::vec2(1.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, 1.0f, -1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, 1.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec2(1.0f, 1.0f), glm::vec3(0.0f, 1.0f, 0.0f)),

		//Bottom
		Vertex(glm::vec3(-1.0f, -1.0f, 1.0f), glm::vec2(1.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, 1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, -1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, -1.0f), glm::vec2(1.0f, 1.0f), glm::vec3(0.0f, -1.0f, 0.0f))
	};

	glm::vec3 cVerts[] = 
	{
		//Front
		glm::vec3(1.0f, 1.0f, -1.0f),
		glm::vec3(-1.0f, 1.0f, -1.0f),
		glm::vec3(-1.0f, -1.0f, -1.0f),
		glm::vec3(1.0f, -1.0f, -1.0f),

		//Back
		glm::vec3(1.0f, 1.0f, 1.0f),
		glm::vec3(-1.0f, 1.0f, 1.0f),
		glm::vec3(-1.0f, -1.0f, 1.0f),
		glm::vec3(1.0f, -1.0f, 1.0f), 
	};

	
	for(int i = 0; i < 24; i++)
	{
		//cubeVerts[i].texCoords.x = 1.0f - cubeVerts[i].texCoords.x;
		cubeVerts[i].texCoords.y = 1.0f - cubeVerts[i].texCoords.y;
		verts.push_back(cubeVerts[i]);
	}

	for(int i = 0; i < 8; i++)
		collisionVerts.push_back(cVerts[i]);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVerts), cubeVerts, GL_STATIC_DRAW);
}

void BasicCube::CreateIndexBuffer(void)
{
	

	

	unsigned int indices[] = 
	{
		// Front
		0, 1, 3,
		3, 1, 2,
		// Back
		5, 4, 6,
		6, 4, 7,
		// Right
		8, 9, 11,
		11, 9, 10,
		// Left
		13, 12, 14,
		14, 12, 15,
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

void BasicCube::Update(float t, float dt)
{
	rigidBody->Update(t, dt);
	TransformVertices(rigidBody->GetCurrentState()->worldMat);
	AABB->GenerateAABB(transformedVerts);
}

void BasicCube::Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix, double interpolatedAlpha) //viewMatrix, glm::mat4 &projMatrix)
{
	GLint shaderID;
	glGetIntegerv(GL_CURRENT_PROGRAM, &shaderID);

	glm::mat4 scaleMat = glm::scale(scale);
	glm::mat4 transformationMat = projMatrix * viewMatrix * rigidBody->GetInterpolatedState(interpolatedAlpha).worldMat;// * scaleMat;

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

void BasicCube::RenderAABB(void)
{
	AABB->RenderAABB();
}

//void BasicCube::Render(glm::mat4 &transformationMat) //viewMatrix, glm::mat4 &projMatrix)
//{
//	GLint shaderID;
//	glGetIntegerv(GL_CURRENT_PROGRAM, &shaderID);
//
//    glUniformMatrix4fv(glGetUniformLocation(shaderID, "gWVP"), 1, GL_FALSE, &transformationMat[0][0]);
//
//	glEnableVertexAttribArray(0);
//	glEnableVertexAttribArray(1);
//
//    glBindBuffer(GL_ARRAY_BUFFER, VBO);
//    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), 0);
//	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)12);
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
//
//	cubeTexture->Bind(GL_TEXTURE0);
//    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
//
//    glDisableVertexAttribArray(0);
//	glDisableVertexAttribArray(1);
//}