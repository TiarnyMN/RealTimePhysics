#include "BasicCube.h"

BasicCube::BasicCube(void)
{
	CreateVertexBuffer();
	CreateIndexBuffer();

	polyTexture = new Texture(GL_TEXTURE_2D, "Textures\\futurecrate.png");
	if (!polyTexture->Load()) 
		assert(polyTexture == NULL);

	rigidBody = new RigidBody();
	AABB = new AxisAlignedBoundingBox();

	rigidBody->SetGravity(true);
	rigidBody->SetDamping(true);
}

BasicCube::BasicCube(glm::vec3 initPosition, glm::vec3 initMomentum, glm::vec3 initAngularMomentum, float m, bool gravityEnabled)
{
	CreateVertexBuffer();
	CreateIndexBuffer();

	polyTexture = new Texture(GL_TEXTURE_2D, "Textures\\futurecrate.png");
	if (!polyTexture->Load()) 
		assert(polyTexture == NULL);

	float height, depth, width, size, mass;
	height = depth = width = size = 2.0f;
	mass = m;

	glm::mat3 approxCubeTensor = glm::mat3(
		(1.0f /12.0f) * mass * (height * height + depth * depth), 0.0f, 0.0f,
		0.0f, (1.0f / 12.0f) * mass * (width * width + depth * depth), 0.0f,
		0.0f, 0.0f, (1.0f / 12.0f) * mass * (width * width + height * height));

	rigidBody = new RigidBody(mass, initPosition, initMomentum, glm::quat(), initAngularMomentum, approxCubeTensor);
	AABB = new AxisAlignedBoundingBox();

	rigidBody->SetGravity(gravityEnabled);
	rigidBody->SetDamping(gravityEnabled);
}

void BasicCube::CreateVertexBuffer(void)
{
	Vertex cubeVerts[] = 
	{
		//Front
		Vertex(glm::vec3(1.0f, 1.0f, -1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec2(0.5f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, -1.0f), glm::vec2(0.5f, 0.5f), glm::vec3(0.0f, 0.0f, -1.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, -1.0f), glm::vec2(0.0f, 0.5f), glm::vec3(0.0f, 0.0f, -1.0f)),

		//Back
		Vertex(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec2(0.5f, 0.5f), glm::vec3(0.0f, 0.0f, 1.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, 1.0f), glm::vec2(0.0f, 0.5f), glm::vec3(0.0f, 0.0f, 1.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, 1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, 0.0f, 1.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, 1.0f), glm::vec2(0.5f, 1.0f), glm::vec3(0.0f, 0.0f, 1.0f)),

		//Right
		Vertex(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec2(0.5f, 0.5f), glm::vec3(-1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, 1.0f), glm::vec2(0.0f, 0.5f), glm::vec3(-1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, 1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(-1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, -1.0f, -1.0f), glm::vec2(0.5f, 1.0f), glm::vec3(-1.0f, 0.0f, 0.0f)),	

		//Left
		Vertex(glm::vec3(1.0f, 1.0f, -1.0f), glm::vec2(0.5f, 0.5f), glm::vec3(1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec2(0.0f, 0.5f), glm::vec3(1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, 1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(1.0f, 0.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, -1.0f), glm::vec2(0.5f, 1.0f), glm::vec3(1.0f, 0.0f, 0.0f)),		

		//Top
		Vertex(glm::vec3(-1.0f, 1.0f, 1.0f), glm::vec2(1.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec2(0.5f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, 1.0f, -1.0f), glm::vec2(0.5f, 0.5f), glm::vec3(0.0f, 1.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec2(1.0f, 0.5f), glm::vec3(0.0f, 1.0f, 0.0f)),

		//Bottom
		Vertex(glm::vec3(-1.0f, -1.0f, 1.0f), glm::vec2(1.0f, 0.5f), glm::vec3(0.0f, -1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, 1.0f), glm::vec2(0.5f, 0.5f), glm::vec3(0.0f, -1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -1.0f, -1.0f), glm::vec2(0.5f, 1.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
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

	renderVertsCount = 24;
	colliderVertsCount = 8;
	faceIndicesCount = 30;
	
	verts.reserve(renderVertsCount);
	for(int i = 0; i < renderVertsCount; i++)
		verts.push_back(cubeVerts[i]);

	collisionVerts.reserve(colliderVertsCount);
	for(int i = 0; i < colliderVertsCount; i++)
		collisionVerts.push_back(cVerts[i]);

	const int arr[] = 
	{
		4, 0,3,2,1,	//Front
		4, 4,5,6,7,	//Back
		4, 1,2,6,5,	//Right
		4, 0,4,7,3,	//Left
		4, 5,4,0,1,	//Top
		4, 6,2,3,7	//Bottom
	};

	faceIndices.reserve(faceIndicesCount);
	for(int i = 0; i < faceIndicesCount; i++)
		faceIndices.push_back(arr[i]);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVerts), cubeVerts, GL_STATIC_DRAW);
}

void BasicCube::CreateIndexBuffer(void)
{
	unsigned int oIndices[] = 
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
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(oIndices), oIndices, GL_STATIC_DRAW);

	indices.reserve(36);
	for(int i = 0; i < 36; i++)
		indices.push_back(oIndices[i]);
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

	polyTexture->Bind(GL_TEXTURE0);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}