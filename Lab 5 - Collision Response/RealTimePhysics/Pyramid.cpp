#include "Pyramid.h"

Pyramid::Pyramid(void)
{
	CreateVertexBuffer();
	CreateIndexBuffer();

	polyTexture = new Texture(GL_TEXTURE_2D, "Textures\\metal.png");
	if (!polyTexture->Load()) 
		assert(polyTexture == NULL);

	float height, width, mass;
	height = 1.0f;
	width = 1.0f;
	mass = 1.0f;

	glm::mat3 pyramidTensor = glm::mat3(
		0.36111111111111111111f, 0.0f, 0.0f,
		0.0f, 0.53333333333333333f, 0.0f,
		0.0f, 0.0f, 0.36111111111111111111f);

	rigidBody = new RigidBody(1.0f, glm::vec3(), glm::vec3(), glm::quat(), glm::vec3(), pyramidTensor);
	AABB = new AxisAlignedBoundingBox();

	rigidBody->SetGravity(true);
	rigidBody->SetDamping(true);
}

Pyramid::Pyramid(glm::vec3 initPosition, glm::vec3 initMomentum, glm::vec3 initAngularMomentum, glm::vec3 initRotation, float m, bool gravityEnabled)
{
	CreateVertexBuffer();
	CreateIndexBuffer();

	polyTexture = new Texture(GL_TEXTURE_2D, "Textures\\metal.png");
	if (!polyTexture->Load()) 
		assert(polyTexture == NULL);

	float height, width, mass;
	height = 1.0f;
	width = 1.0f;
	mass = m;

	glm::mat3 approxPyramidTensor = glm::mat3(
		(3.0f / 5.0f) * mass * (height * height) + (3.0f / 20.0f) * mass * (width * width), 0.0f, 0.0f,
		0.0f, (3.0f / 5.0f) * mass * (height * height) + (3.0f / 20.0f) * mass * (width * width), 0.0f,
		0.0f, 0.0f, (3.0f / 10.0f) * mass * (width * width));

	glm::quat initialOrientation = glm::quat(initRotation);

	rigidBody = new RigidBody(mass, initPosition, initMomentum, initialOrientation, initAngularMomentum, approxPyramidTensor);
	AABB = new AxisAlignedBoundingBox();

	rigidBody->SetGravity(gravityEnabled);
	rigidBody->SetDamping(gravityEnabled);
}

void Pyramid::CreateVertexBuffer(void)
{
	Vertex pyramidVerts[] = 
	{
		//Front
		Vertex(glm::vec3(0.0f, 0.75f, 0.0f), glm::vec2(0.5f, 0.5f), glm::vec3(0.0f, 0.5f, 0.5f)),
		Vertex(glm::vec3(1.0f, -0.25f, 1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, 0.5f, 0.5f)),
		Vertex(glm::vec3(-1.0f, -0.25f, 1.0f), glm::vec2(1.0f, 1.0f), glm::vec3(0.0f, 0.5f, 0.5f)),

		//Back
		Vertex(glm::vec3(0.0f, 0.75, 0.0f), glm::vec2(0.5f, 0.5f), glm::vec3(0.0f, 0.5f, -0.5f)),
		Vertex(glm::vec3(1.0f, -0.25f, -1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, 0.5f, -0.5f)),
		Vertex(glm::vec3(-1.0f, -0.25f, -1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(0.0f, 0.5f, -0.5f)),

		//Right
		Vertex(glm::vec3(0.0f, 0.75f, 0.0f), glm::vec2(0.5f, 0.5f), glm::vec3(-0.5f, 0.5f, 0.0f)),
		Vertex(glm::vec3(-1.0f, -0.25f, 1.0f), glm::vec2(1.0f, 1.0f), glm::vec3(-0.5f, 0.5f, 0.0f)),
		Vertex(glm::vec3(-1.0f, -0.25f, -1.0f), glm::vec2(1.0f, 0.0f), glm::vec3(-0.5f, 0.5f, 0.0f)),

		//Left
		Vertex(glm::vec3(0.0f, 0.75f, 0.0f), glm::vec2(0.5f, 0.5f), glm::vec3(0.5f, 0.5f, 0.0f)),
		Vertex(glm::vec3(1.0f, -0.25f, 1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.5f, 0.5f, 0.0f)),
		Vertex(glm::vec3(1.0f, -0.25f, -1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(0.5f, 0.5f, 0.0f)),

		//Bottom
		Vertex(glm::vec3(-1.0f, -0.25f, 1.0f), glm::vec2(1.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -0.25f, 1.0f), glm::vec2(0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
		Vertex(glm::vec3(1.0f, -0.25f, -1.0f), glm::vec2(0.0f, 1.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
		Vertex(glm::vec3(-1.0f, -0.25f, -1.0f), glm::vec2(1.0f, 1.0f), glm::vec3(0.0f, -1.0f, 0.0f))
	};

	glm::vec3 cVerts[] = 
	{
		//Top
		glm::vec3(0.0f, 0.75f, 0.0f),

		//Bottom
		glm::vec3(1.0f, -0.25f, 1.0f),
		glm::vec3(-1.0f, -0.25f, 1.0f),
		glm::vec3(-1.0f, -0.25f, -1.0f),
		glm::vec3(1.0f, -0.25f, -1.0f), 
	};

	renderVertsCount = 16;
	colliderVertsCount = 5;
	faceIndicesCount = 21;
	
	verts.reserve(renderVertsCount);
	for(int i = 0; i < renderVertsCount; i++)
		verts.push_back(pyramidVerts[i]);

	collisionVerts.reserve(colliderVertsCount);
	for(int i = 0; i < colliderVertsCount; i++)
		collisionVerts.push_back(cVerts[i]);

	const int arr[] = 
	{
		3, 0,2,1,	//Front
		3, 0,4,3,	//Back
		3, 0,3,2,	//Right
		3, 0,1,4,	//Left
		4, 1,2,3,4	//Bottom
	};

	faceIndices.reserve(faceIndicesCount);
	for(int i = 0; i < faceIndicesCount; i++)
		faceIndices.push_back(arr[i]);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(pyramidVerts), pyramidVerts, GL_STATIC_DRAW);
}

void Pyramid::CreateIndexBuffer(void)
{
	unsigned int oIndices[] = 
	{
		// Front
		0, 1, 2,
		// Back
		3, 4, 5,
		// Right
		6, 7, 8,
		// Left
		9, 10, 11,
		// Top
		12, 13, 14,
		12, 14, 15,
	};

	glGenBuffers(1, &IBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(oIndices), oIndices, GL_STATIC_DRAW);

	indices.reserve(18);
	for(int i = 0; i < 18; i++)
		indices.push_back(oIndices[i]);
}

void Pyramid::Update(float t, float dt)
{
	rigidBody->Update(t, dt);
	TransformVertices(rigidBody->GetCurrentState()->worldMat);
	AABB->GenerateAABB(transformedVerts);
}

void Pyramid::Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix, double interpolatedAlpha) //viewMatrix, glm::mat4 &projMatrix)
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
    glDrawElements(GL_TRIANGLES, 18, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}