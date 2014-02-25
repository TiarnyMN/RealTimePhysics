#include "Plane.h"

Plane::Plane(void)
{
	CreateVertexBuffer();
	CreateIndexBuffer();

	polyTexture = new Texture(GL_TEXTURE_2D, "Textures\\floor.png");
	if (!polyTexture->Load()) 
		assert(polyTexture == NULL);

	rigidBody = new RigidBody();
	AABB = new AxisAlignedBoundingBox();

	rigidBody->SetGravity(false);
	rigidBody->SetDamping(false);
}

Plane::Plane(glm::vec3 initialPosition, const float size)
{
	CreateVertexBuffer(size);
	CreateIndexBuffer();

	polyTexture = new Texture(GL_TEXTURE_2D, "Textures\\floor.png");
	if (!polyTexture->Load()) 
		assert(polyTexture == NULL);

	rigidBody = new RigidBody(0.0f, initialPosition, glm::vec3(), glm::quat(), glm::vec3(), glm::mat3());
	AABB = new AxisAlignedBoundingBox();

	rigidBody->SetGravity(false);
	rigidBody->SetDamping(false);
}

void Plane::CreateVertexBuffer(void)
{
	CreateVertexBuffer(50.0f);
}

void Plane::CreateVertexBuffer(const float size)
{
	Vertex planeVerts[] = 
	{
		//Plane
		Vertex(glm::vec3(size, 0.0f, size), glm::vec2(10.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)),
		Vertex(glm::vec3(-size, 0.0f, size), glm::vec2(0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)),
		Vertex(glm::vec3(-size, 0.0f, -size), glm::vec2(0.0f, 10.0f), glm::vec3(0.0f, 1.0f, 0.0f)),
		Vertex(glm::vec3(size, 0.0f, -size), glm::vec2(10.0f, 10.0f), glm::vec3(0.0f, 1.0f, 0.0f)),
	};

	glm::vec3 cVerts[] = 
	{
		//Plane
		glm::vec3(size, 0.0f, size),
		glm::vec3(-size, 0.0f, size),
		glm::vec3(-size, 0.0f, -size),
		glm::vec3(size, 0.0f, -size)
	};

	renderVertsCount = 4;
	colliderVertsCount = 4;
	faceIndicesCount = 5;
	
	verts.reserve(renderVertsCount);
	for(int i = 0; i < renderVertsCount; i++)
	{
		planeVerts[i].texCoords.y = 1.0f - planeVerts[i].texCoords.y;
		verts.push_back(planeVerts[i]);
	}

	collisionVerts.reserve(colliderVertsCount);
	for(int i = 0; i < colliderVertsCount; i++)
		collisionVerts.push_back(cVerts[i]);

	const int arr[] = 
	{
		4, 0,3,2,1,	//Plane
	};

	faceIndices.reserve(faceIndicesCount);
	for(int i = 0; i < faceIndicesCount; i++)
		faceIndices.push_back(arr[i]);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(planeVerts), planeVerts, GL_STATIC_DRAW);
}

void Plane::CreateIndexBuffer(void)
{
	unsigned int indices[] = 
	{
		// Plane
		0, 1, 3,
		3, 1, 2,
	};

	glGenBuffers(1, &IBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
}

void Plane::Update(float t, float dt)
{
	rigidBody->Update(t, dt);

	TransformVertices(rigidBody->GetCurrentState()->worldMat);
	AABB->GenerateAABB(transformedVerts);
}

void Plane::Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix, double interpolatedAlpha) //viewMatrix, glm::mat4 &projMatrix)
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
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}
