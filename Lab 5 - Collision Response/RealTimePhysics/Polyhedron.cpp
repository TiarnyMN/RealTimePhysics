#include "Polyhedron.h"

BodyState* Polyhedron::GetCurrentState(void)
{
	return rigidBody->GetCurrentState();
}

AxisAlignedBoundingBox* Polyhedron::GetBoundingBox(void)
{
	return AABB;
}

vector<int> Polyhedron::GetFaceIndices(void)
{
	return faceIndices;
}

vector<Vertex> Polyhedron::GetTransformedVertices(void)
{
	return transformedVerts;
}

vector<glm::vec3> Polyhedron::GetTransformedCollisionVertices(void)
{
	return transformedCollisionVerts;
}

void Polyhedron::RenderAABB(void)
{
	AABB->RenderAABB();
}

void Polyhedron::SetGravity(const bool enabled)
{
	rigidBody->SetGravity(enabled);
}

void Polyhedron::SetDamping(const bool enabled)
{
	rigidBody->SetDamping(enabled);
}

void Polyhedron::TransformVertices(glm::mat4 transformationMat)
{
	transformedVerts.clear();

	for(int i = 0; i < renderVertsCount; i++)
		transformedVerts.push_back(Vertex(glm::vec3(transformationMat * glm::vec4(verts[i].position, 1.0f)), verts[i].texCoords, glm::normalize(glm::vec3(transformationMat * glm::vec4(verts[i].normal, 0.0f)))));

	transformedCollisionVerts.clear();

	for(int i = 0; i < colliderVertsCount; i++)
		transformedCollisionVerts.push_back(glm::vec3(transformationMat * glm::vec4(collisionVerts[i], 1.0f)));
}

void Polyhedron::RenderNormals(void)
{
	if(transformedCollisionVerts.size() == 0)
		return;

	const float normalSize = 0.75f;

	int faceID = 0;
	while(faceID < faceIndices.size())
	{
		int vertsInFace = faceIndices[faceID];
		++faceID;

		glm::vec3 centrePoint = transformedCollisionVerts[faceIndices[faceID]] + transformedCollisionVerts[faceIndices[faceID + 1]] + transformedCollisionVerts[faceIndices[faceID + 2]];
		centrePoint = centrePoint / 3.0f;

		glm::vec3 lhs = transformedCollisionVerts[faceIndices[faceID]] - transformedCollisionVerts[faceIndices[faceID + 1]];
		glm::vec3 rhs = transformedCollisionVerts[faceIndices[faceID]] - transformedCollisionVerts[faceIndices[faceID + 2]];

		glm::vec3 normal = glm::cross(lhs, rhs);

		if(glm::length(normal) > 0.0f)
			normal = glm::normalize(normal);

		glBegin(GL_LINES);
		glVertex3f(centrePoint.x, centrePoint.y, centrePoint.z);
		glVertex3f(centrePoint.x + (normal.x * normalSize), centrePoint.y + (normal.y * normalSize), centrePoint.z + (normal.z * normalSize));
		glEnd();

		faceID += vertsInFace;
	}
}
