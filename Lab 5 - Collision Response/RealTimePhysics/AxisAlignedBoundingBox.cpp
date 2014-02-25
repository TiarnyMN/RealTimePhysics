#include "AxisAlignedBoundingBox.h"

int AxisAlignedBoundingBox::boxes = 0;

AxisAlignedBoundingBox::AxisAlignedBoundingBox(void)
{
	boundingID = boxes++;

	for(int i = 0; i < 3; i++)
	{
		minPoints[i] = new SweepPoint();
		minPoints[i]->isMinVal = true;
	}

	for(int i = 0; i < 3; i++)
	{
		maxPoints[i] = new SweepPoint();
		maxPoints[i]->isMinVal = false;
	}
}

void AxisAlignedBoundingBox::GenerateAABB(vector<Vertex> transformedVertices)
{
	minPoint = glm::vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	maxPoint = glm::vec3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

	for(int i = 0; i < transformedVertices.size(); i++)
	{
		Vertex curVert = transformedVertices[i];

		minPoint.x = glm::min(minPoint.x, curVert.position.x);
		minPoint.y = glm::min(minPoint.y, curVert.position.y);
		minPoint.z = glm::min(minPoint.z, curVert.position.z);

		maxPoint.x = glm::max(maxPoint.x, curVert.position.x);
		maxPoint.y = glm::max(maxPoint.y, curVert.position.y);
		maxPoint.z = glm::max(maxPoint.z, curVert.position.z);
	}

	minPoints[0]->mValue = minPoint.x;
	minPoints[1]->mValue = minPoint.y;
	minPoints[2]->mValue = minPoint.z;

	maxPoints[0]->mValue = maxPoint.x;
	maxPoints[1]->mValue = maxPoint.y;
	maxPoints[2]->mValue = maxPoint.z;
}

void AxisAlignedBoundingBox::RenderAABB(void)
{
	if(isColliding)
		glColor3f(1.0f, 0.0f, 0.0f);
	else
		glColor3f(0.0f, 1.0f, 0.0f);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_QUADS);

		glVertex3f(maxPoint.x, maxPoint.y, maxPoint.z);
		glVertex3f(minPoint.x, maxPoint.y, maxPoint.z);
		glVertex3f(minPoint.x, minPoint.y, maxPoint.z);
		glVertex3f(maxPoint.x, minPoint.y, maxPoint.z);

		glVertex3f(maxPoint.x, maxPoint.y, minPoint.z);
		glVertex3f(minPoint.x, maxPoint.y, minPoint.z);
		glVertex3f(minPoint.x, minPoint.y, minPoint.z);
		glVertex3f(maxPoint.x, minPoint.y, minPoint.z);

		glVertex3f(maxPoint.x, maxPoint.y, maxPoint.z);
		glVertex3f(maxPoint.x, minPoint.y, maxPoint.z);
		glVertex3f(maxPoint.x, minPoint.y, minPoint.z);
		glVertex3f(maxPoint.x, maxPoint.y, minPoint.z);

		glVertex3f(minPoint.x, maxPoint.y, maxPoint.z);
		glVertex3f(minPoint.x, minPoint.y, maxPoint.z);
		glVertex3f(minPoint.x, minPoint.y, minPoint.z);
		glVertex3f(minPoint.x, maxPoint.y, minPoint.z);

		glVertex3f(maxPoint.x, maxPoint.y, maxPoint.z);
		glVertex3f(minPoint.x, maxPoint.y, maxPoint.z);
		glVertex3f(minPoint.x, maxPoint.y, minPoint.z);
		glVertex3f(maxPoint.x, maxPoint.y, minPoint.z);

		glVertex3f(maxPoint.x, minPoint.y, maxPoint.z);
		glVertex3f(minPoint.x, minPoint.y, maxPoint.z);
		glVertex3f(minPoint.x, minPoint.y, minPoint.z);
		glVertex3f(maxPoint.x, minPoint.y, minPoint.z);

	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}
