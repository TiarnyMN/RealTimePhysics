#pragma once
#include "shared.h"

struct SweepPoint
{
	int boundID;		//Owner box ID | Min/Max Flag
	float mValue;		//Min or max axis value for the current box.
	bool isMinVal;
};

class AxisAlignedBoundingBox
{
public:
	AxisAlignedBoundingBox(void);

	void GenerateAABB(vector<Vertex> transformedVertices);
	void RenderAABB(void);
	void ResetID(void);

	glm::vec3 minPoint;
	glm::vec3 maxPoint;

	SweepPoint* minPoints[3];
	SweepPoint* maxPoints[3];

	bool isColliding;
	int boundingID;

	static int boxes;
};

