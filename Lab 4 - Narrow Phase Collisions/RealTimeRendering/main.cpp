#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <GL/glew.h>
#include <GL/freeglut.h>

#include "shared.h"
#include "shader.h"
#include "BasicCube.h"
#include "SimpleCamera.h"
#include "RigidBody.h"
#include "AxisAlignedBoundingBox.h"

#include <algorithm>
#include <vector>

#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 720

struct Projection
{
	Projection(float mi, float ma)
	{
		min = mi;
		max = ma;
	}

	float min;
	float max;
};

struct ContactManifold
{
	vector<glm::vec3> verts;
};

struct CollisionFace
{
	vector<int> vertexIDs;
	glm::vec3 faceNormal;
};

struct FaceCollection
{
	int collisionFaceID;
	vector<CollisionFace> vertexFaces;
	vector<glm::vec3> verts;
};

enum DemoMode
{
	World,
	Collision,
	Paused,
};

SimpleCamera* simpleCamera = NULL;

int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;

Shader* bodyShader;
Shader* lineShader;

vector<BasicCube*> sceneBodies;
BasicCube* curSelectedCube;
ContactManifold* contactM;

double accumulator = 0.0f;
double interpolationAlpha = 0.0f;

const int rigidBodyCount = 2;

vector<std::list<AxisAlignedBoundingBox*>::iterator> activeColliders;
vector<vector<int>> collisionTable;

DemoMode curDemoMode = DemoMode::World;

vector<AxisAlignedBoundingBox*> xAxisStart;
vector<AxisAlignedBoundingBox*> yAxisStart;
vector<AxisAlignedBoundingBox*> zAxisStart;

vector<AxisAlignedBoundingBox*> xAxisEnd;
vector<AxisAlignedBoundingBox*> yAxisEnd;
vector<AxisAlignedBoundingBox*> zAxisEnd;

bool collisionFound = false;
bool broadCollisionFound = false;
glm::vec3 collisionNormal;

bool planeCollision = false;
float planePos = -4.0f;
int satCount = 0;

float GetRandomValue(float minVal, float maxVal)
{
	float rndVal = ((float)rand()) / (float)RAND_MAX;
	return (rndVal * (maxVal - minVal)) + minVal;
}

void HandleRegularInput(unsigned char curKey, int x, int y)
{
	switch (curKey )
	{
		case ' ':
			simpleCamera->spinEnabled = !simpleCamera->spinEnabled;
			break;

		case '1':
			curDemoMode = DemoMode::World;
			simpleCamera->spinEnabled = false;
			simpleCamera->distance = 10.0f;
			simpleCamera->verticalOffset = 0.0f;
			simpleCamera->camSpin = 0.0f;
			break;

		case '2':
			curDemoMode = DemoMode::Collision;
			break;

		case '3':
			curDemoMode = DemoMode::Paused;
			break;

		case 'w':
			curSelectedCube->GetCurrentState()->momentum += glm::vec3(0.0f, 1.0f, 0.0f);
			break;

		case 's':
			curSelectedCube->GetCurrentState()->momentum -= glm::vec3(0.0f, 1.0f, 0.0f);
			break;

		case 'a':
			curSelectedCube->GetCurrentState()->momentum += glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'd':
			curSelectedCube->GetCurrentState()->momentum -= glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'i':
			curSelectedCube->GetCurrentState()->angularMomentum += glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'k':
			curSelectedCube->GetCurrentState()->angularMomentum -= glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'j':
			curSelectedCube->GetCurrentState()->angularMomentum -= glm::vec3(0.0f, 0.0f, 1.0f);
			break;

		case 'l':
			curSelectedCube->GetCurrentState()->angularMomentum += glm::vec3(0.0f, 0.0f, 1.0f);
			break;

		case 'b':
			sceneBodies[1]->GetCurrentState()->angularMomentum += glm::vec3(1.0f, 1.0f, 1.0f);
			sceneBodies[1]->GetCurrentState()->momentum += glm::vec3(-1.0f, 0.0f, 0.0f);
			break;
    }
}

void AddToSweepPrune(AxisAlignedBoundingBox* newAABB)
{
	xAxisStart.push_back(newAABB);
	yAxisStart.push_back(newAABB);
	zAxisStart.push_back(newAABB);

	xAxisEnd.push_back(newAABB);
	yAxisEnd.push_back(newAABB);
	zAxisEnd.push_back(newAABB);
}

void SortList(vector<AxisAlignedBoundingBox*> &unsortedList, int axis, bool sortMinValues)
{
	for(int i = 0; i < unsortedList.size(); i++)
	{
		AxisAlignedBoundingBox* curAABB = unsortedList[i];
		float j = i;

		while(j > 0 && 
			(sortMinValues ? unsortedList[j - 1]->minPoints[axis]->mValue : unsortedList[j - 1]->maxPoints[axis]->mValue) 
			> (sortMinValues ? curAABB->minPoints[axis]->mValue  : curAABB->maxPoints[axis]->mValue))
		{
			unsortedList[j] = unsortedList[j-1];
			j = j - 1;
		}

		unsortedList[j] = curAABB;
	}
}

void SortSweepPrune(void)
{
	SortList(xAxisStart, 0, true);
	SortList(yAxisStart, 1, true);
	SortList(zAxisStart, 2, true);

	SortList(xAxisEnd, 0, false);
	SortList(yAxisEnd, 1, false);
	SortList(zAxisEnd, 2, false);
}

void SweepAxis(vector<AxisAlignedBoundingBox*> &startAxis, vector<AxisAlignedBoundingBox*> &endAxis, vector<vector<int>> &colTable, int axisID)
{
	int comparisonCount = startAxis.size();
	int startID = 0;
	int endID = 0;

	list<AxisAlignedBoundingBox*> potentialColliders = list<AxisAlignedBoundingBox*>();

	while(startID < comparisonCount)
	{
		float curStartVal = startAxis[startID]->minPoints[axisID]->mValue;
		float curEndVal = endAxis[endID]->maxPoints[axisID]->mValue;

 		if(curStartVal < curEndVal)
		{
			potentialColliders.push_back(startAxis[startID]);
			activeColliders[startAxis[startID]->boundingID] = --potentialColliders.end();
			startID++;
		}
		else
		{
			potentialColliders.erase(activeColliders[endAxis[endID]->boundingID]);
			for(auto it = potentialColliders.begin(); it != potentialColliders.end(); ++it)
				++colTable[std::min((*it)->boundingID, endAxis[endID]->boundingID)][std::max((*it)->boundingID, endAxis[endID]->boundingID)];

			endID++;
		}
	}

	while(endID < comparisonCount)
	{
		potentialColliders.erase(activeColliders[endAxis[endID]->boundingID]);
		for(auto it = potentialColliders.begin(); it != potentialColliders.end(); ++it)
			++colTable[std::min((*it)->boundingID, endAxis[endID]->boundingID)][std::max((*it)->boundingID, endAxis[endID]->boundingID)];
		
		endID++;
	}
}

glm::vec3 GetPerpendicularVector(const glm::vec3 &initialVector)
{
	glm::vec3 perpVector = glm::vec3(initialVector.y, -initialVector.x, initialVector.z);

	float dot = glm::dot(perpVector, initialVector);
	// 0 = Perpendicular, 1 = Parallel

	if(dot < 0.0f && dot > 0.1f)
	{
		perpVector = glm::vec3(initialVector.x, -initialVector.z, initialVector.y);
		dot = glm::dot(perpVector, initialVector);

		if(dot < 0.0f && dot > 0.1f)
		{
			perpVector = glm::vec3(-initialVector.z, -initialVector.y, initialVector.x);
			dot = glm::dot(perpVector, initialVector);

			if(dot < 0.0f && dot > 0.1f)
			{
				perpVector = glm::vec3(1.0f, 1.0f, 1.0f);
			}
		}
	}
		
	return glm::normalize(perpVector);
}

vector<glm::vec3> GetAxes(const vector<Vertex> &verts1, const vector<Vertex> verts2)
{
	vector<glm::vec3> axes = vector<glm::vec3>();

	glm::vec3 vec1, vec2, vec3, vec4, edge1, edge2;
	for(int i = 0; i < verts1.size(); i++)
	{
		vec1 = verts1[i].position;
		vec2 = verts1[i % 4 == 3 ? i - 3 : i + 1].position;
		
		edge1 = vec1 - vec2;

		for(int j = 0; j < verts2.size(); j++)
		{
			vec3 = verts2[j].position;
			vec4 = verts2[j % 4 == 3 ? j - 3 : j + 1].position;
			edge2 = vec3 - vec4;

			glm::vec3 crossed = glm::cross(edge1, edge2);

			if(glm::length(crossed) != 0.0f)
			{
				crossed = glm::normalize(crossed);
				axes.push_back(crossed);
			}
		}
	}

	return axes;
}

void GetNormalAxes(vector<glm::vec3> &axes, const vector<Vertex> &vert1, const vector<Vertex> vert2)
{
	//This can be further optimised.
	for(int i = 0; i < 6; i++)
		axes.push_back(glm::normalize(vert1[4 * i].normal));

	for(int i = 0; i < 6; i++)
		axes.push_back(glm::normalize(vert2[4 * i].normal));
}

Projection ProjectOnAxis(glm::vec3 axis, const vector<Vertex> &verts)
{
	double min = glm::dot(axis, verts[0].position);
	double max = min;

	for(int i = 1; i < verts.size(); i++)
	{
		double p = glm::dot(axis, verts[i].position);

		if(p < min)
			min = p;
		else if(p > max)
			max = p;
	}

	return Projection(min, max);
}

float CheckForAxisOverlap(Projection p1, Projection p2)
{
	if(p2.min < p1.max && p1.min < p2.max)
		return min(abs(p2.max - p1.min), abs(p1.max - p2.min));

	if(p2.min > p1.min && p2.max < p1.max)
		return abs(p2.max - p2.min);

	if(p1.min > p2.min && p1.max < p2.max)
		return abs(p1.max - p1.min);

	return 0.0f;
}

bool vertexSort(glm::vec3 a, glm::vec3 b)
{
	if(abs(a.x) > abs(b.x))
		return true;

	if(abs(b.x) > abs(a.x))
		return false;

	if(abs(a.y) > abs(b.y))
		return true;

	if(abs(b.y) > abs(a.y))
		return false;

	if(abs(a.z) > abs(b.z))
		return true;

	if(abs(b.z) > abs(a.z))
		return false;

	return false;
}

bool vertexCompare(glm::vec3 a, glm::vec3 b)
{
	static double EPSILON = 0.001f;
	bool retVal = abs(a.x) - abs(b.x) < EPSILON;
	retVal &= abs(a.y) - abs(b.y) < EPSILON;
	retVal &= abs(a.z) - abs(b.z) < EPSILON;

	return retVal;

	//return glm::vec3(abs(a.x), abs(a.y), abs(a.z)) == glm::vec3(abs(b.x), abs(b.y), abs(b.z));
}

glm::vec3 GenerateFaceNormal(const glm::vec3 &vec1, const glm::vec3 &vec2, const glm::vec3 &vec3)
{
	glm::vec3 lhs = vec1 - vec2;
	glm::vec3 rhs = vec1 - vec3;

	glm::vec3 norm = glm::cross(lhs, rhs);

	glm::vec3 norm2 = glm::vec3(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x);

	if(glm::length(norm) > 0.0f)
		norm = glm::normalize(norm);

	if(glm::length(norm2) > 0.0f)
		norm2 = glm::normalize(norm2);

	return norm;
}

FaceCollection FindClosestFace(const glm::vec3 separationAxis, vector<glm::vec3> &verts)
{
	//STEP 1 - We want to find the closest face.
	float maxDist = -std::numeric_limits<float>::max();
	vector<int> closestFaceVerts = vector<int>();

	int closestVert = 0;
	//cout << separationAxis.x << " " << separationAxis.y << " " << separationAxis.z << "\n";
	for(int i = 0; i < verts.size(); i++)
	{
		//cout << verts[i].x << " " << verts[i].y << " " << verts[i].z << " - " << glm::dot(separationAxis, verts[i]) << "\n";
		float projection = glm::dot(separationAxis, verts[i]);

		if(projection == maxDist)
			closestFaceVerts.push_back(i);
		else if(projection > maxDist)
		{
			closestFaceVerts.clear();
			maxDist = projection;
			closestFaceVerts.push_back(i);
		}
	}

	const int arr[] = 
	{
		4, 0,3,2,1,	//Front
		4, 4,5,6,7,	//Back
		4, 1,2,6,5,	//Right
		4, 0,4,7,3,	//Left
		4, 5,4,0,1,	//Top
		4, 6,2,3,7	//Bottom
	};

	vector<int> objectFaces(arr, arr + sizeof(arr) / sizeof(arr[0]));
	int faceID = 0;

	FaceCollection collisionFaces = FaceCollection();

	while(faceID < objectFaces.size())
	{
		bool validFace = false;
		CollisionFace curFace = CollisionFace();
		int faceVertexCount = objectFaces[faceID];
		++faceID;
		for(int i = 0; i < faceVertexCount; i++)
		{
			curFace.vertexIDs.push_back(objectFaces[faceID + i]);

			for(int j = 0; j < closestFaceVerts.size(); j++)
			{
				if(objectFaces[faceID + i] == closestFaceVerts[j])
				{
					validFace = true;
					break;
				}
			}
		}

		if(!validFace)
		{
			faceID += faceVertexCount;
			continue;
		}

		curFace.faceNormal = GenerateFaceNormal(verts[objectFaces[faceID]], verts[objectFaces[faceID + 1]], verts[objectFaces[faceID + 2]]);
		
		collisionFaces.vertexFaces.push_back(curFace);
		faceID += faceVertexCount;
	}

	int collisionFaceID = 0;
	for(int i = 1; i < collisionFaces.vertexFaces.size(); i++)
	{
		float sDot1 = glm::dot(collisionFaces.vertexFaces[i].faceNormal, separationAxis);
		float sDot2 = glm::dot(collisionFaces.vertexFaces[collisionFaceID].faceNormal, separationAxis);

		//cout << sDot1 << " " << sDot2 << "\n";
		
		if(abs(glm::dot(collisionFaces.vertexFaces[i].faceNormal, separationAxis)) > abs(glm::dot(collisionFaces.vertexFaces[collisionFaceID].faceNormal, separationAxis)))
			collisionFaceID = i;
	}

	collisionFaces.collisionFaceID = collisionFaceID;
	collisionFaces.verts = verts;

	return collisionFaces;
}

glm::vec3 FindClosestFace(const glm::vec3 separationAxis, vector<Vertex> &verts) 
{
	//STEP 1 - We want to find the closest face.
	float maxDist = -std::numeric_limits<float>::max();
	int closestVert = 0;
	for(int i = 0; i < verts.size(); i++)
	{
		float projection = glm::dot(separationAxis, verts[i].position);
		if(projection > maxDist)
		{
			maxDist = projection;
			closestVert = i;
		}
	}

	const unsigned int faces[] = 
	{
		4, 0,3,2,1,
		4, 4,5,6,7,
		4, 1,5,6,2,
		4, 0,3,7,4,
		4, 5,4,0,1,
		4, 6,2,3,7
	};

	vector<GLuint> closestVerts = vector<GLuint>();

	for(int i = 0; i < verts.size(); i++)
	{
		if(verts[closestVert].position == verts[i].position)
			closestVerts.push_back(i);
	}

	vector<glm::vec3> faceNormals;

	for(int i = 0; i < closestVerts.size(); i++)
		faceNormals.push_back(verts[closestVerts[i]].normal);

	glm::vec3 closestFace = faceNormals[0];
	//cout << glm::dot(faceNormals[0], separationAxis) << "\n";
	for(int i = 1; i < faceNormals.size(); i++)
	{
		//cout << glm::dot(faceNormals[i], separationAxis) << "\n";
		if(glm::dot(faceNormals[i], separationAxis) >= glm::dot(closestFace, separationAxis))
			closestFace = faceNormals[i];
	}

	return closestFace;
}

void Clip(vector<glm::vec3> &points, glm::vec3 &clipFace, glm::vec3 refPoint)
{
	vector<glm::vec3> newPoints;
	float dotR = glm::dot(refPoint, clipFace);

	for(int k = 0; k < points.size(); k++)
	{
		glm::vec3 a = points[k];
		glm::vec3 b = points[(k + 1) % points.size()];

		float dotA = glm::dot(a, clipFace);
		float dotB = glm::dot(b, clipFace);

		//Four Possibilities
		//1: Both points are fine, add the second point to the contact manifold.
		if(dotR >= dotA && dotR >= dotB)
		{
			newPoints.push_back(b);
			continue;
		}

		//2: The first point is inside, the second point is outside. Add the point of intersection to the manifold.
		if(dotR >= dotA && dotR < dotB)
		{
			float intVal = (dotR - dotA) / (dotB - dotA);
			glm::vec3 intPoint = a + (b - a) * intVal;
			newPoints.push_back(intPoint);
			continue;
		}

		//3: The first point is outside, the second point is inside. Add the point of intersection and the second point.
		if(dotR < dotA && dotR >= dotB)
		{
			float intVal = (dotR - dotA) / (dotB - dotA);
			glm::vec3 intPoint = a + (b - a) * intVal;
			newPoints.push_back(intPoint);
			newPoints.push_back(b);
			continue;
		}

		//4: Both points are outside, so we don't add either.
		continue;
	}

	points = newPoints;
}

bool NarrowPhaseCollision(BasicCube* firstCube, BasicCube* secondCube, ContactManifold* contactManifold)
{
	vector<Vertex> fVerts = firstCube->GetTransformedVertices();
	vector<Vertex> sVerts = secondCube->GetTransformedVertices();

	vector<glm::vec3> axes = GetAxes(fVerts, sVerts);
	GetNormalAxes(axes, sVerts, fVerts);

	vector<glm::vec3> culledAxes = vector<glm::vec3>();

	sort(axes.begin(), axes.end(), vertexSort);
	std::vector<glm::vec3>::iterator it;
	it = std::unique(axes.begin(), axes.end(), vertexCompare);
	axes.resize(std::distance(axes.begin(), it));
	
	satCount = axes.size();

	float minOverlapDistance = std::numeric_limits<float>::max();
	float curOverlap;
	glm::vec3 separationAxis;

	for(int i = 0; i < axes.size(); i++)
	{
		Projection p1 = ProjectOnAxis(axes[i], fVerts);
		Projection p2 = ProjectOnAxis(axes[i], sVerts);

		curOverlap = CheckForAxisOverlap(p1, p2);

		if(curOverlap > 0)
		{
			if(curOverlap < minOverlapDistance)
			{
				minOverlapDistance = curOverlap;
				separationAxis = axes[i];
			}
		}
		else
			return false;
	}

	if(glm::dot(firstCube->GetCurrentState()->position, separationAxis) > glm::dot(secondCube->GetCurrentState()->position, separationAxis))
		separationAxis*= -1;

	collisionNormal = separationAxis;

	vector<glm::vec3> firstCollisionVerts = firstCube->GetTransformedCollisionVertices();
	vector<glm::vec3> secondCollisionVerts = secondCube->GetTransformedCollisionVertices();

	FaceCollection firstCollection = FindClosestFace(separationAxis, firstCollisionVerts);
	FaceCollection secondCollection = FindClosestFace(-separationAxis, secondCollisionVerts);
	//Most parallel is the one to be clipped against

	FaceCollection *clipCollection, *refCollection;

	//cout << abs(glm::dot(separationAxis, firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal)) << " " << abs(glm::dot(separationAxis, secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal)) << "\n";

	bool flippedClip = false;
	//if(abs(glm::dot(firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal, separationAxis)) <= abs(glm::dot(secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal, separationAxis)))
	
	//cout << "A+ " << glm::dot(firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal, separationAxis) << "\nB+ " << glm::dot(secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal, separationAxis) << "\nA- " << glm::dot(firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal, -separationAxis) << "\nB- " << glm::dot(secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal, -separationAxis);
	if(abs(glm::dot(firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal, separationAxis)) <= abs(glm::dot(secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal, separationAxis)))
	{
		clipCollection = &firstCollection;
		refCollection = &secondCollection;
	}
	else
	{
		clipCollection = &secondCollection;
		refCollection = &firstCollection;
		flippedClip = true;				//We set a flag to indicate that the ref and clip plane were flipped so that we do the final clip operation using the right face normal.
	}

	contactManifold->verts.clear();
	for(int i = 0; i < clipCollection->vertexFaces[clipCollection->collisionFaceID].vertexIDs.size(); i++)
		contactManifold->verts.push_back(clipCollection->verts[clipCollection->vertexFaces[clipCollection->collisionFaceID].vertexIDs[i]]);
	
	//Now we've got the COLLISION FACES, the EDGE FACES, and the EDGE FACE NORMALS.
	for(int i = 0; i < refCollection->vertexFaces.size(); i++)
	{
		if(i == refCollection->collisionFaceID)
			continue;

		Clip(contactManifold->verts, refCollection->vertexFaces[i].faceNormal, refCollection->verts[refCollection->vertexFaces[i].vertexIDs[0]]);		
	}

	vector<glm::vec3> newPoints;
	glm::vec3 clipFace = refCollection->vertexFaces[refCollection->collisionFaceID].faceNormal;

	float dotR = glm::dot(refCollection->verts[refCollection->vertexFaces[refCollection->collisionFaceID].vertexIDs[0]], clipFace);
	for(int i = 0; i < contactManifold->verts.size(); i++)
	{
		float pDot = glm::dot(contactManifold->verts[i], clipFace);
		if(dotR > pDot)
			newPoints.push_back(contactManifold->verts[i]);
	}

	contactManifold->verts = newPoints;

	return true;
}

//USE THE FACE INFORMATION TO FIND THE EDGES OF THE UCBE NECCESSARY FOR TESTING.
//PICK ANY VERTEX, THEN GET THE CONNECTING EDGES FROM THE FACE.
//THESE THREE EDGES ARE UNIQUE FOR THE CUBE, AND CROSS THEM AGAINST THE THREE UNIQUE FACES FROM THE OTHER CUBE.
//FOR A QUICK CHECK!

void NiceBroadPhase(void)
{
	list<AxisAlignedBoundingBox*> potentialColliders = list<AxisAlignedBoundingBox*>();
	collisionFound = false;
	broadCollisionFound = false;

	for(int i = 0; i < sceneBodies.size(); i++)
		sceneBodies[i]->GetBoundingBox()->isColliding = false;

	SortSweepPrune();
	SweepAxis(xAxisStart, xAxisEnd, collisionTable, 0);
	SweepAxis(yAxisStart, yAxisEnd, collisionTable, 1);
	SweepAxis(zAxisStart, zAxisEnd, collisionTable, 2);

	for(int i = 0; i < sceneBodies.size(); i++)
	{
		for(int j = 0; j < sceneBodies.size(); j++)
		{
			if(collisionTable[i][j] == 3)
			{
				//We've got a collision!
				sceneBodies[i]->GetBoundingBox()->isColliding = true;
				sceneBodies[j]->GetBoundingBox()->isColliding = true;

				broadCollisionFound = true;

				if(NarrowPhaseCollision(sceneBodies[i], sceneBodies[j], contactM))
				{
					cout << "A collision has been found between " << i << " and " << j << ".\n";
					collisionFound = true;
				}
			}

			collisionTable[i][j] = 0;	//No matter what we reset it to zero for the next frame.
		}
	}

	if(!collisionFound)
		contactM->verts.clear();
	else
		curDemoMode = DemoMode::Collision;

}

void UpdateScene(void)
{
	float t = 0.0f;
	const float dt = 1.0f / 120.0f;

	int timeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	int deltaTime = timeSinceStart - prevFrameTimeSinceStart;
	prevFrameTimeSinceStart = timeSinceStart;

	elapsedTimeStep = ((float)deltaTime / 1000.0f);
	simpleCamera->Update(elapsedTimeStep);

	if(curDemoMode == DemoMode::World)
	{
		if(elapsedTimeStep > 0.25)
			elapsedTimeStep = 0.25;				//Maximum frame time allowed if we want to avoid the spiral of death.

		accumulator += elapsedTimeStep;			//This is wonderful.
		simpleCamera->Update(elapsedTimeStep);

		while(accumulator >= dt)	//so, so wonderful.
		{
			//t = time, dt = deltaTime.
			for(int i = 0; i < sceneBodies.size(); i++)
				sceneBodies[i]->Update(t, dt);
		
			accumulator -= dt;
			t += dt;
		}

		planeCollision = false;
		contactM->verts.clear();
		vector<glm::vec3> playerVerts = sceneBodies[0]->GetTransformedCollisionVertices();
		for(int i = 0; i < playerVerts.size(); i++)
		{
			if(playerVerts[i].y <= planePos)
			{
				contactM->verts.push_back(playerVerts[i]);
				planeCollision = true;
				collisionNormal = glm::vec3(0.0f, 1.0f, 0.0f);
			}
		}

		if(planeCollision)
			curDemoMode = DemoMode::Collision;
		else
			NiceBroadPhase();

		//Remaining "time" in the current frame, to interpolate our current physics step by.
		interpolationAlpha = accumulator / dt;
	}

	glutPostRedisplay();
}

void RenderScene(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(curDemoMode == DemoMode::Collision)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glUseProgram(bodyShader->GetShaderID());

	for(int i = 0; i < sceneBodies.size(); i++)
		sceneBodies[i]->Render(simpleCamera->GetViewMatrix(), simpleCamera->GetProjMatrix(), interpolationAlpha);

	glUseProgram(lineShader->GetShaderID());
	glm::mat4 identityMat = simpleCamera->GetProjMatrix() * simpleCamera->GetViewMatrix();
	glUniformMatrix4fv(glGetUniformLocation(lineShader->GetShaderID(), "gWVP"), 1, GL_FALSE, &identityMat[0][0]);

	glLineWidth(1.0f);
	for(int i = 0; i < sceneBodies.size(); i++)
	{
		if(sceneBodies[i]->GetBoundingBox()->isColliding)
			glColor3f(1.0f, 0.0f, 0.0f);
		else
			glColor3f(0.0f, 1.0f, 0.0f);

		sceneBodies[i]->RenderAABB();
	}

	glColor3f(0.6f, 0.6f, 0.6f);
	glBegin(GL_QUADS);
		glVertex3f(-1000.0f, planePos, -1000.0f);
		glVertex3f(-1000.0f, planePos, 1000.0f);
		glVertex3f(1000.0f, planePos, 1000.0f);
		glVertex3f(1000.0f, planePos, -1000.0f);	
	glEnd();	

	if(curDemoMode == DemoMode::Collision)
	{
		glPointSize(10);
		glColor3f(1.0f, 0.0f, 0.0f);

		glBegin(GL_POINTS);	

		glm::vec3 coreCollisionPoint = glm::vec3();

		for(int i = 0; i < contactM->verts.size(); i++)
		{
			coreCollisionPoint += contactM->verts[i];
			glVertex3f(contactM->verts[i].x, contactM->verts[i].y, contactM->verts[i].z);
		}

		coreCollisionPoint = glm::vec3(coreCollisionPoint.x / contactM->verts.size(), coreCollisionPoint.y / contactM->verts.size(), coreCollisionPoint.z / contactM->verts.size());

		glEnd();

		glLineWidth(2);
		glColor3f(0.0f, 1.0f, 0.0f);

		glBegin(GL_LINES);
		glVertex3f(coreCollisionPoint.x, coreCollisionPoint.y, coreCollisionPoint.z);
		glVertex3f(coreCollisionPoint.x + (collisionNormal.x * 1.0f), coreCollisionPoint.y + (collisionNormal.y * 1.0f), coreCollisionPoint.z + (collisionNormal.z * 1.0f));
		glEnd();


		glColor3f(0.0f, 0.0f, 1.0f);
		glBegin(GL_LINES);
		glVertex3f(coreCollisionPoint.x, coreCollisionPoint.y, coreCollisionPoint.z);
		glVertex3f(coreCollisionPoint.x - (collisionNormal.x * 1.0f), coreCollisionPoint.y - (collisionNormal.y * 1.0f), coreCollisionPoint.z - (collisionNormal.z * 1.0f));
		glEnd();

		glLineWidth(1);

		/*for(int i = 0; i < sceneBodies.size(); i++)
		{
			glLineWidth(4);
			glColor3f(collisionNormal.x, collisionNormal.y, collisionNormal.z);
			glBegin(GL_LINES);

			glVertex3f(sceneBodies[i]->GetCurrentState()->position.x, sceneBodies[i]->GetCurrentState()->position.y, sceneBodies[i]->GetCurrentState()->position.z);
			glVertex3f(sceneBodies[i]->GetCurrentState()->position.x + (collisionNormal.x * 2.0f), sceneBodies[i]->GetCurrentState()->position.y + (collisionNormal.y * 2.0f), sceneBodies[i]->GetCurrentState()->position.z + (collisionNormal.z * 2.0f));

			glEnd();
		}*/
	}

	stringstream ss;
	ss << "Narrow Phase Collision - Separating Axis Theorem (SAT) & Culling\n";
	if(broadCollisionFound)
	{
		if(collisionFound)
			ss << "A narrowphase collision has been found.\n";
		else
			ss << "A broadphase collision has been found, but narrow phase revealed no collision.\n";

		ss << "Separating Axes - " << satCount << "\n";
	}
	else
	{
		ss << "No collisions found.\n";
	}

	if(contactM->verts.size() > 0)
	{
		ss << "\nCollision Information\n";
		ss << "Collision Normal - (" << collisionNormal.x << ", " << collisionNormal.y << ", " << collisionNormal.z << ")\n";
	
		if(contactM->verts.size() == 1)
			ss << "Collision Point - (" << contactM->verts[0].x << ", " << contactM->verts[0].y << ", " << contactM->verts[0].z << ")\n";
		else
		{
			ss << "Collision Points\n";
			for(int i = 0; i < contactM->verts.size(); i++)
			{
				ss << "(" << contactM->verts[i].x << ", " << contactM->verts[i].y << ", " << contactM->verts[i].z << ")\n";
			}
		}
	}

	string text = ss.str();
	glUseProgram(0);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glRasterPos2f(-1.0f, 0.95f);	
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char*)text.c_str());

    glutSwapBuffers();
}

void InitialiseScene(void)
{
	//Code for the rigid bodies is now nice and tidy, but look at cleaning up Sweep and Prune / Better visualising it!
	contactM = new ContactManifold();

	prevFrameTimeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	
	bodyShader = new Shader("vTexture.shader", "fTexture.shader");
	lineShader = new Shader("vLine.shader", "fLine.shader");
	simpleCamera = new SimpleCamera();

	sceneBodies.push_back(new BasicCube(glm::vec3(-5.0f, 0.0f, 0.0f), glm::vec3(), glm::vec3()));

	for(int i = 1; i < rigidBodyCount; i++)
		sceneBodies.push_back(new BasicCube());

	curSelectedCube = sceneBodies[0];

	for(int i = 0; i < sceneBodies.size(); i++)
		AddToSweepPrune(sceneBodies[i]->GetBoundingBox());

	activeColliders.resize(sceneBodies.size());

	collisionTable = vector<vector<int>>();
	collisionTable.resize(sceneBodies.size());
	for(int i = 0; i < sceneBodies.size(); i++)
	{
		collisionTable[i].resize(sceneBodies.size());
		for(int j = 0; j < sceneBodies.size(); j++)
			collisionTable[i][j] = 0;
	}

	SortSweepPrune();
}

void InitialiseCallbacks()
{
	glutKeyboardFunc(HandleRegularInput);
	//glutPassiveMotionFunc(HandleMouseMovement);
	//glutSpecialFunc(HandleSpecialInput);
	//glutMouseFunc(HandleMouseInput);
	glutIdleFunc(UpdateScene);							// Tell glut where the update function is, to execute when system is idle (after drawing completed).
	glutDisplayFunc(RenderScene);						// Tell glut where the display function is
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH)-1280)/2, (glutGet(GLUT_SCREEN_HEIGHT)-720)/2);
	glutCreateWindow("Real Time Physics - Lab 4 - Narrow Phase Collisions - SAT and Culling");
	//glutGameModeString("1280x720@32");
    //glutEnterGameMode();

	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_LINE_SMOOTH);
	//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	//glEnable(GL_POINT_SMOOTH);
	//glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//glEnable(GL_MULTISAMPLE);
	//glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);

	InitialiseCallbacks();
	
	GLenum res = glewInit();
    if (res != GLEW_OK)
	{
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		return 1;
    }

	glClearColor(0.3f, 0.3f, 0.3f, 1.0f);

	ilInit();
	iluInit();
	ilutInit();
	ilutRenderer(ILUT_OPENGL);

	InitialiseScene();
	glutMainLoop();
	return 0;
}

