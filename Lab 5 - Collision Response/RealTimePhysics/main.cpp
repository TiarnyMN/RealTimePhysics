#pragma region Includes and Defines
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <GL/glew.h>
#include <GL/freeglut.h>

#include "shared.h"
#include "shader.h"
#include "Plane.h"
#include "BasicCube.h"
#include "Cuboid.h"
#include "Pyramid.h"
#include "SimpleCamera.h"
#include "RigidBody.h"
#include "AxisAlignedBoundingBox.h"

#include <algorithm>
#include <vector>

#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 720
#pragma endregion

#pragma region Structs and Enums
//The min and max points of a rigid body projected on to a separation axis.
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

//Holds the points in a Contact Manifold.
struct ContactManifold
{
	vector<glm::vec3> verts;
};

//Holds the IDs of each vertex, and the normal, of each face needed for Clipping.
struct CollisionFace
{
	vector<int> vertexIDs;
	glm::vec3 faceNormal;
};

//Holds all the information needed for clipping - the collision face, the reference faces, and all the verts.
struct FaceCollection
{
	int collisionFaceID;
	vector<CollisionFace> vertexFaces;
	vector<glm::vec3> verts;
};

//Contains all the contact information needed for collision response to generate the correct impulse response.
struct ContactInfo
{
	ContactInfo()
	{
		penetrationDepth = 0.0f;
		denominator = 1.0f;

		totalImpulse = 0.0f;
		penetrationImpulse = 0.0f;
	}

	glm::vec3 rA;
	glm::vec3 rB;

	glm::vec3 position;

	float penetrationDepth;
	float denominator;

	float totalImpulse;
	float penetrationImpulse;
	float restitutionImpulse;
};

//Holds all the information related to a collision between two rigid bodies.
struct CollisionInfo
{
	CollisionInfo()
	{
		vector<ContactInfo> contactPoints = vector<ContactInfo>();

		largestImpulseDifference = std::numeric_limits<float>::max();
		converged = false;
	}

	int aID;
	int bID;

	Polyhedron* a;
	Polyhedron* b;

	BodyState* A;
	BodyState* B;

	glm::vec3 collisionNormal;
	ContactManifold contactManifold;
	vector<ContactInfo> contactPoints;

	float largestImpulseDifference;

	bool converged;	//If they have sufficiently convered, we can stop iterating over them.
};

//Enums which control how the demo works.
enum DemoMode
{
	Barrel,
	Stacking,
	Game,
};

enum WorldView
{
	Regular,
	Full,
	Simple,
	Wireframe,
	Collision
};

enum GameState
{
	Running,
	Paused
};
#pragma endregion

#pragma region Variables
//Camera
SimpleCamera* simpleCamera = NULL;

//Shaders
Shader* bodyShader;
Shader* lineShader;

//Timestep
float relativeGameSpeed = 1.0f;		//For slow motion collisions.
int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;
double accumulator = 0.0f;
double interpolationAlpha = 0.0f;

//Demo manipulators / Debug counters
DemoMode curDemoMode = DemoMode::Barrel;
WorldView curWorldView = WorldView::Regular;
GameState curGameState = GameState::Running;

int broadPhaseCollisionsFound = 0;
int narrowPhaseCollisionsFound = 0;
int highlightedCollision = 1000;

bool gravityEnabled = false;
bool pauseOnCollisions = false;

//Rigid Body Info
vector<Polyhedron*> sceneBodies;
Polyhedron* curSelectedPoly;

//Narrow Phase Collision Info
vector<CollisionInfo> collisionPairs = vector<CollisionInfo>();

//Broad Phase Collision Info
vector<vector<int>> collisionTable;

//Sweep and Prune
vector<AxisAlignedBoundingBox*> xAxisStart;
vector<AxisAlignedBoundingBox*> yAxisStart;
vector<AxisAlignedBoundingBox*> zAxisStart;

vector<AxisAlignedBoundingBox*> xAxisEnd;
vector<AxisAlignedBoundingBox*> yAxisEnd;
vector<AxisAlignedBoundingBox*> zAxisEnd;
#pragma endregion

//Generates a random value between minVal and maxVal.
float GetRandomValue(float minVal, float maxVal)
{
	float rndVal = ((float)rand()) / (float)RAND_MAX;
	return (rndVal * (maxVal - minVal)) + minVal;
}

#pragma region Sweep & Prune
//Adds the AABB to our Sweep and Prune considerations.
void AddToSweepPrune(AxisAlignedBoundingBox* newAABB)
{
	xAxisStart.push_back(newAABB);
	yAxisStart.push_back(newAABB);
	zAxisStart.push_back(newAABB);

	xAxisEnd.push_back(newAABB);
	yAxisEnd.push_back(newAABB);
	zAxisEnd.push_back(newAABB);
}

//Sorts the list, using an Insertion Sort, for all points on the given axis.
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

//Sorts all our axes for a Sweep and Prune update.
void SortSweepPrune(void)
{
	SortList(xAxisStart, 0, true);
	SortList(yAxisStart, 1, true);
	SortList(zAxisStart, 2, true);

	SortList(xAxisEnd, 0, false);
	SortList(yAxisEnd, 1, false);
	SortList(zAxisEnd, 2, false);
}

//Performs a Sweep and Prune operation, searching a given axis for overlaps and recording them in the collision table.
bool SweepAxis(vector<AxisAlignedBoundingBox*> &startAxis, vector<AxisAlignedBoundingBox*> &endAxis, vector<vector<int>> &colTable, int axisID)
{
	int comparisonCount = startAxis.size();
	int startID = 0;
	int endID = 0;

	bool retVal = false;	//Allows us to quit out early if no collisions are found.

	list<AxisAlignedBoundingBox*> potentialColliders = list<AxisAlignedBoundingBox*>();
	
	vector<std::list<AxisAlignedBoundingBox*>::iterator> activeColliders;
	activeColliders.resize(sceneBodies.size());

	while(startID < comparisonCount)
	{
		float curStartVal = startAxis[startID]->minPoints[axisID]->mValue;
		float curEndVal = endAxis[endID]->maxPoints[axisID]->mValue;

 		if(curStartVal <= curEndVal)
		{
			potentialColliders.push_back(startAxis[startID]);
			activeColliders[startAxis[startID]->boundingID] = --potentialColliders.end();
			startID++;
		}
		else
		{
			potentialColliders.erase(activeColliders[endAxis[endID]->boundingID]);
			
			if(potentialColliders.size() > 0)
				retVal = true;

			for(auto it = potentialColliders.begin(); it != potentialColliders.end(); ++it)
				++colTable[std::min((*it)->boundingID, endAxis[endID]->boundingID)][std::max((*it)->boundingID, endAxis[endID]->boundingID)];

			endID++;
		}
	}

	while(endID < comparisonCount)
	{
		potentialColliders.erase(activeColliders[endAxis[endID]->boundingID]);

		if(potentialColliders.size() > 0)
				retVal = true;

		for(auto it = potentialColliders.begin(); it != potentialColliders.end(); ++it)
			++colTable[std::min((*it)->boundingID, endAxis[endID]->boundingID)][std::max((*it)->boundingID, endAxis[endID]->boundingID)];
		
		endID++;
	}

	return retVal;
}
#pragma endregion

#pragma region Separating Axis Theroem (SAT) & Clipping
//Gets a vector perpendicular to the initial vector.
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

//Gets all the potential separation axes between two objects.
vector<glm::vec3> GetAxes(const vector<glm::vec3> &fVerts, const vector<glm::vec3> sVerts, const vector<int> &fFaces, const vector<int> &sFaces)
{
	vector<glm::vec3> axes = vector<glm::vec3>();

	glm::vec3 fVert, sVert, tVert, foVert, fEdge, sEdge;
	int firstFaceID = 0;
	int secondFaceID = 0;
	while(firstFaceID < fFaces.size())
	{
		int firstFaceVertexCount = fFaces[firstFaceID];
		++firstFaceID;
		for(int i = 0; i < firstFaceVertexCount; i++)
		{
			fVert = fVerts[fFaces[firstFaceID + i]];
			sVert = fVerts[fFaces[firstFaceID + (i + 1) % firstFaceVertexCount]];

			fEdge = fVert - sVert;

			secondFaceID = 0;
			while(secondFaceID < sFaces.size())
			{
				int secondFaceVertexCount = sFaces[secondFaceID];
				++secondFaceID;

				for(int j = 0; j < secondFaceVertexCount; j++)
				{
					tVert = sVerts[sFaces[secondFaceID + j]];
					foVert = sVerts[sFaces[secondFaceID + (j + 1) % secondFaceVertexCount]];

					sEdge = tVert - foVert;

					glm::vec3 crossedAxes = glm::cross(fEdge, sEdge);

					if(glm::length(crossedAxes) != 0.0f)
					{
						axes.push_back(glm::normalize(crossedAxes));
					}
				}

				secondFaceID += secondFaceVertexCount;
			}
		}

		firstFaceID += firstFaceVertexCount;		
	}
	return axes;
}

//Adds the normals of each object to the list of separation axes.
void GetNormalAxes(vector<glm::vec3> &axes, const vector<Vertex> &vert1, const vector<Vertex> vert2, const vector<int> &fFaces, const vector<int> &sFaces)
{
	int faceID = 0;
	int normalID = 0;
	while(faceID < fFaces.size())
	{
		axes.push_back(glm::normalize(vert1[normalID].normal));
		normalID += fFaces[faceID];
		faceID += fFaces[faceID] + 1;
	}

	faceID = 0;
	normalID = 0;
	while(faceID < sFaces.size())
	{
		axes.push_back(glm::normalize(vert2[normalID].normal));
		normalID += sFaces[faceID];
		faceID += sFaces[faceID] + 1;
	}
}

//Projects the object on to the given axes, and returns the projection results.
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

//Checks to see if two objects overlap on a given axes.
float CheckForAxisOverlap(Projection p1, Projection p2)
{
	//Checks for overlap.
	if(p2.min < p1.max && p1.min < p2.max)
		return min(abs(p2.max - p1.min), abs(p1.max - p2.min));

	//Checks for containment of p2 inside p1.
	if(p2.min > p1.min && p2.max < p1.max)
		return abs(p2.max - p2.min);

	//Checks for containment of p1 inside p2.
	if(p1.min > p2.min && p1.max < p2.max)
		return abs(p1.max - p1.min);

	return 0.0f;
}

//Given three verts on a face, generate its normal.
glm::vec3 GenerateFaceNormal(const glm::vec3 &vec1, const glm::vec3 &vec2, const glm::vec3 &vec3)
{
	glm::vec3 lhs = vec1 - vec2;
	glm::vec3 rhs = vec1 - vec3;

	glm::vec3 norm = glm::cross(lhs, rhs);

	if(glm::length(norm) > 0.0f)
		norm = glm::normalize(norm);

	return norm;
}

//Find the closest face involved in a collision.
FaceCollection FindClosestFace(const glm::vec3 separationAxis, vector<glm::vec3> &verts, vector<int> &faceIndices)
{
	float maxDist = -std::numeric_limits<float>::max();
	vector<int> closestFaceVerts = vector<int>();

	int closestVert = 0;
	for(int i = 0; i < verts.size(); i++)
	{
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

	int faceID = 0;

	FaceCollection collisionFaces = FaceCollection();

	while(faceID < faceIndices.size())
	{
		bool validFace = false;
		CollisionFace curFace = CollisionFace();
		int faceVertexCount = faceIndices[faceID];
		++faceID;
		for(int i = 0; i < faceVertexCount; i++)
		{
			curFace.vertexIDs.push_back(faceIndices[faceID + i]);

			for(int j = 0; j < closestFaceVerts.size(); j++)
			{
				if(faceIndices[faceID + i] == closestFaceVerts[j])
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

		curFace.faceNormal = GenerateFaceNormal(verts[faceIndices[faceID]], verts[faceIndices[faceID + 1]], verts[faceIndices[faceID + 2]]);
		
		collisionFaces.vertexFaces.push_back(curFace);
		faceID += faceVertexCount;
	}

	int collisionFaceID = 0;
	for(int i = 1; i < collisionFaces.vertexFaces.size(); i++)
	{	
		if(abs(glm::dot(collisionFaces.vertexFaces[i].faceNormal, separationAxis)) > abs(glm::dot(collisionFaces.vertexFaces[collisionFaceID].faceNormal, separationAxis)))
			collisionFaceID = i;
	}

	collisionFaces.collisionFaceID = collisionFaceID;
	collisionFaces.verts = verts;

	return collisionFaces;
}

//Perform a clipping operation. 2D Source - http://www.codezealot.org/archives/394, 3D Algorithm - http://en.wikipedia.org/wiki/Sutherland-Hodgman_algorithm
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

		
		if(dotR >= dotA && dotR >= dotB)	//Both points are fine, add the second point to the contact manifold.
		{
			newPoints.push_back(b);
			continue;
		}
		else if(dotR >= dotA && dotR < dotB)	//The first point is inside, the second point is outside. Add the point of intersection to the manifold.
		{
			float intVal = (dotR - dotA) / (dotB - dotA);
			glm::vec3 intPoint = a + (b - a) * intVal;
			newPoints.push_back(intPoint);
			continue;
		}
		else if(dotR < dotA && dotR >= dotB)	//The first point is outside, the second point is inside. Add the point of intersection and the second point.
		{
			float intVal = (dotR - dotA) / (dotB - dotA);
			glm::vec3 intPoint = a + (b - a) * intVal;
			newPoints.push_back(intPoint);
			newPoints.push_back(b);
			continue;
		}

		continue;	//Both points are outside, so toss both.
	}

	points = newPoints;
}
#pragma endregion

//A sorting method to prepare the duplicated separation axes for pruning.
bool SortAxes(glm::vec3 a, glm::vec3 b)
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

//Determines if two axes can be considered duplicate and discared. An EPSILON value accounts for floating point precision error.
bool CompareAxes(glm::vec3 a, glm::vec3 b)
{
	static double EPSILON = 0.00001f;

	if(abs(glm::dot(a, b)) > (1 - EPSILON))
		return true;

	return false;
}

//Carries out a Narrow Phase Collision check using SAT & Clipping. SAT 2D Source - http://www.codezealot.org/archives/55 
bool NarrowPhaseCollision(Polyhedron* firstPoly, Polyhedron* secondPoly, CollisionInfo* collision)
{
	//Add the information to our collision container.
	collision->a = firstPoly;
	collision->b = secondPoly;
	collision->A = firstPoly->GetCurrentState();
	collision->B = secondPoly->GetCurrentState();

	//Get our transformed vertices, to be used for the projection.
	vector<Vertex> fVerts = firstPoly->GetTransformedVertices();
	vector<Vertex> sVerts = secondPoly->GetTransformedVertices();

	vector<glm::vec3> axes = GetAxes(firstPoly->GetTransformedCollisionVertices(), secondPoly->GetTransformedCollisionVertices(), firstPoly->GetFaceIndices(), secondPoly->GetFaceIndices());	//Get all the separation axes.
	GetNormalAxes(axes, fVerts, sVerts, firstPoly->GetFaceIndices(), secondPoly->GetFaceIndices());	//And add the normals to our separation axes too.

	//Sort and prune any duplicate separation axes.
	sort(axes.begin(), axes.end(), SortAxes);
	std::vector<glm::vec3>::iterator it;
	it = std::unique(axes.begin(), axes.end(), CompareAxes);
	axes.resize(std::distance(axes.begin(), it));	//Resize it, now that duplicates have been removed.
	
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
				separationAxis = axes[i];			//The minimum overlap is also the collision normal, so we want to store that for later.
			}
		}
		else
			return false;	//If no overlap, then not a collision!
	}

	//We want to make sure the collision normal always points the same way (from A > B)
	if(glm::dot(firstPoly->GetCurrentState()->position, separationAxis) > glm::dot(secondPoly->GetCurrentState()->position, separationAxis))
		separationAxis*= -1;

	//Prepare for the clipping by setting up our face collections, collision normal and collision vertices.
	collision->collisionNormal = separationAxis;

	vector<glm::vec3> firstCollisionVerts = firstPoly->GetTransformedCollisionVertices();
	vector<glm::vec3> secondCollisionVerts = secondPoly->GetTransformedCollisionVertices();

	FaceCollection firstCollection = FindClosestFace(separationAxis, firstCollisionVerts, firstPoly->GetFaceIndices());
	FaceCollection secondCollection = FindClosestFace(-separationAxis, secondCollisionVerts, secondPoly->GetFaceIndices());
	//Most parallel is the one to be clipped against

	FaceCollection *clipCollection, *refCollection;	//These will point to the faces, depending on which one is being clipped against.	
	if(abs(glm::dot(firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal, separationAxis)) <= abs(glm::dot(secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal, separationAxis)))
	{
		clipCollection = &firstCollection;
		refCollection = &secondCollection;
	}
	else
	{
		clipCollection = &secondCollection;
		refCollection = &firstCollection;
	}

	for(int i = 0; i < clipCollection->vertexFaces[clipCollection->collisionFaceID].vertexIDs.size(); i++)
		collision->contactManifold.verts.push_back(clipCollection->verts[clipCollection->vertexFaces[clipCollection->collisionFaceID].vertexIDs[i]]);
	
	//Now we've got the COLLISION FACES, the EDGE FACES, and the EDGE FACE NORMALS. Now we perform the clipping.
	for(int i = 0; i < refCollection->vertexFaces.size(); i++)
	{
		if(i == refCollection->collisionFaceID)
			continue;

		Clip(collision->contactManifold.verts, refCollection->vertexFaces[i].faceNormal, refCollection->verts[refCollection->vertexFaces[i].vertexIDs[0]]);		
	}

	if(curDemoMode == DemoMode::Barrel)
	{
		//Calculate "Frictional Forces" based on the angle and speed of the collision.
		glm::vec3 fA = glm::dot(secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal, collision->A->momentum) * secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal * 0.2f;
		glm::vec3 fB = glm::dot(firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal, collision->B->momentum) * firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal * 0.2f;

		glm::vec3 fmA = glm::dot(secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal, collision->A->angularMomentum) * secondCollection.vertexFaces[secondCollection.collisionFaceID].faceNormal * 0.2f;
		glm::vec3 fmB = glm::dot(firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal, collision->B->angularMomentum) * firstCollection.vertexFaces[firstCollection.collisionFaceID].faceNormal * 0.2f;
	
		collision->A->momentum -= glm::vec3(fA.x, 0.0f, fA.z);
		collision->B->momentum -= glm::vec3(fB.x, 0.0f, fB.z);

		collision->A->angularMomentum -= glm::vec3(fmA.x, fmA.y, fmA.z);
		collision->B->angularMomentum -= glm::vec3(fmB.z, fmB.y, fmB.z);

		glm::mat3 rotationMat = glm::mat3_cast(collision->A->orientation);
		glm::mat3 inverseTensor = glm::transpose(rotationMat) * collision->A->inverseInertialTensor * rotationMat;
		collision->A->velocity = collision->A->momentum * collision->A->inverseMass;
		collision->A->angularVelocity = collision->A->angularMomentum * inverseTensor;

		rotationMat = glm::mat3_cast(collision->B->orientation);
		inverseTensor = glm::transpose(rotationMat) * collision->B->inverseInertialTensor * rotationMat;
		collision->B->velocity = collision->B->momentum * collision->B->inverseMass;
		collision->B->angularVelocity = collision->B->angularMomentum * inverseTensor;
	}


	vector<glm::vec3> newPoints;
	glm::vec3 clipFace = refCollection->vertexFaces[refCollection->collisionFaceID].faceNormal;

	//We need to perform one final clip, this tells us which of the points on the manifold are actually penetrating (e.g. truly involved in the collision).
	//If we find one, we can store the static information to save computing it for every iteration of our impulses.
	float dotR = glm::dot(refCollection->verts[refCollection->vertexFaces[refCollection->collisionFaceID].vertexIDs[0]], clipFace);
	for(int i = 0; i < collision->contactManifold.verts.size(); i++)
	{
		float pDot = glm::dot(collision->contactManifold.verts[i], clipFace);

		if(dotR > pDot)
		{
			newPoints.push_back(collision->contactManifold.verts[i]);
			
			ContactInfo contactInfo = ContactInfo();
			contactInfo.position = collision->contactManifold.verts[i];
			contactInfo.penetrationDepth = abs(pDot - dotR);

			contactInfo.rA = contactInfo.position - collision->A->position;
			contactInfo.rB = contactInfo.position - collision->B->position;

			contactInfo.denominator = (collision->A->inverseMass + collision->B->inverseMass + 
								glm::dot(collision->collisionNormal, collision->A->inverseInertialTensor * glm::cross((glm::cross(contactInfo.rA, collision->collisionNormal)), contactInfo.rA)) +
								glm::dot(collision->collisionNormal, collision->B->inverseInertialTensor * glm::cross((glm::cross(contactInfo.rB, collision->collisionNormal)), contactInfo.rB)));
			
			glm::vec3 pA = collision->A->velocity + glm::cross(collision->A->angularVelocity, contactInfo.rA);
			glm::vec3 pB = collision->B->velocity + glm::cross(collision->B->angularVelocity, contactInfo.rB);

			//Mess with this and the slop value to control rests.
			if(glm::dot(collision->collisionNormal, (pA - pB)) < 0.5f)
				contactInfo.restitutionImpulse = 0.0f;
			else
				contactInfo.restitutionImpulse = (glm::dot(collision->collisionNormal, (pA - pB)) * 0.2f) / contactInfo.denominator;

			if(contactInfo.penetrationDepth < 0.002)	//Allowed "overlap" on resting objects.
				contactInfo.penetrationImpulse = 0.0f;
			else
				contactInfo.penetrationImpulse = (contactInfo.penetrationDepth * 120.0f * 0.09f) / contactInfo.denominator;

			collision->contactPoints.push_back(contactInfo);
		}
	}

	collision->contactManifold.verts = newPoints;	//Save the final contact manifold.

	return true;	//The narrow phase is complete, so let it know that the collision info can be saved for later.
}

//Carries out a Broad Phase Collision check using AABBs & Sweep and Prune. AABBs and Sweep and Prune from lecture notes.
void BroadPhaseCollision(void)
{
	list<AxisAlignedBoundingBox*> potentialColliders = list<AxisAlignedBoundingBox*>();
	collisionPairs.clear();
	broadPhaseCollisionsFound = 0;
	narrowPhaseCollisionsFound = 0;

	for(int i = 0; i < sceneBodies.size(); i++)
		sceneBodies[i]->GetBoundingBox()->isColliding = false;

	SortSweepPrune();

	//If no collisions are found on a given axis, then quit out early. We must overlap on ALL THREE for a collision.
	/*if(!SweepAxis(xAxisStart, xAxisEnd, collisionTable, 0) || !SweepAxis(yAxisStart, yAxisEnd, collisionTable, 1) || !SweepAxis(zAxisStart, zAxisEnd, collisionTable, 2));
		return;*/

	SweepAxis(xAxisStart, xAxisEnd, collisionTable, 0);
	SweepAxis(yAxisStart, yAxisEnd, collisionTable, 1);
	SweepAxis(zAxisStart, zAxisEnd, collisionTable, 2);

	for(int i = 0; i < sceneBodies.size()-1; i++)
	{
		for(int j = i+1; j < sceneBodies.size(); j++)
		{
			if(collisionTable[i][j] == 3)
			{
				//We've got a collision!
				sceneBodies[i]->GetBoundingBox()->isColliding = true;
				sceneBodies[j]->GetBoundingBox()->isColliding = true;

				CollisionInfo collision = CollisionInfo();
				collision.aID = i;
				collision.bID = j;

				++broadPhaseCollisionsFound;

				if(NarrowPhaseCollision(sceneBodies[i], sceneBodies[j], &collision))
				{
					++narrowPhaseCollisionsFound;
					collisionPairs.push_back(collision);
				}
			}

			collisionTable[i][j] = 0;	//No matter what we reset it to zero for the next frame.
		}
	}

	const float CONVERGED_EPSILON = 0.0001f;
	for(int i = 0; i < 100; i++)
	{
		for(int j = 0; j < collisionPairs.size(); j++)
		{
			CollisionInfo* collisionInfo = &collisionPairs[j];
			
			//if(collisionInfo->converged)	//If we've converged early, we don't need to waste time calculating all the values.
			//	continue;

			collisionInfo->largestImpulseDifference = -std::numeric_limits<float>::max();
			for(int k = 0; k < collisionPairs[j].contactManifold.verts.size(); k++)
			{
				ContactInfo* contactInfo = &collisionInfo->contactPoints[k];
				
				glm::vec3 pA = collisionInfo->A->velocity + glm::cross(collisionInfo->A->angularVelocity, contactInfo->rA);
				glm::vec3 pB = collisionInfo->B->velocity + glm::cross(collisionInfo->B->angularVelocity, contactInfo->rB);
				
				float impulseMagnitude = -glm::dot(collisionInfo->collisionNormal, (pA - pB)); 	//The negative of the relative velocity.
		
				impulseMagnitude = impulseMagnitude / contactInfo->denominator;
				impulseMagnitude -= contactInfo->restitutionImpulse;
				impulseMagnitude -= contactInfo->penetrationImpulse;  

				float oldImpulse = contactInfo->totalImpulse;
				contactInfo->totalImpulse += impulseMagnitude;
				contactInfo->totalImpulse = glm::min(contactInfo->totalImpulse, 0.0f);

				float impulseDiff = contactInfo->totalImpulse - oldImpulse;
				glm::vec3 impulseVector = impulseDiff * collisionInfo->collisionNormal;

				if(collisionInfo->largestImpulseDifference < impulseDiff)
					collisionInfo->largestImpulseDifference = impulseDiff;

				collisionInfo->A->momentum += impulseVector;
				collisionInfo->B->momentum -= impulseVector;

				collisionInfo->A->angularMomentum += glm::cross(contactInfo->rA, impulseVector);
				collisionInfo->B->angularMomentum -= glm::cross(contactInfo->rB, impulseVector);

				glm::mat3 rotationMat = glm::mat3_cast(collisionInfo->A->orientation);
				glm::mat3 inverseTensor = glm::transpose(rotationMat) * collisionInfo->A->inverseInertialTensor * rotationMat;
				collisionInfo->A->velocity = collisionInfo->A->momentum * collisionInfo->A->inverseMass;
				collisionInfo->A->angularVelocity = collisionInfo->A->angularMomentum * inverseTensor;

				rotationMat = glm::mat3_cast(collisionInfo->B->orientation);
				inverseTensor = glm::transpose(rotationMat) * collisionInfo->B->inverseInertialTensor * rotationMat;
				collisionInfo->B->velocity = collisionInfo->B->momentum * collisionInfo->B->inverseMass;
				collisionInfo->B->angularVelocity = collisionInfo->B->angularMomentum * inverseTensor;
			}

			if(collisionInfo->largestImpulseDifference <= CONVERGED_EPSILON)
				collisionInfo->converged = true;
		}
	}

	//If detecting collisions for debug viewing, then trigger it if one is found.
	if(pauseOnCollisions && collisionPairs.size() > 0)
	{
		curGameState = GameState::Paused;
		curWorldView = WorldView::Collision;
	}

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

	if(curWorldView != WorldView::Collision && curGameState != GameState::Paused)
	{
		if(elapsedTimeStep > 0.25)
			elapsedTimeStep = 0.25;				//Maximum frame time allowed if we want to avoid the spiral of death.

		accumulator += elapsedTimeStep;
		simpleCamera->Update(elapsedTimeStep);

		while(accumulator >= dt)
		{
			//t = time, dt = deltaTime.
			for(int i = 0; i < sceneBodies.size(); i++)
				sceneBodies[i]->Update(t, dt);

			BroadPhaseCollision();

			#pragma region Barrel Mode Boundary Collisions
			if(curDemoMode == DemoMode::Barrel)
			{
				for(int i = 0; i < sceneBodies.size(); i++)
				{
					BodyState* state = sceneBodies[i]->GetCurrentState();
					if((state->position.x < -10.0f && state->velocity.x < 0.0f) || (state->position.x > 10.0f && state->velocity.x > 0.0f))
						state->momentum.x *= -1;
	
					if((state->position.y < -10.0f  && state->velocity.y < 0.0f) || (state->position.y > 10.0f && state->velocity.y > 0.0f))
						state->momentum.y *= -1;

					if((state->position.z < -10.0f  && state->velocity.z < 0.0f) || (state->position.z > 10.0f && state->velocity.z > 0.0f))
						state->momentum.z *= -1;
				}
			}
			#pragma endregion
		
			accumulator -= dt;
			t += dt;
		}	

		//Remaining "time" in the current frame, to interpolate our current physics step by.
		interpolationAlpha = accumulator / dt;
	}

	glutPostRedisplay();
}

void RenderScene(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(curWorldView == WorldView::Collision || curWorldView == WorldView::Wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	#pragma region Render Scene Objects
	glUseProgram(bodyShader->GetShaderID());
	for(int i = 0; i < sceneBodies.size(); i++)
		sceneBodies[i]->Render(simpleCamera->GetViewMatrix(), simpleCamera->GetProjMatrix(), interpolationAlpha);
	#pragma endregion

	glUseProgram(lineShader->GetShaderID());
	glm::mat4 identityMat = simpleCamera->GetProjMatrix() * simpleCamera->GetViewMatrix();
	glUniformMatrix4fv(glGetUniformLocation(lineShader->GetShaderID(), "gWVP"), 1, GL_FALSE, &identityMat[0][0]);
	glLineWidth(1.0f);
	glPointSize(10);

	#pragma region Render AABBs
	if(curWorldView != WorldView::Simple)	//We don't render AABBs in simple mode.
	{
		for(int i = 0; i < sceneBodies.size(); i++)
		{
			if(sceneBodies[i]->GetBoundingBox()->isColliding)
				glColor3f(1.0f, 0.0f, 0.0f);
			else
				glColor3f(0.0f, 1.0f, 0.0f);

			sceneBodies[i]->RenderAABB();
		}
	}
	#pragma endregion

	#pragma region Render Barrel Boundary
	if(curDemoMode == DemoMode::Barrel)
	{
		glColor3f(0.0f, 0.3f, 0.3f);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glBegin(GL_QUADS);
			glVertex3f(-10.0f, 10.0f, -10.0f);
			glVertex3f(10.0f, 10.0f, -10.0f);
			glVertex3f(10.0f, 10.0f, 10.0f);
			glVertex3f(-10.0f, 10.0f, 10.0f);

			glVertex3f(-10.0f, -10.0f, -10.0f);
			glVertex3f(10.0f, -10.0f, -10.0f);
			glVertex3f(10.0f, -10.0f, 10.0f);
			glVertex3f(-10.0f, -10.0f, 10.0f);

			glVertex3f(-10.0f, -10.0f, -10.0f);
			glVertex3f(-10.0f, -10.0f, 10.0f);
			glVertex3f(-10.0f, 10.0f, 10.0f);
			glVertex3f(-10.0f, 10.0f, -10.0f);

			glVertex3f(10.0f, -10.0f, -10.0f);
			glVertex3f(10.0f, -10.0f, 10.0f);
			glVertex3f(10.0f, 10.0f, 10.0f);
			glVertex3f(10.0f, 10.0f, -10.0f);

			glVertex3f(-10.0f, -10.0f, 10.0f);
			glVertex3f(10.0f, -10.0f, 10.0f);
			glVertex3f(10.0f, 10.0f, 10.0f);
			glVertex3f(-10.0f, 10.0f, 10.0f);

			glVertex3f(-10.0f, -10.0f, -10.0f);
			glVertex3f(10.0f, -10.0f, -10.0f);
			glVertex3f(10.0f, 10.0f, -10.0f);
			glVertex3f(-10.0f, 10.0f, -10.0f);

		glEnd();
	}
	#pragma endregion

	if(curWorldView == WorldView::Collision)
	{
		const float collisionNormalSize = 0.75f;
		glm::vec3 coreCollisionPoint = glm::vec3();

		for(int j = 0; j < collisionPairs.size(); j++)
		{
			if(highlightedCollision != -1 && highlightedCollision != j)
				continue;

			#pragma region Render Contact Points
			glColor3f(1.0f, 0.0f, 0.0f);
			glBegin(GL_POINTS);	

			for(int i = 0; i < collisionPairs[j].contactManifold.verts.size(); i++)
			{
				coreCollisionPoint += collisionPairs[j].contactManifold.verts[i];
				glVertex3f(collisionPairs[j].contactManifold.verts[i].x, collisionPairs[j].contactManifold.verts[i].y, collisionPairs[j].contactManifold.verts[i].z);
			}

			coreCollisionPoint = glm::vec3(coreCollisionPoint.x / collisionPairs[j].contactManifold.verts.size(), coreCollisionPoint.y / collisionPairs[j].contactManifold.verts.size(), coreCollisionPoint.z / collisionPairs[j].contactManifold.verts.size());
			glEnd();
			#pragma endregion

			#pragma region Render Collision Normals
			glLineWidth(2);
			glColor3f(0.0f, 1.0f, 0.0f);

			glBegin(GL_LINES);
			glVertex3f(coreCollisionPoint.x, coreCollisionPoint.y, coreCollisionPoint.z);
			glVertex3f(coreCollisionPoint.x + (collisionPairs[j].collisionNormal.x * collisionNormalSize), coreCollisionPoint.y + (collisionPairs[j].collisionNormal.y * collisionNormalSize), coreCollisionPoint.z + (collisionPairs[j].collisionNormal.z * collisionNormalSize));
			glEnd();

			glColor3f(0.0f, 0.0f, 1.0f);
			glBegin(GL_LINES);
			glVertex3f(coreCollisionPoint.x, coreCollisionPoint.y, coreCollisionPoint.z);
			glVertex3f(coreCollisionPoint.x - (collisionPairs[j].collisionNormal.x * collisionNormalSize), coreCollisionPoint.y - (collisionPairs[j].collisionNormal.y * collisionNormalSize), coreCollisionPoint.z - (collisionPairs[j].collisionNormal.z * collisionNormalSize));
			glEnd();
			#pragma endregion

		}
	}

	if(curWorldView == WorldView::Full)
	{
		glColor3f(0.0f, 1.0f, 0.0f);
		#pragma region Render Collision Normals
		for(int i = 0; i < sceneBodies.size(); i++)
			sceneBodies[i]->RenderNormals();
		#pragma endregion
	}

	#pragma region Print Collision Information
	stringstream ss;
	ss << "Collision Response - Collision Resolution using Cumulative Impulses\n";
	if(broadPhaseCollisionsFound > 0)
	{
		if(narrowPhaseCollisionsFound > 0)
			ss << broadPhaseCollisionsFound << " Broad Phase collisions have been found, of which " << narrowPhaseCollisionsFound << " were also Narrow Phase collisions." << endl;
		else
			ss << broadPhaseCollisionsFound << " Broad Phase collisions have been found, but none of them were also Narrow Phase collisions." << endl;
	}
	else
		ss << "No collisions found." << endl;

	ss << "Pause on Collisions: " << (pauseOnCollisions ? "ON" : "OFF") << endl;
	ss << "World View: ";
	switch(curWorldView)
	{
		case WorldView::Collision:
			ss << "Collision";
			break;
		case WorldView::Full:
			ss << "Show AABBs & Collision Normals";
			break;
		case WorldView::Regular:
			ss << "Show AABBs";
			break;
		case WorldView::Simple:
				ss << "Basic";
			break;
		case WorldView::Wireframe:
			ss << "Wireframe";
			break;
	}
	ss << endl << "Game State: " << (curGameState == GameState::Paused || curWorldView == WorldView::Collision ? "PAUSED" : "RUNNING") << endl;

	if(curDemoMode == DemoMode::Barrel)
		ss << "Gravity: " << (gravityEnabled ? "ON" : "OFF") << endl;

	ss << endl;
	for(int j = 0; j < collisionPairs.size(); j++)
	{
		if(highlightedCollision != j)
			continue;

		if(collisionPairs[j].contactManifold.verts.size() > 0)
		{
			ss << "\nCollision Information for " << collisionPairs[j].aID << " and " << collisionPairs[j].bID << ".\n";
			ss << "Collision Normal - (" << collisionPairs[j].collisionNormal.x << ", " << collisionPairs[j].collisionNormal.y << ", " << collisionPairs[j].collisionNormal.z << ")\n";
	
			if(collisionPairs[j].contactManifold.verts.size() == 1)
				ss << "Collision Point - (" << collisionPairs[j].contactManifold.verts[0].x << ", " << collisionPairs[j].contactManifold.verts[0].y << ", " << collisionPairs[j].contactManifold.verts[0].z << ")\n";
			else
			{
				ss << "Collision Points\n";
				for(int i = 0; i < collisionPairs[j].contactManifold.verts.size(); i++)
				{
					ss << "(" << collisionPairs[j].contactManifold.verts[i].x << ", " << collisionPairs[j].contactManifold.verts[i].y << ", " << collisionPairs[j].contactManifold.verts[i].z << ")\n";
				}
			}
		}

		cout << endl << endl;
	}

	glActiveTexture(GL_TEXTURE0);
	glDisable(GL_TEXTURE_2D);

	string text = ss.str();
	glUseProgram(0);
	glColor3f(1.0f, 1.0f, 1.0f);
	glRasterPos2f(-1.0f, 0.95f);	
	glutBitmapString(GLUT_BITMAP_HELVETICA_18, (const unsigned char*)text.c_str());
	
	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_2D);
	#pragma endregion

    glutSwapBuffers();
}

void InitialiseBarrelDemo(void)
{
	AxisAlignedBoundingBox::boxes = 0;

	for(int i = 0; i < sceneBodies.size(); i++)
		delete(sceneBodies[i]);

	sceneBodies.clear();
	const int objectCount = 8;

	for(int i = 0; i < objectCount; i++)
	{
		sceneBodies.push_back(new BasicCube(
			glm::vec3(GetRandomValue(-10.0f, 10.0f), GetRandomValue(-10.0f, 10.0f), GetRandomValue(-10.0f, 10.0f)),
			glm::vec3(GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f)),
			glm::vec3(GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f)),
			GetRandomValue(1.0f, 10.0f), false));
	}

	for(int i = 0; i < objectCount; i++)
	{
		sceneBodies.push_back(new Pyramid(
			glm::vec3(GetRandomValue(-10.0f, 10.0f), GetRandomValue(-10.0f, 10.0f), GetRandomValue(-10.0f, 10.0f)),
			glm::vec3(GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f)),
			glm::vec3(GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f)),
			glm::vec3(), GetRandomValue(1.0f, 10.0f), false));
	}

	for(int i = 0; i < objectCount; i++)
	{
		sceneBodies.push_back(new Cuboid(
			glm::vec3(GetRandomValue(-10.0f, 10.0f), GetRandomValue(-10.0f, 10.0f), GetRandomValue(-10.0f, 10.0f)),
			glm::vec3(GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f)),
			glm::vec3(GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f), GetRandomValue(-100.0f, 100.0f)),
			GetRandomValue(1.0f, 10.0f), GetRandomValue(1.0f, 3.0f), GetRandomValue(1.0f, 5.0f), GetRandomValue(1.0f, 3.0f), false));
	}

	sceneBodies.push_back(new Plane(glm::vec3(0.0f, -9.5f, 0.0f), 10.0f));

	simpleCamera->distance = 30.0f;
	simpleCamera->verticalOffset = 0.0f;
	simpleCamera->spinEnabled = true;

	curSelectedPoly = sceneBodies[0];
	highlightedCollision = -1;

	#pragma region Initialise Sweep and Prune
	xAxisStart.clear();
	yAxisStart.clear();
	zAxisStart.clear();

	xAxisEnd.clear();
	yAxisEnd.clear();
	zAxisEnd.clear();

	for(int i = 0; i < sceneBodies.size(); i++)
		AddToSweepPrune(sceneBodies[i]->GetBoundingBox());

	collisionTable = vector<vector<int>>();
	collisionTable.resize(sceneBodies.size());
	for(int i = 0; i < sceneBodies.size(); i++)
	{
		collisionTable[i].resize(sceneBodies.size());
		for(int j = 0; j < sceneBodies.size(); j++)
			collisionTable[i][j] = 0;
	}

	SortSweepPrune();
	#pragma endregion

	curDemoMode = DemoMode::Barrel;
}

void InitialiseStackingDemo(void)
{
	AxisAlignedBoundingBox::boxes = 0;

	for(int i = 0; i < sceneBodies.size(); i++)
		delete(sceneBodies[i]);

	sceneBodies.clear();

	sceneBodies.push_back(new BasicCube(glm::vec3(10.0f, 0.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, false));

	sceneBodies.push_back(new BasicCube());
	sceneBodies.push_back(new BasicCube(glm::vec3(0.0f, 2.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new BasicCube(glm::vec3(0.0f, 4.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new BasicCube(glm::vec3(0.0f, 6.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new BasicCube(glm::vec3(0.0f, 8.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new Pyramid(glm::vec3(0.0f, 10.0f, 0.0f), glm::vec3(), glm::vec3(), glm::vec3(), 1.0f, true));

	sceneBodies.push_back(new BasicCube(glm::vec3(4.0f, 0.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new BasicCube(glm::vec3(4.0f, 2.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new BasicCube(glm::vec3(4.0f, 4.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new BasicCube(glm::vec3(4.0f, 6.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new Pyramid(glm::vec3(4.0f, 8.0f, 0.0f), glm::vec3(), glm::vec3(), glm::vec3(), 1.0f, true));

	sceneBodies.push_back(new BasicCube(glm::vec3(-4.0f, 0.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new BasicCube(glm::vec3(-4.0f, 2.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new BasicCube(glm::vec3(-4.0f, 4.0f, 0.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new Pyramid(glm::vec3(-4.0f, 6.0f, 0.0f), glm::vec3(), glm::vec3(), glm::vec3(), 1.0f, true));

	sceneBodies.push_back(new BasicCube(glm::vec3(0.0f, 0.0f, -4.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new BasicCube(glm::vec3(0.0f, 2.0f, -4.0f), glm::vec3(), glm::vec3(), 1.0f, true));
	sceneBodies.push_back(new Pyramid(glm::vec3(0.0f, 4.0f, -4.0f), glm::vec3(), glm::vec3(), glm::vec3(), 1.0f, true));

	sceneBodies.push_back(new Plane(glm::vec3(0.0f, -3.0f, 0.0f), 50.0f));

	simpleCamera->distance = 15.0f;
	simpleCamera->verticalOffset = 10.0f;
	simpleCamera->camSpin = 0.0f;
	simpleCamera->spinEnabled = false;

	curSelectedPoly = sceneBodies[0];
	highlightedCollision = -1;

	#pragma region Initialise Sweep and Prune
	xAxisStart.clear();
	yAxisStart.clear();
	zAxisStart.clear();

	xAxisEnd.clear();
	yAxisEnd.clear();
	zAxisEnd.clear();

	for(int i = 0; i < sceneBodies.size(); i++)
		AddToSweepPrune(sceneBodies[i]->GetBoundingBox());

	collisionTable = vector<vector<int>>();
	collisionTable.resize(sceneBodies.size());
	for(int i = 0; i < sceneBodies.size(); i++)
	{
		collisionTable[i].resize(sceneBodies.size());
		for(int j = 0; j < sceneBodies.size(); j++)
			collisionTable[i][j] = 0;
	}

	SortSweepPrune();
	#pragma endregion

	curDemoMode = DemoMode::Stacking;
}

void InitialiseScene(void)
{
	prevFrameTimeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	
	bodyShader = new Shader("vTexture.shader", "fTexture.shader");
	lineShader = new Shader("vLine.shader", "fLine.shader");
	simpleCamera = new SimpleCamera();

	InitialiseBarrelDemo();
}

void HandleRegularInput(unsigned char curKey, int x, int y)
{
	glm::vec3 focusPoint;
	switch (curKey)
	{
		case ' ':
			simpleCamera->spinEnabled = !simpleCamera->spinEnabled;
			break;

		case '-':
			highlightedCollision = glm::max(highlightedCollision-1, -1);
			if(highlightedCollision >= 0)
			{
				simpleCamera->distance = 5.0f;
				simpleCamera->verticalOffset = 1.0f;

				for(int i = 0; i < collisionPairs[highlightedCollision].contactPoints.size(); i++)
					focusPoint += collisionPairs[highlightedCollision].contactPoints[i].position;

				focusPoint /= collisionPairs[highlightedCollision].contactPoints.size();
				simpleCamera->SetNewTarget(focusPoint);
			}
			break;

		case '=':
			++highlightedCollision;
			if(highlightedCollision >= collisionPairs.size())
				highlightedCollision = collisionPairs.size() - 1;
			
			simpleCamera->distance = 5.0f;
			simpleCamera->verticalOffset = 1.0f;

			for(int i = 0; i < collisionPairs[highlightedCollision].contactPoints.size(); i++)
				focusPoint += collisionPairs[highlightedCollision].contactPoints[i].position;

			focusPoint /= collisionPairs[highlightedCollision].contactPoints.size();
			simpleCamera->SetNewTarget(focusPoint);
			break;

		case '0':
			if(curDemoMode == DemoMode::Barrel)
			{
				gravityEnabled = !gravityEnabled;
				for(int i = 0; i < sceneBodies.size(); i++)
				{
					sceneBodies[i]->SetGravity(gravityEnabled);
					sceneBodies[i]->SetDamping(gravityEnabled);
				}
			}
			break;

		case '1':
			if(curWorldView != WorldView::Full && curWorldView != WorldView::Simple)
			{
				simpleCamera->spinEnabled = false;

				if(curDemoMode == DemoMode::Barrel)
				{
					simpleCamera->distance = 30.0f;
					simpleCamera->verticalOffset = 0.0f;
				}
				else
				{
					simpleCamera->distance = 15.0f;
					simpleCamera->verticalOffset = 10.0f;
				}

				simpleCamera->camSpin = 0.0f;
				simpleCamera->SetNewTarget(glm::vec3());
			}
			curWorldView = WorldView::Regular;
			break;

		case '2':
			if(curWorldView != WorldView::Regular && curWorldView != WorldView::Simple)
			{
				simpleCamera->spinEnabled = false;
				
				if(curDemoMode == DemoMode::Barrel)
				{
					simpleCamera->distance = 30.0f;
					simpleCamera->verticalOffset = 0.0f;
				}
				else
				{
					simpleCamera->distance = 15.0f;
					simpleCamera->verticalOffset = 10.0f;
				}

				simpleCamera->camSpin = 0.0f;
				simpleCamera->SetNewTarget(glm::vec3());
			}
			curWorldView = WorldView::Full;
			break;

		case '3':
			if(curWorldView != WorldView::Regular && curWorldView != WorldView::Full)
			{
				simpleCamera->spinEnabled = false;
				
				if(curDemoMode == DemoMode::Barrel)
				{
					simpleCamera->distance = 30.0f;
					simpleCamera->verticalOffset = 0.0f;
				}
				else
				{
					simpleCamera->distance = 15.0f;
					simpleCamera->verticalOffset = 10.0f;
				}

				simpleCamera->camSpin = 0.0f;
				simpleCamera->SetNewTarget(glm::vec3());
			}
			curWorldView = WorldView::Simple;
			break;

		case '4':
			pauseOnCollisions = !pauseOnCollisions;
			break;

		case '5':
			if(curGameState == GameState::Paused)
				curGameState = GameState::Running;
			else
				curGameState = GameState::Paused;
			break;

		case '7':
			if(curDemoMode != DemoMode::Barrel)
				InitialiseBarrelDemo();
			break;

		case '8':
			if(curDemoMode != DemoMode::Stacking)
				InitialiseStackingDemo();
			break;

		case 'w':
			curSelectedPoly->GetCurrentState()->momentum += glm::vec3(0.0f, 1.0f, 0.0f);
			break;

		case 's':
			curSelectedPoly->GetCurrentState()->momentum -= glm::vec3(0.0f, 1.0f, 0.0f);
			break;

		case 'a':
			curSelectedPoly->GetCurrentState()->momentum += glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'd':
			curSelectedPoly->GetCurrentState()->momentum -= glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'i':
			curSelectedPoly->GetCurrentState()->angularMomentum += glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'k':
			curSelectedPoly->GetCurrentState()->angularMomentum -= glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'j':
			curSelectedPoly->GetCurrentState()->angularMomentum -= glm::vec3(0.0f, 0.0f, 1.0f);
			break;

		case 'l':
			curSelectedPoly->GetCurrentState()->angularMomentum += glm::vec3(0.0f, 0.0f, 1.0f);
			break;
    }
}

void InitialiseCallbacks()
{
	glutKeyboardFunc(HandleRegularInput);
	glutIdleFunc(UpdateScene);							// Tell glut where the update function is, to execute when system is idle (after drawing completed).
	glutDisplayFunc(RenderScene);						// Tell glut where the display function is
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH)-1280)/2, (glutGet(GLUT_SCREEN_HEIGHT)-720)/2);
	glutCreateWindow("Real Time Physics - Lab 5 - Collision Response");

	glEnable(GL_DEPTH_TEST);
	InitialiseCallbacks();
	
	GLenum res = glewInit();
    if (res != GLEW_OK)
	{
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		return 1;
    }

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	ilInit();
	iluInit();
	ilutInit();
	ilutRenderer(ILUT_OPENGL);

	InitialiseScene();
	glutMainLoop();
	return 0;
}

