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

enum DemoMode
{
	Basic,
	Springs,
	Free
};

SimpleCamera* simpleCamera = NULL;

int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;

Shader* bodyShader;
Shader* lineShader;

vector<BasicCube*> sceneBodies;
BasicCube* curSelectedCube;

double accumulator = 0.0f;
double interpolationAlpha = 0.0f;

int rigidBodyCount = 32;
int collidingCount;

vector<std::list<AxisAlignedBoundingBox*>::iterator> activeColliders;
vector<vector<int>> collisionTable;

DemoMode curDemoMode = DemoMode::Basic;

vector<AxisAlignedBoundingBox*> xAxisStart;
vector<AxisAlignedBoundingBox*> yAxisStart;
vector<AxisAlignedBoundingBox*> zAxisStart;

vector<AxisAlignedBoundingBox*> xAxisEnd;
vector<AxisAlignedBoundingBox*> yAxisEnd;
vector<AxisAlignedBoundingBox*> zAxisEnd;

bool simulationRunning = true;

float GetRandomValue(float minVal, float maxVal)
{
	float rndVal = ((float)rand()) / (float)RAND_MAX;
	return (rndVal * (maxVal - minVal)) + minVal;
}

void HandleRegularInput(unsigned char curKey, int x, int y)
{
	switch (curKey) 
	{
        case '1':
			simulationRunning = !simulationRunning;
			break;

		case '2':
			simpleCamera->spinEnabled = !simpleCamera->spinEnabled;
			break;

		case '3':
			simpleCamera->distance += 1.0f;
			break;

		case '4':
			simpleCamera->distance -= 1.0f;
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

		case 'q':
			curSelectedCube->GetCurrentState()->momentum += glm::vec3(0.0f, 0.0f, 1.0f);
			break;

		case 'e':
			curSelectedCube->GetCurrentState()->momentum -= glm::vec3(0.0f, 0.0f, 1.0f);
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

void NiceBroadPhase(void)
{
	collidingCount = 0;
	list<AxisAlignedBoundingBox*> potentialColliders = list<AxisAlignedBoundingBox*>();
	bool collisionFound = false;

	for(int i = 0; i < sceneBodies.size(); i++)
		sceneBodies[i]->GetBoundingBox()->isColliding = false;

	SortSweepPrune();
	SweepAxis(xAxisStart, xAxisEnd, collisionTable, 0);
	SweepAxis(yAxisStart, yAxisEnd, collisionTable, 1);
	SweepAxis(zAxisStart, zAxisEnd, collisionTable, 2);

	for(int i = 0; i < sceneBodies.size(); i++)
	{
		for(int j = i; j < sceneBodies.size(); j++)
		{
			if(collisionTable[i][j] == 3)
			{
				//We've got a collision!
				sceneBodies[i]->GetBoundingBox()->isColliding = true;
				sceneBodies[j]->GetBoundingBox()->isColliding = true;
			}

			collisionTable[i][j] = 0;	//No matter what we reset it to zero for the next frame.
		}
	}

	for(int i = 0; i < sceneBodies.size(); i++)
	{
		if(sceneBodies[i]->GetBoundingBox()->isColliding)
			collidingCount++;
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

	if(simulationRunning)
	{
		if(elapsedTimeStep > 0.25)
			elapsedTimeStep = 0.25;		//Maximum frame time allowed if we want to avoid the spiral of death.

		accumulator += elapsedTimeStep;		//This is wonderful.
		simpleCamera->Update(elapsedTimeStep);

		while(accumulator >= dt)	//so, so wonderful.
		{
			//t = time, dt = deltaTime.
			for(int i = 0; i < sceneBodies.size(); i++)
				sceneBodies[i]->Update(t, dt);
		
			accumulator -= dt;
			t += dt;
		}

		NiceBroadPhase();

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

		//Remaining "time" in the current frame, to interpolate our current physics step by.
		interpolationAlpha = accumulator / dt;
	}
	glutPostRedisplay();
}

void RenderScene(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glLoadIdentity();

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

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	stringstream ss;
	ss << "Broad Phase Demo - Sweep & Prune\n";
	ss << collidingCount << " of " << rigidBodyCount << " rigid bodies are currently colliding.\n";

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
	prevFrameTimeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	
	bodyShader = new Shader("vTexture.shader", "fTexture.shader");
	lineShader = new Shader("vLine.shader", "fLine.shader");
	simpleCamera = new SimpleCamera();

	sceneBodies.push_back(new BasicCube());
	for(int i = 1; i < rigidBodyCount; i++)
		sceneBodies.push_back(new BasicCube(glm::vec3(GetRandomValue(-1.0f, 1.0f), GetRandomValue(-1.0f, 1.0f), GetRandomValue(-1.0f, 1.0f)),
											glm::vec3(GetRandomValue(-1.0f, 1.0f), GetRandomValue(-1.0f, 1.0f), GetRandomValue(-1.0f, 1.0f)),
											glm::vec3(GetRandomValue(-1.0f, 1.0f), GetRandomValue(-1.0f, 1.0f), GetRandomValue(-1.0f, 1.0f))));

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
		{
			collisionTable[i][j] = 0;
		}
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
	glutInitWindowSize(1280, 720);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH)-1280)/2, (glutGet(GLUT_SCREEN_HEIGHT)-720)/2);
	glutCreateWindow("Real Time Physics - Lab 3 - Broad Phase Collisions - AABB Sweep and Prune");
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

