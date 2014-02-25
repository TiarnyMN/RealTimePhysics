#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <GL/glew.h>
#include <GL/freeglut.h>

#include "shared.h"
#include "shader.h"
#include "SimpleCamera.h"

#include <algorithm>
#include <vector>

#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 720

SimpleCamera* simpleCamera = NULL;

int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;

Shader* lineShader;
vector<glm::vec3> triVerts;

glm::vec3 pointPos = glm::vec3(5.0f, 0.0f, 0.0f);

int voronoiVert = -1;
pair<int, int> voronoiEdge = pair<int, int>(-1, -1);
glm::vec3 edgePoint;
bool voronoiFace = false;

stringstream ss;

float GetRandomValue(float minVal, float maxVal)
{
	float rndVal = ((float)rand()) / (float)RAND_MAX;
	return (rndVal * (maxVal - minVal)) + minVal;
}

void HandleRegularInput(unsigned char curKey, int x, int y)
{
	switch (curKey) 
	{
		case 'w':
			pointPos.y += 0.1f;
			break;

		case 's':
			pointPos.y -= 0.1f;
			break;

		case 'a':
			pointPos.x -= 0.1f;
			break;

		case 'd':
			pointPos.x += 0.1f;
			break;
    }
}

void CreateVoronoiTriangle(void)
{
	glm::vec3 verts[] = 
	{
		glm::vec3(0.0f, 3.0f, 0.0f),
		glm::vec3(-6.0f, -3.0f, 0.0f),
		glm::vec3(6.0f, -3.0f, 0.0f)
	};

	triVerts = vector<glm::vec3>();
	triVerts.reserve(3);

	for(unsigned int i = 0; i < 3; i++)
		triVerts.push_back(verts[i]);
}

void PerformVoronoiRegionCheck(void)
{
	voronoiVert = -1;
	voronoiEdge = pair<int, int>(-1, -1);
	voronoiFace = false;

	//First check if it is in any of the vertex regions.
	for(unsigned int i = 0; i < triVerts.size(); i++)
	{
		glm::vec3 ax = triVerts[i] - pointPos;
		glm::vec3 ab = triVerts[i] - triVerts[(i + 1) % triVerts.size()];
		glm::vec3 ac = triVerts[i] - triVerts[(i + 2) % triVerts.size()];

		if(glm::dot(ax, ab) <= 0.0f && glm::dot(ax, ac) <= 0.0f)
		{
			voronoiVert = i;

			ss << "Region: Vertex - " << i << "\n";
			ss << "Closest Point: (" << triVerts[i].x << ", " << triVerts[i].y << ", " << triVerts[i].z << ")\n";
			return;
		}
	}

	//Now check if it is in any of the edge regions.
	for(unsigned int i = 0; i < triVerts.size(); i++)
	{
		glm::vec3 bc = triVerts[(i + 1) % triVerts.size()] - triVerts[(i + 2) % triVerts.size()];
		glm::vec3 ba = triVerts[(i + 1) % triVerts.size()] - triVerts[i];
		glm::vec3 bx = triVerts[(i + 1) % triVerts.size()] - pointPos;
		glm::vec3 ax = triVerts[i] - pointPos;
		glm::vec3 ab = triVerts[i] - triVerts[(i + 1) % triVerts.size()];

		if(glm::dot(glm::cross(glm::cross(bc, ba), ba), bx) >= 0.0f && glm::dot(ax, ab) >= 0.0f && glm::dot(bx, ba) >= 0.0f)
		{
			glm::vec3 startPoint = triVerts[i];
			glm::vec3 direction = triVerts[(i + 1) % triVerts.size()] - triVerts[i];
			direction = glm::normalize(direction);

			voronoiEdge = pair<int, int>(i, (i + 1) % triVerts.size());
			edgePoint = startPoint + (glm::dot((pointPos - startPoint), glm::normalize(direction)) * direction);

			ss << "Region: Edge - " << i << " ~ " << (i + 1) % triVerts.size() << "\n";
			ss << "Closest Point: (" << edgePoint.x << ", " << edgePoint.y << ", " << edgePoint.z << ")\n";
			return;
		}
	}

	//Then we are somewhere above the face!
	voronoiFace = true;

	ss << "Region: Face\n";
	ss << "Closest Point: (" << pointPos.x << ", " << pointPos.y << ", " << pointPos.z << ")\n";
	return;
}

void UpdateScene(void)
{
	int timeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	int deltaTime = timeSinceStart - prevFrameTimeSinceStart;
	prevFrameTimeSinceStart = timeSinceStart;
	elapsedTimeStep = ((float)deltaTime / 1000.0f);

	ss << "Voronoi Region Test\n";
	PerformVoronoiRegionCheck();

	glutPostRedisplay();
}

void RenderScene(void)
{
	glUseProgram(lineShader->GetShaderID());
	glm::mat4 identityMat = simpleCamera->GetProjMatrix() * simpleCamera->GetViewMatrix();
	glUniformMatrix4fv(glGetUniformLocation(lineShader->GetShaderID(), "gWVP"), 1, GL_FALSE, &identityMat[0][0]);

	glPointSize(8);
	glLineWidth(4);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glColor3f(1.0f, 1.0f, 1.0f);
	glBegin(GL_POINTS);
		glVertex3f(pointPos.x, pointPos.y, pointPos.z);
	glEnd();
	
	if(voronoiVert != -1)
	{
		glColor3f(0.0f, 1.0f, 0.0f);
		glBegin(GL_POINTS);
			glVertex3f(triVerts[voronoiVert].x, triVerts[voronoiVert].y, triVerts[voronoiVert].z);
		glEnd();

		glColor3f(1.0f, 0.0f, 0.0f);
		glLineWidth(3);
		glBegin(GL_LINES);
			glVertex3f(pointPos.x, pointPos.y, pointPos.z);
			glVertex3f(triVerts[voronoiVert].x, triVerts[voronoiVert].y, triVerts[voronoiVert].z);
		glEnd();
	}
	else if(voronoiEdge != pair<int, int>(-1, -1))
	{
		glColor3f(0.0f, 1.0f, 0.0f);
		glBegin(GL_LINES);
			glVertex3f(triVerts[voronoiEdge.first].x, triVerts[voronoiEdge.first].y, triVerts[voronoiEdge.first].z);
			glVertex3f(triVerts[voronoiEdge.second].x, triVerts[voronoiEdge.second].y, triVerts[voronoiEdge.second].z);
		glEnd();

		glColor3f(1.0f, 0.0f, 0.0f);
		glLineWidth(3);
		glBegin(GL_LINES);
			glVertex3f(pointPos.x, pointPos.y, pointPos.z);
			glVertex3f(edgePoint.x, edgePoint.y, edgePoint.z);
		glEnd();
	}

	if(voronoiFace)
		glColor3f(0.0f, 1.0f, 0.0f);
	else
		glColor3f(0.0f, 0.3f, 0.7f);

	glBegin(GL_TRIANGLES);
		for(unsigned int i = 0; i < triVerts.size(); i++)
			glVertex3f(triVerts[i].x, triVerts[i].y, triVerts[i].z);
	glEnd();

	string screenText = ss.str();
	glUseProgram(0);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	//glRasterPos2f(0.0f, 0.0f);
	glRasterPos2f(-1.0f, 0.95f);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char*)screenText.c_str());

	ss.str(std::string());

    glutSwapBuffers();
}

void InitialiseScene(void)
{
	prevFrameTimeSinceStart = glutGet(GLUT_ELAPSED_TIME);

	lineShader = new Shader("vLine.shader", "fLine.shader");
	simpleCamera = new SimpleCamera();

	glUseProgram(lineShader->GetShaderID());
	glm::mat4 identityMat = simpleCamera->GetProjMatrix() * simpleCamera->GetViewMatrix();
	glUniformMatrix4fv(glGetUniformLocation(lineShader->GetShaderID(), "gWVP"), 1, GL_FALSE, &identityMat[0][0]);

	CreateVoronoiTriangle();
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
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH)-1280)/2, (glutGet(GLUT_SCREEN_HEIGHT)-720)/2);
	glutCreateWindow("Real Time Physics - Lab 4 - Narrow Phase - Voronoi Region");

	glEnable(GL_DEPTH_TEST);

	InitialiseCallbacks();
	
	GLenum res = glewInit();
    if (res != GLEW_OK)
	{
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		return 1;
    }

	glClearColor(0.3f, 0.3f, 0.3f, 1.0f);

	InitialiseScene();
	glutMainLoop();
	return 0;
}

