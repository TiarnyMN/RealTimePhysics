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

#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 720

SimpleCamera* simpleCamera = NULL;

int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;

Shader* basicShader;
RigidBody* basicRigidBody;

RigidBody* basicSpringBody;
RigidBody* heavySpringBody;
RigidBody* tenseSpringBody;

RigidBody* freeSpringBody;

double accumulator = 0.0f;

glm::vec2 mousePos;
glm::vec3 mouseWorldPos;

enum DemoMode
{
	Basic,
	Springs,
	Free
};

DemoMode curDemoMode = DemoMode::Basic;

bool movingCube = false;

void HandleRegularInput(unsigned char curKey, int x, int y)
{
	switch (curKey) 
	{
		case 'w':
			basicRigidBody->curState.momentum += glm::vec3(0.0f, 0.0f, 1.0f);
			break;

		case 's':
			basicRigidBody->curState.momentum -= glm::vec3(0.0f, 0.0f, 1.0f);
			break;

		case 'a':
			basicRigidBody->curState.momentum += glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'd':
			basicRigidBody->curState.momentum -= glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'i':
			basicRigidBody->curState.angularMomentum += glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'k':
			basicRigidBody->curState.angularMomentum -= glm::vec3(1.0f, 0.0f, 0.0f);
			break;

		case 'j':
			basicRigidBody->curState.angularMomentum -= glm::vec3(0.0f, 0.0f, 1.0f);
			break;

		case 'l':
			basicRigidBody->curState.angularMomentum += glm::vec3(0.0f, 0.0f, 1.0f);
			break;

		case ' ':
			if(curDemoMode == DemoMode::Basic)
				curDemoMode = DemoMode::Springs;
			else if(curDemoMode == DemoMode::Springs)
				curDemoMode = DemoMode::Free;
			else
				curDemoMode = DemoMode::Basic;
			break;
    }
}

void HandleMouseInput(int button, int state, int x, int y)
{
	if(button == GLUT_LEFT_BUTTON)
	{
		if(state == GLUT_DOWN)
			movingCube = true;
		else
			movingCube = false;
	}
}

void HandleMouseMovement(int x, int y)
{
	mousePos = glm::vec2(x, y);
}

void UpdateScene(void)
{
	double t = 0.0f;
	const double dt = 1.0f / 120.0f;

	int timeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	int deltaTime = timeSinceStart - prevFrameTimeSinceStart;
	prevFrameTimeSinceStart = timeSinceStart;

	elapsedTimeStep = ((float)deltaTime / 1000.0f);

	if(elapsedTimeStep > 0.25)
		elapsedTimeStep = 0.25;		//Maximum frame time allowed if we want to avoid the spiral of death.

	//Casts a ray from the mouse into the world
	glm::vec3 rayNDS = glm::vec3();
	rayNDS.x = mousePos.x / glutGet(GLUT_WINDOW_WIDTH);
	rayNDS.x = -1.0f + (rayNDS.x * 2.0f);
	 
	rayNDS.y = mousePos.y / glutGet(GLUT_WINDOW_HEIGHT);
	rayNDS.y = 1.0f - (rayNDS.y * 2.0f);

	rayNDS.z = 1.0f;

	glm::vec4 rayClip = glm::vec4(rayNDS.x, rayNDS.y, -1.0f, 1.0f);
	glm::vec4 rayEye = glm::inverse(simpleCamera->GetProjMatrix()) * rayClip;
	rayEye = glm::vec4(rayEye.x, rayEye.y, -1.0f, 0.0f);

	rayEye = glm::inverse(simpleCamera->GetViewMatrix()) * rayEye;

	glm::vec3 rayWorld = glm::vec3(rayEye.x, rayEye.y, rayEye.z);
	rayWorld = glm::normalize(rayWorld);	//A ray from the mouse position, through the camera, and into the world.

	mouseWorldPos = simpleCamera->GetPos() + rayWorld * 50.0f;	//We set the "world position" of the cursor to be in front of the camera, along the projected ray.

	accumulator += elapsedTimeStep;						//We add the current time delta to the accumulator, to know how far to advance our physics time step.

	if(movingCube && curDemoMode == DemoMode::Free)
		freeSpringBody->curState.position = mouseWorldPos;

	while(accumulator >= dt)	//We update our physics more often than we render, so we need to update as much as many steps as the accumulator has gained.
	{
		//Integrate our physics using RK4
		if(curDemoMode == DemoMode::Basic)
		{
			basicRigidBody->Update(t, dt, glm::vec3(-100.0f, -10.0f, 15.0f), 5.0f, 1.0f);
		}
		else if (curDemoMode == DemoMode::Springs)
		{
			basicSpringBody->Update(t, dt, glm::vec3(-10.0f, 10.0f, 15.0f), 2.0f, 1.0f);
			heavySpringBody->Update(t, dt, glm::vec3(0.0f, 10.0f, 15.0f), 1.0f, 1.0f);
			tenseSpringBody->Update(t, dt, glm::vec3(10.0f, 10.0f, 15.0f), 5.0f, 1.0f);
		}
		else
		{
			freeSpringBody->Update(t, dt, glm::vec3(0.0f, 0.0f, 15.0f), 5.0f, 1.0f);
		}
		
		accumulator -= dt;	//Take our timestep (1/120) away from the accumulator to see if we need another update.
		t += dt;
	}

	//To smooth out the rendering, we'll interpolate between our current and previous state based on the time remaining until the next physics step update.
	const double alpha = accumulator / dt;
	if(curDemoMode == DemoMode::Basic)
	{
		basicRigidBody->Catchup(alpha);
	}
	else if (curDemoMode == DemoMode::Springs)
	{
		basicSpringBody->Catchup(alpha);
		heavySpringBody->Catchup(alpha);
		tenseSpringBody->Catchup(alpha);
	}
	else
	{
		freeSpringBody->Catchup(alpha);
	}

	glutPostRedisplay();
}

void RenderScene(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glUseProgram(basicShader->shaderID);

	glLoadIdentity();
	glLineWidth(2.0f);	

	//Render our cubes.
	if(curDemoMode == DemoMode::Basic)
	{
		basicRigidBody->Render(simpleCamera->GetViewMatrix(), simpleCamera->GetProjMatrix());
	}
	else if (curDemoMode == DemoMode::Springs)
	{
		basicSpringBody->Render(simpleCamera->GetViewMatrix(), simpleCamera->GetProjMatrix());
		heavySpringBody->Render(simpleCamera->GetViewMatrix(), simpleCamera->GetProjMatrix());
		tenseSpringBody->Render(simpleCamera->GetViewMatrix(), simpleCamera->GetProjMatrix());
	}
	else
	{
		freeSpringBody->Render(simpleCamera->GetViewMatrix(), simpleCamera->GetProjMatrix());
	}

	GLint shaderID;
	glGetIntegerv(GL_CURRENT_PROGRAM, &shaderID);

	glm::mat4 identityMat = simpleCamera->GetProjMatrix() * simpleCamera->GetViewMatrix();
	glUniformMatrix4fv(glGetUniformLocation(shaderID, "gWVP"), 1, GL_FALSE, &identityMat[0][0]);

	//Render the springs.
	if(curDemoMode != DemoMode::Basic)
	{
		glBegin(GL_LINES);
		if(curDemoMode == DemoMode::Springs)
		{
			glVertex3f(-10.0f, 10.0f, 15.0f);
			glVertex3f(basicSpringBody->cube->GetTransformedVertex(0).position.x, basicSpringBody->cube->GetTransformedVertex(0).position.y, basicSpringBody->cube->GetTransformedVertex(0).position.z);

			glVertex3f(0.0f, 10.0f, 15.0f);
			glVertex3f(heavySpringBody->cube->GetTransformedVertex(0).position.x, heavySpringBody->cube->GetTransformedVertex(0).position.y, heavySpringBody->cube->GetTransformedVertex(0).position.z);

			glVertex3f(10.0f, 10.0f, 15.0f);
			glVertex3f(tenseSpringBody->cube->GetTransformedVertex(0).position.x, tenseSpringBody->cube->GetTransformedVertex(0).position.y, tenseSpringBody->cube->GetTransformedVertex(0).position.z);
		}
		else
		{
			glVertex3f(0.0f, 0.0f, 15.0f);
			glVertex3f(freeSpringBody->cube->GetTransformedVertex(0).position.x, freeSpringBody->cube->GetTransformedVertex(0).position.y, freeSpringBody->cube->GetTransformedVertex(0).position.z);
		}
		glEnd();
	}

	//Gather and then print out the debug information to the screen.
	stringstream ss;
	ss << "Rigid Bodies - Impulses & Springs\n";
	ss << "Current Mode - ";

	if(curDemoMode == DemoMode::Basic)
		ss << "Simple Impulses" << endl;
	else if(curDemoMode == DemoMode::Springs)
		ss << "Springs" << endl;
	else
		ss << "Free Control" << endl;

	ss << endl;

	if(curDemoMode == DemoMode::Basic)
	{
		ss << "Rigid Body 1" << endl;
		ss << "Mass: " << basicRigidBody->curState.mass << endl;
		ss << "Momentum: (" << basicRigidBody->curState.momentum.x << ", " << basicRigidBody->curState.momentum.y << ", " << basicRigidBody->curState.momentum.z << ")" << endl;
		ss << "Velocity: (" << basicRigidBody->curState.velocity.x << ", " << basicRigidBody->curState.velocity.y << ", " << basicRigidBody->curState.velocity.z << ")" << endl;

		ss << "Angular Momentum: (" << basicRigidBody->curState.angularMomentum.x << ", " << basicRigidBody->curState.angularMomentum.y << ", " << basicRigidBody->curState.angularMomentum.z << ")" << endl;
		ss << "Angular Velocity: (" << basicRigidBody->curState.angularVelocity.x << ", " << basicRigidBody->curState.angularVelocity.y << ", " << basicRigidBody->curState.angularVelocity.z << ")" << endl;
	}
	else if(curDemoMode == DemoMode::Springs)
	{
		ss << "Rigid Body 1" << endl;
		ss << "Mass: " << basicSpringBody->curState.mass << endl;
		ss << "Momentum: (" << basicSpringBody->curState.momentum.x << ", " << basicSpringBody->curState.momentum.y << ", " << basicSpringBody->curState.momentum.z << ")" << endl;
		ss << "Velocity: (" << basicSpringBody->curState.velocity.x << ", " << basicSpringBody->curState.velocity.y << ", " << basicSpringBody->curState.velocity.z << ")" << endl;

		ss << "Angular Momentum: (" << basicSpringBody->curState.angularMomentum.x << ", " << basicSpringBody->curState.angularMomentum.y << ", " << basicSpringBody->curState.angularMomentum.z << ")" << endl;
		ss << "Angular Velocity: (" << basicSpringBody->curState.angularVelocity.x << ", " << basicSpringBody->curState.angularVelocity.y << ", " << basicSpringBody->curState.angularVelocity.z << ")" << endl;

		ss << "Spring Stiffness: 2.0f\nSpring Damping Coefficient: 1.0f\n" << endl << endl;

		ss << "Rigid Body 2" << endl;
		ss << "Mass: " << heavySpringBody->curState.mass << endl;
		ss << "Momentum: (" << heavySpringBody->curState.momentum.x << ", " << heavySpringBody->curState.momentum.y << ", " << heavySpringBody->curState.momentum.z << ")" << endl;
		ss << "Velocity: (" << heavySpringBody->curState.velocity.x << ", " << heavySpringBody->curState.velocity.y << ", " << heavySpringBody->curState.velocity.z << ")" << endl;

		ss << "Angular Momentum: (" << heavySpringBody->curState.angularMomentum.x << ", " << heavySpringBody->curState.angularMomentum.y << ", " << heavySpringBody->curState.angularMomentum.z << ")" << endl;
		ss << "Angular Velocity: (" << heavySpringBody->curState.angularVelocity.x << ", " << heavySpringBody->curState.angularVelocity.y << ", " << heavySpringBody->curState.angularVelocity.z << ")" << endl;

		ss << "Spring Stiffness: 1.0f\nSpring Damping Coefficient: 1.0f\n" << endl;

		ss << "Rigid Body 3" << endl;
		ss << "Mass: " << tenseSpringBody->curState.mass << endl;
		ss << "Momentum: (" << tenseSpringBody->curState.momentum.x << ", " << tenseSpringBody->curState.momentum.y << ", " << tenseSpringBody->curState.momentum.z << ")" << endl;
		ss << "Velocity: (" << tenseSpringBody->curState.velocity.x << ", " << tenseSpringBody->curState.velocity.y << ", " << tenseSpringBody->curState.velocity.z << ")" << endl;

		ss << "Angular Momentum: (" << tenseSpringBody->curState.angularMomentum.x << ", " << tenseSpringBody->curState.angularMomentum.y << ", " << tenseSpringBody->curState.angularMomentum.z << ")" << endl;
		ss << "Angular Velocity: (" << tenseSpringBody->curState.angularVelocity.x << ", " << tenseSpringBody->curState.angularVelocity.y << ", " << tenseSpringBody->curState.angularVelocity.z << ")" << endl;

		ss << "Spring Stiffness: 5.0f\nSpring Damping Coefficient: 1.0f\n" << endl;
	}
	else
	{
		ss << "Rigid Body 1" << endl;
		ss << "Mass: " << freeSpringBody->curState.mass << endl;
		ss << "Momentum: (" << freeSpringBody->curState.momentum.x << ", " << freeSpringBody->curState.momentum.y << ", " << freeSpringBody->curState.momentum.z << ")" << endl;
		ss << "Velocity: (" << freeSpringBody->curState.velocity.x << ", " << freeSpringBody->curState.velocity.y << ", " << freeSpringBody->curState.velocity.z << ")" << endl;

		ss << "Angular Momentum: (" << freeSpringBody->curState.angularMomentum.x << ", " << freeSpringBody->curState.angularMomentum.y << ", " << freeSpringBody->curState.angularMomentum.z << ")" << endl;
		ss << "Angular Velocity: (" << freeSpringBody->curState.angularVelocity.x << ", " << freeSpringBody->curState.angularVelocity.y << ", " << freeSpringBody->curState.angularVelocity.z << ")" << endl;

		ss << "Spring Stiffness: 5.0f\nSpring Damping Coefficient: 1.0f\n" << endl << endl;
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

    glutSwapBuffers();
}

//Creates our shader, camera and rigid bodies.
void InitialiseScene(void)
{
	prevFrameTimeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	
	basicShader = new Shader("vTexture.shader", "fTexture.shader");
	simpleCamera = new SimpleCamera();

	basicRigidBody = new RigidBody();
	basicSpringBody = new RigidBody();

	heavySpringBody = new RigidBody();
	heavySpringBody->curState.mass = 3.0f;
	heavySpringBody->curState.inverseMass = 1.0f / heavySpringBody->curState.mass;

	tenseSpringBody = new RigidBody();
	freeSpringBody = new RigidBody();
}

void InitialiseCallbacks()
{
	glutKeyboardFunc(HandleRegularInput);
	glutPassiveMotionFunc(HandleMouseMovement);
	glutMouseFunc(HandleMouseInput);
	glutIdleFunc(UpdateScene);							// Tell glut where the update function is, to execute when system is idle (after drawing completed).
	glutDisplayFunc(RenderScene);						// Tell glut where the display function is
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(1280, 720);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH)-1280)/2, (glutGet(GLUT_SCREEN_HEIGHT)-720)/2);
	glutCreateWindow("Real Time Physics - Lab 2 - Rigid Bodies & Springs");

	glEnable(GL_DEPTH_TEST);

	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_MULTISAMPLE);
	glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);

	InitialiseCallbacks();
	
	GLenum res = glewInit();
    if (res != GLEW_OK)
	{
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		return 1;
    }

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	ilInit();
	iluInit();
	ilutInit();
	ilutRenderer(ILUT_OPENGL);

	InitialiseScene();
	glutMainLoop();
	return 0;
}

