#include "shared.h"
#include "WaterSim.h"

int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;
vector<Plane*> planes;

Shader particleShader;
Shader planeShader;
Shader sphereShader;
Shader lineShader;

Camera particleCamera;

WaterSim* fluidSimulation;
GLUquadricObj* refSphere;

void HandleRegularInput(unsigned char curKey, int x, int y)
{
	curKey = tolower(curKey);
	//fluidSimulation.HandleRegularInput(curKey);
	//switch(curKey)
	//{
	//	case ' ':
	//		wallEnabled = false;
	//		return;

	//	case 'a':
	//		debuggingEnabled = !debuggingEnabled;
	//		return;
	//}
}

void update(void)
{
	//Get the time elapsed since the last frame, in fractions of a second.
	int timeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	int deltaTime = timeSinceStart - prevFrameTimeSinceStart;
	prevFrameTimeSinceStart = timeSinceStart;

	elapsedTimeStep = ((float)deltaTime / 1000.0f);
	particleCamera.Update(elapsedTimeStep);

	fluidSimulation->Update();    
	glutPostRedisplay();
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glm::mat4 viewMat = particleCamera.GetViewMatrix();
	glm::mat4 projMat = particleCamera.GetProjMatrix();

	//Draw the planes.
	int aID = planeShader.GetShaderID();
	glUseProgram(aID);

	glUniformMatrix4fv(glGetUniformLocation(aID, "camMat"), 1, GL_FALSE, glm::value_ptr(viewMat));
	glUniformMatrix4fv(glGetUniformLocation(aID, "perMat"), 1, GL_FALSE, glm::value_ptr(projMat));

	glBegin(GL_QUADS);
	glColor4f(0.4f, 0.4f, 0.4f, 1.0f);

	for(int i = 0; i < planes.size(); i++)
		planes[i]->Render();

	glEnd();

	GLuint shaderID = sphereShader.GetShaderID();
	glUseProgram(shaderID);

	glm::mat4 finalTransform;
	int transformLoc = glGetUniformLocation(shaderID, "mTransform");
	int colourLoc = glGetUniformLocation(shaderID, "vColor");

	glUniformMatrix4fv(glGetUniformLocation(shaderID, "mView"), 1, GL_FALSE, glm::value_ptr(viewMat));
	glUniformMatrix4fv(glGetUniformLocation(shaderID, "mProjection"), 1, GL_FALSE, glm::value_ptr(projMat));
	
	glUniformMatrix4fv(glGetUniformLocation(shaderID, "gWVP"), 1, GL_FALSE, glm::value_ptr(projMat * viewMat));

	fluidSimulation->Render(shaderID);
	glutSwapBuffers();
}

//Load our shaders, particles, planes and forces.
void initialise(void)
{
	srand(0);

	prevFrameTimeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	
	particleShader = Shader("vParticle.shader", "fParticle.shader");
	sphereShader = Shader("vSphere.shader", "fSphere.shader");
	planeShader = Shader("vPlane.shader", "fPlane.shader");
	lineShader = Shader("vLine.shader", "fLine.shader");
	Plane* flatPlane = new Plane();
	flatPlane->GenerateVertices(glm::vec3(0.0f, -20.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	planes.push_back(flatPlane);

	fluidSimulation = new WaterSim();

	refSphere = gluNewQuadric();
	gluQuadricTexture(refSphere, TRUE);
	gluQuadricNormals(refSphere, TRUE);
}


void InitialiseCallbacks()
{
	glutKeyboardFunc(HandleRegularInput);
	//glutIdleFunc(UpdateScene);						// Tell glut where the update function is, to execute when system is idle (after drawing completed).
	//glutDisplayFunc(Render);						// Tell glut where the display function is
	//glutMouseFunc(MouseCB);
	//glutMotionFunc(MouseMotionCB);
	//glutPassiveMotionFunc(MousePassiveCB);
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
	glutInitWindowSize(1280, 720);

	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH)-1280)/2, (glutGet(GLUT_SCREEN_HEIGHT)-720)/2);
	glutCreateWindow("Particle System Window");

	glEnable(GL_DEPTH_TEST);						// Enable the depth buffer to display the hand correctly.
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glutIdleFunc(update);							// Tell glut where the update function is, to execute when system is idle (after drawing completed).
	glutDisplayFunc(display);						// Tell glut where the display function is
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	InitialiseCallbacks();

    GLenum res = glewInit();						// A call to glewInit() must be done after glut is initialized!
    if (res != GLEW_OK) 
	{
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		return 1;
    }

	initialise();
	glutMainLoop();

	return 0;
}