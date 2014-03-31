#include "shared.h"
#include "WaterSim.h"
#include "MucusSim.h"

int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;
vector<Plane*> planes;

Shader particleShader;
Shader planeShader;
Shader sphereShader;
Shader lineShader;

Camera particleCamera;

FluidSim* fluidSimulation;
FluidSim* mucusSim;
FluidSim* waterSim;
GLUquadricObj* refSphere;

bool showDetailedDebug = true;

void HandleRegularInput(unsigned char curKey, int x, int y)
{
	curKey = tolower(curKey);
	if(!fluidSimulation->HandleRegularInput(curKey))
	{
		switch(curKey)
		{
			case '#':
				particleCamera.spinEnabled = !particleCamera.spinEnabled;
				break;

			case '/':
				if(fluidSimulation == waterSim)
					fluidSimulation = mucusSim;
				else
					fluidSimulation = waterSim;
				break;

			case ']':
				showDetailedDebug = !showDetailedDebug;
				break;
		}
	}
}

void update(void)
{
	//Get the time elapsed since the last frame, in fractions of a second.
	int timeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	int deltaTime = timeSinceStart - prevFrameTimeSinceStart;
	prevFrameTimeSinceStart = timeSinceStart;

	elapsedTimeStep = ((float)deltaTime / 1000.0f);
	particleCamera.Update(elapsedTimeStep, fluidSimulation->container.centre);

	fluidSimulation->Update();    
	glutPostRedisplay();
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	//Draw the planes.
	int shaderID = planeShader.GetShaderID();
	glUseProgram(shaderID);

	glm::mat4 viewMat = particleCamera.GetViewMatrix();
	glm::mat4 projMat = particleCamera.GetProjMatrix();
	glm::mat4 transformMat = glm::translate(glm::vec3());

	glUniformMatrix4fv(glGetUniformLocation(shaderID, "mView"), 1, GL_FALSE, glm::value_ptr(viewMat));
	glUniformMatrix4fv(glGetUniformLocation(shaderID, "mProjection"), 1, GL_FALSE, glm::value_ptr(projMat));
	glUniformMatrix4fv(glGetUniformLocation(shaderID, "mTransform"), 1, GL_FALSE, glm::value_ptr(transformMat));
	glUniform3f(glGetUniformLocation(shaderID, "vColor"), 0.3f, 0.3f, 0.3f);	

	glBegin(GL_QUADS);
	glColor4f(0.4f, 0.4f, 0.4f, 1.0f);

	for(int i = 0; i < planes.size(); i++)
		planes[i]->Render();

	glEnd();

	//Draw everything to do with the fluid simulation.
	fluidSimulation->Render(sphereShader.GetShaderID(), planeShader.GetShaderID(), viewMat, projMat);

	#pragma region Print Debug Information
	stringstream ss;
	ss << "Fluid Simulation - Smoothed Particle Hydrodynamics" << endl;

	if(fluidSimulation == waterSim)
		ss << "Fluid Type: Water" << endl;
	else
		ss << "Fluid Type: Mucus" << endl;

	if(showDetailedDebug)
	{
		fluidSimulation->PrintDebug(ss);

		ss << endl << "CONTROLS" << endl;
		ss << "# ~ Spin Camera | / ~ Swap Fluid Type" << endl;
		ss << "WASD ~ Move Capsule | IK ~ Rotate Capsule" << endl;
		ss << "SPACE ~ Drop Capsule" << endl;
		ss << "0 ~ Default Mode" << endl;
		ss << "1 ~ Show | 2 ~ Hide Surface Particles" << endl;
		ss << "3 ~ Highlight Neighbours | - + ~ Change Highlight" << endl;
		ss << "4 ~ Fat | 5 ~ Thin Inner Capsule" << endl;
		ss << "6 ~ Tall | 7 ~ Short Inner Capsule" << endl;
	}

	ss << endl << "] ~ Toggle Debug Information" << endl;

	glActiveTexture(GL_TEXTURE0);
	glDisable(GL_TEXTURE_2D);

	string text = ss.str();
	glUseProgram(0);
	glColor3f(1.0f, 1.0f, 1.0f);
	glRasterPos2f(-0.99f, 0.95f);	
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char*)text.c_str());
	
	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_2D);
	#pragma endregion

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

	mucusSim = new MucusSim();
	waterSim = new WaterSim();
	fluidSimulation = waterSim;

	refSphere = gluNewQuadric();
	gluQuadricTexture(refSphere, TRUE);
	gluQuadricNormals(refSphere, TRUE);
}


void InitialiseCallbacks()
{
	glutKeyboardFunc(HandleRegularInput);
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
	glClearColor(0.33f, 0.66f, 0.95f, 1.0f);

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