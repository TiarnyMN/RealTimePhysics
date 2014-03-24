#include "shared.h"
#include <unordered_map>

int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;

vector<Particle*> particles;
vector<Plane*> planes;

Shader particleShader;
Shader planeShader;
Shader sphereShader;
Shader lineShader;

Camera particleCamera;

GLUquadricObj* refSphere;

class Hasher
{
public:
	int operator()(glm::vec3 val)
	{
		int x = val.x * 73856093;
		int y = val.y * 19349663;
		int z = val.z * 83492791;
		
		return (x ^ y ^ z);
	}

	bool operator()(glm::vec3 first, glm::vec3 second)
	{
		return first == second;
	}
};

std::unordered_multimap<glm::vec3, Particle*, Hasher, Hasher> hashMap;

bool wallEnabled = true;
bool debuggingEnabled = false;

//int particleCount = 6;
//int maxParticles;
//float kernelSize;

void HandleRegularInput(unsigned char curKey, int x, int y)
{
	curKey = tolower(curKey);
	switch(curKey)
	{
		case ' ':
			wallEnabled = false;
			return;

		case 'a':
			debuggingEnabled = !debuggingEnabled;
			return;
	}
}

//Generates a random float between minVal and maxVal
float GetRandomValue(float minVal, float maxVal)
{
	float rndVal = ((float)rand()) / (float)RAND_MAX;
	return (rndVal * (maxVal - minVal)) + minVal;
}

float CalculateKernal(Particle* curParticle, Particle* compareParticle, float d)
{
	float r2 = glm::distance2(curParticle->position, compareParticle->position);
	float d2 = glm::pow(d, 2.0f);
	float pi = 3.14f;

	return 315.0f / (64.0f * pi * glm::pow(d, 9.0f)) * glm::pow((d2 - r2), 3.0f);
}

glm::vec3 CalculateNormalKernal(glm::vec3 r, float h)
{
	if(glm::length(r) == 0.0f)
		return glm::vec3();

	return (-945.0f / (32.0f * glm::pi<float>() * glm::pow(h, 9.0f))) * (r / glm::length(r)) * glm::pow(glm::pow(h, 2.0f) - glm::pow(glm::length(r), 2.0f), 2.0f);
}

float CalculateSurfaceKernal(glm::vec3 r, float h)
{
	if(glm::length(r) == 0.0f)
		return 0.0f;

	return (-945.0f / (32.0f * glm::pi<float>() * glm::pow(h, 9.0f))) * (glm::pow(h, 2.0f) - glm::pow(glm::length(r), 2.0f)) * (3.0f * glm::pow(h, 2.0f) - 7.0f * glm::pow(glm::length(r), 2.0f));
}

glm::vec3 CalculatePressureKernal(glm::vec3 r, float h)
{
	if(glm::length(r) == 0.0f)
		return glm::vec3();

	return (-45.0f  / (glm::pi<float>() * glm::pow(h, 6.0f))) * ((r / glm::length(r)) * glm::pow((h - glm::length(r)), 2.0f));
}

float CalculateViscosityKernal(glm::vec3 r, float h)
{
	return	(45.0f / (glm::pi<float>() * glm::pow(h, 6.0f))) * (h - glm::length(r));
}

void update(void)
{
	//Get the time elapsed since the last frame, in fractions of a second.
	int timeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	int deltaTime = timeSinceStart - prevFrameTimeSinceStart;
	prevFrameTimeSinceStart = timeSinceStart;

	elapsedTimeStep = ((float)deltaTime / 1000.0f);

	particleCamera.Update(elapsedTimeStep);

	//https://www8.cs.umu.se/kurser/5DV058/VT10/lectures/Lecture8.pdf
	//http://image.diku.dk/projects/media/kelager.06.pdf

	hashMap.clear();
	for(int i = 0; i < PARTICLE_COUNT; i++)
	{
		hashMap.insert(std::pair<glm::vec3, Particle*>(glm::floor(particles[i]->position / SMOOTHING_KERNEL_SIZE), particles[i]));
		Particle* curParticle = particles[i];
		curParticle->neighbours.clear();

		for(int x = -1; x < 2; x++)
			for(int y = -1; y < 2; y++)
				for(int z = -1; z < 2; z++)
				{
					auto iterator = hashMap.equal_range(glm::vec3(x, y, z) + glm::floor(particles[i]->position / SMOOTHING_KERNEL_SIZE));
					for(auto it = iterator.first; it != iterator.second; it++)
					{
						Particle* compParticle = (*it).second;

						if(glm::distance(curParticle->position, compParticle->position) <= SMOOTHING_KERNEL_SIZE)
						{
							curParticle->neighbours.push_back((*it).second);

							if(compParticle != curParticle)
								compParticle->neighbours.push_back(curParticle);
						}
					}
				}
	}

	//5.6.2 (ii + iii) Compute Density and Pressure
	//for(int i = 0; i < PARTICLE_COUNT; i++)
	//{
	//	Particle* curParticle = particles[i];
	//	curParticle->neighbours.clear();

	//	for(int x = -1; x < 2; x++)
	//		for(int y = -1; y < 2; y++)
	//			for(int z = -1; z < 2; z++)
	//			{
	//				auto iterator = hashMap.equal_range(glm::vec3(x, y, z) + glm::floor(particles[i]->position / SMOOTHING_KERNEL_SIZE));
	//				for(auto it = iterator.first; it != iterator.second; it++)
	//				{
	//					curParticle->neighbours.push_back((*it).second);
	//				}
	//			}
	//}


	for(int i = 0; i < PARTICLE_COUNT; i++)
	{
		Particle* curParticle = particles[i];
		curParticle->density = 0.0f;
		curParticle->pressure = 0.0f;

		for(int j = 0; j < curParticle->neighbours.size(); j++)
		{
			Particle* compareParticle = curParticle->neighbours[j];
			if(debuggingEnabled)
				cout << i << " (" << j << ") - " << compareParticle->mass * CalculateKernal(curParticle, compareParticle, SMOOTHING_KERNEL_SIZE) << endl;
			curParticle->density += compareParticle->mass * CalculateKernal(curParticle, compareParticle, SMOOTHING_KERNEL_SIZE);
		}
		if(debuggingEnabled)
			cout << i << " - " << curParticle->density << endl << endl;
		curParticle->pressure = STIFFNESS_CONSTANT * (curParticle->density);	//Compute pressure p(i).
	}

	//cout << endl;

	//5.6.3 (ii) Compute the pressure force density
	for(int i = 0; i < PARTICLE_COUNT; i++)
	{
		Particle* curParticle = particles[i];
		curParticle->forcePressure = glm::vec3();
		curParticle->forceViscosity = glm::vec3();
		curParticle->forceSurface = glm::vec3();
		curParticle->color = glm::vec3(0.0f, 0.0f, 1.0f);

		for(int j = 0; j < curParticle->neighbours.size(); j++)
		{
			Particle* compareParticle = curParticle->neighbours[j];

			if(curParticle == compareParticle)
				continue;

			//Calculating the pressure gradient.			
			float firstPart = (curParticle->mass / curParticle->density) * ((curParticle->pressure + compareParticle->pressure) * 0.5f);
			//Looking at using the pressure multipler (mass / density) * ((pressure[i] + pressure[j]) / 2)
			curParticle->forcePressure -= firstPart * CalculatePressureKernal(curParticle->position - compareParticle->position, SMOOTHING_KERNEL_SIZE);
			
			glm::vec3 firstHalf = (compareParticle->velocity - curParticle->velocity) * (compareParticle->mass / compareParticle->density);
			curParticle->forceViscosity += firstHalf * CalculateViscosityKernal(curParticle->position - compareParticle->position, SMOOTHING_KERNEL_SIZE);			//Should this be -= ???
		}

		curParticle->forceViscosity *= VISCOSITY_CONSTANT;
		curParticle->forceInternal = curParticle->forcePressure + curParticle->forceViscosity;

		glm::vec3 surfaceNormal = glm::vec3();	//Inward Surface Normal

		for(int j = 0; j < curParticle->neighbours.size(); j++)
		{
			Particle* compareParticle = curParticle->neighbours[j];
			surfaceNormal += (compareParticle->mass / compareParticle->density) * CalculateNormalKernal(curParticle->position - compareParticle->position, SMOOTHING_KERNEL_SIZE);
		}

		if(glm::length(surfaceNormal) > 0.5f)
		{
			//Colour the particle to indicate it is on the surface.
			//curParticle->color = glm::vec3(1.0f, 1.0f, 1.0f);

			float laplacianC = 0.0f;
			for(int j = 0; j < curParticle->neighbours.size(); j++)
			{
				Particle* compareParticle = curParticle->neighbours[j];
				laplacianC += CalculateSurfaceKernal(curParticle->position - compareParticle->position, SMOOTHING_KERNEL_SIZE);
			}

			curParticle->forceSurface = -SIGMA * laplacianC * glm::normalize(surfaceNormal);
		}
		
	}

	//5.6.4 (i) + 5.6.5 (ii) Add gravitational force and update particle acceleration. 
	for(int i = 0; i < PARTICLE_COUNT; i++)
	{
		Particle* curParticle = particles[i];
		curParticle->forceExternal = glm::vec3(0.0f, -9.81f, 0.0f) * curParticle->density;	//Fix this later.

		curParticle->acceleration = (curParticle->forceInternal + curParticle->forceExternal + curParticle->forceSurface) / curParticle->density;
	}

	//5.6.5 (iiu) Use leapfrog integration to update particle position.
	for(int i = 0; i < PARTICLE_COUNT; i++)
	{
		Particle* curParticle = particles[i];
		glm::vec3 newVel = curParticle->oldVelocity + FLUID_TIME_STEP * curParticle->acceleration;
		//glm::vec3 halfVel = (curParticle->oldVelocity + newVel) / 2.0f;

		glm::vec3 oldPos = curParticle->position;
		curParticle->position = curParticle->position + newVel * FLUID_TIME_STEP;

		curParticle->velocity = newVel + (FLUID_TIME_STEP / 2.0f) * curParticle->acceleration;
		curParticle->oldVelocity = newVel;
	}

	for(int i = 0; i < PARTICLE_COUNT; i++)
	{
		Particle* curParticle = particles[i];
		if(curParticle->position.y < -HORIZONTAL_BOUNDS && curParticle->velocity.y < 0.0f)
		{
			curParticle->position.y = -HORIZONTAL_BOUNDS;
			curParticle->velocity.y = -curParticle->velocity.y * DAMPENING_STRENGTH;
			curParticle->oldVelocity.y = -curParticle->oldVelocity.y * DAMPENING_STRENGTH;
		}

		if(wallEnabled && curParticle->position.x < -HORIZONTAL_BOUNDS)
		{
			curParticle->position.x = -HORIZONTAL_BOUNDS;
			curParticle->velocity.x = -curParticle->velocity.x * DAMPENING_STRENGTH;
			curParticle->oldVelocity.x = -curParticle->oldVelocity.x * DAMPENING_STRENGTH;
		}
		else if(curParticle->position.x < -HORIZONTAL_BOUNDS * 3.0f)
		{
			curParticle->position.x = -HORIZONTAL_BOUNDS * 3.0f;
			curParticle->velocity.x = -curParticle->velocity.x * DAMPENING_STRENGTH;
			curParticle->oldVelocity.x = -curParticle->oldVelocity.x * DAMPENING_STRENGTH;
		}
		else if(curParticle->position.x > HORIZONTAL_BOUNDS)
		{
			curParticle->position.x = HORIZONTAL_BOUNDS;
			curParticle->velocity.x = -curParticle->velocity.x * DAMPENING_STRENGTH;
			curParticle->oldVelocity.x = -curParticle->oldVelocity.x * DAMPENING_STRENGTH;
		}

		if(curParticle->position.z < -HORIZONTAL_BOUNDS)
		{
			curParticle->position.z = -HORIZONTAL_BOUNDS;
			curParticle->velocity.z = -curParticle->velocity.z * DAMPENING_STRENGTH;
			curParticle->oldVelocity.z = -curParticle->oldVelocity.z * DAMPENING_STRENGTH;
		}
		//else if(curParticle->position.z < -HORIZONTAL_BOUNDS * 3.0f)
		//{
		//	curParticle->position.z = -HORIZONTAL_BOUNDS * 3.0f;
		//	curParticle->velocity.z = -curParticle->velocity.z * DAMPENING_STRENGTH;
		//	curParticle->oldVelocity.z = -curParticle->oldVelocity.z * DAMPENING_STRENGTH;
		//}
		else if(curParticle->position.z > HORIZONTAL_BOUNDS)
		{
			curParticle->position.z = HORIZONTAL_BOUNDS;
			curParticle->velocity.z = -curParticle->velocity.z * DAMPENING_STRENGTH;
			curParticle->oldVelocity.z = -curParticle->oldVelocity.z * DAMPENING_STRENGTH;
		}
	}
    
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

	int shaderID = sphereShader.GetShaderID();
	glUseProgram(shaderID);

	glm::mat4 finalTransform;
	int transformLoc = glGetUniformLocation(shaderID, "mTransform");
	int colourLoc = glGetUniformLocation(shaderID, "vColor");

	glUniformMatrix4fv(glGetUniformLocation(shaderID, "mView"), 1, GL_FALSE, glm::value_ptr(viewMat));
	glUniformMatrix4fv(glGetUniformLocation(shaderID, "mProjection"), 1, GL_FALSE, glm::value_ptr(projMat));
	
	glUniformMatrix4fv(glGetUniformLocation(shaderID, "gWVP"), 1, GL_FALSE, glm::value_ptr(projMat * viewMat));

	//glPointSize(10.0f);
	//glBegin(GL_POINTS);

	for(int i = 0; i < particles.size(); i++)
	{
		//if(particles[i]->color != glm::vec3(1.0f, 1.0f, 1.0f))
		//	continue;

		finalTransform = glm::translate(particles[i]->GetPosition());
		//finalTransform = glm::translate(glm::vec3());
		glUniformMatrix4fv(transformLoc, 1, GL_FALSE, &finalTransform[0][0]);
		glUniform3f(colourLoc, particles[i]->color.x, particles[i]->color.y, particles[i]->color.z);

		float particleRenderSize = glm::pow((3.0f * particles[i]->mass) / (4.0f * glm::pi<float>() * particles[i]->density), 1.0f / 3.0f);
		//glColor3f(particles[i]->color.x, particles[i]->color.y, particles[i]->color.z);
		//glVertex3f(particles[i]->GetPosition().x, particles[i]->GetPosition().y, particles[i]->GetPosition().z);		
		gluSphere(refSphere, particleRenderSize, 5, 5);
	}
	//glEnd();

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

	particles.reserve(PARTICLE_COUNT);

	for(int x = 0; x < PARTICLES_PER_DIMENSION; x++)
		for(int y = 0; y < PARTICLES_PER_DIMENSION; y++)
			for(int z = 0; z < PARTICLES_PER_DIMENSION; z++)
			{
				Particle* particle = new Particle();
				particle->SetPosition(GetRandomValue(-HORIZONTAL_BOUNDS, HORIZONTAL_BOUNDS), GetRandomValue(0.0f, HORIZONTAL_BOUNDS), GetRandomValue(-HORIZONTAL_BOUNDS, HORIZONTAL_BOUNDS));
				particle->SetVelocity(0.0f, 0.0f, 0.0f);
				particles.push_back(particle);
			}

	Plane* flatPlane = new Plane();
	flatPlane->GenerateVertices(glm::vec3(0.0f, -20.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	planes.push_back(flatPlane);

	refSphere = gluNewQuadric();
	gluQuadricTexture(refSphere, TRUE);
	gluQuadricNormals(refSphere, TRUE);

	cout << "REST DENSITY - " << REST_DENSITY << endl;
	cout << "PARTICLES PER DIMENSION - " << PARTICLES_PER_DIMENSION << endl;
	cout << "PARTICLE COUNT - " << PARTICLE_COUNT << endl;

	cout << "STIFFNESS CONSTANT - " << STIFFNESS_CONSTANT << endl;
	cout << "VISCOSITY CONSTANT - " << VISCOSITY_CONSTANT << endl;
	cout << "AVERAGE KERNEL PARTICLES - " << AVERAGE_KERNEL_PARTICLES << endl;
	cout << "SURFACE THRESHOLD - " << SURFACE_THRESHOLD << endl;

	cout << "SIGMA - " << SIGMA << endl;

	cout << "FLUID TIME STEP - " << FLUID_TIME_STEP << endl;
	cout << "VERTICAL BOUNDS - " << VERTICAL_BOUNDS << endl;
	cout << "HORIZONTAL BOUNDS - " << HORIZONTAL_BOUNDS << endl;
	cout << "BOX VOLUME - " << BOX_VOLUME << endl;

	cout << "DAMPENING STRENGTH - " << DAMPENING_STRENGTH << endl;
	cout << "PARTICLE MASS - " << PARTICLE_MASS << endl;

	cout << "SMOOTHING KERNEL SIZE - " << SMOOTHING_KERNEL_SIZE << endl;
	cout << "PARTICLE SIZE - " << PARTICLE_SIZE << endl;
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

	//glutKeyboardFunc(HandleRegularInput);
	//glutMotionFunc(HandleMouseMovement);
	//glutPassiveMotionFunc(HandleMouseMovement);
	//glutMouseFunc(HandleMouseClick);

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