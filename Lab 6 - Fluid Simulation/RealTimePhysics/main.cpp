#include "shared.h"
#include "Fluid.h"

int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;

enum DemoMode
{
	Celestial,
	Snow,
	Plume
};

vector<Particle*> particles;
vector<Plane*> planes;

Shader particleShader;
Shader planeShader;

Camera particleCamera;

float epsilon = 0.1f;	//A collision smoothing / buffer value, how close a particle can actually get to a plane before we consider it a collision.
float dampeningStrength = 0.25f;	//We apply a slight dampening, e.g. wind resistance, to reduce the forces applied to the particles over time.
float elasticity = 0.8f;	//Elasticity of surface, lower values = stronger frictional forces / less rebound.

bool spawningParticles = false;

int maxParticles = 500;
int oldestParticle = 0;
int recycledParticles = 0;

const float particleLifetime = 10.0f;
bool particlesShouldAge = false;

DemoMode currentDemoMode = DemoMode::Celestial;

stringstream ss;
glm::vec2 mousePos;

glm::vec3 windVector = glm::vec3();

FluidCube* fluidCube;
int fluidSize = 50;
float particleSpacing = 0.5f;
float particleSize = 5.0f;

//Generates a random float between minVal and maxVal
float GetRandomValue(float minVal, float maxVal)
{
	float rndVal = ((float)rand()) / (float)RAND_MAX;
	return (rndVal * (maxVal - minVal)) + minVal;
}

void update(void)
{
	
	//Get the time elapsed since the last frame, in fractions of a second.
	int timeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	int deltaTime = timeSinceStart - prevFrameTimeSinceStart;
	prevFrameTimeSinceStart = timeSinceStart;

	elapsedTimeStep = ((float)deltaTime / 1000.0f);

	particleCamera.Update(elapsedTimeStep);

	//FluidCubeAddDensity(fluidCube, glm::round(fluidSize / 2.0f), glm::round(fluidSize / 2.0f), glm::round(fluidSize / 2.0f), 100.0f);
	FluidCubeAddDensity(fluidCube, glm::round(fluidSize / 2.0f), glm::round(fluidSize / 2.0f), glm::round(fluidSize / 2.0f), 5000.0f);
	FluidCubeStep(fluidCube);

	////Casts a ray from the mouse into the world
	//glm::vec3 rayNDS = glm::vec3();
	//rayNDS.x = mousePos.x / glutGet(GLUT_WINDOW_WIDTH);
	//rayNDS.x = -1.0f + (rayNDS.x * 2.0f);
	// 
	//rayNDS.y = mousePos.y / glutGet(GLUT_WINDOW_HEIGHT);
	//rayNDS.y = 1.0f - (rayNDS.y * 2.0f);

	//rayNDS.z = 1.0f;

	//glm::vec4 rayClip = glm::vec4(rayNDS.x, rayNDS.y, -1.0f, 1.0f);
	//glm::vec4 rayEye = glm::inverse(particleCamera.GetProjMatrix()) * rayClip;
	//rayEye = glm::vec4(rayEye.x, rayEye.y, -1.0f, 0.0f);

	//rayEye = glm::inverse(particleCamera.GetViewMatrix()) * rayEye;

	//glm::vec3 rayWorld = glm::vec3(rayEye.x, rayEye.y, rayEye.z);
	//rayWorld = glm::normalize(rayWorld);	//A ray from the mouse position, through the camera, and into the world.

	////Finds out if the mouse is hovering over a plane and, if so, which one (if overlapping two we take the closest)
	//int planeID = -1;
	//float shortestDistance = std::numeric_limits<float>::max();
	//float distanceAlongRay = 0.0f;

	////Find out if the ray is intersecting any of the planes, and if it hits multiple find the closest one.
	//for(int i = 0; i < planes.size(); i++)
	//{
	//	distanceAlongRay = glm::dot((planes[i]->position - particleCamera.position), planes[i]->normal) / glm::dot(rayWorld, planes[i]->normal);

	//	if(distanceAlongRay > 0 && distanceAlongRay < shortestDistance)
	//	{
	//		shortestDistance = distanceAlongRay;
	//		planeID = i;
	//	}
	//}

	////If particles should age, then age them. If they have "died", then recycle them!
	//if(particlesShouldAge)
	//{
	//	for(int i = 0; i < particles.size(); i++)
	//	{
	//		particles[i]->age += elapsedTimeStep;

	//		if(particles[i]->age >= particleLifetime)
	//		{
	//			++recycledParticles;
	//			particles[i]->SetPosition(GetRandomValue(-50.0f, 50.0f), GetRandomValue(0.0f, 100.0f), GetRandomValue(-50.0f, 50.0f));
	//			particles[i]->SetVelocity(windVector.x + (-10.0f, 10.0f), windVector.y + 0.0f, windVector.z + (-10.0f, 10.0f));
	//			particles[i]->age = GetRandomValue(-30.0f, 0.0f);
	//		}
	//	}
	//}

	////If we are hovering over a plane, update the force position.
	//if(planeID != -1)
	//{
	//	glm::vec3 pointAtPlane = particleCamera.position + (rayWorld * shortestDistance);
	//	if(curSelectedForce != NULL)
	//	{
	//		curSelectedForce->position = pointAtPlane + glm::vec3(0.0f, 5.0f, 0.0f);

	//		if(curSelectedForce->position.y < -50.0f)
	//			curSelectedForce->position.y = -50.0f;
	//	}

	//	//If active, then add some more particles to the scene!
	//	if(spawningParticles)
	//	{
	//		for(int i = 0; i < 5; i++)
	//		{
	//			particles[oldestParticle]->SetPosition(pointAtPlane.x, pointAtPlane.y + 15.0f, pointAtPlane.z);
	//			particles[oldestParticle]->SetVelocity(GetRandomValue(-20.0f, 20.0f), 0.0f, GetRandomValue(-20.0f, 20.0f));	//Give them a random horizontal velocity so they appear to shoot out of the mouse cursor.
	//			particles[oldestParticle]->age = 0.0f;

	//			oldestParticle = (oldestParticle + 1) % particles.size();
	//			++recycledParticles;
	//		}
	//	}
	//}

	////Mid Point Method - Update once, but using the value as it was in the middle.
	//float midPointTimeStep = elapsedTimeStep / 2;
	//for(int i = 0; i < particles.size(); i++)
	//{
	//	Particle* curParticle = particles[i];
	//	particles[i]->position += particles[i]->velocity * elapsedTimeStep;

	//	particles[i]->velocity.y -= 9.81f * midPointTimeStep;

	//	//Collision handling from lecture slides.
	//	for(int j = 0; j < planes.size(); j++)
	//	{
	//		Plane* curPlane = planes[j];
	//		glm::vec3 relPos = particles[i]->position - curPlane->position;
	//		float dotProd = glm::dot(relPos, curPlane->normal);		//(x-p).n

	//		float dotProd2 = glm::dot(curPlane->normal, particles[i]->velocity);		//(n.v)

	//		if(dotProd < epsilon && dotProd2 < epsilon)
	//		{
	//			glm::vec3 vN = glm::dot(curPlane->normal, particles[i]->velocity) * curPlane->normal;
	//			glm::vec3 vT = particles[i]->velocity - vN;
	//			vT *= 0.75f;	//"Frictional" force

	//			particles[i]->velocity = vT - elasticity * vN;

	//			if(currentDemoMode == DemoMode::Snow)
	//				particles[i]->velocity = glm::vec3(0.0f, 0.0f, 0.0f);
	//		}
	//	}		
	//	particles[i]->velocity -= particles[i]->velocity * dampeningStrength * midPointTimeStep;	//Damp the value slightly.

		////If the scene has active forces, then apply the forces to the particles to pull / push them about!
		//if(currentDemoMode == DemoMode::Celestial || currentDemoMode == DemoMode::Plume)
		//{
		//	for(int j = 0; j < forces.size(); j++)
		//	{
		//		if(forces[j]->forceEnabled)
		//		{
		//			glm::vec3 distFromForce = forces[j]->position - particles[i]->position;
		//			float distSquared = glm::dot(distFromForce, distFromForce);

		//			if(forces[j]->forcePushing)
		//				distSquared *= -1;

		//			float relativeForceStrength = 1 / distSquared;
		//			particles[i]->velocity += (distFromForce * relativeForceStrength * forces[j]->strength * midPointTimeStep); 
		//		}
		//	}
		//}
	//}
    
	glutPostRedisplay();
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glm::mat4 viewMat = particleCamera.GetViewMatrix();
	glm::mat4 projMat = particleCamera.GetProjMatrix();

	//Draw the planes.
	//int aID = planeShader.GetShaderID();
	//glUseProgram(aID);

	//glUniformMatrix4fv(glGetUniformLocation(aID, "camMat"), 1, GL_FALSE, glm::value_ptr(viewMat));
	//glUniformMatrix4fv(glGetUniformLocation(aID, "perMat"), 1, GL_FALSE, glm::value_ptr(projMat));

	//glBegin(GL_QUADS);
	//glColor4f(0.4f, 0.4f, 0.4f, 1.0f);

	//for(int i = 0; i < planes.size(); i++)
	//	planes[i]->Render();

	//glEnd();

	//Draw the particles.
	int bID = particleShader.GetShaderID();
	glUseProgram(bID);

	glUniformMatrix4fv(glGetUniformLocation(bID, "camMat"), 1, GL_FALSE, glm::value_ptr(viewMat));
	glUniformMatrix4fv(glGetUniformLocation(bID, "perMat"), 1, GL_FALSE, glm::value_ptr(projMat));

	glEnable (GL_BLEND);
	glDepthMask (GL_FALSE);

	glPointSize(particleSize);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glBegin(GL_POINTS);

	for(int x = 0; x < fluidSize; x++)
		for(int y = 0; y < fluidSize; y++)
			for(int z = 0; z < fluidSize; z++)
			{
				float density = FluidCubeGetDensity(fluidCube, x, y, z);
				/*if(density > 0.0f)
					cout << density << endl;*/

				glColor4f(density, 0.0f, 0.0f, density);//, 1.0f, 1.0f, 1.0f);
				glVertex3f(x * particleSpacing - (fluidSize * particleSpacing / 2), y * particleSpacing - (fluidSize * particleSpacing / 2), z * particleSpacing - (fluidSize * particleSpacing / 2));
			}

	glEnd();

	glDepthMask (GL_TRUE);
	glDisable (GL_BLEND);

	//for(float x = 0.0f; x < 100.0f; x+=1.0f)
	//{
	//	for(float y = 0.0f; y < 100.0f; y+=1.0f)
	//	{
	//		for(float z = 0.0f; z < 100.0f; z+= 1.0f)
	//		{
	//			glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	//			glm::vec3 pos = particles[x]->GetPosition();
	//			glVertex3f(pos.x, pos.y, pos.z);
	//		}
	//	}
	//}

	

	////Draw the point at which the forces are acting, if in the right mode, as white squares.
	//if(currentDemoMode == DemoMode::Celestial || currentDemoMode == DemoMode::Plume)
	//{
	//	glPointSize(10);
	//	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	//	glBegin(GL_POINTS);	

	//	for(int i = 0; i < forces.size(); i++)
	//	{
	//		if(forces[i]->forceEnabled)
	//			glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	//		else
	//			glColor4f(0.1f, 0.1f, 0.1f, 1.0f);

	//		glVertex3f(forces[i]->position.x, forces[i]->position.y, forces[i]->position.z);
	//	}

	//	glEnd();
	//}

	////We use big "particles" for the snow and plume modes, otherwise smaller ones.
	//if(currentDemoMode != DemoMode::Celestial)
	//	glPointSize(10);
	//else
	//	glPointSize(5);

	//glBegin(GL_POINTS);

	////Draw and colour all the particles.
	//for(int i = 0; i < particles.size(); i++)
	//{
	//	glm::vec3 pos = particles[i]->GetPosition();
	//	glVertex3f(pos.x, pos.y, pos.z);

	//	if(currentDemoMode == DemoMode::Snow)
	//	{
	//		glColor4f(1.0f, 1.0f, 1.0f, (particlesShouldAge ? 1.0f - glm::min(1.0f, (particles[i]->age / particleLifetime)) : 1.0f));
	//	}
	//	else
	//	{
	//		glm::vec3 colorStrength = glm::normalize(pos);
	//		colorStrength *= particles[i]->velocity;

	//		colorStrength.x = glm::max(0.0f, glm::min(colorStrength.x, 1.0f));
	//		colorStrength.y = glm::max(0.0f, glm::min(colorStrength.y, 1.0f));
	//		colorStrength.z = glm::max(0.0f, glm::min(colorStrength.z, 1.0f));

	//		glColor4f(colorStrength.x, colorStrength.y, colorStrength.z, (particlesShouldAge ? 1.0f - glm::min(1.0f, (particles[i]->age / particleLifetime)) : 1.0f));
	//	}
	//}

	//glEnd();

	////Prints out debug information about the current scene state to the screen.
	//ss.str(std::string());
	//ss << "Particle Systems - Particle / Plane Collisions\n";
	//ss << "Mode: ";

	//if(currentDemoMode == DemoMode::Celestial)
	//	ss << "Celestial" << endl;
	//else if(currentDemoMode == DemoMode::Snow)
	//	ss << "Snow Fall" << endl;
	//else
	//	ss << "High Friction Plume" << endl;

	//ss << "Maximum Particles: " << maxParticles << endl;
	//ss << "Recycled Particles: " << recycledParticles << endl;

	//ss << endl;
	//if(currentDemoMode != DemoMode::Snow)
	//{
	//	for(int i = 0; i < forces.size(); i++)
	//	{
	//		ss << "Force " << i+1 << " - " << (forces[i]->forceEnabled ? "ON" : "OFF") << (curSelectedForce == forces[i] ? " - SELECTED" : "") << endl;
	//		ss << "Pos: (" << forces[i]->position.x << ", " << forces[i]->position.y << ", " << forces[i]->position.z << ")" << endl;
	//		if(forces[i]->forceEnabled)
	//		{
	//			ss << "Strength: " << forces[i]->strength << endl;
	//			ss << "Force Mode: " << (forces[i]->forcePushing ? "PUSH" : "PULL") << endl;
	//		}
	//		ss << endl;
	//	}
	//}

	////We need to disable textures if we want the text to render correctly.
	//glActiveTexture(GL_TEXTURE0);
	//glDisable(GL_TEXTURE_2D);

	//string text = ss.str();
	//glUseProgram(0);

	//if(currentDemoMode == DemoMode::Snow)
	//	glColor3f(1.0f, 0.0f, 0.0f);
	//else
	//	glColor3f(1.0f, 1.0f, 1.0f);

	//glRasterPos2f(-1.0f, 0.95f);
	//glutBitmapString(GLUT_BITMAP_HELVETICA_18, (const unsigned char*)text.c_str());

	//glActiveTexture(GL_TEXTURE0);
	//glEnable(GL_TEXTURE_2D);	

	glutSwapBuffers();
}

//Load our shaders, particles, planes and forces.
void initialise(void)
{
	prevFrameTimeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	
	particleShader = Shader("vParticle.shader", "fParticle.shader");
	planeShader = Shader("vPlane.shader", "fPlane.shader");

	particles.reserve(maxParticles);
	for(int i = 0; i < maxParticles; i++)
	{
		Particle* particle = new Particle();
		particle->SetPosition(GetRandomValue(-50.0f, 50.0f), GetRandomValue(0.0f, 100.0f), GetRandomValue(-50.0f, 50.0f));
		particle->SetVelocity(0.0f, 0.0f, 0.0f);
		particles.push_back(particle);
	}

	Plane* flatPlane = new Plane();
	flatPlane->GenerateVertices(glm::vec3(0.0f, -75.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	planes.push_back(flatPlane);

	Plane* rightPlane = new Plane();
	rightPlane->GenerateVertices(glm::vec3(25.0f, 0.0f, 0.0f), glm::vec3(0.5f, 1.0f, 0.0f));
	planes.push_back(rightPlane);

	Plane* leftPlane = new Plane();
	leftPlane->GenerateVertices(glm::vec3(-25.0f, -100.0f, 0.0f), glm::vec3(-0.5f, 1.0f, 0.0f));
	planes.push_back(leftPlane);

	fluidCube = FluidCubeCreate(fluidSize, 0.1f, 100.0f, 0.0001f);
	FluidCubeAddDensity(fluidCube, glm::round(fluidSize / 2.0f), glm::round(fluidSize / 2.0f), glm::round(fluidSize / 2.0f), 5000.0f);

	for(int x = 0; x < fluidSize; x++)
		for(int y = 0; y < fluidSize; y++)
			for(int z = 0; z < fluidSize; z++)
				FluidCubeAddVelocity(fluidCube, x, y, z, 0.0f, -100.0f, 0.0f);
}

//Toggle between pushing and pulling forces, or spawn new particles.
void HandleMouseClick(int button, int state, int x, int y)
{
	if(button == GLUT_LEFT_BUTTON)
	{
		if(state == GLUT_DOWN)
			FluidCubeAddDensity(fluidCube, glm::round(fluidSize / 2.0f), glm::round(fluidSize / 2.0f), glm::round(fluidSize / 2.0f), 5000.0f);
	}

	//if(button == GLUT_LEFT_BUTTON && curSelectedForce != NULL)
	//{
	//	if(state == GLUT_DOWN)
	//		curSelectedForce->forcePushing = true;
	//	else
	//		curSelectedForce->forcePushing = false;
	//}

	//if(button == GLUT_RIGHT_BUTTON)
	//{
	//	if(state == GLUT_DOWN)
	//		spawningParticles = true;
	//	else
	//		spawningParticles = false;
	//}
}

void SwitchDemoMode(DemoMode demoMode)
{
	//if(demoMode == DemoMode::Plume)
	//	elasticity = 0.1f;
	//else
	//	elasticity = 0.8f;

	////Set up the particles for snow fall.
	//if(demoMode == DemoMode::Snow)
	//{
	//	particleCamera.position = glm::vec3(0.0f, -60.0f, 100.0f);
	//	particleCamera.forwardVec = glm::normalize(particleCamera.target - particleCamera.position);
	//	particlesShouldAge = true;

	//	for(int i = 0; i < maxParticles; i++)
	//	{
	//		particles[i]->SetPosition(GetRandomValue(-100.0f, 100.0f), GetRandomValue(0.0f, 100.0f), GetRandomValue(-100.0f, 100.0f));
	//		particles[i]->SetVelocity(GetRandomValue(-10.0f, 10.0f), 0.0f, 0.0f);
	//		particles[i]->age = GetRandomValue(-30.0f, 0.0f);
	//	}
	//}
	//else
	//{
	//	//Set the scene up for plume mode.
	//	if(demoMode == DemoMode::Plume)
	//	{
	//		particleCamera.position = glm::vec3(0.0f, -50.0f, 100.0f);
	//		particleCamera.forwardVec = glm::normalize(particleCamera.target - particleCamera.position);
	//		particlesShouldAge = false;

	//		for(int i = 0; i < forces.size(); i++)
	//			forces[i]->forceEnabled = false;

	//		curSelectedForce = forces[0];
	//	}
	//	else
	//	{
	//		particleCamera.position = glm::vec3(0.0f, 0.0f, 100.0f);
	//		particleCamera.forwardVec = glm::normalize(particleCamera.target - particleCamera.position);
	//		particlesShouldAge = true;

	//		for(int i = 0; i < forces.size(); i++)
	//			forces[i]->forceEnabled = true;

	//		curSelectedForce = NULL;
	//	}

	//	//Either way, we want to randomise their spawn position.
	//	for(int i = 0; i < maxParticles; i++)
	//	{
	//		particles[i]->SetPosition(GetRandomValue(-50.0f, 50.0f), GetRandomValue(0.0f, 100.0f), GetRandomValue(-50.0f, 50.0f));
	//		particles[i]->SetVelocity(0.0f, 0.0f, 0.0f);
	//		particles[i]->age = GetRandomValue(-30.0f, 0.0f);
	//	}
	//}
}

void HandleRegularInput(unsigned char key, int x, int y)
{
	/*switch(key)
	{
		case ' ':
			if(curSelectedForce != NULL)
				curSelectedForce->forceEnabled = !curSelectedForce->forceEnabled;
			break;

		case 'w':
			particleCamera.MoveForwards(elapsedTimeStep);
			break;

		case 's':
			particleCamera.MoveBackwards(elapsedTimeStep);
			break;

		case 'a':
			particleCamera.RotateX(elapsedTimeStep);
			break;

		case 'd':
			particleCamera.RotateX(-elapsedTimeStep);
			break;

		case 'q':
			particleCamera.RotateY(elapsedTimeStep);
			break;

		case 'e':
			particleCamera.RotateY(-elapsedTimeStep);
			break;

		case '0':
			curSelectedForce = NULL;
			break;

		case '1':
			curSelectedForce = forces[0];
			break;

		case '2':
			curSelectedForce = forces[1];
			break;

		case '3':
			curSelectedForce = forces[2];
			break;

		case '4':
			curSelectedForce = forces[3];
			break;

		case '-':
			curSelectedForce->strength -= 100.0f;
			break;

		case '=':
			curSelectedForce->strength += 100.0f;
			break;

		case 'p':
			if(currentDemoMode == DemoMode::Celestial)
				currentDemoMode = DemoMode::Snow;
			else if(currentDemoMode == DemoMode::Snow)
				currentDemoMode = DemoMode::Plume;
			else
				currentDemoMode = DemoMode::Celestial;

			SwitchDemoMode(currentDemoMode);
			break;
	}*/
}

void HandleMouseMovement(int x, int y)
{
	mousePos = glm::vec2(x, y);
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

	glutKeyboardFunc(HandleRegularInput);
	glutMotionFunc(HandleMouseMovement);
	glutPassiveMotionFunc(HandleMouseMovement);
	glutMouseFunc(HandleMouseClick);

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