#include "shared.h"

int prevFrameTimeSinceStart = 0;
GLfloat elapsedTimeStep = 0.0f;

vector<Particle*> particles;
vector<Force*> forces;
Force* curSelectedForce;
vector<Plane*> planes;

Shader particleShader;
Shader planeShader;

Camera particleCamera;

float epsilon = 0.1f;	//A collision smoothing / buffer value, how close a particle can actually get to a plane before we consider it a collision.
float dampeningStrength = 0.25f;	//We apply a slight dampening, e.g. wind resistance, to reduce the forces applied to the particles over time.

bool spawningParticles = false;

//bool snowing = false;

glm::vec2 mousePos;

int maxParticles = 5000;
int oldestParticle = 0;

const float particleLifetime = 10.0f;
bool particlesShouldAge = false;

enum DemoMode
{
	Celestial,
	Snow,
	Plume
};

DemoMode currentDemoMode = DemoMode::Celestial;

float elasticity = 0.8f;	//Elasticity of surface, lower values = stronger frictional forces / less rebound.

//TO DO: Give each Snow particle it's own slightly random velocity...?

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

	/*for(int i = 0; i < forces.size(); i++)
		cout << i << ". " << forces[i]->position.x << " : " << forces[i]->position.y << " : " << forces[i]->position.z << "\n";*/

	//curSelectedForce->strength = glm::sin((float)timeSinceStart) * 1000.0f + 500.0f;
	//cout << curSelectedForce->strength << "\n";

	//Casts a ray from the mouse into the world
	glm::vec3 rayNDS = glm::vec3();
	rayNDS.x = mousePos.x / glutGet(GLUT_WINDOW_WIDTH);
	rayNDS.x = -1.0f + (rayNDS.x * 2.0f);
	 
	rayNDS.y = mousePos.y / glutGet(GLUT_WINDOW_HEIGHT);
	rayNDS.y = 1.0f - (rayNDS.y * 2.0f);

	rayNDS.z = 1.0f;

	glm::vec4 rayClip = glm::vec4(rayNDS.x, rayNDS.y, -1.0f, 1.0f);
	glm::vec4 rayEye = glm::inverse(particleCamera.GenerateProjMatrix()) * rayClip;
	rayEye = glm::vec4(rayEye.x, rayEye.y, -1.0f, 0.0f);

	rayEye = glm::inverse(particleCamera.GenerateViewMatrix()) * rayEye;

	glm::vec3 rayWorld = glm::vec3(rayEye.x, rayEye.y, rayEye.z);
	rayWorld = glm::normalize(rayWorld);	//A ray from the mouse position, through the camera, and into the world.

	//Finds out if the mouse is hovering over a plane and, if so, which one (if overlapping two we take the closest)
	int planeID = -1;
	float shortestDistance = std::numeric_limits<float>::max();
	float distanceAlongRay = 0.0f;

	for(int i = 0; i < planes.size(); i++)
	{
		distanceAlongRay = glm::dot((planes[i]->position - particleCamera.position), planes[i]->normal) / glm::dot(rayWorld, planes[i]->normal);

		if(distanceAlongRay > 0 && distanceAlongRay < shortestDistance)
		{
			shortestDistance = distanceAlongRay;
			planeID = i;
		}
	}

	if(particlesShouldAge)
	{
		for(int i = 0; i < particles.size(); i++)
		{
			particles[i]->age += elapsedTimeStep;

			if(particles[i]->age >= particleLifetime)
			{
				particles[i]->SetPosition(GetRandomValue(-50.0f, 50.0f), GetRandomValue(0.0f, 100.0f), GetRandomValue(-50.0f, 50.0f));
				particles[i]->SetVelocity(GetRandomValue(-10.0f, 10.0f), 0.0f, GetRandomValue(-10.0f, 10.0f));
				particles[i]->age = GetRandomValue(-30.0f, 0.0f);
			}
		}
	}

	//If we are hovering over a plane, update the force position.
	if(planeID != -1)
	{
		glm::vec3 pointAtPlane = particleCamera.position + (rayWorld * shortestDistance);
		if(curSelectedForce != NULL)
		{
			curSelectedForce->position = pointAtPlane + glm::vec3(0.0f, 5.0f, 0.0f);

			if(curSelectedForce->position.y < -50.0f)
				curSelectedForce->position.y = -50.0f;
		}

		//If active, then add some more particles to the scene!
		if(spawningParticles)
		{
			for(int i = 0; i < 5; i++)
			{
				particles[oldestParticle]->SetPosition(pointAtPlane.x, pointAtPlane.y + 15.0f, pointAtPlane.z);
				particles[oldestParticle]->SetVelocity(GetRandomValue(-20.0f, 20.0f), 0.0f, GetRandomValue(-20.0f, 20.0f));	//Give them a random horizontal velocity so they appear to shoot out of the mouse cursor.
				particles[oldestParticle]->age = 0.0f;

				oldestParticle = (oldestParticle + 1) % particles.size();
			}
		}
	}

	//Mid Point Method - Update once, but using the value as it was in the middle.
	float midPointTimeStep = elapsedTimeStep / 2;

	//To implement the "update twice, using halve values" approach uncomment the lines below.
	
	/*for(int j = 0; j < 2; j++)
	{*/
	for(int i = 0; i < particles.size(); i++)
	{
		Particle* curParticle = particles[i];
		particles[i]->position += particles[i]->velocity * elapsedTimeStep;	//ALT TWO UPDATE APPROACH - particles[i]->position += particles[i]->velocity * midPointTimeStep;

		particles[i]->velocity.y -= 9.81f * midPointTimeStep;

		for(int j = 0; j < planes.size(); j++)
		{
			Plane* curPlane = planes[j];
			glm::vec3 relPos = particles[i]->position - curPlane->position;
			float dotProd = glm::dot(relPos, curPlane->normal);		//(x-p).n

			float dotProd2 = glm::dot(curPlane->normal, particles[i]->velocity);		//(n.v)

			if(dotProd < epsilon && dotProd2 < epsilon)
			{
				glm::vec3 vN = glm::dot(curPlane->normal, particles[i]->velocity) * curPlane->normal;
				glm::vec3 vT = particles[i]->velocity - vN;
				vT *= 0.75f;	//Restitution...? (Check slides for term)

				particles[i]->velocity = vT - elasticity * vN;

				if(currentDemoMode == DemoMode::Snow)
					particles[i]->velocity = glm::vec3(0.0f, 0.0f, 0.0f);
			}
		}		
		particles[i]->velocity -= particles[i]->velocity * dampeningStrength * midPointTimeStep;

		if(currentDemoMode == DemoMode::Celestial || currentDemoMode == DemoMode::Plume)
		{
			for(int j = 0; j < forces.size(); j++)
			{
				if(forces[j]->forceEnabled)
				{
					glm::vec3 distFromForce = forces[j]->position - particles[i]->position;
					float distSquared = glm::dot(distFromForce, distFromForce);

					if(forces[j]->forcePushing)
						distSquared *= -1;

					float relativeForceStrength = 1 / distSquared;
					particles[i]->velocity += (distFromForce * relativeForceStrength * forces[j]->strength * midPointTimeStep); 
				}
			}
		}
	}
	//}

	glutPostRedisplay();
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	int aID = planeShader.GetShaderID();
	glUseProgram(aID);

	glUniformMatrix4fv(glGetUniformLocation(aID, "camMat"), 1, GL_FALSE, glm::value_ptr(particleCamera.GenerateViewMatrix()));
	glUniformMatrix4fv(glGetUniformLocation(aID, "perMat"), 1, GL_FALSE, glm::value_ptr(particleCamera.GenerateProjMatrix()));

	glBegin(GL_QUADS);
	glColor4f(0.4f, 0.4f, 0.4f, 1.0f);

	for(int i = 0; i < planes.size(); i++)
		planes[i]->DrawPlane();

	glEnd();

	int bID = particleShader.GetShaderID();
	glUseProgram(bID);

	glUniformMatrix4fv(glGetUniformLocation(bID, "camMat"), 1, GL_FALSE, glm::value_ptr(particleCamera.GenerateViewMatrix()));
	glUniformMatrix4fv(glGetUniformLocation(bID, "perMat"), 1, GL_FALSE, glm::value_ptr(particleCamera.GenerateProjMatrix()));

	if(currentDemoMode == DemoMode::Celestial || currentDemoMode == DemoMode::Plume)
	{
		glPointSize(10);
		glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		glBegin(GL_POINTS);	

		for(int i = 0; i < forces.size(); i++)
		{
			if(forces[i]->forceEnabled)
				glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
			else
				glColor4f(0.1f, 0.1f, 0.1f, 1.0f);

			glVertex3f(forces[i]->position.x, forces[i]->position.y, forces[i]->position.z);
		}

		glEnd();
	}

	if(currentDemoMode != DemoMode::Celestial)
		glPointSize(10);
	else
		glPointSize(5);

	glBegin(GL_POINTS);

	for(int i = 0; i < particles.size(); i++)
	{
		glm::vec3 pos = particles[i]->GetPosition();
		glVertex3f(pos.x, pos.y, pos.z);

		if(currentDemoMode == DemoMode::Snow)
		{
			glColor4f(1.0f, 1.0f, 1.0f, (particlesShouldAge ? 1.0f - glm::min(1.0f, (particles[i]->age / particleLifetime)) : 1.0f));
		}
		else
		{
			glm::vec3 colorStrength = glm::normalize(pos);
			colorStrength *= particles[i]->velocity;

			colorStrength.x = glm::max(0.0f, glm::min(colorStrength.x, 1.0f));
			colorStrength.y = glm::max(0.0f, glm::min(colorStrength.y, 1.0f));
			colorStrength.z = glm::max(0.0f, glm::min(colorStrength.z, 1.0f));

			glColor4f(colorStrength.x, colorStrength.y, colorStrength.z, (particlesShouldAge ? 1.0f - glm::min(1.0f, (particles[i]->age / particleLifetime)) : 1.0f));
		}
	}

	glEnd();

	glutSwapBuffers();	//Swap buffers, showing our new image.
}



void initialise(void)
{
	prevFrameTimeSinceStart = glutGet(GLUT_ELAPSED_TIME);
	
	particleShader = Shader();
	particleShader.CreateShaderFromFile("vParticle.vShader", "pParticle.pShader");

	planeShader = Shader();
	planeShader.CreateShaderFromFile("vPlane.vShader", "pPlane.pShader");

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
	rightPlane->GenerateVertices(glm::vec3(25.0f, -100.0f, 0.0f), glm::vec3(0.5f, 1.0f, 0.0f));
	planes.push_back(rightPlane);

	Plane* leftPlane = new Plane();
	leftPlane->GenerateVertices(glm::vec3(-25.0f, -100.0f, 0.0f), glm::vec3(-0.5f, 1.0f, 0.0f));
	planes.push_back(leftPlane);	

	Force* firstForce = new Force(glm::vec3(-25.0f, -50.0f, -25.0f), 500.0f);
	forces.push_back(firstForce);

	Force* secondForce = new Force(glm::vec3(25.0f, -50.0f, 25.0f), 500.0f);
	forces.push_back(secondForce);

	Force* thirdForce = new Force(glm::vec3(-25.0f, -50.0f, 25.0f), 500.0f);
	forces.push_back(thirdForce);

	Force* fourthForce = new Force(glm::vec3(25.0f, -50.0f, -25.0f), 500.0f);
	forces.push_back(fourthForce);
}

void HandleMouseClick(int button, int state, int x, int y)
{
	if(button == GLUT_LEFT_BUTTON && curSelectedForce != NULL)
	{
		if(state == GLUT_DOWN)
			curSelectedForce->forcePushing = true;
		else
			curSelectedForce->forcePushing = false;
	}

	if(button == GLUT_RIGHT_BUTTON)
	{
		if(state == GLUT_DOWN)
			spawningParticles = true;
		else
			spawningParticles = false;
	}
}

void HandleRegularRelease(unsigned char key, int x, int y)
{
	switch(key)
	{
	}
}

void SwitchDemoMode(DemoMode demoMode)
{
	if(demoMode == DemoMode::Plume)
		elasticity = 0.1f;
	else
		elasticity = 0.8f;

	if(demoMode == DemoMode::Snow)
	{
		particleCamera.position = glm::vec3(0.0f, -60.0f, 100.0f);
		particleCamera.forwardVec = glm::normalize(particleCamera.target - particleCamera.position);
		particlesShouldAge = true;

		for(int i = 0; i < maxParticles; i++)
		{
			particles[i]->SetPosition(GetRandomValue(-100.0f, 100.0f), GetRandomValue(0.0f, 100.0f), GetRandomValue(-100.0f, 100.0f));
			particles[i]->SetVelocity(GetRandomValue(-10.0f, 10.0f), 0.0f, 0.0f);
			particles[i]->age = GetRandomValue(-30.0f, 0.0f);
		}
	}
	else
	{
		if(demoMode == DemoMode::Plume)
		{
			particleCamera.position = glm::vec3(0.0f, -50.0f, 100.0f);
			particleCamera.forwardVec = glm::normalize(particleCamera.target - particleCamera.position);
			particlesShouldAge = false;

			for(int i = 0; i < forces.size(); i++)
				forces[i]->forceEnabled = false;

			curSelectedForce = forces[0];
		}
		else
		{
			particleCamera.position = glm::vec3(0.0f, 0.0f, 100.0f);
			particleCamera.forwardVec = glm::normalize(particleCamera.target - particleCamera.position);
			particlesShouldAge = true;

			for(int i = 0; i < forces.size(); i++)
				forces[i]->forceEnabled = true;

			curSelectedForce = NULL;
		}

		for(int i = 0; i < maxParticles; i++)
		{
			particles[i]->SetPosition(GetRandomValue(-50.0f, 50.0f), GetRandomValue(0.0f, 100.0f), GetRandomValue(-50.0f, 50.0f));
			particles[i]->SetVelocity(0.0f, 0.0f, 0.0f);
			particles[i]->age = GetRandomValue(-30.0f, 0.0f);
		}
	}
}

void HandleRegularInput(unsigned char key, int x, int y)
{
	switch(key)
	{
		case ' ':
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
	}
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
	glutKeyboardUpFunc(HandleRegularRelease);
	glutMotionFunc(HandleMouseMovement);
	glutPassiveMotionFunc(HandleMouseMovement);
	glutMouseFunc(HandleMouseClick);
	//glutSpecialFunc(HandleSpecialInput);

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

