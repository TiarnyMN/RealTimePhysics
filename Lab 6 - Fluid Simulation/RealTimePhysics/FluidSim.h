#pragma once
#include "shared.h";
#include <map>
#include <unordered_map>
#include <gtx\fast_exponential.hpp>
#include "Hasher.h"
#include "glm/gtx/euler_angles.hpp"

class FluidSim
{
public:
	virtual void Update(void);
	virtual void Render(GLuint capsuleShaderID, GLuint particleShaderID, glm::mat4 &viewMat, glm::mat4 &projMat);

	vector<Particle*> GetParticles(void);

	bool HandleRegularInput(unsigned char curKey);
	void PrintDebug(std::stringstream& ss);
	Capsule container;

private:
	float CalculateKernal(Particle* curParticle, Particle* compareParticle);
	glm::vec3 CalculateNormalKernal(glm::vec3 &r, float &rLength);
	float CalculateSurfaceKernal(glm::vec3 &r, float &rLength);
	glm::vec3 CalculatePressureKernal(glm::vec3 &r, float &rLength);
	float CalculateViscosityKernal(glm::vec3 &r, float &rLength);

	bool CalculateInternalCapsuleCollision(Capsule* container, Particle* particle);
	bool CalculateExternalCapsuleCollision(Capsule* container, Particle* particle);
	float CalculateCapsulePoint(Capsule* container, Particle* particle);
	float CalculateCapsuleFunc(Capsule* container, Particle* particle);

	//std::multimap<glm::vec3, Particle*, Hasher> hashMap;
	std::unordered_multimap<glm::vec3, Particle*, Hasher, Hasher> hashMap;

	vector<Particle*> particles;

protected:
	void InitialiseParticles(void);

	float RestDensity;
	float ParticlesPerDimension;
	float ParticleCount;
	float SurfaceThreshold;

	float StiffnessConstant;
	float ViscosityConstant;
	float AverageParticlesPerKernel;
	float FluidStep;

	float VerticalBounds;
	float HorizontalBounds;

	float BoxVolume;

	float DampeningStrength, FrictionStrength;

	float ParticleMass;

	float SmoothingKernelSize;
	float ParticleSize;

	float Sigma;

	glm::vec3 Color;

	bool wallEnabled;

	float KernelPow9, KernelPow2, KernelPow6;
	float DefaultKernelConstant, NormalKernelConstant, SurfaceKernelConstant, PressureKernelConstant, ViscosityKernelConstant;

	glm::vec3 containerVelocity;
	bool crashCapsule, brokenCapsule, showSurfaceParticles, hideSurfaceParticles, showNeighbours;
	float innerRadius, innerLength;
	int particleHighlightID;
	
	Capsule innerContainer;
};

