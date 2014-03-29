#pragma once
#include "shared.h";
#include <map>
#include <unordered_map>
#include "Hasher.h"



class FluidSim
{
public:
	virtual void Update(void);
	virtual void Render(GLuint shaderID);

	vector<Particle*> GetParticles(void);

private:
	float CalculateKernal(Particle* curParticle, Particle* compareParticle);
	glm::vec3 CalculateNormalKernal(glm::vec3 &r, float &rLength);
	float CalculateSurfaceKernal(glm::vec3 &r, float &rLength);
	glm::vec3 CalculatePressureKernal(glm::vec3 &r, float &rLength);
	float CalculateViscosityKernal(glm::vec3 &r, float &rLength);

	bool CalculateCapsuleCollision(Capsule* container, Particle* particle);
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

	float DampeningStrength;

	float ParticleMass;

	float SmoothingKernelSize;
	float ParticleSize;

	float Sigma;

	bool wallEnabled;

	float KernelPow9, KernelPow2, KernelPow6;
	float DefaultKernelConstant, NormalKernelConstant, SurfaceKernelConstant, PressureKernelConstant, ViscosityKernelConstant;

	Capsule container;
};

