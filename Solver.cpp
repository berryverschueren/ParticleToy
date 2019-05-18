#include "Particle.h"
#include "Force.h"
#include <vector>
#include "ParticleSystem.h"

#define DAMP 0.98f

void skeletonSolver(ParticleSystem *particleSystem, float dt) {
	int ii, pSize = particleSystem->getParticles().size(), fSize = particleSystem->getForces().size();

	// clear force accumulators
	for (ii=0; ii<pSize; ii++) {
		particleSystem->getParticles()[ii]->reset_force();
	}

	// sum forces for all particles
	for (ii=0; ii<fSize; ii++) {
		particleSystem->getForces()[ii]->apply();
	}

	// gather step
	for(ii=0; ii<pSize; ii++)
	{
		particleSystem->getParticles()[ii]->m_Position += dt*particleSystem->getParticles()[ii]->m_Velocity;
		particleSystem->getParticles()[ii]->m_Velocity += dt*particleSystem->getParticles()[ii]->m_Force/particleSystem->getParticles()[ii]->m_Mass;
	}
}

void explicitEulerSolver(ParticleSystem *particleSystem, float dt) {
}

void midpointSolver(ParticleSystem *particleSystem, float dt) {
}

void rungeKuttaSolver(ParticleSystem *particleSystem, float dt) {
}

void simulation_step(ParticleSystem *particleSystem, float dt, int solverVersion)
{	
	switch (solverVersion) {
		case 1:
			explicitEulerSolver(particleSystem, dt);
		break;
		case 2:
			midpointSolver(particleSystem, dt);
		break;
		case 3: 
			rungeKuttaSolver(particleSystem, dt);
		break;
		case 0:
		default:
			skeletonSolver(particleSystem, dt);
		break;
	}
}