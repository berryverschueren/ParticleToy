#include "Particle.h"
#include "Force.h"
#include <vector>
#include "ParticleSystem.h"

#define DAMP 0.98f

std::vector<float> vectorTimesScalar(std::vector<float> vec, float scalar) {
	int ii, size = vec.size();
	for (ii=0; ii<size; ii++) {
		vec[ii] *= scalar;
	}
	return vec;
}

std::vector<float> vectorAddition(std::vector<float> vec1, std::vector<float> vec2) {
	int ii, size = vec1.size();
	for (ii=0; ii<size; ii++) {
		vec1[ii] = vec1[ii] + vec2[ii];
	}
	return vec1;
}

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
	auto state = particleSystem->getState();
	auto derivative = particleSystem->derivEval();
	derivative = vectorTimesScalar(derivative, dt);
	auto newState = vectorAddition(state, derivative);

	newState = particleSystem->simpleCollision(newState);

	particleSystem->setState(newState);
}

void midpointSolver(ParticleSystem *particleSystem, float dt) {
	auto state = particleSystem->getState();
	auto derivative1 = particleSystem->derivEval();
	derivative1 = vectorTimesScalar(derivative1, (dt * 0.5f));
	auto midpoint = vectorAddition(state, derivative1);
	particleSystem->setState(midpoint);
	auto derivative2 = particleSystem->derivEval();
	derivative2 = vectorTimesScalar(derivative2, dt);
	auto newState = vectorAddition(state, derivative2);
	particleSystem->setState(newState);
}

// TODO::
// not working properly
void rungeKuttaSolver(ParticleSystem *particleSystem, float dt) {
	auto state = particleSystem->getState();
	auto k1 = vectorTimesScalar(particleSystem->derivEval(), dt);
	auto y1 = vectorAddition(state, vectorTimesScalar(k1, 0.5f));
	particleSystem->setState(y1);
	auto k2 = vectorTimesScalar(particleSystem->derivEval(), dt);
	auto y2 = vectorAddition(state, vectorTimesScalar(k2, 0.5f));
	particleSystem->setState(y2);
	auto k3 = vectorTimesScalar(particleSystem->derivEval(), dt);
	auto y3 = vectorAddition(state, k3);
	particleSystem->setState(y3);
	auto k4 = vectorTimesScalar(particleSystem->derivEval(), dt);
	auto wK1 = vectorTimesScalar(k2, 1.0f/6.0f);
	auto wK2 = vectorTimesScalar(k2, 1.0f/3.0f);
	auto wK3 = vectorTimesScalar(k3, 1.0f/3.0f);
	auto wK4 = vectorTimesScalar(k2, 1.0f/6.0f);
	auto combined = vectorAddition(vectorAddition(wK1, wK2), vectorAddition(wK3, wK4));
	auto newState = vectorAddition(state, combined);
	particleSystem->setState(newState);
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