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
	auto derivative1 = particleSystem->derivEval();
	derivative1 = vectorTimesScalar(derivative1, dt);
	auto state1 = derivative1;
	auto newState1 = vectorAddition(state, vectorTimesScalar(state1, 0.5f));
	particleSystem->setState(newState1);
	auto derivative2 = particleSystem->derivEval();
	derivative2 = vectorTimesScalar(derivative2, dt);
	auto state2 = derivative2;
	auto newState2 = vectorAddition(state, vectorTimesScalar(state2, 0.5f));
	particleSystem->setState(newState2);
	auto derivative3 = particleSystem->derivEval();
	derivative3 = vectorTimesScalar(derivative3, dt);
	auto state3 = derivative3;
	auto newState3 = vectorAddition(state, vectorTimesScalar(state3, 0.5f));
	particleSystem->setState(newState3);
	auto derivative4 = particleSystem->derivEval();
	derivative4 = vectorTimesScalar(derivative4, dt);
	auto state4 = derivative4;
	auto newState4 = vectorAddition(vectorAddition(vectorAddition(state, 
		vectorTimesScalar(newState1, 1.0f/6.0f)), vectorTimesScalar(newState2, 1.0f/3.0f)), 
		vectorAddition(vectorTimesScalar(newState3, 1.0f/3.0f), vectorTimesScalar(derivative4, 1.0f/6.0f)));
	particleSystem->setState(newState4);
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