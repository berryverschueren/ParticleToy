#include "Particle.h"
#include "Force.h"

#include <vector>

#define DAMP 0.98f
#define RAND (((rand()%2000)/1000.f)-1.f)
void simulation_step( std::vector<Particle*> pVector, std::vector<Force*> fVector, float dt )
{
	int ii, pSize = pVector.size(), fSize = fVector.size();
	
	// clear force accumulators
	for (ii=0; ii<pSize; ii++) {
		pVector[ii]->reset_force();
	}

	// sum forces for all particles
	for (ii=0; ii<fSize; ii++) {
		fVector[ii]->apply();
	}

	// gather step
	for(ii=0; ii<pSize; ii++)
	{
		pVector[ii]->m_Position += dt*pVector[ii]->m_Velocity;
		pVector[ii]->m_Velocity += dt*pVector[ii]->m_Force/pVector[ii]->m_Mass;
	}

}

