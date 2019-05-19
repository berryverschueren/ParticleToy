#include "ParticleSystem.h"
#include <cmath>
#include <vector>

ParticleSystem::ParticleSystem() {};

int ParticleSystem::getDim() {
    return pVector.size() * 4;
}

std::vector<float> ParticleSystem::getState() {
    std::vector<float> state;
    state.resize(this->getDim());
    int ii, size = pVector.size();
    for (ii = 0; ii<size; ii++) {
        Particle *p = pVector[ii];
        state[ii*4] = p->m_Position[0];
        state[ii*4+1] = p->m_Position[1];
        state[ii*4+2] = p->m_Velocity[0];
        state[ii*4+3] = p->m_Velocity[1];
    }
    return state;
}

std::vector<float> ParticleSystem::derivEval() {
    this->resetForces();
    this->applyForces();
    this->applyConstraints();

    std::vector<float> state;
    state.resize(this->getDim());
    int ii, size = pVector.size();
    for (ii = 0; ii<size; ii++) {
        Particle *p = pVector[ii];
        state[ii*4] = p->m_Velocity[0];
        state[ii*4+1] = p->m_Velocity[1];
        state[ii*4+2] = p->m_Force[0] / p->m_Mass;
        state[ii*4+3] = p->m_Force[1] / p->m_Mass;
    }
    return state;
}

void ParticleSystem::setState(std::vector<float> state) {
    int ii, size = pVector.size();
    for (ii=0; ii<size; ii++) {
        pVector[ii]->m_Position[0] = state[ii*4];
        pVector[ii]->m_Position[1] = state[ii*4+1];
        pVector[ii]->m_Velocity[0] = state[ii*4+2];
        pVector[ii]->m_Velocity[1] = state[ii*4+3];
    }
}

void ParticleSystem::resetForces() {
    int ii, size = pVector.size();
    for (ii=0; ii<size; ii++) {
        pVector[ii]->reset_force();
    }
}

void ParticleSystem::applyForces() {
    int ii, size = fVector.size();
    for (ii=0; ii<size; ii++) {
        fVector[ii]->apply();
    }
}
// NOTE need particles, constraints and force
void ParticleSystem::applyConstraints() {
    // we got : C, Cder, J, Jder
    // need to make : M (diagonal mass matrix 3n*3n), W (M^-1), Q (3n long force field), q (3n state vector)
    // calculate lambda^T = (JWJ^T)^-1 (-Jder qder^T - J W Q^T)
    // next calculate ^Q = lambda J
}

std::vector<Particle*> ParticleSystem::getParticles() {
    return pVector;
}

std::vector<Force*> ParticleSystem::getForces() {
    return fVector;
}

std::vector<Constraint*> ParticleSystem::getConstraints() {
    return cVector;
}

void ParticleSystem::addParticle(Particle *p) {
    pVector.push_back(p);
}

void ParticleSystem::addForce(Force *f) {
    fVector.push_back(f);
}

void ParticleSystem::addConstraint(Constraint *c) {
    cVector.push_back(c);
}