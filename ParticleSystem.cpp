#include "ParticleSystem.h"

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

void ParticleSystem::applyConstraints() {
    int ii, size = cVector.size();
    for (ii=0; ii<size; ii++) {
        // TODO::
        //cVector[ii]->apply();
    }
}

std::vector<float> ParticleSystem::simpleCollision(std::vector<float> state){

    float floor = -1.0; float ceiling = 1.0;
    float wallL = -1.0; float wallR = 1.0;
    float k_r = 0.8;

    //floor
    for (int i = 0; i < pVector.size() ; ++i) {
        if(state[i*4+1] < floor){
            float vN = state[i*4+3];
            if(vN >= 0) {
                state[i * 4 + 3] = -k_r * vN;
            } else{
                state[i*4+3] = 0.001;
            }
        }
    }

    //ceiling
    for (int i = 0; i < pVector.size() ; ++i) {
        if(state[i*4+1] > ceiling){
            float vN = state[i*4+3];
            if(vN <= 0) {
                state[i * 4 + 3] = -k_r * vN;
            } else{
                state[i*4+3] = -0.001;
            }

        }
    }

    //walls
    for (int i = 0; i < pVector.size() ; ++i) {
        float vN = state[i*4+2];
        float reversedForce = -k_r*vN;
        if(state[i*4] < wallL){
            if(reversedForce >= 0){
                state[i*4+2] = reversedForce;
            } else{
                state[i*4+2] = 0.001;
            }
        }
        if(state[i*4] > wallR){
            if(reversedForce <= 0) {
                state[i*4+2] = reversedForce;
            } else{
                state[i*4+2] = -0.001;
            }
        }
    }

    return state;
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

void ParticleSystem::removeLastParticle() {
    pVector.pop_back();
}

void ParticleSystem::addForce(Force *f) {
    fVector.push_back(f);
}

void ParticleSystem::removeLastForce() {
    fVector.pop_back();
}

void ParticleSystem::addConstraint(Constraint *c) {
    cVector.push_back(c);
}

void ParticleSystem::deleteAll() {
    pVector.clear();
    fVector.clear();
    cVector.clear();
}
