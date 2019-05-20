#pragma once

#include "Particle.h"
#include "Force.h"
#include "Constraint.h"
#include "linearSolver.h"

class ParticleSystem {
    public: 
    ParticleSystem();

    std::vector<Particle*> pVector;
    std::vector<Force*> fVector;
    std::vector<Constraint*> cVector;

    int getDim();
    std::vector<float> getState();
    std::vector<float> derivEval();
    void setState(std::vector<float> state);
    void resetForces();
    void applyForces();
    void applyConstraints();
    std::vector<Particle*> getParticles();
    std::vector<Force*> getForces();
    std::vector<Constraint*> getConstraints();
    void addParticle(Particle *p);
    void addForce(Force *f);
    void addConstraint(Constraint *c);
};