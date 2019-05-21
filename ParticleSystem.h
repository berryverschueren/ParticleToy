#pragma once

#include "Particle.h"
#include "Force.h"
#include "Constraint.h"
#include "linearSolver.h"


#include "include/Eigen/Dense"
#include "include/Eigen/IterativeLinearSolvers"
using namespace Eigen;

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
    std::vector<float> simpleCollision(std::vector<float> state);
    void resetForces();
    void applyForces();
    void applyConstraints();
    std::vector<Particle*> getParticles();
    std::vector<Force*> getForces();
    std::vector<Constraint*> getConstraints();
    void addParticle(Particle *p);
    void removeLastParticle();
    void addForce(Force *f);
    void removeLastForce();
    void addConstraint(Constraint *c);
    void deleteAll();
    int getPosition(Particle *p);
};