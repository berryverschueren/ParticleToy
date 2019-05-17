#include "GravityForce.h"

GravityForce::GravityForce(std::vector<Particle*> particles) : {
   this->particles = particles;
};

void GravityForce::apply() {
    for (Particle* p : particles) {
        p->m_Force += p->mass * 10;
    }
};

void GravityForce::draw() {};
