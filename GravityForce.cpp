#include "GravityForce.h"

GravityForce::GravityForce(std::vector<Particle*> p) {
   this->particles = p;
};

void GravityForce::apply() {
    for (Particle* p : particles) {
        p->m_Force += Vec2f(0.0, p->m_Mass * -0.001f);
    }
};

void GravityForce::draw() {};
