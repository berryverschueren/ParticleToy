#include "GravityForce.h"
#include <vector>

GravityForce::GravityForce(std::vector<Particle*> p, Vec2f v): m_v(v) {
   this->particles = p;
};

void GravityForce::apply() {
    for (Particle* p : particles) {
       p->m_Force += Vec2f(p->m_Mass * m_v[0] ,p->m_Mass * m_v[1]);
    }
};

void GravityForce::draw() {};
