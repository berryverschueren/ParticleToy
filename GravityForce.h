#pragma once

#include "Particle.h"
#include "Force.h"

class GravityForce: public Force {
 public:
  GravityForce(std::vector<Particle*> p, Vec2f v);

  void apply() override;
  void draw() override;
 private:

  Vec2f const m_v;
};