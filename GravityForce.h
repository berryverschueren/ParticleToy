#pragma once

#include "Particle.h"
#include "Force.h"

class GravityForce: public Force {
 public:
  GravityForce(std::vector<Particle*> particles);

  void apply() override;
  void draw() override;
};