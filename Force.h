#pragma once

#include "Particle.h"
#include <vector>

class Force {
 public:
  std::vector<Particle*> particles;
  virtual void apply() = 0;
  virtual void draw() = 0;
};
