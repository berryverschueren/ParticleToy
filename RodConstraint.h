#pragma once

#include "Particle.h"
#include <vector>

class RodConstraint {
 public:
  RodConstraint(Particle *p1, Particle * p2, double dist);

  float constraint_value();
  float constraint_derivative_value();
  std::vector<Vec2f> jacobian_value();
  std::vector<Vec2f> jacobian_derivative_value();
  void draw();

 private:

  Particle * const m_p1;
  Particle * const m_p2;
  double const m_dist;
};
