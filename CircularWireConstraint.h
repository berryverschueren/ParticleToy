#pragma once

#include "Particle.h"
#include <vector>

class CircularWireConstraint {
 public:
  CircularWireConstraint(Particle *p, const Vec2f & center, const double radius);
  float constraint_value();
  float constraint_derivative_value();
  std::vector<Vec2f> jacobian_value();
  std::vector<Vec2f> jacobian_derivative_value();

  void draw();

 private:

  Particle * const m_p;
  Vec2f const m_center;
  double const m_radius;
};
