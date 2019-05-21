#pragma once

#include "Particle.h"
#include "Constraint.h"
#include <vector>

class CircularWireConstraint: public Constraint {
 public:
  CircularWireConstraint(Particle *p, const Vec2f & center, const double radius);
  float constraint_value() override ;
  float constraint_derivative_value() override ;
  std::vector<Vec2f> jacobian_value() override ;
  std::vector<Vec2f> jacobian_derivative_value() override ;
  void draw() override ;

  Particle * const m_p;
  Vec2f const m_center;
  double const m_radius;
};
