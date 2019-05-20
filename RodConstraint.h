#pragma once

#include "Particle.h"
#include "Constraint.h"
#include <vector>

class RodConstraint: public Constraint {
 public:
  RodConstraint(Particle *p1, Particle * p2, double dist);

  float constraint_value() override ;
  float constraint_derivative_value() override ;
  std::vector<Vec2f> jacobian_value() override ;
  std::vector<Vec2f> jacobian_derivative_value() override ;
  void draw() override ;

  Particle * m_p1;
  Particle * m_p2;
  double m_dist;
};
