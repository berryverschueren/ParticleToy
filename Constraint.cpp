#include "Constraint.h"
#include <vector>

Constraint::Constraint(Particle *p1, Particle *p2, double dist) {};

float constraint_value();
float constraint_derivative_value();
std::vector<Vec2f> jacobian_value();
std::vector<Vec2f> jacobian_derivative_value();

void Constraint::draw() {};
