#include "Constraint.h"
#include <vector>

float Constraint::constraint_value();
float Constraint::constraint_derivative_value();
std::vector<Vec2f> Constraint::jacobian_value();
std::vector<Vec2f> Constraint::jacobian_derivative_value();
void Constraint::draw() {};
