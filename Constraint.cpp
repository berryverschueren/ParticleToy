#include "Constraint.h"
#include <vector>

float constraint_value();
float constraint_derivative_value();
std::vector<Vec2f> jacobian_value();
std::vector<Vec2f> jacobian_derivative_value();

void Constraint::draw() {};
