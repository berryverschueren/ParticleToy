#pragma once

#include "Particle.h"
#include <vector>

class Constraint {
public:
    Constraint();
    float constraint_value();
    float constraint_derivative_value();
    std::vector<Vec2f> jacobian_value();
    std::vector<Vec2f> jacobian_derivative_value();
    void draw();
};