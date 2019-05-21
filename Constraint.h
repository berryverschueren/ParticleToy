#pragma once

#include "Particle.h"
#include <vector>

class Constraint {
public:
    virtual float constraint_value() = 0;
    virtual float constraint_derivative_value() = 0;
    virtual std::vector<Vec2f> jacobian_value() = 0;
    virtual std::vector<Vec2f> jacobian_derivative_value() = 0;
    virtual void draw() = 0;
};