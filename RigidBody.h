#pragma once
#include <gfx/vec2.h>
#include "include/Eigen/Dense"
#include "include/Eigen/IterativeLinearSolvers"
using namespace Eigen;

class RigidBody
{
    public:
        RigidBody(const Vector2f & ConstructPos);
        virtual ~RigidBody(void);

        void ResetForce();
        void Reset();
        void Draw();

        Vector2f _original;
        Vector2f _center;
        Vector2f _force;
        Vector2f _velocity;
        float _mass;
        float _width;
        float _height;

        MatrixXf _orientation;      // rotation --> derivative == _orientation * _angularVelocity
        VectorXf _angularVelocity;  // omega    --> derivative == relative position * force       == torque
};
