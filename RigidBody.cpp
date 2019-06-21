#include "RigidBody.h"
#include <GL/glut.h>

RigidBody::RigidBody(const Vector2f & ConstructPos) :
	_original(ConstructPos), _center(ConstructPos)
	, _velocity(Vector2f(0.0f, 0.0f)), _force(Vector2f(0.0f, 0.0f)), 
    _mass(1.0f), _width(10.0f/64), _height(10.0f/64), _orientation(MatrixXf::Zero(3,3)),
	_angularVelocity(VectorXf::Zero(3)) { 

		_orientation << 0.f,-1.f,0.f,
						1.f,0.f,0.f,
						0.f,0.f,1.f;

	}

RigidBody::~RigidBody(void) { }

void RigidBody::ResetForce() {
	_force = Vector2f(0.0f, 0.0f);
}

void RigidBody::Reset() {
	_center = _original;
	_velocity = Vector2f(0.0f, 0.0f);
	_force = Vector2f(0.0f, 0.0f);
	_mass = 1.0f;
	_width = 10.0f/64;
	_height = 10.0f/64;

	_orientation = MatrixXf::Zero(3,3);
		_orientation << 0.f,-1.f,0.f,
						1.f,0.f,0.f,
						0.f,0.f,1.f;
	_angularVelocity = VectorXf::Zero(3);
}

void RigidBody::Draw() {
    auto h = 1.0f / (64);
	glColor3f(1.f, 1.f, 1.f); 
	glBegin(GL_QUADS);
	glVertex2f((_center[0]-_width/2.0)*h, (_center[1]-_height/2.0)*h);
	glVertex2f((_center[0]+_width/2.0)*h, (_center[1]-_height/2.0)*h);
	glVertex2f((_center[0]+_width/2.0)*h, (_center[1]+_height/2.0)*h);
	glVertex2f((_center[0]-_width/2.0)*h, (_center[1]+_height/2.0)*h);

	//printf ("Current center: %g, %g\n", (_center[0]-_width/2.0)*h, (_center[1]-_height/2.0)*h);
	glEnd();
}
