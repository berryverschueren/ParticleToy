#include "RodConstraint.h"
#include <GL/glut.h>
#include <cmath>
#include <vector>

RodConstraint::RodConstraint(Particle *p1, Particle *p2, double dist) {
  this->m_p1 = p1;
  this->m_p2 = p2;
  this->m_dist = dist;
}

float RodConstraint::constraint_value() {
  float xDiff = pow(m_p1->m_Position[0] - m_p2->m_Position[0], 2);
  float yDiff = pow(m_p1->m_Position[1] - m_p2->m_Position[1], 2);
  return (float) (xDiff + yDiff - pow(m_dist, 2));
}

float RodConstraint::constraint_derivative_value() {
  float xDiff = 2 * (m_p1->m_Position[0] - m_p2->m_Position[0]);
  float yDiff = 2 * (m_p1->m_Position[1] - m_p2->m_Position[1]);
  float xVelocityDiff = 2 * (m_p1->m_Velocity[0] - m_p2->m_Velocity[0]);
  float yVelocityDiff = 2 * (m_p1->m_Velocity[1] - m_p2->m_Velocity[1]);
  return (float) (xDiff * xVelocityDiff * yDiff * yVelocityDiff);
}

std::vector<Vec2f> RodConstraint::jacobian_value() {
  float xDiffNormal = 2 * (m_p1->m_Position[0] - m_p2->m_Position[0]);
  float yDiffNormal = 2 * (m_p1->m_Position[1] - m_p2->m_Position[1]);
  float xDiffReverse = 2 * (m_p2->m_Position[0] - m_p1->m_Position[0]);
  float yDiffReverse = 2 * (m_p2->m_Position[1] - m_p1->m_Position[1]);
  return std::vector<Vec2f> { Vec2f(xDiffNormal, yDiffNormal), Vec2f(xDiffReverse, yDiffReverse) };
}

std::vector<Vec2f> RodConstraint::jacobian_derivative_value() {
  float xVelocityDiffNormal = 2 * (m_p1->m_Velocity[0] - m_p2->m_Velocity[0]);
  float yVelocityDiffNormal = 2 * (m_p1->m_Velocity[1] - m_p2->m_Velocity[1]);
  float xVelocityDiffReverse = 2 * (m_p2->m_Velocity[0] - m_p1->m_Velocity[0]);
  float yVelocityDiffReverse = 2 * (m_p2->m_Velocity[1] - m_p1->m_Velocity[1]);
  return std::vector<Vec2f> { Vec2f(xVelocityDiffNormal, yVelocityDiffNormal), Vec2f(xVelocityDiffReverse, yVelocityDiffReverse) };
}
  
void RodConstraint::draw()
{
  glBegin( GL_LINES );
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f( m_p1->m_Position[0], m_p1->m_Position[1] );
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f( m_p2->m_Position[0], m_p2->m_Position[1] );
  glEnd();
}
