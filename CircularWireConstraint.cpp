#include "CircularWireConstraint.h"
#include <GL/glut.h>
#include <vector>

#define PI 3.1415926535897932384626433832795
static void draw_circle(const Vec2f & vect, float radius)
{
    glBegin(GL_LINE_LOOP);
    glColor3f(0.0,1.0,0.0);
    for (int i=0; i<360; i=i+18)
    {
        float degInRad = i*PI/180;
        glVertex2f(vect[0]+cos(degInRad)*radius,vect[1]+sin(degInRad)*radius);
    }
    glEnd();
}

CircularWireConstraint::CircularWireConstraint(Particle *p, const Vec2f & center, const double radius) :
    m_p(p), m_center(center), m_radius(radius) {}

float CircularWireConstraint::constraint_value() {
    float xDiff = pow(m_p->m_Position[0] - m_center[0], 2);
    float yDiff = pow(m_p->m_Position[1] - m_center[1], 2);
    return (float) (xDiff + yDiff - pow(m_radius, 2));
}

float CircularWireConstraint::constraint_derivative_value() {
    float xDiff = 2 * (m_p->m_Position[0] - m_center[0]);
    float yDiff = 2 * (m_p->m_Position[1] - m_center[1]);
    float xVelocityDiff = 2 * (m_p->m_Velocity[0] - 0);
    float yVelocityDiff = 2 * (m_p->m_Velocity[1] - 0);
    return (float) (xDiff * xVelocityDiff * yDiff * yVelocityDiff);
}

std::vector<Vec2f> CircularWireConstraint::jacobian_value() {
    float xDiffNormal = 2 * (m_p->m_Position[0] - m_center[0]);
    float yDiffNormal = 2 * (m_p->m_Position[1] - m_center[1]);
    float xDiffReverse = 2 * (0 - m_p->m_Position[0]);
    float yDiffReverse = 2 * (0 - m_p->m_Position[1]);
    return std::vector<Vec2f> { Vec2f(xDiffNormal, yDiffNormal), Vec2f(xDiffReverse, yDiffReverse) };
}

std::vector<Vec2f> CircularWireConstraint::jacobian_derivative_value() {
    float xVelocityDiffNormal = 2 * (m_p->m_Velocity[0] - 0);
    float yVelocityDiffNormal = 2 * (m_p->m_Velocity[1] - 0);
    float xVelocityDiffReverse = 2 * (0 - m_p->m_Velocity[0]);
    float yVelocityDiffReverse = 2 * (0 - m_p->m_Velocity[1]);
    return std::vector<Vec2f> { Vec2f(xVelocityDiffNormal, yVelocityDiffNormal), Vec2f(xVelocityDiffReverse, yVelocityDiffReverse) };
}






void CircularWireConstraint::draw()
{
	draw_circle(m_center, m_radius);
}
