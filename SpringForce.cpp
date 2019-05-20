#include "SpringForce.h"
#include <GL/glut.h>

SpringForce::SpringForce(Particle* p1, Particle* p2, double dist, double ks, double kd) :
  m_p1(p1), m_p2(p2), m_dist(dist), m_ks(ks), m_kd(kd) {}

void SpringForce::draw()
{
  glBegin( GL_LINES );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( m_p1->m_Position[0], m_p1->m_Position[1] );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( m_p2->m_Position[0], m_p2->m_Position[1] );
  glEnd();
}

void SpringForce::apply(){
    Vec2f l = m_p1->m_Position - m_p2->m_Position;
    Vec2f i = m_p1->m_Velocity - m_p2->m_Velocity;

    float normL = std::sqrt(std::abs(std::pow(l[0],2))+std::abs(std::pow(l[1],2)));
    float dot = i[0]*l[0]+i[1]*l[1];

    Vec2f f1 = -(m_ks*(normL-m_dist)+m_kd*(dot/normL))*(l/normL);
    Vec2f f2 = -f1;

    m_p1->m_Force += f1;
    m_p2->m_Force += f2;
}
