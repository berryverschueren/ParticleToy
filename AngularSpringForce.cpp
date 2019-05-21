#include "AngularSpringForce.h"
#include <GL/glut.h>

AngularSpringForce::AngularSpringForce(Particle *p1, Particle *p2, Particle *p3, double dist, double ks, double kd):
        m_p1(p1), m_p2(p2), m_p3(p3), m_dist(dist), m_ks(ks), m_kd(kd) {}

void AngularSpringForce::draw()
{
    glBegin( GL_LINES );
    glColor3f(0.6, 0.7, 0.9);
    glVertex2f( m_p1->m_Position[0], m_p1->m_Position[1] );
    glColor3f(0.6, 0.7, 0.9);
    glVertex2f( m_p2->m_Position[0], m_p2->m_Position[1] );
    glColor3f(0.6, 0.7, 0.9);
    glVertex2f( m_p3->m_Position[0], m_p3->m_Position[1] );
    glEnd();
}

void AngularSpringForce::apply(){
    ////not sure
    Vec2f v1 = m_p1->m_Position - m_p2->m_Position;
    float lengthV1 = std::sqrt(std::abs(std::pow(v1[0],2))+std::abs(std::pow(v1[1],2)));
    Vec2f v2 = m_p2->m_Position - m_p3->m_Position;
    float lengthV2 = std::sqrt(std::abs(std::pow(v2[0],2))+std::abs(std::pow(v2[1],2)));

    Vec2f normV1 = v1/lengthV1;
    Vec2f normV2 = v2/lengthV2;

    float dotProduct = normV1[0]*normV2[0]+normV1[1]*normV2[1];

    float angle = std::acos(dotProduct);

    //m_p1->m_Force +=


}
