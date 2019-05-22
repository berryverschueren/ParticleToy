#include "AngularSpringForce.h"
#include <GL/glut.h>

AngularSpringForce::AngularSpringForce(Particle *p1, Particle *p2, Particle *p3, double alpha, double ks, double kd):
        m_p1(p1), m_p2(p2), m_p3(p3), m_alpha(alpha), m_ks(ks), m_kd(kd) {}

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
    Vec2f v1 = m_p1->m_Position - m_p2->m_Position;
    float lengthV1 = std::sqrt(std::abs(std::pow(v1[0],2))+std::abs(std::pow(v1[1],2)));
    Vec2f v2 = m_p2->m_Position - m_p3->m_Position;
    float lengthV2 = std::sqrt(std::abs(std::pow(v2[0],2))+std::abs(std::pow(v2[1],2)));

    Vec2f normV1 = v1/lengthV1;
    Vec2f normV2 = v2/lengthV2;

    float dotProduct = normV1[0]*normV2[0]+normV1[1]*normV2[1];

    float C = std::acos(dotProduct)-m_alpha;

    /*Vec2f l = m_p1->m_Position - m_p2->m_Position;
    Vec2f vi1 = m_p1->m_Velocity - m_p2->m_Velocity;
    Vec2f vi1 = m_p1->m_Velocity - m_p2->m_Velocity;

    float lengthL = std::sqrt(std::abs(std::pow(l[0],2))+std::abs(std::pow(l[1],2)));
    float dot = i[0]*l[0]+i[1]*l[1];

    Vec2f f1 = -(m_ks*(lengthL-m_dist)+m_kd*(dot/lengthL))*(l/lengthL);
    Vec2f f2 = -f1;

    m_p1->m_Force += f1;
    m_p2->m_Force += f2;
*/
    //m_p1->m_Force +=


}
