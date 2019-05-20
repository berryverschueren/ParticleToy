#include "ParticleSystem.h"
#include <cmath>
#include <vector>
#include "Matrices.h"
#include "linearSolver.h"

ParticleSystem::ParticleSystem() {};

int ParticleSystem::getDim() {
    return pVector.size() * 4;
}
int getPosition(Particle *p) {
    int position = std::find(pVector.begin(), pVector.end(), p) - pVector.begin();
}
std::vector<float> ParticleSystem::getState() {
    std::vector<float> state;
    state.resize(this->getDim());
    int ii, size = pVector.size();
    for (ii = 0; ii<size; ii++) {
        Particle *p = pVector[ii];
        state[ii*4] = p->m_Position[0];
        state[ii*4+1] = p->m_Position[1];
        state[ii*4+2] = p->m_Velocity[0];
        state[ii*4+3] = p->m_Velocity[1];
    }
    return state;
}

std::vector<float> ParticleSystem::derivEval() {
    this->resetForces();
    this->applyForces();
    this->applyConstraints();

    std::vector<float> state;
    state.resize(this->getDim());
    int ii, size = pVector.size();
    for (ii = 0; ii<size; ii++) {
        Particle *p = pVector[ii];
        state[ii*4] = p->m_Velocity[0];
        state[ii*4+1] = p->m_Velocity[1];
        state[ii*4+2] = p->m_Force[0] / p->m_Mass;
        state[ii*4+3] = p->m_Force[1] / p->m_Mass;
    }
    return state;
}

void ParticleSystem::setState(std::vector<float> state) {
    int ii, size = pVector.size();
    for (ii=0; ii<size; ii++) {
        pVector[ii]->m_Position[0] = state[ii*4];
        pVector[ii]->m_Position[1] = state[ii*4+1];
        pVector[ii]->m_Velocity[0] = state[ii*4+2];
        pVector[ii]->m_Velocity[1] = state[ii*4+3];
    }
}

void ParticleSystem::resetForces() {
    int ii, size = pVector.size();
    for (ii=0; ii<size; ii++) {
        pVector[ii]->reset_force();
    }
}

void ParticleSystem::applyForces() {
    int ii, size = fVector.size();
    for (ii=0; ii<size; ii++) {
        fVector[ii]->apply();
    }
}
void ParticleSystem::applyConstraints() {
    // we got : C, Cder, J, Jder
    // need to make : M (diagonal mass matrix 3n*3n), W (M^-1), Q (3n long force field), q (3n state vector)
    // calculate J W JT lambda  = -Jd qd - J W Q - ks C - kd Cd
    // use conj gradient to calculate lambda
    // next calculate ^Q = lambda J

    // Define forces/particles/constraints
    std::vector<Particle*> p = pVector;
    std::vector<Constraint*> c = cVector;
    std::vector<Force*> q = fVector;

    // Define sizes
    int noParticles = pVector.size();
    int dim = 2;
    int sizePart = noParticles * dim;
    int noConstraints = cVector.size();

    // define M, W, Q and qder ++ C, Cder, J, Jd and  J^T
    std::vector<std::vector<float>> M(sizePart, std::vector<float>(sizePart));
    std::vector<std::vector<float>> W(sizePart, std::vector<float>(sizePart));

    std::vector<float> Q;
    Q.resize(sizePart);
    std::vector<float> qd;
    qd.resize(sizePart);

    std::vector<float> C;
    C.resize(noConstraints);
    std::vector<float> Cd;
    Cd.resize(noConstraints);
    std::vector<std::vector<float>> J(noConstraints, std::vector<float>(sizePart));
    std::vector<std::vector<float>> Jd(noConstraints, std::vector<float>(sizePart));
    std::vector<std::vector<float>> JTranspose(sizePart, std::vector<float>(noConstraints));


    // compute mass along main diagonal for the particles
    for (int ii = 0; ii<noParticles; ii++) {
        Particle *p = pVector[ii];
        for (int d = 0; d < dim; d++) {
            M[dim * ii + d][dim*ii + d] = p->m_Mass; // only mass of the particles along the diagonal
            W[dim * ii + d][dim*ii + d] = 1 / p->mass; // inverse diagonal matrix is 1/..
            }
    }

    // compute Q and qd for the particles
    for (int i =0; i<noParticles; i++) {
        Particle *p = pVector[i];
        for (int d = 0; d<dim; d++) {
            Q[dim*i + d] = p->m_Force[d]; // list of all forces
            qd[dim*i + d] = p->m_Velocity[d]; // derivative q is velocity
        }
    }

    for (int ii=0; ii<noConstraints; ii++) {
        Constraint* c = cVector[ii];
        C[ii] = c->constraint_value();
        Cd[ii] = c->constraint_derivative_value();
        std::vector<Vec2f> jtemp = c->jacobian_value();
        std::vector<Vec2f> jdtemp = c->jacobian_derivative_value();
        std::vector<Particle *> currentParticles = c->cVector; //TODO true? WHich particles?

        for (int j = 0; j < currentParticles.size(); j++)
        {
            int position = getPosition(currentParticles[j]);
            if (position != -1)
            {
                int index = position * dim;
                for (int d = 0; d < dim; d++)
                {
                    Jd[ii][index + d] = jdtemp[j][d];
                    J[ii][index + d] = jtemp[j][d];
                    JTranspose[index + d][ii] = jtemp[j][d];
                }
            }
            else
            {
                std::cout << "Error position -1";
            }
        }
    }

    // J W JT lambda  = -Jd qd - J W Q - ks C - kd Cd
    double ks = 0.05;
    double kd =0.5;
    std::vector<std::vector<float>> JW = Matrices::matrixMultiplication(J, W);
    std::vector<std::vector<float>> JWJTranspose = Matrices::matrixMultiplication(JW ,JTranspose);
    std::vector<float>  Jdqd =  Matrices::matrixMultiplicationScalar(matrixMultiplication(Jd, qd), -1);
    std::vector<float> JWQ = Matrices::matrixMultiplicationScalar(matrixMultiplication(JW, Q), -1);
    std::vector<float> ksC = Matrices::matrixMultiplicationScalar(C, ks);
    std::vector<float> kdCd = Matrices::matrixMultiplicationScalar(Cd, kd);

    std::vector<std::vector<float>> right = Jdqd - JWQ - ksC - kdCd;

    //TODO
    double ConjGrad(int n, implicitMatrix *A, double x[], double b[],
                    double epsilon,	// how low should we go?
                    int    *steps);

    ConjugateGradient<std::vector<std::vector<float>>, Lower|Upper> conjugategradient;
    auto lambda = conjugategradient.compute(JWJTranspose).solve(right);

    std::vector<float> Qhat = Matrices::matrixMultiplicationScalar(JTranspose, lambda);

    for (int i = 0; i < noParticles; i++) {
        Particle *p = pVector[i];
        int index = i * dim;
        p->m_Force[0] += Qhat[index];
        p->m_Force[1] += Qhat[index + 1];
    }

}

std::vector<Particle*> ParticleSystem::getParticles() {
    return pVector;
}

std::vector<Force*> ParticleSystem::getForces() {
    return fVector;
}

std::vector<Constraint*> ParticleSystem::getConstraints() {
    return cVector;
}

void ParticleSystem::addParticle(Particle *p) {
    pVector.push_back(p);
}

void ParticleSystem::addForce(Force *f) {
    fVector.push_back(f);
}

void ParticleSystem::addConstraint(Constraint *c) {
    cVector.push_back(c);
}