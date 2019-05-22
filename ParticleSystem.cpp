#include "ParticleSystem.h"
#include <cmath>
#include <vector>
#include "Matrices.h"
#include "linearSolver.h"
#include <typeinfo>
#include "CircularWireConstraint.h"
#include "RodConstraint.h"
#include <algorithm>

ParticleSystem::ParticleSystem() {};

int ParticleSystem::getDim() {
    return pVector.size() * 4;
}

int ParticleSystem::getPosition(Particle *p) {
	int pos = std::find(pVector.begin(), pVector.end(), p) - pVector.begin();
	if (pos < pVector.size())
	{

		return pos;
	}
	return -1;
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
    const int dimensions = 2;
	int vectorSize = pVector.size() * dimensions;
	int constraintsSize = cVector.size();
	float ks = 0.005f;
	float kd = 0.05f;

	VectorXf q = VectorXf::Zero(vectorSize);
	VectorXf Q = VectorXf::Zero(vectorSize);
	MatrixXf W = MatrixXf::Zero(vectorSize, vectorSize);
  
	for (int i = 0; i < pVector.size(); i ++) {
		Particle *p = pVector[i];
		for (int d = 0; d < dimensions; d++) {
			W(dimensions * i + d,dimensions*i + d) = 1 / p->m_Mass;
			Q[dimensions*i + d] = p->m_Force[d];
			q[dimensions*i + d] = p->m_Velocity[d];
		}
	}

	VectorXf C = VectorXf::Zero(constraintsSize);
	VectorXf Cder = VectorXf::Zero(constraintsSize);

	for (int i = 0; i < constraintsSize; i++) {
		Constraint *c = cVector[i];
        std::vector<Vec2f> j;
        std::vector<Vec2f> jd;
        std::vector<Particle *> currentParticles;

        if(CircularWireConstraint* con = dynamic_cast<CircularWireConstraint*>(c)) {
            C[i] = con->constraint_value();
		    Cder[i] = con->constraint_derivative_value();
		    j = con->jacobian_value();
    		jd = con->jacobian_derivative_value();
            currentParticles.push_back(con->m_p);
        } else if (RodConstraint* con = dynamic_cast<RodConstraint*>(c)) {
            C[i] = con->constraint_value();
		    Cder[i] = con->constraint_derivative_value();
		    j = con->jacobian_value();
    		jd = con->jacobian_derivative_value();
            currentParticles.push_back(con->m_p1);
            currentParticles.push_back(con->m_p2);
        }
        else {
            std::cout << "Couldnt cast constraint";
        }

	    MatrixXf J = MatrixXf::Zero(constraintsSize, vectorSize);
	    MatrixXf Jt = MatrixXf::Zero(vectorSize, constraintsSize);
        MatrixXf Jder = MatrixXf::Zero(constraintsSize, vectorSize);
		
        for (int k = 0; k < currentParticles.size(); k++) {
			int currentPos = getPosition(currentParticles[k]);
			if (currentPos != -1) {
				int pIndex = currentPos * dimensions;
				for (int d = 0; d < dimensions; d++) {
					Jder(i,pIndex + d) = jd[k][d];
					J(i,pIndex + d) = j[k][d];
					Jt(pIndex + d,i) = j[k][d];
				}
			}
		}
	}

	MatrixXf JW = J * W;
	MatrixXf JWJt = JW * Jt;
	VectorXf Jderq = -1 * Jder * q;
	VectorXf JWQ = JW * Q;
	VectorXf KsC = ks * C;
	VectorXf KdCd = kd * Cder;
	VectorXf rhs = Jderq - JWQ - KsC - KdCd;

	ConjugateGradient<MatrixXf, Lower|Upper> cg;
	auto lambda = cg.compute(JWJt).solve(rhs);
	VectorXf Qhat = Jt * lambda;

	for (int i = 0; i < pVector.size(); i++) {
		Particle *p = pVector[i];
		int index = i * dimensions;
		p->m_Force[0] += Qhat[index];
		p->m_Force[1] += Qhat[index + 1];
	}
}

std::vector<float> ParticleSystem::simpleCollision(std::vector<float> state){

    float floor = -1.0; float ceiling = 1.0;
    float wallL = -1.0; float wallR = 1.0;
    float k_r = 0.8;

    //floor
    for (int i = 0; i < pVector.size() ; ++i) {
        if(state[i*4+1] <= floor){
            float vN = state[i*4+3];
            if(vN >= 0) {
                state[i * 4 + 3] = -k_r * vN;
            } else{
                state[i*4+3] = 0.001;
            }
        }
    }

    //ceiling
    for (int i = 0; i < pVector.size() ; ++i) {
        if(state[i*4+1] >= ceiling){
            float vN = state[i*4+3];
            if(vN <= 0) {
                state[i * 4 + 3] = -k_r * vN;
            } else{
                state[i*4+3] = -0.001;
            }

        }
    }

    //walls
    for (int i = 0; i < pVector.size() ; ++i) {
        float vN = state[i*4+2];
        float reversedForce = -k_r*vN;
        if(state[i*4] < wallL){
            if(reversedForce >= 0){
                state[i*4+2] = reversedForce;
            } else{
                state[i*4+2] = 0.001;
            }
        }
        if(state[i*4] > wallR){
            if(reversedForce <= 0) {
                state[i*4+2] = reversedForce;
            } else{
                state[i*4+2] = -0.001;
            }
        }
    }

    return state;
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

void ParticleSystem::removeLastParticle() {
    pVector.pop_back();
}

void ParticleSystem::addForce(Force *f) {
    fVector.push_back(f);
}

void ParticleSystem::removeLastForce() {
    fVector.pop_back();
}

void ParticleSystem::addConstraint(Constraint *c) {
    cVector.push_back(c);
}

void ParticleSystem::removeLastConstraint() {
    cVector.pop_back();
}

void ParticleSystem::deleteAll() {
    pVector.clear();
    fVector.clear();
    cVector.clear();
}
