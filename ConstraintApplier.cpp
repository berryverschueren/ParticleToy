#include "ConstraintApplier.h"
#include <vector>
#include <cmath>

ConstraintApplier::ConstraintApplier(ParticleSystem *s : s) {
    // we got : C, Cder, J, Jder
    // need to make : M (diagonal mass matrix 3n*3n), W (M^-1), Q (3n long force field), q (3n state vector)
    // calculate lambda^T = (JWJ^T)^-1 (-Jder qder^T - J W Q^T)
    // next calculate ^Q = lambda J

    //std::vector<std::vector<float>> matrix(row, std::vector<float>(column));

    // Define forces/particles/constraints
    // TODO fix
    vector<Particle*> p = s->pVector;
    vector<Constraint*> c = s->cVector;
    vector<Force*> q = s->fVector;


    // Define sizes
    int noParticles = p->pVector.size();
    int dim = 2;
    int sizePart = noParticles * dim;
    int noConstraints = c->cVector.size();

    // define M, W, Q and qder
    std::vector<std::vector<float>> M(sizePart, std::vector<float>(sizePart));
    std::vector<std::vector<float>> W(sizePart, std::vector<float>(sizePart));
    //float M[totSizePart][totSizePart];
    //float W[totSizePart][totSizePart];
    std::vector<float> Q;
    Q.resize(sizePart);
    std::vector<float> qd;
    qd.resize(sizePart);

    // compute mass along main diagonal for the particles
    for (int i = 0; i<sizePart; i+= dim) {
        Particle *p = p[floor(i / dim)];
        M[i][i]= p->m_Mass; // only mass of the particles along the diagonal
        W[i][i]= 1 / p->m_Mass; // inverse of a matrix with only a diagonal
    }

    // compute Q and qd for the particles
    for (int i =0; i<sizePart; i+=dim) {
        Particle *p = p[i/dim];
        for (int j =0; j<dim; j++) {
            Q[i+j]=p->m_Force[j];  // list of all forces
            qd[i+j]= p->m_Velocity[j]; // derivative position is velocity
        }
    }

    // define C, Cder, J, Jder and  J^T
    std::vector<float> C;
    C.resize(noConstraints);
    std::vector<float> Cd;
    Cd.resize(noConstraints);
    std::vector<std::vector<float>> J(noConstraints, std::vector<float>(sizePart));
    std::vector<std::vector<float>> Jd(noConstraints, std::vector<float>(sizePart));
    std::vector<std::vector<float>> JTranspose(sizePart, std::vector<float>(noConstraints));

    //float J[sizeConst][totSizePart];
    //float Jd[sizeConst][totSizePart];
    //float JTranspose[totSizePart][sizeConst];

    for (int ii=0; ii<noConstraints; ii++) {
        Constraint* c = s->cVector[ii];
        C[ii] = c->constraint_value();
        Cd[ii] = c->constraint_derivative_value();
        std::vector<Vec2f> jtemp = c->jacobian_value();
        std::vector<Vec2f> tdtemp = c->jacobian_derivative_value();
        // TODO fill matrix
    }

    // lambda T = (J W JT)^-1 (-Jd qdT - J W QT)
    // Qhat = lambda J



    // calculate left hand side
    //JWJTransposeLambda = std::vector<float> {- Jd * qd - J * W * Q};

    // calculate conjugate gradient of left hand side
    //ConjugateGradient<MatrixXf, Lower|Upper> conjGrad;

    //conjGrad.compute(J*W*JTranspose);
    //lambda = std::vector<float> {cg.solve(JWJTransposeLambda)};

}


