#include "FluidSystem.h"
#include <cmath>
#include <vector>
#include <typeinfo>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>

FluidSystem::FluidSystem(int _N, float _dt, float _diff, float _visc, float _force, float _source, int _dvel) {
    this->N = _N;
    this->dt = _dt;
    this->diff = _diff;
    this->visc = _visc;
    this->force = _force;
    this->source = _source;
    this->dvel = _dvel;
};

void FluidSystem::free_data() {
    uForces.clear();
    vForces.clear();
    uForcesPrev.clear();
    vForcesPrev.clear();
    densities.clear();
    densitiesPrev.clear();
}

void FluidSystem::clear_data() {
    int i, size = (N+2) * (N+2);
    for (i = 0; i < size; i++) {
        uForces[i] = vForces[i] = uForcesPrev[i] = vForcesPrev[i] = densities[i] = densitiesPrev[i] = 0.0f;
    }
}

int FluidSystem::allocate_data() {
    int size = (N+2) * (N+2);
    uForces.resize(size);
    vForces.resize(size);
    uForcesPrev.resize(size);
    vForcesPrev.resize(size);
    densities.resize(size);
    densitiesPrev.resize(size);
    return 1;
}

void FluidSystem::draw_velocity() {
    int i, j;
	float x, y, h;
	h = 1.0f/N; 

	glColor3f ( 1.0f, 1.0f, 1.0f );
	glLineWidth ( 1.0f );

	glBegin ( GL_LINES );

		for ( i=1 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=1 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;

				glVertex2f ( x, y );
				glVertex2f ( x+uForces[((i)+(N+2)*(j))], y+vForces[((i)+(N+2)*(j))] );
			}
		}

	glEnd ();
}

void FluidSystem::draw_density() {
	int i, j;
	float x, y, h, d00, d01, d10, d11;
	h = 1.0f/N;

	glBegin ( GL_QUADS );

		for ( i=0 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=0 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;

				d00 = densities[((i)+(N+2)*(j))];
				d01 = densities[((i)+(N+2)*(j+1))];
				d10 = densities[((i+1)+(N+2)*(j))];
				d11 = densities[((i+1)+(N+2)*(j+1))];

				glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
				glColor3f ( d10, d10, d10 ); glVertex2f ( x+h, y );
				glColor3f ( d11, d11, d11 ); glVertex2f ( x+h, y+h );
				glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+h );
			}
		}

	glEnd ();
}

void FluidSystem::add_source(std::vector<float> &x, std::vector<float> &s) {
    int i, size = (N+2) * (N+2);
    for (i = 0; i < size; i++) {
        x[i] += dt * s[i];
    }
}

void FluidSystem::linear_solve(int b, std::vector<float> &x, std::vector<float> &x0, float a, float c) {
    int i,j,k;
    for (k = 0; k < 20; k++) {
        for (i = 1; i <= N; i++ ) { 
            for (j = 1; j <= N; j++ ) {
                x[((i)+(N+2)*(j))] = (x0[((i)+(N+2)*(j))] + a*(x[((i-1)+(N+2)*(j))]+x[((i+1)+(N+2)*(j))]+x[((i)+(N+2)*(j-1))]+x[((i)+(N+2)*(j+1))]))/c;
            }
        }
        set_boundary(b, x);
    }    
}

void FluidSystem::diffuse(int b, std::vector<float> &x, std::vector<float> &x0) {
    float a = dt * diff * N * N;
    linear_solve(b, x, x0, a, (1 + 4 * a));
}

void FluidSystem::advect(int b, std::vector<float> &d, std::vector<float> &d0, std::vector<float> &u, std::vector<float> &v) {
    int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;
    dt0 = dt * N;
    
    for (i = 1; i <= N; i++ ) { 
        for (j = 1; j <= N; j++ ) {
            x = i - dt0 * u[((i)+(N+2)*(j))]; 
            y = j - dt0 * v[((i)+(N+2)*(j))];
            
            if (x < 0.5f) {
                x = 0.5f;
            } 
            
            if (x > N + 0.5f) {
                x = N + 0.5f; 
            }
            
            i0 = (int)x; 
            i1 = i0 + 1;
            
            if (y < 0.5f) {
                y = 0.5f;
            } 
            
            if (y > N + 0.5f) {
                y = N + 0.5f;
            }
             
            j0 = (int)y; 
            j1 = j0 + 1;

            s1 = x - i0; 
            s0 = 1 - s1; 
            t1 = y - j0; 
            t0 = 1 - t1;

            d[((i)+(N+2)*(j))] = s0 * (t0 * d0[((i0)+(N+2)*(j0))]
                + t1 * d0[((i0)+(N+2)*(j1))]) 
                + s1 * (t0 * d0[((i1)+(N+2)*(j0))] 
                + t1 * d0[((i1)+(N+2)*(j1))]);
        }
    }

    set_boundary(b, d);
}

void FluidSystem::project(std::vector<float> &u, std::vector<float> &v, std::vector<float> &p, std::vector<float> &div) {
    int i,j;

    for (i = 1; i <= N; i++ ) { 
        for (j = 1; j <= N; j++ ) {
            div[((i)+(N+2)*(j))] = -0.5f*(u[((i+1)+(N+2)*(j))]-u[((i-1)+(N+2)*(j))]+v[((i)+(N+2)*(j+1))]-v[((i)+(N+2)*(j-1))])/N;
		    p[((i)+(N+2)*(j))] = 0;
        }
    }
    
    set_boundary(0, div);
    set_boundary(0, p);
    linear_solve(0, p, div, 1, 4);

    for (i = 1; i <= N; i++ ) { 
        for (j = 1; j <= N; j++ ) {
            u[((i)+(N+2)*(j))] -= 0.5f*N*(p[((i+1)+(N+2)*(j))]-p[((i-1)+(N+2)*(j))]);
            v[((i)+(N+2)*(j))] -= 0.5f*N*(p[((i)+(N+2)*(j+1))]-p[((i)+(N+2)*(j-1))]);
        }
    }

    set_boundary(1, u);
    set_boundary(2, v);
}

void FluidSystem::density_step() {
    add_source(densities, densitiesPrev);
    densitiesPrev.swap(densities);
    diffuse(0, densities, densitiesPrev);
    densitiesPrev.swap(densities);
    advect(0, densities, densitiesPrev, uForces, vForces);
}

void FluidSystem::velocity_step() {
    add_source(uForces, uForcesPrev);
    add_source(vForces, vForcesPrev);
    uForcesPrev.swap(uForces);
    diffuse(1, uForces, uForcesPrev);
    vForcesPrev.swap(vForces);
    diffuse(2, vForces, vForcesPrev);
    project(uForces, vForces, uForcesPrev, vForcesPrev);
    uForcesPrev.swap(uForces);
    vForcesPrev.swap(vForces);
    advect(1, uForces, uForcesPrev, uForcesPrev, vForcesPrev);
    advect(2, vForces, vForcesPrev, uForcesPrev, vForcesPrev);
    project(uForces, vForces, uForcesPrev, vForcesPrev);
}

void FluidSystem::set_boundary(int b, std::vector<float> &x) {
    int i;
    for (i = 0; i <= N; i++) {
        x[((0)+(N+2)*(i))] = b==1 ? -x[((1)+(N+2)*(i))] : x[((1)+(N+2)*(i))];
		x[((N+1)+(N+2)*(i))] = b==1 ? -x[((N)+(N+2)*(i))] : x[((N)+(N+2)*(i))];
		x[((i)+(N+2)*(0))] = b==2 ? -x[((i)+(N+2)*(1))] : x[((i)+(N+2)*(1))];
		x[((i)+(N+2)*(N+1))] = b==2 ? -x[((i)+(N+2)*(N))] : x[((i)+(N+2)*(N))];
    }
    x[((0)+(N+2)*(0))] = 0.5f*(x[((1)+(N+2)*(0))]+x[((0)+(N+2)*(1))]);
	x[((0)+(N+2)*(N+1))] = 0.5f*(x[((1)+(N+2)*(N+1))]+x[((0)+(N+2)*(N))]);
	x[((N+1)+(N+2)*(0))] = 0.5f*(x[((N)+(N+2)*(0))]+x[((N+1)+(N+2)*(1))]);
	x[((N+1)+(N+2)*(N+1))] = 0.5f*(x[((N)+(N+2)*(N+1))]+x[((N+1)+(N+2)*(N))]);
}