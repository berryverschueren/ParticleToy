#define IX(i,j) ((i)+(N+2)*(j)) // macro to get cell of matrix (converted to 1d array instead of 2d for efficiency)
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;} // swap two pointer arrays
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) { // loop over entire grid
#define END_FOR }}
#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include <string> 
#include <math.h> 
#include "include/Eigen/Dense"
#include "include/Eigen/IterativeLinearSolvers"
using namespace Eigen;

void add_source ( int N, float * x, float * s, float dt )
{
	int i, size=(N+2)*(N+2);
	for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}

void set_bnd ( int N, int b, float * x, float * grid )
{
	int i;

	for ( i=1 ; i<=N ; i++ ) {
		x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
		x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
		x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
		x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
	}
	x[IX(0  ,0  )] = 0.5f*(x[IX(1,0  )]+x[IX(0  ,1)]);
	x[IX(0  ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0  ,N)]);
	x[IX(N+1,0  )] = 0.5f*(x[IX(N,0  )]+x[IX(N+1,1)]);
	x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N)]);

	int j;

	FOR_EACH_CELL
		if (grid[IX(i,j)]==1) {			
			float div = 1.0f;
			float count = 0.0f;
			float force_left = 0.0f,force_right  = 0.0f,force_up= 0.0f,force_down= 0.0f;

			if (grid[IX(i-1,j)]==0) {
				force_left = b==1 ? -x[IX(i-1,j)] : x[IX(i-1,j)];
				x[IX(i,j)] = force_left;
				count++;
			}
			
			if (grid[IX(i+1,j)]==0) {
				force_right = b==1 ? -x[IX(i+1,j)] : x[IX(i+1,j)];
				x[IX(i,j)] = force_right;
				count++;
			}
			if (grid[IX(i,j-1)]==0) {
				force_down = b==2 ? -x[IX(i,j-1)] : x[IX(i,j-1)];
				x[IX(i,j)] = force_down;
				count++;
			}
			if (grid[IX(i,j+1)]==0) {
				force_up = b==2 ? -x[IX(i,j+1)] : x[IX(i,j+1)];
				x[IX(i,j)] = force_up;
				count++;
			}			
			
			if (count > 1) {
				x[IX(i,j)]=div/count*(force_left+force_right+force_up+force_down);
			}
		}
	END_FOR
}

void lin_solve ( int N, int b, float * x, float * x0, float a, float c, float * grid )
{
	int i, j, k;

	for ( k=0 ; k<20 ; k++ ) {
		FOR_EACH_CELL
			x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;		
		END_FOR
		set_bnd ( N, b, x, grid );
	}
}

void diffuse ( int N, int b, float * x, float * x0, float diff, float dt, float * grid)
{
	float a=dt*diff*N*N;
	lin_solve ( N, b, x, x0, a, 1+4*a, grid );
}

void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt, float * grid, Vector2f &genForce)
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt*N;
	
	FOR_EACH_CELL
		if (grid[IX(i,j)] == 2 || grid[IX(i,j)] == 1) {
			// ignore cell inside object
		} else {
			x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];

			/* todo:: implement path clipping
			if x,y == inside object -> find closest cell outside
			set x,y to coordinates of that cell */

			if (x<0.5f) x=0.5f; 
			if (x>N+0.5f) x=N+0.5f; 
			i0=(int)x; i1=i0+1;
			if (y<0.5f) y=0.5f;
			if (y>N+0.5f) y=N+0.5f; 
			j0=(int)y; j1=j0+1;
			s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;

			// check if interpolation cells are inside the object
			// if so, include forces from rigid body in the new
			// velocity field
			auto a = grid[IX(i0,j0)] == 1 || grid[IX(i0,j0)] == 2;
			auto c = grid[IX(i0,j1)] == 1 || grid[IX(i0,j1)] == 2;
			auto e = grid[IX(i1,j0)] == 1 || grid[IX(i1,j0)] == 2;
			auto f = grid[IX(i1,j1)] == 1 || grid[IX(i1,j1)] == 2;
			
			// don't take special care for advect in dens_step (yet...)
			if (b != 1 && b != 2) { a = c = e = f = false; }

			// initialize participation values for all interpolation cells
			auto aVal = 0.0f, cVal = 0.0f, eVal = 0.0f, fVal = 0.0f;

			// compute force from object (for either u or v)
			auto additionFromGrid = b == 1 ? genForce[0] : genForce[1];

			// if interpolation cell is inside, use the participation value
			if (a) aVal = additionFromGrid;
			if (c) cVal = additionFromGrid;
			if (e) eVal = additionFromGrid;
			if (f) fVal = additionFromGrid;

			// if none of the interpolation cells are inside
			// the object, just do normal advection
			if (!a && !c && !e && !f) {
				d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
						s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);	
			} 
			// if one of them is inside, take special care
			else {
				d[IX(i,j)] = (aVal+cVal+eVal+fVal)/4;	
				// d[IX(i,j)] = (aVal+cVal+eVal+fVal);	
				// d[IX(i,j)] = additionFromGrid;	
			}
		}
	END_FOR
	set_bnd ( N, b, d, grid );
}

void project ( int N, float * u, float * v, float * p, float * div, float * grid )
{
	int i, j;

	FOR_EACH_CELL
		div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
		p[IX(i,j)] = 0;
	END_FOR	
	set_bnd ( N, 0, div, grid ); set_bnd ( N, 0, p, grid );

	lin_solve ( N, 0, p, div, 1, 4, grid);

	FOR_EACH_CELL
		u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
		v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
	END_FOR
	set_bnd ( N, 1, u, grid ); set_bnd ( N, 2, v, grid );
}

void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt, float* grid, Vector2f &genForce)
{
	add_source ( N, x, x0, dt );
	SWAP ( x0, x ); 
	diffuse ( N, 0, x, x0, diff, dt, grid);
	SWAP ( x0, x ); 
	advect ( N, 0, x, x0, u, v, dt, grid, genForce);
}

void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc,
 	float dt, int xx, int yy, float * grid, Vector2f &genForce)
{
	add_source ( N, u, u0, dt ); 
	add_source ( N, v, v0, dt );
	SWAP ( u0, u ); 
	diffuse ( N, 1, u, u0, visc, dt, grid);
	SWAP ( v0, v ); 
	diffuse ( N, 2, v, v0, visc, dt, grid);
	project ( N, u, v, u0, v0, grid);
	SWAP ( u0, u ); 
	SWAP ( v0, v );
	advect ( N, 1, u, u0, u0, v0, dt, grid, genForce); 
	advect ( N, 2, v, v0, u0, v0, dt, grid, genForce);
	project ( N, u, v, u0, v0, grid);
}

