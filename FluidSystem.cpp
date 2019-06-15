#define IX(i,j) ((i)+(N+2)*(j)) // macro to get cell of matrix (converted to 1d array instead of 2d for efficiency)
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;} // swap two pointer arrays
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) { // loop over entire grid
#define END_FOR }}

#include <math.h>

// add density sources to the grid based on user interaction
// s is an array containing the density values for specific cells
void add_source ( int N, float * x, float * s, float dt )
{
	int i, size=(N+2)*(N+2); // size of all arrays ( +2 to simplify bounding box handling)
	for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}

// enforce horizontal velocity = 0 on vertical walls
// and vertical velocity = 0 on horizontal walls
// for density and other fields we assume continuity
void set_bnd ( int N, int b, float * x )
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
}

// Gauss-Seidel relaxation to maintain stable system for large timestep / diff rate etc..
void lin_solve ( int N, int b, float * x, float * x0, float a, float c )
{
	int i, j, k;

	for ( k=0 ; k<20 ; k++ ) {
		FOR_EACH_CELL
			// exchange density with 4 direct neighbors of the cell (net difference)
			x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
		END_FOR
		set_bnd ( N, b, x );
	}
}

void diffuse ( int N, int b, float * x, float * x0, float diff, float dt )
{
	// compute diffusion rate a
	float a=dt*diff*N*N;
	// Gauss-Seidel relaxation used as a stable way to solve diffusion for large timestep / diffusion rate
	lin_solve ( N, b, x, x0, a, 1+4*a );
}

// linear back-trace 
// instead of moving cell centers forward in time
// we look for particles that end up exactly at the cell centers
// by tracing backwards in time from the cell centers
void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	// Start with two grids: one that contains the density values from the previous time step and one
	// that will contain the new values. For each grid cell of the latter we trace the cell’s center
	// position backwards through the velocity field. We then linearly interpolate from the grid of
	// previous density values and assign this value to the current grid cell.

	dt0 = dt*N;
	FOR_EACH_CELL
		// amount of density that the particle carries is obtained by 
		// linearly interpolating the density at their starting location
		// from the four closest neighbors
		x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
		if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
		if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
		s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
		d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
					 s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
	END_FOR
	set_bnd ( N, b, d );
}

// make the fluid mass conserving
void project ( int N, float * u, float * v, float * p, float * div )
{
	int i, j;
	
	// Hodge decomposition:
	// every velocity field is a sum of a mass conserving field and a gradient field
	// computing the gradient is equivalent to computing a height field (Poisson equation)
	// then subtract the height field from the velocity field to obtain the mass conserving field

	FOR_EACH_CELL
		div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
		p[IX(i,j)] = 0;
	END_FOR	
	set_bnd ( N, 0, div ); set_bnd ( N, 0, p );

	lin_solve ( N, 0, p, div, 1, 4 );

	FOR_EACH_CELL
		u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
		v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
	END_FOR
	set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
}

// h and epsilon are constants control spatial discretion and small scale detail
void vorticity_force ( int N, float * u, float * v, float dt, double eps, double h )
{
    //double h = 1.;
    //double eps = 1.;
    int i,j;
    FOR_EACH_CELL
            // calculate force for current position and velocity
            float f = eps * h * pow(2, 0.5) * ( u[IX(i,j)] - v[IX(i,j)] );
            // get current velocity += dt * force
            u[IX(i,j)] += dt * f;
            v[IX(i,j)] += dt * -f;
    END_FOR
}

void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt )
{
	// source density initially contained in the array x0
	add_source ( N, x, x0, dt ); // new density sources may be added by user interaction
	SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt ); // diffuse 
	SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, dt ); // advect
}

void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt, double eps, double h )
{
	// add forces stored in the arrays u0 and v0 re-using the add_source function
	add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
	vorticity_force( N, u, v, dt, eps, h); // vorticity force added to u and v
	SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
	SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
	project ( N, u, v, u0, v0 ); // conserve mass (1st time)
	SWAP ( u0, u ); SWAP ( v0, v ); // swap force arrays
	advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt ); // more accurate advection because of project
	project ( N, u, v, u0, v0 ); // conserve mass (2nd time)
}


