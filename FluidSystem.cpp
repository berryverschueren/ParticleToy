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

int width = 20;
int heigt = 10;

float posX = 20.2;
float posY = 30.2;

int startW;
int endW;

int startH;
int endH;

// add density sources to the grid based on user interaction
// s is an array containing the density values for specific cells
void add_source ( int N, float * x, float * s, float dt )
{
	int i, size=(N+2)*(N+2); // size of all arrays ( +2 to simplify bounding box handling)
	for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}


void newPos(){
	startW = std::floor(posX - width/2.0);
	endW = std::ceil(posX + width/2.0);

	startH = std::floor(posY - heigt/2.0);
	endH = std::ceil(posY + heigt/2.0);

	//std::cout<<"startH "<<startH<< " endH "<<endH<<"\n";
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

	for ( i=1 ; i<endW-startW ; i++ ) {
		x[IX(startW+i,startH  )] = b==2 ? -x[IX(startW+i,startH-1)] : x[IX(startW+i,startH-1)];
		x[IX(startW+i,endH)] = b==2 ? -x[IX(startW+i,endH+1)] : x[IX(startW+i,endH+1)];
	}
	for ( i=1 ; i<endH-startH ; i++ ) {
		x[IX(startW ,startH+i)] = b==1 ? -x[IX(startW-1,startH+i)] : x[IX(startW-1,startH+i)];
		x[IX(endW,startH+i)] = b==1 ? -x[IX(endW+1,startH+i)] : x[IX(endW+1,startH+i)];
	}

	x[IX(startW,startH)] = 0.5f*(x[IX(startW-1,startH)]+x[IX(startW,startH-1)]);
	x[IX(startW,endH  )] = 0.5f*(x[IX(startW-1,endH  )]+x[IX(startW,endH+1  )]);
	x[IX(endW  ,startH)] = 0.5f*(x[IX(endW+1,startH  )]+x[IX(endW,startH-1  )]);
	x[IX(endW  ,endH  )] = 0.5f*(x[IX(endW+1,endH    )]+x[IX(endW,endH+1    )]);
}

// Gauss-Seidel relaxation to maintain stable system for large timestep / diff rate etc..
void lin_solve ( int N, int b, float * x, float * x0, float a, float c )
{
	int i, j, k;

	for ( k=0 ; k<20 ; k++ ) {
		FOR_EACH_CELL
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

bool te= false;
void update(int x,int y,int i,int j, float *d, float *d0, int N, float * u, float * v, bool t){
	if(te){
		d[IX(i+x,j+y)] = -0.5;	
	}else{
		d[IX(i-x,j-y)] = -0.5;
	}
	
	/*if(te){
		d[IX(i,j)] = d0[IX(i,j)] + d[IX(i+x,j+y)];	
		d[IX(i+x,j+y)] = 0;
		if(t) d[IX(i-x,j-y)] = 0.1;
	}else{
		d[IX(i,j)] = d0[IX(i,j)] + d[IX(i+x,j+y)];
		d[IX(i+x,j+y)] = 0;
		if(t)d[IX(i-x,j-y)] = -0.1;
	}*/
	
	//std::cout << i<<" "<<j<< "\n"; 
}

float newX=-1.0;
float newY=-1.0;

// linear back-trace 
// instead of moving cell centers forward in time
// we look for particles that end up exactly at the cell centers
// by tracing backwards in time from the cell centers
void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt, bool t)
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;
	//box
	std::vector<std::string> inside;
	std::vector<std::string> border;
	std::vector<std::string> all;

	/*object forces */
	float forceX=0.0, forceY=0.0;
	if(newX != -1.0){
		forceX = newX - posX;
		posX = newX;
	}
	if(newY != -1.0){
		forceY = newY - posY;
		posY = newY;	
	}
	newPos();
	
	for (i = startW; i <= endW; i++)
	{
		for (j = startH; j <= endH; j++)
		{
			if(j==startH ||j==endH ||i==startW||i==endW){
				border.push_back(std::to_string(i)+std::to_string(j));
			}else{
			    inside.push_back(std::to_string(i)+std::to_string(j));
			}
			all.push_back(std::to_string(i)+std::to_string(j));
		}
	}
	
	// Start with two grids: one that contains the density values from the previous time step and one
	// that will contain the new values. For each grid cell of the latter we trace the cell’s center
	// position backwards through the velocity field. We then linearly interpolate from the grid of
	// previous density values and assign this value to the current grid cell.
	//int ii = 0;
	dt0 = dt*N;
	float force = 0.0;
	FOR_EACH_CELL
		// amount of density that the particle carries is obtained by 
		// linearly interpolating the density at their starting location
		// from the four closest neighbors
		auto cell = std::to_string(i)+std::to_string(j);

		if(std::find(border.begin(), border.end(), cell) != border.end()) {
			
			if(t){
				force += d[IX(i,j)];
				//std::cout<< "f " <<d[IX(i,j)]<<"\n";
			}
		}
	
		if(std::find(all.begin(), all.end(), cell) != all.end()) {
			//in object
			//is force
			d[IX(i,j)] = 0;
			if(forceX != 0.0 && forceY != 0.0){
				if(t){
					if(b==1){
						//horizontal
						d[IX(i,j)] = forceX;		
					}else if(b==2){
						//vertical
						d[IX(i,j)] = forceY;
					}
				}
			}
			//}
		}else{
			x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];

			if (x<0.5f) x=0.5f; 
			if (x>N+0.5f) x=N+0.5f; 
			i0=(int)x; i1=i0+1;
			
			if (y<0.5f) y=0.5f;
			if (y>N+0.5f) y=N+0.5f; 
			j0=(int)y; j1=j0+1;

			s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;

			d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
						s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
			
		}
		
		//if(std::find(border.begin(), border.end(), cell) != border.end()) {
			
			//if(t){
				//force += d[IX(i,j)];
				//std::cout<< "f " <<d[IX(i,j)]<<"\n";
			//}
			//ii++;
			//d[IX(i,j)] = 1;
			//find nighboring border
			/*auto cellr = std::to_string(i+1)+std::to_string(j);
			auto celll = std::to_string(i-1)+std::to_string(j);
			auto cellu = std::to_string(i)+std::to_string(j+1);
			auto celld = std::to_string(i)+std::to_string(j-1);

			if(std::find(inside.begin(), inside.end(), cellr) != inside.end()){
				update(1,0,i,j,d,d0,N, u, v, t);
				
			}if(std::find(inside.begin(), inside.end(), celll) != inside.end()){
				update(-1,0,i,j,d,d0,N, u, v, t);
				
			}if(std::find(inside.begin(), inside.end(), cellu) != inside.end()){
				update(0,1,i,j,d,d0,N,u, v, t);
				
			}if(std::find(inside.begin(), inside.end(), celld) != inside.end()){
				update(0,-1,i,j,d,d0,N,u,v, t);
			}*/
		//}
	END_FOR

	if(b==1){
		//horizontal force
		//posX = posX - force*20.0;
		//std::cout<< "h " <<force<<"\n";
	}else if(b==2){
		//posY = posY - force*20.0;
		//std::cout<< "v " <<posY<<"\n";
	}

	//std::cout<<"count "<< std::to_string(ii)<<" num ="<<width*heigt << " border= "<<border.size() << "\n";
	//clip_path( N, b, d);
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

void new_object_position(int N, int xx, int yy){
	int i, j;
	
	for (i = startW; i <= endW; i++)
	{
		for (j = startH; j <= endH; j++)
		{
			
			if(i*8 <=xx && ((i+1)*8-1)> xx){
				//std::cout<<"i "<<i*8<<" x "<<xx<<"\n";
				if(j*8 <=512-yy && ((j+1)*8-1)> 512-yy){
					//std::cout<<"here "<<"\n";
					newX = i+0.5;
					newY = j+0.5;
					//std::cout<<"newx "<<newX<<" newY "<<newY<<"\n";
				}else{
					//newX = -1.0;
					//newY = -1.0;
				}
			}
		}
	}
	
	
}

void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt )
{
	//newPos();
	// source density initially contained in the array x0
	add_source ( N, x, x0, dt ); // new density sources may be added by user interaction
	SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt ); // diffuse 
	SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, dt, false); // advect
}

void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt, int xx, int yy)
{
	newPos();
	new_object_position(N, xx, yy);

	// add forces stored in the arrays u0 and v0 re-using the add_source function
	add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
	// 
	SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
	SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );

	project ( N, u, v, u0, v0 ); // conserve mass (1st time)

	SWAP ( u0, u ); SWAP ( v0, v ); // swap force arrays
	advect ( N, 1, u, u0, u0, v0, dt, true); advect ( N, 2, v, v0, u0, v0, dt, true); // more accurate advection because of project
	project ( N, u, v, u0, v0 ); // conserve mass (2nd time)
}

