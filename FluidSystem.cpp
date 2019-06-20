#define IX(i,j) ((i)+(N+2)*(j)) // macro to get cell of matrix (converted to 1d array instead of 2d for efficiency)
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;} // swap two pointer arrays
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) { // loop over entire grid
#define END_FOR }}
#define _USE_MATH_DEFINES
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
void set_bnd ( int N, int b, float * x, float * grid)
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
	//std::cout<<grid[IX(i,j)]<<"\n";
		if(grid[IX(i,j)]==1){
			
			float div = 1.0f;
			float count = 0.0f;
			float force_left = 0.0f,force_right  = 0.0f,force_up= 0.0f,force_down= 0.0f;

			//horizontal
			//check nighbors
			if(grid[IX(i-1,j)]==0){
				force_left = b==1 ? -x[IX(i-1,j)] : x[IX(i-1,j)];
				x[IX(i,j)] = force_left;
				count++;
			}
			if(grid[IX(i+1,j)]==0){
				force_right = b==1 ? -x[IX(i+1,j)] : x[IX(i+1,j)];
				x[IX(i,j)] = force_right;
				count++;
			}
			if(grid[IX(i,j-1)]==0){
				force_down = b==2 ? -x[IX(i,j-1)] : x[IX(i,j-1)];
				x[IX(i,j)] = force_down;
				count++;
			}
			if(grid[IX(i,j+1)]==0){
				force_up = b==2 ? -x[IX(i,j+1)] : x[IX(i,j+1)];
				x[IX(i,j)] = force_up;
				count++;
			}
			
			
			if(count > 1){
				//std::cout<<(force_left+force_right+force_up+force_down)<<"\n";
				x[IX(i,j)]=div/count*(force_left+force_right+force_up+force_down);
			}
		}
	END_FOR
 
 /*
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
	*/
}

// Gauss-Seidel relaxation to maintain stable system for large timestep / diff rate etc..
void lin_solve ( int N, int b, float * x, float * x0, float a, float c, float * grid )
{
	int i, j, k;

	for ( k=0 ; k<20 ; k++ ) {
		FOR_EACH_CELL
			x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;		
		END_FOR
		set_bnd ( N, b, x, grid);
	}
}

void diffuse ( int N, int b, float * x, float * x0, float diff, float dt, float * grid  )
{
	// compute diffusion rate a
	float a=dt*diff*N*N;
	// Gauss-Seidel relaxation used as a stable way to solve diffusion for large timestep / diffusion rate
	lin_solve ( N, b, x, x0, a, 1+4*a, grid );
}

bool te= false;
bool ttt= true;

// linear back-trace 
// instead of moving cell centers forward in time
// we look for particles that end up exactly at the cell centers
// by tracing backwards in time from the cell centers
void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt, bool t, float * grid, float xx, float yy, int &centerX, int &centerY)
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;
	
	// Start with two grids: one that contains the density values from the previous time step and one
	// that will contain the new values. For each grid cell of the latter we trace the cell’s center
	// position backwards through the velocity field. We then linearly interpolate from the grid of
	// previous density values and assign this value to the current grid cell.
	//int ii = 0;
	dt0 = dt*N;
	float force = 0.0;
	FOR_EACH_CELL

		if(grid[IX(i,j)]==1 || grid[IX(i,j)]==2){
			//in object
			//is force
			int xxx = (int)xx;
			int yyy = (int)yy;
			//if((xxx>0||yyy>0) && grid[IX(i,j)]==2)std::cout<<xxx<<" "<<yyy<<"\n";

			if(!t){
				//std::cout<<"test"<<xxx<<"\n";
				if(grid[IX(i,j)]==2){
					if(xxx!=0||yyy!=0){	
						d0[IX(i+xxx,j+yyy)] += d0[IX(i,j)];
						//std::cout<<"test"<<xxx<<"\n";
						d0[IX(i,j)] = 0;
						d[IX(i,j)] = 0;
					}else{
						bool border = false;
						int c=0;
						if(d0[IX(i,j)] != 0){
							
							while(!border){
								c++;
								if(grid[IX(i+c,j)]==0){
									border = true;
									
									d0[IX(i+c,j)] += d0[IX(i,j)]; 
									d0[IX(i,j)] = 0;
									d[IX(i,j)] = 0;
								}
							}
						}
					}
				}
			}
			
			
			if(xxx != 0.0 && yyy != 0.0){
				if(t){
					if(b==1){
						//horizontal
						d0[IX(i,j)] = xxx;		
					}else if(b==2){
						//vertical
						d0[IX(i,j)] = yyy;
					}
				}
			}
		}

	END_FOR

	FOR_EACH_CELL

		// amount of density that the particle carries is obtained by 
		// linearly interpolating the density at their starting location
		// from the four closest neighbors


		if(grid[IX(i,j)]==1){	
			if(t){
				force -= d0[IX(i,j)];
			}
		}
	

		if(grid[IX(i,j)]==1 || grid[IX(i,j)]==2){
			/*
			x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];

			if (x<0.5f) x=0.5f; 
			if (x>N+0.5f) x=N+0.5f; 
			i0=(int)x; i1=i0+1;
			
			if (y<0.5f) y=0.5f;
			if (y>N+0.5f) y=N+0.5f; 
			j0=(int)y; j1=j0+1;

			s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;

			if(grid[IX(i-1,j)]==0){
				d[IX(i-1,j)] += s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
						s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
			}else if(grid[IX(i+1,j)]==0){
				d[IX(i+1,j)] += s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
						s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
			}else if(grid[IX(i,j-1)]==0){
				d[IX(i,j-1)] += s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
						s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
			}else if(grid[IX(i,j+1)]==0){
				d[IX(i,j+1)] += s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
						s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
			}*/
						
		}else{
			

			x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];

			if(x!=i &&y!=j){ 
				//std::cout<<"i "<<i<<" j "<<j<<"\n";
				//std::cout<<"x "<<x<<" y "<<y<<"\n";
			}

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
	
	END_FOR







	//if(b==1){
		//if((force > 0.5f || force < -0.5f)){
			

		//}

	//}else if(b==2 ){
		/*
		if((force > 0.5f || force < -0.5f)){
			int moveY = (force/std::abs(force));

			new_grid = (float *) malloc (size*sizeof(float));
			FOR_EACH_CELL
				new_grid[IX(i,j)] = grid[IX(i, j-moveY)];
			END_FOR
			FOR_EACH_CELL
				grid[IX(i,j)] = new_grid[IX(i, j)];
			END_FOR
			centerY = centerY +moveY;
		}*/
	//}

	

	//std::cout<<"count "<< std::to_string(ii)<<" num ="<<width*heigt << " border= "<<border.size() << "\n";
	//clip_path( N, b, d);
	set_bnd ( N, b, d, grid);

}

void angle(int &x, int &y, float i, float j, int centerX, int centerY, float zz){
	x = (int)((j - centerY) *std::sin(M_PI*(zz)) + (i - centerX)*std::cos(M_PI*(zz))+ centerX);
	y = (int)((j - centerY) *std::cos(M_PI*(zz)) - (i - centerX)*std::sin(M_PI*(zz))+ centerY);
}

float prev_zz;

void add(int N, float * grid, int centerX, int centerY, float zz, float * prev_grid){
	int size = (N+2)*(N+2), i,j;
	float* new_grid;
	

	if(ttt && prev_zz!=zz){
		//ttt=false;
				
		new_grid = (float *) malloc (size*sizeof(float));
		//int moveX = (force/std::abs(force));
		//std::cout<<"centerX "<<centerX<<" centerY "<<centerY<<"\n";
		FOR_EACH_CELL
			new_grid[IX(i,j)] = 0.0f;
		END_FOR
		FOR_EACH_CELL
		if(prev_grid[IX(i,j)] != 0){

			//int x = (int)((j - centerY) *std::sin(3.14*(1.0/4.0)) + (i - centerX)*std::cos(3.14*(1.0/4.0))+ centerX);
			//int y = (int)((j - centerY) *std::cos(3.14*(1.0/4.0)) - (i - centerX)*std::sin(3.14*(1.0/4.0))+ centerY);
			int x,y;
			angle(x,y,i,j,centerX,centerY, zz);

			if(new_grid[IX(x,y)] != 1){
				new_grid[IX(x,y)] = grid[IX(i, j)];
			}

			angle(x,y,i+0.3,j+0.3,centerX,centerY, zz);
			
			if(new_grid[IX(x,y)] != 1){
				new_grid[IX(x,y)] = grid[IX(i, j)];
			}
			//int x1,y1;
			//angle(x,y,i+1,j+1,centerX,centerY);

			//if(new_grid[IX(x1,y1)] != 1){
				//new_grid[IX(x1,y1)] = grid[IX(i, j)];
			//}
/*
			angle(x,y,i,j+1,centerX,centerY);

			if(new_grid[IX(x,y)] != 1){
				new_grid[IX(x,y)] = grid[IX(i, j)];
			}

			angle(x,y,i+1,j+1,centerX,centerY);

			if(new_grid[IX(x,y)] != 1){
				new_grid[IX(x,y)] = grid[IX(i, j)];
			}*/
			
		}
		END_FOR
		FOR_EACH_CELL
			grid[IX(i,j)] = new_grid[IX(i, j)];
		END_FOR
		/*
		FOR_EACH_CELL
			if(grid[IX(i,j)]==0){
				if(grid[IX(i+1,j)] == (2 || 1)){
					if(grid[IX(i-1,j)] == (2 || 1)){
						grid[IX(i,j)]=2;
					}
				}

				if(grid[IX(i,j+1)] == (2||1) && grid[IX(i,j-1)] == 2){
					grid[IX(i,j)]=2;
				}

			}
			
		END_FOR*/
		//centerX = centerX +moveX;
	}
	prev_zz = zz;
}



// make the fluid mass conserving
void project ( int N, float * u, float * v, float * p, float * div, float * grid)
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
	set_bnd ( N, 0, div, grid ); set_bnd ( N, 0, p, grid );

	lin_solve ( N, 0, p, div, 1, 4, grid );

	FOR_EACH_CELL
		u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
		v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
	END_FOR
	set_bnd ( N, 1, u, grid ); set_bnd ( N, 2, v, grid);
}

void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt, float * grid, int xx, int yy,  int &centerX, int &centerY)
{
	// source density initially contained in the array x0
	add_source ( N, x, x0, dt ); // new density sources may be added by user interaction
	SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt, grid ); // diffuse 
	SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, dt, false, grid, xx, yy, centerX, centerY); // advect
}

void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt, int xx, int yy, float zz, float * grid, int &centerX, int &centerY, float* prev_grid)
{
	add(N, grid, centerX, centerY, zz, prev_grid);
	// add forces stored in the arrays u0 and v0 re-using the add_source function
	add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
	// 
	SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt, grid );
	SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt, grid );

	project ( N, u, v, u0, v0, grid); // conserve mass (1st time)

	SWAP ( u0, u ); SWAP ( v0, v ); // swap force arrays
	advect ( N, 1, u, u0, u0, v0, dt, true, grid, xx, yy, centerX, centerY); advect ( N, 2, v, v0, u0, v0, dt, true, grid, xx, yy,centerX, centerY); // more accurate advection because of project
	project ( N, u, v, u0, v0, grid); // conserve mass (2nd time)
}

