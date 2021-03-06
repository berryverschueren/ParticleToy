// project 2
#include "imageio.h"
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#include <iostream>
#include <string> 
#include "RigidBody.h"
#include "include/Eigen/Dense"
#include "include/Eigen/IterativeLinearSolvers"

using namespace Eigen;

#define IX(i,j) ((i)+(N+2)*(j))
extern void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, 
	float dt, float* grid, float* grid_prev, Vector2f &genForce, float eps);

extern void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, 
	float dt, float* grid, Vector2f &genForce);

extern void acc_step (int N, float * u, float * v, float* grid, float dt, Vector2f &genForce, Vector2f center, Vector3f &torq, int fluidSolid);

static int N, dvel, ovel, lvel, rvel, fluidSolid, solidFluid;
static float dt, diff, visc, force, source;
static float * u, * v, * u_prev, * v_prev; 
static float * dens, * dens_prev;
static RigidBody * rb = new RigidBody(Vector2f(0.5f, 0.5f));
static Vector2f genForce = Vector2f(0.0f, 0.0f);
static Vector3f torq = Vector3f(0.0f, 0.0f, 0.0f);

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;
static int dump_frames;
static int frame_number;

static float * grid, *grid_prev;
static float eps;



/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/


static void free_data ( void )
{
	if ( u ) free ( u );
	if ( v ) free ( v );
	if ( u_prev ) free ( u_prev );
	if ( v_prev ) free ( v_prev );
	if ( dens ) free ( dens );
	if ( dens_prev ) free ( dens_prev );
	if (grid) free(grid);
	if (grid_prev) free(grid_prev);

	rb->Reset();
}

static void clear_data ( void )
{
	int i, size=(N+2)*(N+2); // size of all arrays ( +2 to simplify bounding box handling)

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
		grid[i] = grid_prev[i] = 0.0f;
	}

	rb->Reset();
}

static int allocate_data ( void )
{
	int size = (N+2)*(N+2); // size of all arrays ( +2 to simplify bounding box handling)

	u			= (float *) malloc ( size*sizeof(float) );
	v			= (float *) malloc ( size*sizeof(float) );
	u_prev		= (float *) malloc ( size*sizeof(float) );
	v_prev		= (float *) malloc ( size*sizeof(float) );
	dens		= (float *) malloc ( size*sizeof(float) );	
	dens_prev	= (float *) malloc ( size*sizeof(float) );
	grid 		= (float *) malloc (size*sizeof(float));
	grid_prev 		= (float *) malloc (size*sizeof(float));

	rb->Reset();

	if ( !u || !v || !u_prev || !v_prev || !dens || !dens_prev || !grid || !grid_prev ) {
		fprintf ( stderr, "cannot allocate data\n" );
		return ( 0 );
	}

	return ( 1 );
}

static void body_step(RigidBody * rb, float dt) {
	rb->ResetForce();

	rb->_force = genForce;
	rb->_center = rb->_center + ((rb->_velocity) * dt);
	rb->_velocity = 0.9*rb->_velocity + ((rb->_force / rb->_mass) * dt);

	// 0    -az   ay
	// az   0     -ax
	// -ay  ax    0
	MatrixXf w_star = MatrixXf(3,3);
	w_star(0,0) = 0.0f;
	w_star(0,1) = -rb->_angularVelocity[2];
	w_star(0,2) = rb->_angularVelocity[1];
	w_star(1,0) = rb->_angularVelocity[2];
	w_star(1,1) = 0.0f;
	w_star(1,2) = -rb->_angularVelocity[0];
	w_star(2,0) = -rb->_angularVelocity[1];
	w_star(2,1) = rb->_angularVelocity[0];
	w_star(2,2) = 0.0f;

	auto test = w_star * rb->_orientation;
	auto newOrientation = test * dt;
	rb->_orientation += newOrientation;
	if(rvel){
		rb->_angularVelocity = Vector3f(torq[0] * dt, torq[1] * dt, torq[2] * dt);
	}else{
		rb->_angularVelocity = Vector3f(0.0f,0.0f,0.0f);
	}
}

static void VoxelizeRigidBody(RigidBody * rb, float * grid, float * grid_prev) {
	// convert coordinates of rb to grid cells
	int width = (int)(rb->_width*N);
	int height = (int)(rb->_height*N);
	int startX = (int)(rb->_center[0]*N) - (width / 2);
	int startY = (int)(rb->_center[1]*N) - (height / 2);

	// init grid
	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			grid_prev[IX(i,j)] = grid[IX(i, j)];
			grid[IX(i, j)] = 0.0f;
		}
	}

	// content cells
	for (float i = startX; i < (startX + width); i+=0.5f) {
		for (float j = startY; j < (startY + height); j+=0.5f) {

			//printf("Computing for: i=%g, j=%g\n", i, j);

			auto x = (i * 1.f/N) - rb->_center[0];
			auto y = (j * 1.f/N) - rb->_center[1];

			//printf("Without 0-1 grid: x=%g, y=%g\n", x, y);
			//printf("R[0]: x=%g, y=%g\n", rb->_orientation(0,0), rb->_orientation(0,1));
			//printf("R[1]: x=%g, y=%g\n", rb->_orientation(1,0), rb->_orientation(1,1));

			auto newPos = Vector2f(rb->_orientation(0,0)*x + rb->_orientation(0,1)*y + rb->_center[0],
								rb->_orientation(1,0)*x + rb->_orientation(1,1)*y + rb->_center[1]);

			//printf("Rotated to: x=%g, y=%g\n", newPos[0], newPos[1]);

			// convert newPos to cell
			auto newI = (int)(newPos[0]*N);
			auto newJ = (int)(newPos[1]*N);

			//printf("In cells: i=%d, j=%d\n", newI, newJ);

			grid[IX(newI, newJ)] = 2.0f;
			// grid[IX(i, j)] = 2.0f;
		}
	}

	// top border
	for (float i = 0; i <= width; i+=0.5f) {		
			int j = startY;			
			auto x = ((i+startX) * 1.f/N) - rb->_center[0];
			auto y = (j * 1.f/N) - rb->_center[1];
			auto newPos = Vector2f(rb->_orientation(0,0)*x + rb->_orientation(0,1)*y + rb->_center[0],
								rb->_orientation(1,0)*x + rb->_orientation(1,1)*y + rb->_center[1]);

			// convert newPos to cell
			auto newI = (int)(newPos[0]*N);
			auto newJ = (int)(newPos[1]*N);

			grid[IX(newI, newJ)] = 1.0f;
		// grid[IX(startX + i, startY)] = 1.0f;
	}
	
	// bottom border
	for (float i = 0; i <= width; i+=0.5f) {
			int j = startY + height;
			auto x = ((i+startX) * 1.f/N) - rb->_center[0];
			auto y = (j * 1.f/N) - rb->_center[1];
			auto newPos = Vector2f(rb->_orientation(0,0)*x + rb->_orientation(0,1)*y + rb->_center[0],
								rb->_orientation(1,0)*x + rb->_orientation(1,1)*y + rb->_center[1]);

			// convert newPos to cell
			auto newI = (int)(newPos[0]*N);
			auto newJ = (int)(newPos[1]*N);

			grid[IX(newI, newJ)] = 1.0f;
		// grid[IX(startX + i, startY + height)] = 1.0f;
	}
	
	// left border
	for (float j = 0; j <= height; j+=0.5f) {
			int i = startX;
			auto x = (i * 1.f/N) - rb->_center[0];
			auto y = ((j + startY) * 1.f/N) - rb->_center[1];
			auto newPos = Vector2f(rb->_orientation(0,0)*x + rb->_orientation(0,1)*y + rb->_center[0],
								rb->_orientation(1,0)*x + rb->_orientation(1,1)*y + rb->_center[1]);

			// convert newPos to cell
			auto newI = (int)(newPos[0]*N);
			auto newJ = (int)(newPos[1]*N);

			grid[IX(newI, newJ)] = 1.0f;
		// grid[IX(startX, startY + i)] = 1.0f;
	}
	
	// right border
	for (float j = 0; j <= height; j+=0.5f) {
			int i = startX + width;
			auto x = (i * 1.f/N) - rb->_center[0];
			auto y = ((j + startY) * 1.f/N) - rb->_center[1];
			auto newPos = Vector2f(rb->_orientation(0,0)*x + rb->_orientation(0,1)*y + rb->_center[0],
								rb->_orientation(1,0)*x + rb->_orientation(1,1)*y + rb->_center[1]);

			// convert newPos to cell
			auto newI = (int)(newPos[0]*N);
			auto newJ = (int)(newPos[1]*N);

			grid[IX(newI, newJ)] = 1.0f;
		// grid[IX(startX + width, startY + i)] = 1.0f;
	}
}


/*
  ----------------------------------------------------------------------
   OpenGL specific drawing routines
  ----------------------------------------------------------------------
*/

static void pre_display ( void )
{
	glViewport ( 0, 0, win_x, win_y );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ();
	gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
}

static void post_display ( void )
{
	// Write frames if necessary.
	if (dump_frames) {
		const int FRAME_INTERVAL = 4;
		if ((frame_number % FRAME_INTERVAL) == 0) {
			const unsigned int w = glutGet(GLUT_WINDOW_WIDTH);
			const unsigned int h = glutGet(GLUT_WINDOW_HEIGHT);
			unsigned char * buffer = (unsigned char *) malloc(w * h * 4 * sizeof(unsigned char));
			if (!buffer)
				exit(-1);
			// glRasterPos2i(0, 0);
			glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
			static char filename[80];
			sprintf(filename, "../snapshots/img%.5i.png", frame_number / FRAME_INTERVAL);
			printf("Dumped %s.\n", filename);
			saveImageRGBA(filename, buffer, w, h);
			
			free(buffer);
		}
	}
	frame_number++;
	glutSwapBuffers ();
}

static void draw_velocity ( void )
{
	int i, j;
	float x, y, h;

	h = 1.0f/N; // grid spacing (assume cells have physical size 1.0)

	glColor3f ( 1.0f, 1.0f, 1.0f );
	glLineWidth ( 1.0f );

	glBegin ( GL_LINES );

		for ( i=1 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=1 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;

				glVertex2f ( x, y );
				glVertex2f ( x+u[IX(i,j)], y+v[IX(i,j)] );
			}
		}

	glEnd ();

	if ( ovel ) {  	
		glBegin( GL_QUADS );

			for ( i=0 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=0 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;

				glColor3f(0.3f * grid[IX(i,j)], 0.2f * grid[IX(i,j)], 0.1f * grid[IX(i,j)]);
				glVertex2f ( x, y );
				glVertex2f ( x+h, y );
				glVertex2f ( x+h, y+h );
				glVertex2f ( x, y+h );
			}
		}

		glEnd();
	}
}

static void draw_density ( void )
{
	int i, j;
	float x, y, h, d00, d01, d10, d11;

	h = 1.0f/N; // grid spacing (assume cells have physical size 1.0)

	glBegin ( GL_QUADS );

		for ( i=0 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=0 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;

				d00 = dens[IX(i,j)];
				d01 = dens[IX(i,j+1)];
				d10 = dens[IX(i+1,j)];
				d11 = dens[IX(i+1,j+1)];

				glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
				glColor3f ( d10, d10, d10 ); glVertex2f ( x+h, y );
				glColor3f ( d11, d11, d11 ); glVertex2f ( x+h, y+h );
				glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+h );
			}
		}

	glEnd ();

	if ( ovel ) {  	
		glBegin( GL_QUADS );

			for ( i=0 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=0 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;

				glColor3f(0.3f * grid[IX(i,j)], 0.2f * grid[IX(i,j)], 0.1f * grid[IX(i,j)]);
				glVertex2f ( x, y );
				glVertex2f ( x+h, y );
				glVertex2f ( x+h, y+h );
				glVertex2f ( x, y+h );
			}
		}

		glEnd();
	}
}

/*
  ----------------------------------------------------------------------
   relates mouse movements to forces sources
  ----------------------------------------------------------------------
*/

int xx = 0;
int yy = 0;

static void get_from_UI ( float * d, float * u, float * v, float * grid )
{
	int i, j, size = (N+2)*(N+2);

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = d[i] = 0.0f;
	}
	genForce = Vector2f(0.0f,0.0f);

	if ( !mouse_down[0] && !mouse_down[2] ) return;
	
	i = (int)((       mx /(float)win_x)*N+1);
	j = (int)(((win_y-my)/(float)win_y)*N+1);

	if ( i<1 || i>N || j<1 || j>N ) return;

	if ( mouse_down[0] ) {		
		if ( lvel ) {
			auto mouseX = mx/(float)win_x;
			auto mouseY = (win_y-my)/(float)win_y;
			auto disX = mouseX - rb->_center[0];
			auto disY = mouseY - rb->_center[1];
			auto length = std::sqrt(std::pow(disX,2.0) + std::pow(disY,2.0));

			if(solidFluid){
				genForce[0] = disX/length*0.01f;
				genForce[1] = disY/length*0.01f;
			}else{
				genForce[0] = disX/length*0.01f *4;
				genForce[1] = disY/length*0.01f *4;
			}
			
		} else {
			u[IX(i,j)] = force * (mx-omx);
			v[IX(i,j)] = force * (omy-my);	
		}
	}

	if ( mouse_down[2] ) {
		d[IX(i,j)] = source;
	}

	omx = mx;
	omy = my;

	return;
}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/

static void key_func ( unsigned char key, int x, int y )
{
	switch ( key )
	{
		case 'c':
		case 'C':
			clear_data ();
			break;

		case 'd':
		case 'D':
		dump_frames = !dump_frames;
		break;

		case 'q':
		case 'Q':
			free_data ();
			exit ( 0 );
			break;

		case 'v':
		case 'V':
			dvel = !dvel;
			break;

		case 'b':
		case 'B':
			ovel = !ovel;
			break;

		case 'l':
		case 'L':
			lvel = !lvel;
			break;



		case 'i':
		case 'I':
			if(eps != 2.5){
				eps = 1.5f;		
			}else{
				eps = 0.1f;
			}
			break;

		case 'o':
		case 'O':
			//ob fluid
			solidFluid = !solidFluid;
			break;

		case 'p':
		case 'P':
			//fluid solid
			fluidSolid = !fluidSolid;
			break;

		case 'K':
		case 'k':
			//rot
			rvel = !rvel;
			break;

		case 'm':
		case 'M':
			//big
			if(rb->_width != 40.0f/N){
				rb->_width = 40.0f/N;
				rb->_mass = (rb->_width*64)*(rb->_height*64)/100;
			}else{
				rb->_width = 10.0f/N;
				rb->_mass = (rb->_width*64)*(rb->_height*64)/100;
			}
			
			break;
	}
}

static void mouse_func ( int button, int state, int x, int y )
{
	omx = mx = x;
	omx = my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func ( int x, int y )
{
	mx = x;
	my = y;
}

static void reshape_func ( int width, int height )
{
	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	win_x = width;
	win_y = height;
}

static void idle_func ( void )
{
	// get forces implied on body
	get_from_UI ( dens_prev, u_prev, v_prev, grid );
	// voxelize current rb, without moving it yet
	VoxelizeRigidBody(rb, grid, grid_prev);
	// accumulate forces on body
	acc_step(N, u, v, grid, dt, genForce, rb->_center, torq, fluidSolid); 
	// actually move the body using euler solver
	body_step(rb, dt);
	// use voxelized rb and its implied force in fluid solver
	vel_step ( N, u, v, u_prev, v_prev, visc, dt, grid, rb->_velocity);
	dens_step ( N, dens, dens_prev, u, v, diff, dt , grid, grid_prev, rb->_velocity, eps);

	glutSetWindow ( win_id );
	glutPostRedisplay ();
}

static void display_func ( void )
{
	pre_display ();

		if ( dvel ) { 
			draw_velocity ();
		}
		else {		
			draw_density ();
		}

	post_display ();
}


/*
  ----------------------------------------------------------------------
   open_glut_window --- open a glut compatible window and set callbacks
  ----------------------------------------------------------------------
*/

static void open_glut_window ( void )
{
	glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

	glutInitWindowPosition ( 0, 0 );
	glutInitWindowSize ( win_x, win_y );
	win_id = glutCreateWindow ( "Alias | wavefront" );

	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();

	pre_display ();

	glutKeyboardFunc ( key_func );
	glutMouseFunc ( mouse_func );
	glutMotionFunc ( motion_func );
	glutReshapeFunc ( reshape_func );
	glutIdleFunc ( idle_func );
	glutDisplayFunc ( display_func );
}


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{
	glutInit ( &argc, argv );
    if ( argc == 1 ) {
        N = 64;
        dt = 0.1f;
        diff = 0.0f;
        visc = 0.0f;
        force = 5.0f;
        source = 100.0f;
        eps = 0.1f;
        fprintf ( stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g epsilon=%g\n",
                  N, dt, diff, visc, force, source, eps );
    } else {
        N = atoi(argv[1]);
        dt = atof(argv[2]);
        diff = atof(argv[3]);
        visc = atof(argv[4]);
        force = atof(argv[5]);
        source = atof(argv[6]);
        eps = atof(argv[7]);
    }

	printf ( "\n\nHow to use this demo:\n\n" );
	printf ( "\t Add densities with the right mouse button\n" );
	printf ( "\t Add velocities with the left mouse button and dragging the mouse\n" );
	printf ( "\t Toggle object interaction with the 'l' key, then use the left mouse button and dragging the mouse\n" );
	printf ( "\t Toggle intensity of mouse interaction with the 'o' key\n" );
	printf ( "\t Toggle rigid body rotation with the 'k' key\n" );
	printf ( "\t Toggle density/velocity display with the 'v' key\n" );
	printf ( "\t Toggle object display with the 'b' key\n" );
	printf ( "\t Toggle vorticity intensity with the 'i' key\n" );
	printf ( "\t Toggle fluid-solid interaction with the 'p' key\n" );
	printf ( "\t Toggle large object with the 'm' key\n" );
	printf ( "\t Clear the simulation by pressing the 'c' key\n" );
	printf ( "\t Quit by pressing the 'q' key\n" );

	dvel = 0;
	ovel = 0;
	lvel = 0;
	rvel = 0;
	solidFluid = 1;
	fluidSolid = 0;
	dump_frames = 0;
	frame_number = 0;

	if ( !allocate_data () ) exit ( 1 );
	clear_data ();

	win_x = 512;
	win_y = 512;
	open_glut_window ();

	glutMainLoop ();

	exit ( 0 );
}