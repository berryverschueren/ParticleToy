// ParticleToy.cpp : Defines the entry point for the console application.
//

// #include "Particle.h"
// #include "Force.h"
// #include "GravityForce.h"
// #include "SpringForce.h"
// #include "AngularSpringForce.h"
// #include "RodConstraint.h"
// #include "CircularWireConstraint.h"
// #include "imageio.h"
// #include "ParticleSystem.h"
// #include "FluidSystem.h"

// #include <vector>
// #include <stdlib.h>
// #include <stdio.h>
// #include <GL/glut.h>
// #include <iostream>

// project 2
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#include <iostream>
#include <string> 

#define IX(i,j) ((i)+(N+2)*(j))
extern void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt );
extern void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt, int xx, int yy);
static int N, dvel;
static float dt, diff, visc, force, source;
static float * u, * v, * u_prev, * v_prev; 
static float * dens, * dens_prev;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;

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
}

static void clear_data ( void )
{
	int i, size=(N+2)*(N+2); // size of all arrays ( +2 to simplify bounding box handling)

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
	}
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

	if ( !u || !v || !u_prev || !v_prev || !dens || !dens_prev ) {
		fprintf ( stderr, "cannot allocate data\n" );
		return ( 0 );
	}

	return ( 1 );
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
}

/*
  ----------------------------------------------------------------------
   relates mouse movements to forces sources
  ----------------------------------------------------------------------
*/

int xx = 0;
int yy = 0;

static void get_from_UI ( float * d, float * u, float * v )
{
	int i, j, size = (N+2)*(N+2);

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = d[i] = 0.0f;
	}

	if ( !mouse_down[0] && !mouse_down[2] ) return;

	i = (int)((       mx /(float)win_x)*N+1);
	j = (int)(((win_y-my)/(float)win_y)*N+1);

	if ( i<1 || i>N || j<1 || j>N ) return;

	if ( mouse_down[0] ) {
		u[IX(i,j)] = force * (mx-omx);
		v[IX(i,j)] = force * (omy-my);

		xx = mx;
		yy = my;
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

		case 'q':
		case 'Q':
			free_data ();
			exit ( 0 );
			break;

		case 'v':
		case 'V':
			dvel = !dvel;
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
	get_from_UI ( dens_prev, u_prev, v_prev );
	vel_step ( N, u, v, u_prev, v_prev, visc, dt, xx, yy);
	dens_step ( N, dens, dens_prev, u, v, diff, dt );
	xx=0;
	yy=0;

	glutSetWindow ( win_id );
	glutPostRedisplay ();
}

static void display_func ( void )
{
	pre_display ();

		if ( dvel ) draw_velocity ();
		else		draw_density ();

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

	if ( argc != 1 && argc != 6 ) {
		fprintf ( stderr, "usage : %s N dt diff visc force source\n", argv[0] );
		fprintf ( stderr, "where:\n" );\
		fprintf ( stderr, "\t N      : grid resolution\n" );
		fprintf ( stderr, "\t dt     : time step\n" );
		fprintf ( stderr, "\t diff   : diffusion rate of the density\n" );
		fprintf ( stderr, "\t visc   : viscosity of the fluid\n" );
		fprintf ( stderr, "\t force  : scales the mouse movement that generate a force\n" );
		fprintf ( stderr, "\t source : amount of density that will be deposited\n" );
		exit ( 1 );
	}

	if ( argc == 1 ) {
		N = 64;
		dt = 0.1f;
		diff = 0.0f;
		visc = 0.0f;
		force = 5.0f;
		source = 100.0f;
		fprintf ( stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
			N, dt, diff, visc, force, source );
	} else {
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		diff = atof(argv[3]);
		visc = atof(argv[4]);
		force = atof(argv[5]);
		source = atof(argv[6]);
	}

	printf ( "\n\nHow to use this demo:\n\n" );
	printf ( "\t Add densities with the right mouse button\n" );
	printf ( "\t Add velocities with the left mouse button and dragging the mouse\n" );
	printf ( "\t Toggle density/velocity display with the 'v' key\n" );
	printf ( "\t Clear the simulation by pressing the 'c' key\n" );
	printf ( "\t Quit by pressing the 'q' key\n" );

	dvel = 0;

	if ( !allocate_data () ) exit ( 1 );
	clear_data ();

	win_x = 512;
	win_y = 512;
	open_glut_window ();

	glutMainLoop ();

	exit ( 0 );
}

//









// /* macros */

// /* external definitions (from solver) */
// extern void simulation_step( ParticleSystem *particleSystem, float dt, int solverVersion );

// /* global variables */

// static int N;
// static float dt, d;
// static int dsim;
// static int dump_frames;
// static int frame_number;
// static int solverVersion;
// static bool windForceFlag = false;
// static bool gravityForceFlag = false;
// static bool mouseFlag = false;
// static bool mouseDown = false;

// ParticleSystem *particleSystem = new ParticleSystem();
// FluidSystem *fluidSystem = new FluidSystem(64, 0.1f, 0.0f, 0.0f, 5.0f, 100.0f, 0);

// static int win_id;
// static int win_x, win_y;
// static int mouse_down[3];
// static int mouse_release[3];
// static int mouse_shiftclick[3];
// static int omx, omy, mx, my;
// static int hmx, hmy;

// /*
// ----------------------------------------------------------------------
// free/clear/allocate simulation data
// ----------------------------------------------------------------------
// */

// static void free_data ( void )
// {
// 	particleSystem->deleteAll();
// }

// static void clear_data ( void )
// {
// 	int ii, size = particleSystem->getParticles().size();

// 	for(ii=0; ii<size; ii++){
// 	    particleSystem->getParticles()[ii]->reset();
// 	}
// }

// void static addGravityForce(){
//     if(!gravityForceFlag){
//         std::cout << "Added gravity forces"<< '\n';
//         const Vec2f g(0.0,-0.0015f);
//         auto gravity = new GravityForce(particleSystem->getParticles(), g);
//         particleSystem->addForce(gravity);
//         gravityForceFlag = true;
//     } else{
//         std::cout << "Removed gravity forces"<< '\n';
//         particleSystem->removeLastForce();
//         gravityForceFlag = false;
//     }
// }

// static void createHorizontalContraintCloth(){
//     const Vec2f startingPoint(-0.5f, -0.5f);
//     int ii, jj, maxRow = 10, maxCol = 10;
//     for (ii=0; ii<maxRow; ii++) {
//         for (jj=0; jj<maxCol; jj++) {
//             auto particle = new Particle(startingPoint + Vec2f(1.0f/maxRow*ii, 1.0f/maxCol*jj));
//             particleSystem->addParticle(particle);
//         }
//     }


//     double distance = 0.1, springConstant = 0.1, dampingConstant = 0.5;

//     // springforce particle with particle beneath it
//     for(ii=0; ii<maxRow-1; ii++){
//         for(jj=0; jj<maxCol; jj++){
//             auto spring = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii],
//                                           particleSystem->getParticles()[jj*maxRow+ii+1], distance, springConstant, dampingConstant);
//             particleSystem->addForce(spring);
//         }
//     }
//     //addGravityForce();
//     // springforce particle with particle to the right of it
//     for(jj=0; jj<maxCol-1; jj++){
//         auto rod1 = new RodConstraint(particleSystem->getParticles()[jj*maxRow+0],
//                                       particleSystem->getParticles()[(jj+1)*maxRow+0], 0.1);
//         auto rod2 = new RodConstraint(particleSystem->getParticles()[jj*maxRow+(maxRow-1)],
//                                       particleSystem->getParticles()[(jj+1)*maxRow+(maxRow-1)], 0.1);
//         particleSystem->addConstraint(rod1);
//         particleSystem->addConstraint(rod2);
//     }
//     for(ii=1; ii<maxRow-1; ii++){
//         for(jj=0; jj<maxCol-1; jj++){
//             auto spring = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii],
//                                           particleSystem->getParticles()[(jj+1)*maxRow+ii], distance, springConstant, dampingConstant);
//             particleSystem->addForce(spring);
//         }
//     }

//     // springforce particle with particle to the right and beneath it
//     // springforce particle with particle to the left and beneath it
//     for(ii=0; ii<maxRow-1; ii++){
//         for(jj=0; jj<maxCol-1; jj++){
//             auto spring1 = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii],
//                                            particleSystem->getParticles()[(jj+1)*maxRow+ii+1], distance, springConstant, dampingConstant);
//             auto spring2 = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii+1],
//                                            particleSystem->getParticles()[(jj+1)*maxRow+ii], distance, springConstant, dampingConstant);
//             particleSystem->addForce(spring1);
//             particleSystem->addForce(spring2);
//         }
//     }



//     // springs to hold the cloth in place (and flip it)
//     Particle *p1 = new Particle(Vec2f(-0.9f, 0.9f));
//     Particle *p2 = new Particle(Vec2f(0.9f, 0.9f));
//     SpringForce *sf1 = new SpringForce(p1, particleSystem->getParticles()[maxCol-1], distance, springConstant, dampingConstant);
//     SpringForce *sf2 = new SpringForce(p2, particleSystem->getParticles()[maxRow*(maxCol)-1], distance, springConstant, dampingConstant);
//     particleSystem->addParticle(p1);
//     particleSystem->addParticle(p2);
//     particleSystem->addForce(sf1);
//     particleSystem->addForce(sf2);



//     float circDistance = 0.05f;//sqrt(pow(1.0f/maxRow,2) + pow(1.0f/maxCol,2));
//     const Vec2f allowedOffset(circDistance, 0.0);
//     auto c1 = new CircularWireConstraint(p1, p1->m_ConstructPos + allowedOffset, circDistance);
//     auto c2 = new CircularWireConstraint(p2, p2->m_ConstructPos + allowedOffset, circDistance);
//     particleSystem->addConstraint(c1);
//     particleSystem->addConstraint(c2);
// }

// static void createSpringCloth(){
//     const Vec2f startingPoint(-0.5f, -0.5f);
//     int ii, jj, maxRow = 10, maxCol = 10;
//     for (ii=0; ii<maxRow; ii++) {
//         for (jj=0; jj<maxCol; jj++) {
//             auto particle = new Particle(startingPoint + Vec2f(1.0f/maxRow*ii, 1.0f/maxCol*jj));
//             particleSystem->addParticle(particle);
//         }
//     }


//     double distance = 0.1, springConstant = 0.1, dampingConstant = 0.5;
//     double distanceSpring = 0.25;
//     // springforce particle with particle beneath it
//     for(ii=0; ii<maxRow-1; ii++){
//         for(jj=0; jj<maxCol; jj++){
//             auto spring = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii],
//                                           particleSystem->getParticles()[jj*maxRow+ii+1], distance, springConstant, dampingConstant);
//             particleSystem->addForce(spring);
//         }
//     }
//     //addGravityForce();
//     // springforce particle with particle to the right of it
//     for(ii=0; ii<maxRow; ii++){
//         for(jj=0; jj<maxCol-1; jj++){
//             auto spring = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii],
//                                           particleSystem->getParticles()[(jj+1)*maxRow+ii], distance, springConstant, dampingConstant);
//             particleSystem->addForce(spring);
//         }
//     }

//     // springforce particle with particle to the right and beneath it
//     // springforce particle with particle to the left and beneath it
//     for(ii=0; ii<maxRow-1; ii++){
//         for(jj=0; jj<maxCol-1; jj++){
//             auto spring1 = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii],
//                                            particleSystem->getParticles()[(jj+1)*maxRow+ii+1], distance, springConstant, dampingConstant);
//             auto spring2 = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii+1],
//                                            particleSystem->getParticles()[(jj+1)*maxRow+ii], distance, springConstant, dampingConstant);
//             particleSystem->addForce(spring1);
//             particleSystem->addForce(spring2);
//         }
//     }



//     // springs to hold the cloth in place (and flip it)
//     Particle *p1 = new Particle(Vec2f(-0.9f, 0.9f));
//     Particle *p2 = new Particle(Vec2f(0.9f, 0.9f));
//     SpringForce *sf1 = new SpringForce(p1, particleSystem->getParticles()[maxCol-1], distance, springConstant, dampingConstant);
//     SpringForce *sf2 = new SpringForce(p2, particleSystem->getParticles()[maxRow*(maxCol)-1], distance, springConstant, dampingConstant);
//     SpringForce *sf3 = new SpringForce(p1, particleSystem->getParticles()[0], distanceSpring, springConstant, dampingConstant);
//     SpringForce *sf4 = new SpringForce(p2, particleSystem->getParticles()[maxRow*(maxCol-1)], distanceSpring, springConstant, dampingConstant);
//     particleSystem->addParticle(p1);
//     particleSystem->addParticle(p2);
//     particleSystem->addForce(sf1);
//     particleSystem->addForce(sf2);
//     particleSystem->addForce(sf4);
//     particleSystem->addForce(sf3);

//     float circDistance = 0.05f;//sqrt(pow(1.0f/maxRow,2) + pow(1.0f/maxCol,2));
//     const Vec2f allowedOffset(circDistance, 0.0);
//     auto c1 = new CircularWireConstraint(p1, p1->m_ConstructPos + allowedOffset, circDistance);
//     auto c2 = new CircularWireConstraint(p2, p2->m_ConstructPos + allowedOffset, circDistance);
//     particleSystem->addConstraint(c1);
//     particleSystem->addConstraint(c2);
// }

// static void createCloth() {
// 	const Vec2f startingPoint(-0.5f, -0.5f);
// 	int ii, jj, maxRow = 10, maxCol = 10;
// 	for (ii=0; ii<maxRow; ii++) {
// 		for (jj=0; jj<maxCol; jj++) {
// 			auto particle = new Particle(startingPoint + Vec2f(1.0f/maxRow*ii, 1.0f/maxCol*jj));
// 			particleSystem->addParticle(particle);
// 		}
// 	}


// 	double distance = 0.1, springConstant = 0.1, dampingConstant = 0.5;
//     double distanceSpring = 0.25;
// 	// springforce particle with particle beneath it
//     for(ii=0; ii<maxRow-1; ii++){
//         auto rod1 = new RodConstraint(particleSystem->getParticles()[0*maxRow+ii],
//                                         particleSystem->getParticles()[0*maxRow+ii+1], 0.1);
//         auto rod2 = new RodConstraint(particleSystem->getParticles()[(maxCol-1)*maxRow+ii],
//                                         particleSystem->getParticles()[(maxCol-1)*maxRow+ii+1], 0.1);
//         particleSystem->addConstraint(rod1);
//         particleSystem->addConstraint(rod2);
//     }
//     for(ii=0; ii<maxRow-1; ii++){
//         for(jj=1; jj<maxCol-1; jj++){
//             auto spring = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii],
// 				particleSystem->getParticles()[jj*maxRow+ii+1], distance, springConstant, dampingConstant);
//             particleSystem->addForce(spring);
//         }
//     }
//     //addGravityForce();
// 	// springforce particle with particle to the right of it
//     for(jj=0; jj<maxCol-1; jj++){
//         auto rod1 = new RodConstraint(particleSystem->getParticles()[jj*maxRow+0],
//                                       particleSystem->getParticles()[(jj+1)*maxRow+0], 0.1);
//         auto rod2 = new RodConstraint(particleSystem->getParticles()[jj*maxRow+(maxRow-1)],
//                                       particleSystem->getParticles()[(jj+1)*maxRow+(maxRow-1)], 0.1);
//         particleSystem->addConstraint(rod1);
//         particleSystem->addConstraint(rod2);
//     }
//     for(ii=1; ii<maxRow-1; ii++){
//         for(jj=0; jj<maxCol-1; jj++){
//             auto spring = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii],
// 				particleSystem->getParticles()[(jj+1)*maxRow+ii], distance, springConstant, dampingConstant);
//             particleSystem->addForce(spring);
//         }
//     }

// 	// springforce particle with particle to the right and beneath it
// 	// springforce particle with particle to the left and beneath it
//     for(ii=0; ii<maxRow-1; ii++){
//         for(jj=0; jj<maxCol-1; jj++){
// 			auto spring1 = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii],
// 				particleSystem->getParticles()[(jj+1)*maxRow+ii+1], distance, springConstant, dampingConstant);
// 			auto spring2 = new SpringForce(particleSystem->getParticles()[jj*maxRow+ii+1],
// 				particleSystem->getParticles()[(jj+1)*maxRow+ii], distance, springConstant, dampingConstant);
//             particleSystem->addForce(spring1);
//             particleSystem->addForce(spring2);
//         }
//     }



// 	// springs to hold the cloth in place (and flip it)
// 	Particle *p1 = new Particle(Vec2f(-0.9f, 0.9f));
// 	Particle *p2 = new Particle(Vec2f(0.9f, 0.9f));
// 	SpringForce *sf1 = new SpringForce(p1, particleSystem->getParticles()[maxCol-1], distance, springConstant, dampingConstant);
// 	SpringForce *sf2 = new SpringForce(p2, particleSystem->getParticles()[maxRow*(maxCol)-1], distance, springConstant, dampingConstant);
// 	SpringForce *sf3 = new SpringForce(p1, particleSystem->getParticles()[0], distanceSpring, springConstant*2, dampingConstant);
// 	SpringForce *sf4 = new SpringForce(p2, particleSystem->getParticles()[maxRow*(maxCol-1)], distanceSpring, springConstant*2, dampingConstant);
// 	particleSystem->addParticle(p1);
// 	particleSystem->addParticle(p2);
// 	particleSystem->addForce(sf1);
// 	particleSystem->addForce(sf2);
// 	particleSystem->addForce(sf4);
// 	particleSystem->addForce(sf3);

//     float circDistance = 0.05f;//sqrt(pow(1.0f/maxRow,2) + pow(1.0f/maxCol,2));
// 	const Vec2f allowedOffset(circDistance, 0.0);
// 	auto c1 = new CircularWireConstraint(p1, p1->m_ConstructPos + allowedOffset, circDistance);
// 	auto c2 = new CircularWireConstraint(p2, p2->m_ConstructPos + allowedOffset, circDistance);
// 	particleSystem->addConstraint(c1);
// 	particleSystem->addConstraint(c2);

// }

// static void init_system(void)
// {
// 	// solverVersion = 3;

// 	// const double dist = 0.2;
// 	// const Vec2f center(0.0, 0.0);
// 	// const Vec2f offset(dist, 0.0);

// 	// createCloth();

//     fluidSystem = new FluidSystem(64, 0.1f, 0.0f, 0.0f, 5.0f, 100.0f, 0); 
//     //init(64, 0.1f, 0.0f, 0.0f, 5.0f, 100.0f, 0);
//     fluidSystem->allocate_data();
//     fluidSystem->clear_data();
// }



// void static addHorizontalForce(){
//     if(!windForceFlag){
//         std::cout << "Added horizontal forces"<< '\n';
//         auto gravity = new GravityForce(particleSystem->getParticles(), Vec2f(-0.002f,0.0));
//         particleSystem->addForce(gravity);
//         windForceFlag = true;
//     } else{
//         std::cout << "Removed horizontal forces"<< '\n';
//         particleSystem->removeLastForce();
//         windForceFlag = false;
//     }

// }

// /*
// ----------------------------------------------------------------------
// OpenGL specific drawing routines
// ----------------------------------------------------------------------
// */

// static void pre_display ( void )
// {
// 	// glViewport ( 0, 0, win_x, win_y );
// 	// glMatrixMode ( GL_PROJECTION );
// 	// glLoadIdentity ();
// 	// gluOrtho2D ( -1.0, 1.0, -1.0, 1.0 );
// 	// glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
// 	// glClear ( GL_COLOR_BUFFER_BIT );
// 	glViewport ( 0, 0, win_x, win_y );
// 	glMatrixMode ( GL_PROJECTION );
// 	glLoadIdentity ();
// 	gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
// 	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
// 	glClear ( GL_COLOR_BUFFER_BIT );
// }

// static void post_display ( void )
// {
// 	// // Write frames if necessary.
// 	// if (dump_frames) {
// 	// 	const int FRAME_INTERVAL = 4;
// 	// 	if ((frame_number % FRAME_INTERVAL) == 0) {
// 	// 		const unsigned int w = glutGet(GLUT_WINDOW_WIDTH);
// 	// 		const unsigned int h = glutGet(GLUT_WINDOW_HEIGHT);
// 	// 		unsigned char * buffer = (unsigned char *) malloc(w * h * 4 * sizeof(unsigned char));
// 	// 		if (!buffer)
// 	// 			exit(-1);
// 	// 		// glRasterPos2i(0, 0);
// 	// 		glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
// 	// 		static char filename[80];
// 	// 		sprintf(filename, "../snapshots/img%.5i.png", frame_number / FRAME_INTERVAL);
// 	// 		printf("Dumped %s.\n", filename);
// 	// 		saveImageRGBA(filename, buffer, w, h);
			
// 	// 		free(buffer);
// 	// 	}
// 	// }
// 	// frame_number++;
	
// 	// glutSwapBuffers ();
//     glutSwapBuffers ();
// }

// static void draw_particles ( void )
// {
// 	int size = particleSystem->getParticles().size();

// 	for(int ii=0; ii< size; ii++)
// 	{
// 		particleSystem->getParticles()[ii]->draw();
// 	}
// }

// static void draw_forces ( void )
// {
//     int ii, size = particleSystem->getForces().size();

//     for (ii = 0; ii < size; ii++){
//         particleSystem->getForces()[ii]->draw();
//     }
// }

// static void draw_constraints ( void )
// {
//     int ii, size = particleSystem->getConstraints().size();

//     for (ii = 0; ii < size; ii++) {
//         particleSystem->getConstraints()[ii]->draw();
//     }
// }

// /*
// ----------------------------------------------------------------------
// relates mouse movements to particle toy construction
// ----------------------------------------------------------------------
// */

// // static void get_from_UI ()
// // {
// // 	int i, j;
// // 	// int size, flag;
// // 	int hi, hj;
// // 	// float x, y;
// // 	if ( !mouse_down[0] && !mouse_down[2] && !mouse_release[0] 
// // 	&& !mouse_shiftclick[0] && !mouse_shiftclick[2] ) return;

// // 	i = (int)((       mx /(float)win_x)*N);
// // 	j = (int)(((win_y-my)/(float)win_y)*N);

// // 	if ( i<1 || i>N || j<1 || j>N ) return;

// // 	//left mouse click
// //     if ( mouse_down[0] ) {

// //     }
// //     //right mouse click
// //     if ( mouse_down[2] ) {

// //     }

// //     //start pos
// // 	hi = (int)((       hmx /(float)win_x)*N);
// // 	hj = (int)(((win_y-hmy)/(float)win_y)*N);

// // 	if( mouse_release[0] ) {


// // 	    if(mouseFlag) {
// //             const Vec2f startPos((hi / (N / 2.0)) - 1.0, ((hj / (N / 2.0)) - 1.0));

// //             const Vec2f currentPos((mx / (win_x / 2.0)) - 1.0, -((my / (win_y / 2.0)) - 1.0));

// //             int size = particleSystem->getParticles().size();
// //             Particle *startClosestParticle = NULL;
// //             Particle *currentClosestParticle = NULL;
// //             float sDist = 1000000.0;
// //             float cDist = 1000000.0;
// //             for (int k = 0; k < size; ++k) {
// //                 Vec2f sV = startPos - particleSystem->getParticles()[k]->m_Position;
// //                 Vec2f cV = currentPos - particleSystem->getParticles()[k]->m_Position;

// //                 float sL = std::sqrt(std::abs(std::pow(sV[0], 2)) + std::abs(std::pow(sV[1], 2)));
// //                 float cL = std::sqrt(std::abs(std::pow(cV[0], 2)) + std::abs(std::pow(cV[1], 2)));
// //                 if (sL < sDist) {
// //                     sDist = sL;
// //                     startClosestParticle = particleSystem->getParticles()[k];
// //                 }
// //                 if (cL < cDist) {
// //                     cDist = cL;
// //                     currentClosestParticle = particleSystem->getParticles()[k];
// //                 }
// //             }
// //             if (sDist > 0.1) {
// //                 startClosestParticle = new Particle(startPos);
// //                 particleSystem->addParticle(startClosestParticle);
// //             }
// //             if (cDist > 0.1) {
// //                 currentClosestParticle = new Particle(currentPos);
// //                 particleSystem->addParticle(currentClosestParticle);
// //             }

// //             //Particle* p = new Particle(startPos);
// //             particleSystem->addForce(new SpringForce(startClosestParticle, currentClosestParticle, 1, 0.05, 0.5));
// //             mouse_release[0] = 0;
// //         }

// // 	}

// // 	omx = mx;
// // 	omy = my;
// // }
// static void get_from_UI()
// {
// 	int i, j, size = (fluidSystem->N+2)*(fluidSystem->N+2);

// 	for ( i=0 ; i<size ; i++ ) {
// 		fluidSystem->uForcesPrev[i] = fluidSystem->vForcesPrev[i] = fluidSystem->densitiesPrev[i] = 0.0f;
// 	}

// 	if ( !mouse_down[0] && !mouse_down[2] ) return;

// 	i = (int)((       mx /(float)win_x)*fluidSystem->N+1);
// 	j = (int)(((win_y-my)/(float)win_y)*fluidSystem->N+1);

// 	if ( i<1 || i>fluidSystem->N || j<1 || j>fluidSystem->N ) return;

// 	if ( mouse_down[0] ) {
// 		fluidSystem->uForcesPrev[((i)+(fluidSystem->N+2)*(j))] = fluidSystem->force * (mx-omx);
// 		fluidSystem->vForcesPrev[((i)+(fluidSystem->N+2)*(j))] = fluidSystem->force * (omy-my);
// 	}

// 	if ( mouse_down[2] ) {
//         //printf("\t md2, %d || %d || %f \n", ((i)+(fluidSystem->N+2)*(j)), fluidSystem->densitiesPrev.size(), fluidSystem->source);
// 		fluidSystem->densitiesPrev[((i)+(fluidSystem->N+2)*(j))] = fluidSystem->source;
// 	}

// 	omx = mx;
// 	omy = my;

// 	return;
// }

// static void remap_GUI()
// {
// 	int ii, size = particleSystem->getParticles().size();
// 	for(ii=0; ii<size; ii++)
// 	{
// 		particleSystem->getParticles()[ii]->reset();
// 		//pVector[ii]->m_Position[0] = pVector[ii]->m_ConstructPos[0];
// 		//pVector[ii]->m_Position[1] = pVector[ii]->m_ConstructPos[1];
// 	}
// }

// /*
// ----------------------------------------------------------------------
// GLUT callback routines
// ----------------------------------------------------------------------
// */

// // static void key_func ( unsigned char key, int x, int y )
// // {
// // 	switch ( key )
// // 	{
// // 	case 'c':
// // 	case 'C':
// // 		clear_data ();
// // 		break;

// // 	case 'd':
// // 	case 'D':
// // 		dump_frames = !dump_frames;
// // 		break;

// // 	case 'q':
// // 	case 'Q':
// // 		free_data ();
// // 		exit ( 0 );
// // 		break;

// // 	case 'm':
// //     case 'M':
// //         if(mouseFlag){
// //             mouseFlag = false;
// //         }else{
// //             mouseFlag = true;
// //         }
// //         break;
// // 	case '0':
// // 		solverVersion = 0;
// // 		break;
// // 	case '1':
// // 		solverVersion = 1;
// // 		break;
// // 	case '2':
// // 		solverVersion = 2;
// // 		break;
// // 	case '3':
// // 		solverVersion = 3;
// // 		break;

// // 	case 'g':
// //     case 'G':
// //         addGravityForce();
// //         break;
// //     case 'h':
// //     case 'H':
// //         addHorizontalForce();
// //         break;
// //     case 'j':
// //     case 'J':
// //         particleSystem->deleteAll();

// //         gravityForceFlag = false;
// //         windForceFlag = false;
// //         createSpringCloth();
// //         break;
// //     case 'k':
// //     case 'K':
// //         particleSystem->deleteAll();
// //         gravityForceFlag = false;
// //         windForceFlag = false;
// //         createCloth();
// //         break;
// //     case 'l':
// //     case 'L':
// //         particleSystem->deleteAll();
// //         gravityForceFlag = false;
// //         windForceFlag = false;
// //         createHorizontalContraintCloth();
// //         break;
// // 	case ' ':
// // 		dsim = !dsim;
// // 		break;
// // 	}
// // 	printf("\n Solverversion set to: %u \n", solverVersion);
// // }

// // static void mouse_func ( int button, int state, int x, int y )
// // {
// // 	omx = mx = x;
// // 	omx = my = y;

// // 	if(!mouse_down[0]){
// // 	    hmx=x; hmy=y;
// //         mouseDown = true;
// // 	}
// // 	if(mouse_down[button]) {
// //         mouseDown = false;
// //         mouse_release[button] = state == GLUT_UP;
// // 	}
// // 	if(mouse_down[button]) {
// // 	    mouse_shiftclick[button] = glutGetModifiers()==GLUT_ACTIVE_SHIFT;

// // 	}
// // 	mouse_down[button] = state == GLUT_DOWN;
// // }

// // static void motion_func ( int x, int y )
// // {
// // 	mx = x;
// // 	my = y;
// // }

// // static void reshape_func ( int width, int height )
// // {
// // 	glutSetWindow ( win_id );
// // 	glutReshapeWindow ( width, height );

// // 	win_x = width;
// // 	win_y = height;
// // }

// // static void idle_func ( void )
// // {
// //     if(mouseDown) {
// //         const Vec2f currentPos((mx / (win_x / 2.0)) - 1.0, -((my / (win_y / 2.0)) - 1.0));

// //         //closest particle
// //         int size = particleSystem->getParticles().size();
// //         Particle *closestParticle = NULL;
// //         float dist = 1000000.0;
// //         for (int k = 0; k < size; ++k) {
// //             Vec2f vec = currentPos - particleSystem->getParticles()[k]->m_Position;

// //             float length = std::sqrt(std::abs(std::pow(vec[0], 2)) + std::abs(std::pow(vec[1], 2)));
// //             if (length < dist) {
// //                 dist = length;
// //                 closestParticle = particleSystem->getParticles()[k];
// //             }
// //         }
// //         Particle *p = new Particle(currentPos);
// //         p->m_ConstructPos = currentPos;
// //         particleSystem->addParticle(p);

// //         particleSystem->addForce(new SpringForce(closestParticle, p, 0.5, 0.4, 0.5));
// //     }


// // 	if ( dsim ) simulation_step( particleSystem, dt, solverVersion );
// // 	else        {get_from_UI();remap_GUI();}

// // 	glutSetWindow ( win_id );
// // 	glutPostRedisplay ();

// // 	if(mouseDown) {
// //         particleSystem->removeLastForce();
// //         particleSystem->removeLastParticle();
// //     }


// // }

// // static void display_func ( void )
// // {
// // 	pre_display ();

// // 	draw_forces();
// // 	draw_constraints();
// // 	draw_particles();

// // 	post_display ();
// // }

// static void key_func ( unsigned char key, int x, int y )
// {
// 	switch ( key )
// 	{
// 		case 'c':
// 		case 'C':
// 			clear_data ();
// 			break;

// 		case 'q':
// 		case 'Q':
// 			free_data ();
// 			exit ( 0 );
// 			break;

// 		case 'v':
// 		case 'V':
// 			fluidSystem->dvel = !fluidSystem->dvel;
// 			break;
// 	}
// }

// static void mouse_func ( int button, int state, int x, int y )
// {
// 	omx = mx = x;
// 	omx = my = y;

// 	mouse_down[button] = state == GLUT_DOWN;
// }

// static void motion_func ( int x, int y )
// {
// 	mx = x;
// 	my = y;
// }

// static void reshape_func ( int width, int height )
// {
// 	glutSetWindow ( win_id );
// 	glutReshapeWindow ( width, height );

// 	win_x = width;
// 	win_y = height;
// }

// static void idle_func ( void )
// {
// 	get_from_UI();
//     fluidSystem->velocity_step();
//     fluidSystem->density_step();

// 	glutSetWindow ( win_id );
// 	glutPostRedisplay ();
// }

// static void display_func ( void )
// {
// 	pre_display ();

// 		if ( fluidSystem->dvel ) fluidSystem->draw_velocity();
// 		else		fluidSystem->draw_density();

// 	post_display ();
// }


// /*
// ----------------------------------------------------------------------
// open_glut_window --- open a glut compatible window and set callbacks
// ----------------------------------------------------------------------
// */

// static void open_glut_window ( void )
// {
// 	glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

// 	glutInitWindowPosition ( 0, 0 );
// 	glutInitWindowSize ( win_x, win_y );
// 	win_id = glutCreateWindow ( "Particletoys!" );

// 	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
// 	glClear ( GL_COLOR_BUFFER_BIT );
// 	glutSwapBuffers ();
// 	glClear ( GL_COLOR_BUFFER_BIT );
// 	glutSwapBuffers ();

// 	glEnable(GL_LINE_SMOOTH);
// 	glEnable(GL_POLYGON_SMOOTH);

// 	pre_display ();

// 	glutKeyboardFunc ( key_func );
// 	glutMouseFunc ( mouse_func );
// 	glutMotionFunc ( motion_func );
// 	glutReshapeFunc ( reshape_func );
// 	glutIdleFunc ( idle_func );
// 	glutDisplayFunc ( display_func );
// }


// /*
// ----------------------------------------------------------------------
// main --- main routine
// ----------------------------------------------------------------------
// */

// int main ( int argc, char ** argv )
// {
// 	glutInit ( &argc, argv );

// 	if ( argc == 1 ) {
// 		N = 64;
// 		dt = 0.05f;//0.1f
// 		d = 5.f;
// 		fprintf ( stderr, "Using defaults : N=%d dt=%g d=%g\n",
// 			N, dt, d );
// 	} else {
// 		N = atoi(argv[1]);
// 		dt = atof(argv[2]);
// 		d = atof(argv[3]);
// 	}

// 	printf ( "\n\nHow to use this application:\n\n" );
// 	printf ( "\t Toggle construction/simulation display with the spacebar key\n" );
// 	printf ( "\t Dump frames by pressing the 'd' key\n" );
// 	printf ( "\t Quit by pressing the 'q' key\n" );
// 	printf ( "\t Set solver version 0(skeleton),1(Euler),2(Mid-point),3(Runge-Kutta 4) by pressing the matching number key\n" );
//     printf ( "\t Toggle 'h' for horizontal wind force or 'g' vertical gravity force (global)\n" );
//     printf ( "\t set scene with 'j' full spring cloth\n" );
//     printf ( "\t set scene with 'k' cloth with fully bordered rod constraints\n" );
//     printf ( "\t set scene with 'l' cloth with horizontal bordered rod constraints\n" );
//     printf ( "\t Toggle 'm' for mouse interaction -> spring construction\n" );

// 	dsim = 0;
// 	dump_frames = 0;
// 	frame_number = 0;
	
// 	init_system();
// 	win_x = 512;
// 	win_y = 512;
// 	open_glut_window();

// 	glutMainLoop();

// 	exit(0);
// }

