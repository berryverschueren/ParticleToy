// ParticleToy.cpp : Defines the entry point for the console application.
//

#include "Particle.h"
#include "Force.h"
#include "GravityForce.h"
#include "SpringForce.h"
#include "RodConstraint.h"
#include "CircularWireConstraint.h"
#include "imageio.h"

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#include <iostream>

/* macros */

/* external definitions (from solver) */
extern void simulation_step( std::vector<Particle*> pVector, std::vector<Force*> fVector, float dt );

/* global variables */

static int N;
static float dt, d;
static int dsim;
static int dump_frames;
static int frame_number;

// static Particle *pList;
static std::vector<Particle*> pVector;
// static Force *fList;
static std::vector<Force*> fVector;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int mouse_release[3];
static int mouse_shiftclick[3];
static int omx, omy, mx, my;
static int hmx, hmy;

//static SpringForce * delete_this_dummy_spring = NULL;
static RodConstraint * delete_this_dummy_rod = NULL;
static CircularWireConstraint * delete_this_dummy_wire = NULL;


/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
*/

static void free_data ( void )
{
	pVector.clear();
	if (delete_this_dummy_rod) {
		delete delete_this_dummy_rod;
		delete_this_dummy_rod = NULL;
	}
	fVector.clear();
	//if (delete_this_dummy_spring) {
	//	delete delete_this_dummy_spring;
	//	delete_this_dummy_spring = NULL;
	//}
	if (delete_this_dummy_wire) {
		delete delete_this_dummy_wire;
		delete_this_dummy_wire = NULL;
	}
}

static void clear_data ( void )
{
	int ii, size = pVector.size();

	for(ii=0; ii<size; ii++){
		pVector[ii]->reset();
	}
}

static void init_system(void)
{
	const double dist = 0.2;
	const Vec2f center(0.0, 0.0);
	const Vec2f offset(dist, 0.0);

	// Create three particles, attach them to each other, then add a
	// circular wire constraint to the first.

	pVector.push_back(new Particle(center + offset));
	pVector.push_back(new Particle(center + offset + offset));
	pVector.push_back(new Particle(center + offset + offset + offset));
	
	//fVector.push_back(new GravityForce(pVector));
	fVector.push_back(new SpringForce(pVector[0], pVector[1], dist+dist, 0.05, 0.5));

	// You shoud replace these with a vector generalized forces and one of
	// constraints...
	delete_this_dummy_rod = new RodConstraint(pVector[1], pVector[2], dist);
	delete_this_dummy_wire = new CircularWireConstraint(pVector[0], center, dist);
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
	gluOrtho2D ( -1.0, 1.0, -1.0, 1.0 );
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

static void draw_particles ( void )
{
	int size = pVector.size();

	for(int ii=0; ii< size; ii++)
	{
		pVector[ii]->draw();
	}
}

static void draw_forces ( void )
{
    int ii, size = fVector.size();

    for (ii = 0; ii < size; ii++){
        fVector[ii]->draw();
    }
}

static void draw_constraints ( void )
{
	// change this to iteration over full set
	if (delete_this_dummy_rod)
		delete_this_dummy_rod->draw();
	if (delete_this_dummy_wire)
		delete_this_dummy_wire->draw();
}

/*
----------------------------------------------------------------------
relates mouse movements to particle toy construction
----------------------------------------------------------------------
*/

static Particle* closestParticle(Vec2f pos){
    int size = pVector.size();
    Particle* closestParticle = NULL;
    float dist = 1000000.0;
    for (int k = 0; k < size; ++k) {
        Vec2f l = pos - pVector[k]->m_Position;
        float normL = std::sqrt(std::abs(std::pow(l[0], 2)) + std::abs(std::pow(l[1], 2)));
        if (normL < dist) {
            dist = normL;
            closestParticle = pVector[k];
        }
    }
    if(dist > 1){
        closestParticle = new Particle(pos);
    }
    return closestParticle;
}

static void get_from_UI ()
{
	int i, j;
	// int size, flag;
	int hi, hj;
	// float x, y;
	if ( !mouse_down[0] && !mouse_down[2] && !mouse_release[0] 
	&& !mouse_shiftclick[0] && !mouse_shiftclick[2] ) return;

	i = (int)((       mx /(float)win_x)*N);
	j = (int)(((win_y-my)/(float)win_y)*N);

	if ( i<1 || i>N || j<1 || j>N ) return;

	//left mouse click
    if ( mouse_down[0] ) {

    }
    //right mouse click
    if ( mouse_down[2] ) {

    }

    //start pos
	hi = (int)((       hmx /(float)win_x)*N);
	hj = (int)(((win_y-hmy)/(float)win_y)*N);

	if( mouse_release[0] ) {

	    const Vec2f startPos((hi / (N / 2.0)) - 1.0, ((hj / (N / 2.0)) - 1.0));
	    std::cout << "start=" << startPos[0] << "," << startPos[1] << '\n';

	    const Vec2f currentPos((mx / (win_x / 2.0)) - 1.0, -((my / (win_y / 2.0)) - 1.0));
	    std::cout << "currentPos=" << currentPos[0] << "," << currentPos[1] << '\n';

	    //Particle* startClosestParticle = closestParticle(startPos);
        //Particle* currentClosestParticle = closestParticle(currentPos);

        int size = pVector.size();
        Particle*   startClosestParticle = NULL;
        Particle* currentClosestParticle = NULL;
        float sDist = 1000000.0;
        float cDist = 1000000.0;
        for (int k = 0; k < size; ++k) {
            Vec2f sV = startPos - pVector[k]->m_Position;
            Vec2f cV = currentPos - pVector[k]->m_Position;

            float sL = std::sqrt(std::abs(std::pow(sV[0], 2)) + std::abs(std::pow(sV[1], 2)));
            float cL = std::sqrt(std::abs(std::pow(cV[0], 2)) + std::abs(std::pow(cV[1], 2)));
            if (sL < sDist) {
                sDist = sL;
                startClosestParticle = pVector[k];
            }
            if (cL < cDist) {
                cDist = cL;
                currentClosestParticle = pVector[k];
            }
        }
        if(sDist > 0.1){
            startClosestParticle = new Particle(startPos);
            pVector.push_back(startClosestParticle);
        }
        if(cDist > 0.1){
            currentClosestParticle = new Particle(currentPos);
            pVector.push_back(currentClosestParticle);
        }

        //Particle* p = new Particle(startPos);
        fVector.push_back(new SpringForce(startClosestParticle, currentClosestParticle, 1, 0.05, 0.5));
        mouse_release[0] = 0;

	}

	omx = mx;
	omy = my;
}

static void remap_GUI()
{
	int ii, size = pVector.size();
	for(ii=0; ii<size; ii++)
	{
	    pVector[ii]->reset();
		//pVector[ii]->m_Position[0] = pVector[ii]->m_ConstructPos[0];
		//pVector[ii]->m_Position[1] = pVector[ii]->m_ConstructPos[1];
	}
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

	case ' ':
		dsim = !dsim;
		break;
	}
}

static void mouse_func ( int button, int state, int x, int y )
{
	omx = mx = x;
	omx = my = y;

	if(!mouse_down[0]){hmx=x; hmy=y;}
	if(mouse_down[button]) mouse_release[button] = state == GLUT_UP;
	if(mouse_down[button]) mouse_shiftclick[button] = glutGetModifiers()==GLUT_ACTIVE_SHIFT;
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
	if ( dsim ) simulation_step( pVector, fVector, dt );
	else        {get_from_UI();remap_GUI();}

	glutSetWindow ( win_id );
	glutPostRedisplay ();
}

static void display_func ( void )
{
	pre_display ();

	draw_forces();
	draw_constraints();
	draw_particles();

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
	win_id = glutCreateWindow ( "Particletoys!" );

	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);

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
		d = 5.f;
		fprintf ( stderr, "Using defaults : N=%d dt=%g d=%g\n",
			N, dt, d );
	} else {
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		d = atof(argv[3]);
	}

	printf ( "\n\nHow to use this application:\n\n" );
	printf ( "\t Toggle construction/simulation display with the spacebar key\n" );
	printf ( "\t Dump frames by pressing the 'd' key\n" );
	printf ( "\t Quit by pressing the 'q' key\n" );

	dsim = 0;
	dump_frames = 0;
	frame_number = 0;
	
	init_system();
	
	win_x = 512;
	win_y = 512;
	open_glut_window ();

	glutMainLoop ();

	exit ( 0 );
}

