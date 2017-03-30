#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <time.h>
using namespace std;
//#include <GLUT/glut.h>  // For Macintosh
#include <GL/glut.h>    // For PC or Sun
#define   Macintosh 2   //1 for Macintosh any other number for PC or Sun
#define   halfW         508
#define   halfH         360
/*
For Dev C++
Project->Project options

In the Parameters tab window add the following line in the Linker pane:
	
-lglut32 -lglu32 -lopengl32 -lwinmm -lgdi32

*/

void ChrisGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
int argc;
char **argv;
GLUquadricObj *quadObj;
int i,j;
float drawx[2],PI,angle,r;

int main(void)
{
	cout<<"This program runs forever until 'q' is pressed\n";
PI=4.0*atan(1.0);
j=0;
ChrisGraphics(0.8,0.8,0.8);
return 0;
}
static void Draw(void)
{
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
j+=1;if(j>360)j-=360;
glBegin(GL_LINE_STRIP);
for(i=0;i<=360;i++)
{
glColor3f(1.0*i/360.0,0.0,1.0-1.0*i/360.0);
angle=6.0*PI*(i+j)/360.0;
r=5.0*exp(0.1*angle)*(1.0+0.1*sin(20.0*angle));
drawx[0]=r*cos(angle+0.1*cos(20.0*angle));
drawx[1]=r*sin(angle+0.1*cos(20.0*angle));
glVertex2fv(drawx);
}
glEnd();
glutSwapBuffers();
}
void ChrisGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen)
{
	glutInitWindowSize(2*halfW, 2*halfH);
	if(Macintosh==1)    glutInit(&argc, argv);               // Remove this line for PC or Sun
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Simple OpenGL");
	glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
	glViewport(0, 0, 2*halfW, 2*halfH);
	gluOrtho2D(-halfW,halfW,-halfH,halfH);
	glClear(GL_COLOR_BUFFER_BIT);
	glLineWidth(1.0);
	glPointSize(2.0);
	glutKeyboardFunc(Key);
	glutIdleFunc(Draw);
	glutDisplayFunc(Draw);
	glutMainLoop();
}
static void Key(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'q':gluDeleteQuadric(quadObj);exit(0);
	}
}
