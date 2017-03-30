#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define maxpts 1500
#define LINE_STEP 4

//#include <GLUT/glut.h>  // For Macintosh
#include <GL/glut.h>    // For PC or Sun

void init_graphics(void);
static void Draw(void);
static void Key(unsigned char key, int x, int y);

int argc;
char **argv;

float rv[2][2];	
GLUquadricObj *quadObj;

double PI;

int main(void)
{
    PI=4.0*atan(1.0);
    init_graphics();
  
    glClear(GL_COLOR_BUFFER_BIT);
    glutKeyboardFunc(Key);
    glutIdleFunc(Draw);
    glutDisplayFunc(Draw);
    glutMainLoop();
      
    return 0;
}
    
void init_graphics(void)
{
    int halfwidth =  510;
    int halfheight = 400;

    glutInitWindowSize(2*halfwidth, 2*halfheight);
 //   glutInit(&argc, argv);//Remove for PC
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Random Net");
    glClearColor(1.0, 1.0, 1.0, 0.0);
	glViewport(0, 0, 2*halfwidth, 2*halfheight);
    gluOrtho2D(-halfwidth,halfwidth,-halfheight,halfheight);
}

static void Key(unsigned char key, int x, int y)
{
switch (key) {case 27:gluDeleteQuadric(quadObj);exit(0);}
}

static void Draw(void)
{
double x[3][maxpts],force[3][maxpts],movement[3][maxpts],a,radius;

int    end[2][2*maxpts],weighting[maxpts],
       node,member,lastnode,lastmember,lastcontrolnode,coord,cycle,angle,numangles;


lastnode=maxpts-1;
lastmember=2*maxpts-1;
lastcontrolnode=int(maxpts/10);
a=500.0;
radius=2.0;
numangles=5;

for(node=0;node<=lastnode;node++)
{
for(;;)
{
for(coord=0;coord<=2;coord++)x[coord][node]=a*((2.0*rand())/(1.0*RAND_MAX)-1.0);
if(x[0][node]*x[0][node]+x[1][node]*x[1][node]+x[2][node]*x[2][node]<=a*a)break;
}

for(coord=0;coord<=2;coord++)movement[coord][node]=0.0;

end[0][2*node+0]=node;end[1][2*node+0]=rand()%lastnode;
end[0][2*node+1]=node;end[1][2*node+1]=rand()%lastnode;
weighting[node]=0;
}

for(member=0;member<=lastmember;member++)
{
weighting[end[0][member]]+=1;
weighting[end[1][member]]+=1;
}

for(cycle=0;cycle<=75;cycle++)
{
for(node=0;node<=lastnode;node++)
{
for(coord=0;coord<=2;coord++)force[coord][node]=0.0;
}

for(member=0;member<=lastmember;member++)
{
for(coord=0;coord<=2;coord++)
{
force[coord][end[0][member]]+=x[coord][end[1][member]]-x[coord][end[0][member]];
force[coord][end[1][member]]+=x[coord][end[0][member]]-x[coord][end[1][member]];
}
}

for(node=lastcontrolnode+1;node<=lastnode;node++)
{
for(coord=0;coord<=2;coord++)
{
movement[coord][node]=force[coord][node]/(1.0*weighting[node])+0.7*movement[coord][node];
x[coord][node]+=movement[coord][node];
}
}

glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
glColor3f(0.5,0.5,0.5);
for(member=0;member<=lastmember;member++)
{
rv[0][0]=x[0][end[0][member]];
rv[0][1]=x[1][end[0][member]];

rv[1][0]=x[0][end[1][member]];
rv[1][1]=x[1][end[1][member]];

glBegin(GL_LINE_STRIP);
	glVertex2fv(rv[0]);
	glVertex2fv(rv[1]);
glEnd();
}
for(node=0;node<=lastnode;node++)
{
if(node<=lastcontrolnode)glColor3f(0.0,0.0,1.0);else glColor3f(1.0,0.0,0.0);
for(angle=0;angle<=numangles;angle++)
{
rv[0][0]=x[0][node]+radius*cos((2.0*PI*(angle+0))/(1.0*(numangles+1)));
rv[0][1]=x[1][node]+radius*sin((2.0*PI*(angle+0))/(1.0*(numangles+1)));

rv[1][0]=x[0][node]+radius*cos((2.0*PI*(angle+1))/(1.0*(numangles+1)));
rv[1][1]=x[1][node]+radius*sin((2.0*PI*(angle+1))/(1.0*(numangles+1)));

glBegin(GL_LINE_STRIP);
	glVertex2fv(rv[0]);
	glVertex2fv(rv[1]);
glEnd();
}
}
glutSwapBuffers();
}
}
