#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
using namespace std;
#include <GLUT/glut.h>  // For Macintosh
//#include <GL/glut.h>    // For PC or Sun
#define   Macintosh 1  //1 for Macintosh any other number for PC or Sun
#define   LastNode   100
#define   LastCurve	1000
void Calculation(void);
void ChrisGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
int argc;
char **argv;
GLUquadricObj *quadObj;
int HalfDimension[2],Node,i,Curve,points,lines;
float x[2][LastNode+1][LastCurve+1],Force[2][LastNode+1],Velocity[2][LastNode+1],deltax[2],
L[LastNode],LSq[LastNode],EA[LastNode],
g,deltat,m,mg,ActualLSq,ToverL,mu,OneMinusmu,RedPeak,peak;
int main(void){cout<<"'q' stops the program.\n'p' turns points off and on.\n'l' turns lines off and on.\n";
HalfDimension[0]=630;
HalfDimension[1]=370;
m=1.0;g=0.001;mg=m*g;
for(Node=0;Node<=LastNode-1;Node++){
L[Node]=300.0*(Node+1)/(0.5*(LastNode-1)*LastNode);
LSq[Node]=L[Node]*L[Node];
EA[Node]=(Node+1)/(1.0*LastNode);}
deltat=2.0*sqrt((m*L[0])/(EA[0]/L[0]));
RedPeak=0.1;points=1;lines=1;
for(Curve=0;Curve<=LastCurve;Curve++){
x[0][LastNode][Curve]=0.0;
x[1][LastNode][Curve]=300.0;
for(Node=LastNode-1;Node>=0;Node--){
x[0][Node][Curve]=x[0][Node+1][Curve]+L[Node];
x[1][Node][Curve]=x[1][Node+1][Curve];}
for(i=0;i<=1;i++)Velocity[i][Node]=0.0;}
ChrisGraphics(0.0,0.0,0.0);
return 0;}
static void Draw(void){
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
Calculation();
glutSwapBuffers();}
void Calculation(void){
for(Node=0;Node<=LastNode;Node++){
Force[0][Node]=0.0;Force[1][Node]=-mg*L[Node];}
for(Node=0;Node<=LastNode-1;Node++){
for(i=0;i<=1;i++)
deltax[i]=x[i][Node+1][0]-x[i][Node][0];
ActualLSq=deltax[0]*deltax[0]+deltax[1]*deltax[1];
ToverL=EA[Node]*(ActualLSq-LSq[Node])/(2.0*L[Node]*ActualLSq);
if(ToverL>0.0){
for(i=0;i<=1;i++){
Force[i][Node  ]+=ToverL*deltax[i];
Force[i][Node+1]-=ToverL*deltax[i];}}}
for(Node=0;Node<=LastNode-1;Node++){
for(Curve=LastCurve;Curve>=1;Curve--){
for(i=0;i<=1;i++)x[i][Node][Curve]=x[i][Node][Curve-1];}
for(i=0;i<=1;i++){
Velocity[i][Node]+=deltat*Force[i][Node]/(m*L[Node]);
x[i][Node][0]+=deltat*Velocity[i][Node];}}
glPointSize(5.0);glColor4f(1.0,0.0,0.0,1.0);glBegin(GL_POINTS);
glVertex2f(x[0][LastNode][0],x[1][LastNode][0]);glEnd();
if(points==1){
glPointSize(2.0);glColor4f(1.0,0.0,0.0,1.0);glBegin(GL_POINTS);
for(Node=0;Node<=LastNode-1;Node++)
glVertex2f(x[0][Node][0],x[1][Node][0]);
glEnd();}
if(lines==1){
for(Curve=0;Curve<=LastCurve;Curve++){
glLineWidth(1.5);mu=(1.0*Curve)/(1.0*LastCurve);OneMinusmu=1.0-mu;
peak=exp(-20.0*(mu+RedPeak)*(mu+RedPeak));
glColor4f(0.5+0.5*peak,0.5-0.5*peak,1.0-0.5*mu-0.5*peak,(1.0-exp(-100.0*mu*mu))*OneMinusmu*OneMinusmu*OneMinusmu*OneMinusmu);
glBegin(GL_LINES);
for(Node=0;Node<=LastNode-1;Node++){
glVertex2f(x[0][Node  ][Curve],x[1][Node  ][Curve]);
glVertex2f(x[0][Node+1][Curve],x[1][Node+1][Curve]);}
glEnd();}}}
void ChrisGraphics(float BackGroundRed,float BackGroundGreen,float BackGroundBlue){
glutInitWindowSize(2*HalfDimension[0], 2*HalfDimension[1]);
if(Macintosh==1)    glutInit(&argc, argv);               // Remove this line for PC or Sun
glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
glutCreateWindow("Chain");
glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
glViewport(0, 0, 2*HalfDimension[0], 2*HalfDimension[1]);
gluOrtho2D(-HalfDimension[0],HalfDimension[0],-HalfDimension[1],HalfDimension[1]);
glClear(GL_COLOR_BUFFER_BIT);
glutKeyboardFunc(Key);
glutIdleFunc(Draw);
glutDisplayFunc(Draw);
glEnable(GL_POINT_SMOOTH);
glEnable (GL_BLEND);glBlendFunc (GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
glutMainLoop();}
static void Key(unsigned char key, int x, int y){
switch (key){
case 'p':if(points==0)points=1;else points=0;break;
case 'l':if(lines==0)lines=1;else lines=0;break;
case 'q':gluDeleteQuadric(quadObj);exit(0);}}
