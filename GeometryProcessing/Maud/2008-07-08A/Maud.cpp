//Maud algorithm by Chris J K WIlliams, University of Bath
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
using namespace std;
//The following lines determine from where to get glut.h for Macintosh, Windows etc.
#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "SavePicture.h" //Necessary only if you plan to save images, for example to make a movie. See SavePicture(...,...) below
#define   LastNode      301
#define LastPrevious      5
void Calculation(void);
void Initialise(void);
void Contribution(int,int);
void Invert(int);
void MaudGraphics(float,float,float);
void MouseClick(int, int, int, int);
void MouseDrag(int x, int y);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
int argc;
char **argv;
GLUquadricObj *quadObj;
bool LeftMouseDown=0,RightMouseDown=0;
int  HalfDimension[2],ShowNodes,ShowLines,clickx,clicky,Random;
float x[2][LastNode+1],Force[2][LastNode+1],Velocity[2][LastNode+1],xPrevious[2][LastNode+1][LastPrevious+1],
Stiffness[2][2][LastNode+1],Flex[2][2],deltax[2],HalfPicture[2],ZoneSize[2],xAv[2],xAvChange[2],
PI,a,aSq,L0,L0Sq,MyExp,MyExp0,ToverL,dTbydL,LSq,tempForce,tempStiffness,slow,CarryOver,L,
Zoom,ZoomIncrement,xTrans,xTransInc,yTrans,yTransInc;
int main(int argc, char* argv[]){
glutInit(&argc,argv);
int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
cout<<"\nYour screen resolution is "<<screenWidth<<" x "<<screenHeight<<"\n\n";
HalfDimension[0]=int(float(screenWidth)/2.0)-5;
HalfDimension[1]=int(float(screenHeight)/2.0)-25;
for(int i=0;i<=1;i++)HalfPicture[i]=1.0*HalfDimension[i]-25.0;
cout<<"This program runs for ever until the window is closed or 'ESCAPE' or 'q' is pressed on the keyboard.\n\n";
cout<<"Typing 'n' makes the nodes dissapear or reappear.\n";
cout<<"Typing 'l' makes the lines dissapear or reappear.\n";
cout<<"Typing 'r' produces a random arrangement.\n";
cout<<"Typing 'g' produces a geometrical arrangement.\n";
cout<<"Typing 's' saves an image.\n";
PI=4.0*atan(1.0);srand(time(NULL));
slow=0.001;CarryOver=0.9;
Random=0;
a=50.0;aSq=a*a;
L0=100.0;L0Sq=L0*L0;
MyExp0=exp(-L0Sq/aSq);
Initialise();
MaudGraphics(1.0,1.0,1.0);
return 0;}
static void Draw(void){
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
Calculation();
glutSwapBuffers();}
void Calculation(void){
glColor3f(0.0,0.0,0.0);
glRasterPos2f(-HalfDimension[0]+10,-HalfDimension[1]+20);
/*char *MyText;
int LengthOfString;
MyText="Maud Algorithm - Chris Williams, Bath University";
LengthOfString=strlen(MyText);
for (int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);*/
for(int Node=0;Node<=LastNode;Node++){
for(int i=0;i<=1;i++){
Force[i][Node]=0.0;
for(int j=0;j<=1;j++)Stiffness[i][j][Node]=0.0;}}
for(int Node=0;Node<=LastNode-1;Node++){
for(int OtherNode=Node+1;OtherNode<=LastNode;OtherNode++){
LSq=0.0;
for(int i=0;i<=1;i++){
deltax[i]=x[i][OtherNode]-x[i][Node];
LSq+=deltax[i]*deltax[i];}
if(LSq/aSq>1.0e-24){
Contribution(Node,OtherNode);
if(ShowLines==1){
float MyBlue;
if(MyExp0<MyExp)MyBlue=1.0;else MyBlue=0.0;
glColor4f(0.0,0.0,MyBlue,fabs(MyExp0-MyExp));
glBegin(GL_LINES);
glVertex2f(x[0][Node],x[1][Node]);
glVertex2f(x[0][OtherNode],x[1][OtherNode]);
glEnd();}}}}
for(int Node=0;Node<=LastNode;Node++){
glBegin(GL_LINES);
for(int Previous=0;Previous<=LastPrevious;Previous++){
float factor=float(LastPrevious+1-Previous)/float(LastPrevious+1);
factor=factor*factor;
glColor4f(1.0,0.0,0.0,factor);
if(Previous==0)glVertex2f(x[0][Node],x[1][Node]);
else glVertex2f(xPrevious[0][Node][Previous-1],xPrevious[1][Node][Previous-1]);
glVertex2f(xPrevious[0][Node][Previous],xPrevious[1][Node][Previous]);}
glEnd();
if(ShowNodes==1){
glColor4f(1.0,0.0,0.0,1.0);glBegin(GL_POINTS);
glVertex2f(x[0][Node],x[1][Node]);
glEnd();}
for(int Previous=LastPrevious;Previous>=0;Previous--){
for(int i=0;i<=1;i++){
if(Previous==0)xPrevious[i][Node][0]=x[i][Node];
else xPrevious[i][Node][Previous]=xPrevious[i][Node][Previous-1];}}}
float xMin[2],xMax[2];
for(int i=0;i<=1;i++){xMin[i]=x[i][0];xMax[i]=x[i][0];}
for(int Node=0;Node<=LastNode;Node++){
Invert(Node);
for(int i=0;i<=1;i++){
Velocity[i][Node]=CarryOver*Velocity[i][Node];
for(int j=0;j<=1;j++)Velocity[i][Node]+=slow*Flex[i][j]*Force[j][Node];
x[i][Node]+=Velocity[i][Node];
if(xMin[i]>x[i][Node])xMin[i]=x[i][Node];
if(xMax[i]<x[i][Node])xMax[i]=x[i][Node];}}
for(int i=0;i<=1;i++){
xAvChange[i]=0.9*xAvChange[i]+0.00001*(xAv[i]-(xMax[i]+xMin[i])/2.0);
for(int Node=0;Node<=LastNode;Node++){
for(int i=0;i<=1;i++)x[i][Node]+=xAvChange[i];}}}
void MaudGraphics(float BackGroundRed,float BackGroundGreen,float BackGroundBlue){
glutInitWindowSize(2*HalfDimension[0], 2*HalfDimension[1]);
glutInitWindowPosition(0,0);
glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
glutCreateWindow("Maud");
glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
glViewport(0, 0, 2*HalfDimension[0], 2*HalfDimension[1]);
gluOrtho2D(-HalfDimension[0],HalfDimension[0],-HalfDimension[1],HalfDimension[1]);
glClear(GL_COLOR_BUFFER_BIT);
glLineWidth(2.0);
glPointSize(4.0);
glutKeyboardFunc(Key);
glutIdleFunc(Draw);
glutDisplayFunc(Draw);
glutMouseFunc(MouseClick);
glutMotionFunc(MouseDrag);
glEnable(GL_POINT_SMOOTH);
glEnable(GL_LINE_SMOOTH);
glEnable (GL_BLEND);
glBlendFunc (GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
Zoom=1.0;
glutMainLoop();}
static void Key(unsigned char key, int x, int y){
switch (key){
case 'n':if(ShowNodes==0)ShowNodes=1;else ShowNodes=0;break;
case 'l':if(ShowLines==0)ShowLines=1;else ShowLines=0;break;
case 'g':Random=0;Initialise();break;
case 'r':Random=1;Initialise();break;
case 's':SavePicture(HalfDimension[0],HalfDimension[1]);break;
case 'q':gluDeleteQuadric(quadObj);exit(0);
case '\033':gluDeleteQuadric(quadObj);exit(0);}}
void Contribution(int ThisNode,int ThisOtherNode){
L=sqrt(LSq);
MyExp=exp(-LSq/aSq);
ToverL=(MyExp0-MyExp)/L;
dTbydL=(2.0*L/aSq)*MyExp;
for(int i=0;i<=1;i++){
tempForce=ToverL*deltax[i];
Force[i][ThisNode]+=tempForce;
Force[i][ThisOtherNode]-=tempForce;
Stiffness[i][i][ThisNode]+=ToverL;
Stiffness[i][i][ThisOtherNode]+=ToverL;
for(int j=0;j<=1;j++){
tempStiffness=(dTbydL-ToverL)*deltax[i]*deltax[j]/LSq;
Stiffness[i][j][ThisNode]+=tempStiffness;
Stiffness[i][j][ThisOtherNode]+=tempStiffness;}}}
void Invert(int ThisNode){
float determinant;
determinant=Stiffness[0][0][ThisNode]*Stiffness[1][1][ThisNode]
-Stiffness[0][1][ThisNode]*Stiffness[1][0][ThisNode];
if(fabs(determinant)>1.0e-48){
Flex[0][0]=Stiffness[1][1][ThisNode]/determinant;
Flex[1][1]=Stiffness[0][0][ThisNode]/determinant;
Flex[0][1]=-Stiffness[0][1][ThisNode]/determinant;
Flex[1][0]=-Stiffness[1][0][ThisNode]/determinant;}
else{
Flex[0][0]=0.0;
Flex[1][1]=0.0;
Flex[0][1]=0.0;
Flex[1][0]=0.0;}}

void MouseDrag(int x, int y) {
float ZoomIncrement;

if(RightMouseDown){
ZoomIncrement=(1.0+0.001*float(x-clickx))*(1.0-0.001*float(y-clicky));
glScalef(ZoomIncrement,ZoomIncrement,ZoomIncrement);
Zoom=Zoom*ZoomIncrement;
clickx=x;clicky=y;}
if(LeftMouseDown){
xTransInc=+float(x-clickx)/Zoom;
yTransInc=-float(y-clicky)/Zoom;
glTranslatef(xTransInc,yTransInc,0.0);
clickx=x;clicky=y;}}

void MouseClick(int button, int state, int x, int y) {
if((button==GLUT_LEFT_BUTTON)&&(state==GLUT_DOWN)){LeftMouseDown=1;clickx=x;clicky=y;}
if((button==GLUT_LEFT_BUTTON)&&(state==GLUT_UP))LeftMouseDown=0;
if((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_DOWN)){RightMouseDown=1;clickx=x;clicky=y;}
if((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_UP))RightMouseDown=0;}

void Initialise(void){
ShowNodes=1;ShowLines=0;
for(int i=0;i<=1;i++){xAv[i]=0.0;xAvChange[i]=0.0;}
for(int Node=0;Node<=LastNode;Node++){
if(Random){for(int i=0;i<=1;i++)x[i][Node]=0.75*HalfPicture[i]*((2.0*rand())/(1.0*RAND_MAX)-1.0);}
else{
x[0][Node]=0.8*float(HalfDimension[0]*(2*Node-LastNode))/float(LastNode);
x[1][Node]=0.8*float(HalfDimension[1])*cos(6.0*x[0][Node]/float(HalfDimension[0]));}
for(int i=0;i<=1;i++){
Velocity[i][Node]=0.0;
for(int Previous=0;Previous<=LastPrevious;Previous++)xPrevious[i][Node][Previous]=x[i][Node];}}}

