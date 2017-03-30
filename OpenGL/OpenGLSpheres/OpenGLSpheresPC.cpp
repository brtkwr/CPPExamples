#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

//#include <GLUT/glut.h>  // For Macintosh
#include <GL/glut.h>    // For PC or Sun
#define   Macintosh 0   //1 for Macintosh any other number for PC or Sun

#define AbsoluteLastSphere 135000

#define MaxLastZone0   90
#define MaxLastZone1    6
#define MaxLastZone2   90

#define MaxLastInZone 100

#define HalfWidth     508
#define HalfHeight    360

void ChrisGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen);
void DoTheDrawing(void);
static void Draw(void);
void CalculateMotion(void);
static void Key(unsigned char key, int x, int y);
void MouseClick(int button, int state, int x, int y);
void MouseDrag(int x, int y);
void Rotate(void);
void PlanView(void);
void MoveOrigin(void);
void CreateData(void);

GLUquadricObj *quadObj;

float  Coord[3][AbsoluteLastSphere+1],Radius[AbsoluteLastSphere+1],
       Vely[3][AbsoluteLastSphere+1],Colour[3][AbsoluteLastSphere+1],
	   Boundary[3],
	   PI,
	   	   
	   InitialZoom,Zoom,xTrans,xTransInc,yTrans,yTransInc,zTrans,zTransInc,
	   zRotationIncrement,HorizontalRotationIncrement,zRotation,HorizontalRotation,
	   cosz,cosH,sinz,sinH;

int    Sphere,LastSphere,xyz,
       argc,clickx,clicky,PanOrZoom=1;

char **argv;
bool LeftMouseDown=0,RightMouseDown=0;

int main(void)
{
PI=4.0*atan(1.0);
CreateData();

ChrisGraphics(0.0,0.0,0.0);

return 0;
}
    
static void Draw(void)
{
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
CalculateMotion();
DoTheDrawing();
}

void ChrisGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen)
{
    glutInitWindowSize(2*HalfWidth,2*HalfHeight);
	if(Macintosh==1)    glutInit(&argc, argv);               // Remove this line for PC or Sun
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Spheres");
    glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);

	GLfloat mat_specular[] = {0.8,0.0,0.0, 1.0};
	GLfloat mat_shininess[] = {50};
	
	GLfloat light0_position[] = {-200.0, 10.0, 100.0, 0.0};
	GLfloat light0_ambLevel[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat light0_diffLevel[] = {0.7, 0.7, 0.7, 1.0};
   
	GLfloat light1_position[] = {200.0, 10.0, -80.0, 0.0};
	GLfloat light1_ambLevel[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat light1_diffLevel[] = {1.0, 1.0, 1.0, 1.0};
	
	glShadeModel(GL_SMOOTH);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambLevel);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffLevel);

	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambLevel);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffLevel);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	
	glEnable(GL_COLOR_MATERIAL);
	
	glEnable(GL_DEPTH_TEST);

	glViewport(0,0,2*HalfWidth,2*HalfHeight);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

    glOrtho(-HalfWidth,HalfWidth,-HalfWidth,HalfWidth,-HalfWidth,HalfWidth);
    glClear(GL_COLOR_BUFFER_BIT);
	
    glutKeyboardFunc(Key);
    glutIdleFunc(Draw);
    glutDisplayFunc(Draw);
	
	glutMouseFunc(MouseClick);
	glutMotionFunc(MouseDrag);

	glEnable(GL_POINT_SMOOTH);// glDisable(GL_POINT_SMOOTH) is the reverse. Enable means antialiasing 
	
	glScalef(Zoom,Zoom,Zoom);
	
	zRotation=0.0;
    HorizontalRotation=0.0;
	zRotationIncrement=0.0;
	HorizontalRotationIncrement=-90.0;
	Rotate();
	
    glutMainLoop();	
}

static void Key(unsigned char key, int x, int y)
{
switch (key)
{
case 'p':PlanView();break;
case 'e':PlanView();
         zRotationIncrement=0.0;
		 HorizontalRotationIncrement=-90.0;
		 Rotate();break;
case 'i':PlanView();
		 zRotationIncrement=-45.0;
		 HorizontalRotationIncrement=-45.0;
		 Rotate();break;
case 'o':glScalef(InitialZoom/Zoom,InitialZoom/Zoom,InitialZoom/Zoom);Zoom=InitialZoom;break;
case 's':xTransInc=-xTrans;
		 yTransInc=-yTrans;
		 zTransInc=-zTrans;MoveOrigin();break;
case 'z':if(PanOrZoom==0)PanOrZoom=1;else PanOrZoom=0;break;
case 'q':gluDeleteQuadric(quadObj);exit(0);
case '\033':gluDeleteQuadric(quadObj);exit(0);
}
}

void MouseClick(int button, int state, int x, int y) 
{
	if ((button==GLUT_LEFT_BUTTON)&&(state==GLUT_DOWN)) {
		clickx=x;
		clicky=y;
		LeftMouseDown=1;
	}
	if ((button==GLUT_LEFT_BUTTON)&&(state==GLUT_UP)) {
		LeftMouseDown=0;
	}

	if ((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_DOWN)) {
		clickx=x;
		clicky=y;
		RightMouseDown=1;
	}
	if ((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_UP)) {
		RightMouseDown=0;
	}
}

void MouseDrag(int x, int y) 
{
float ZoomIncrement,TransInc;

	if (LeftMouseDown)
    {
        if(PanOrZoom==1)
		{
		TransInc=+static_cast<float>(x-clickx)/Zoom;
		xTransInc=TransInc*cosz;
		yTransInc=-TransInc*sinz;
		zTransInc=0.0;
		MoveOrigin();
		
		TransInc=-static_cast<float>(y-clicky)/Zoom;
		xTransInc=TransInc*sinz*cosH;
		yTransInc=TransInc*cosz*cosH;
		zTransInc=-TransInc*sinH;		
		MoveOrigin();
    }
		else
		{
		ZoomIncrement=(1.0+0.001*static_cast<float>(x-clickx))*
		              (1.0-0.001*static_cast<float>(y-clicky));
		glScalef(ZoomIncrement,ZoomIncrement,ZoomIncrement);
		Zoom=Zoom*ZoomIncrement;
		}

		clickx=x;
		clicky=y;
	}
	if (RightMouseDown)
	{
	    zRotationIncrement=static_cast<float>(x-clickx)/10;
		HorizontalRotationIncrement=static_cast<float>(y-clicky)/10;
		Rotate();
				
		clickx=x;
		clicky=y;
	}
}

void Rotate(void)
{
        zRotation+=zRotationIncrement;
		HorizontalRotation+=HorizontalRotationIncrement;
		cosz=cos(zRotation*PI/180.0);cosH=cos(HorizontalRotation*PI/180.0);
		sinz=sin(zRotation*PI/180.0);sinH=sin(HorizontalRotation*PI/180.0);
		glRotatef(zRotationIncrement,0.0,0.0,1.0);
		glRotatef(HorizontalRotationIncrement,cosz,-sinz,0.0);
}

void PlanView(void)
{
         zRotationIncrement=-zRotation;
         HorizontalRotationIncrement=-HorizontalRotation;
         Rotate();
}

void MoveOrigin(void)
{
		xTrans+=xTransInc;
		yTrans+=yTransInc;
		zTrans+=zTransInc;
		glTranslatef(xTransInc,yTransInc,zTransInc);
}

void DoTheDrawing(void)
{
for(Sphere=0;Sphere<=LastSphere;Sphere++)
{
glColor3f(Colour[0][Sphere],Colour[1][Sphere],Colour[2][Sphere]);
    glPushMatrix();
	glTranslatef(Coord[0][Sphere],Coord[1][Sphere],Coord[2][Sphere]);
	glutSolidSphere(Radius[Sphere],10,10);
	glPopMatrix();
}

glColor3f(0.0,0.0,0.6);

glBegin(GL_QUADS);
glVertex3f(-HalfHeight,-HalfHeight,0.0);
glVertex3f(+HalfHeight,-HalfHeight,0.0);
glVertex3f(+HalfHeight,+HalfHeight,0.0);
glVertex3f(-HalfHeight,+HalfHeight,0.0);
glEnd();

glutSwapBuffers();
}

void CalculateMotion(void)
{
for(Sphere=0;Sphere<=LastSphere;Sphere++)
{
Vely[2][Sphere]-=0.2;
for(xyz=0;xyz<=2;xyz++)
{
Coord[xyz][Sphere]+=Vely[xyz][Sphere];
if(Coord[xyz][Sphere]>+Boundary[xyz])
{
Coord[xyz][Sphere]=+Boundary[xyz];
if(Vely[xyz][Sphere]>0.0)Vely[xyz][Sphere]=-Vely[xyz][Sphere];
}
if(Coord[xyz][Sphere]<-Boundary[xyz])
{
Coord[xyz][Sphere]=-Boundary[xyz];
if(Vely[xyz][Sphere]<0.0)Vely[xyz][Sphere]=-Vely[xyz][Sphere];
}
}
}
}

void CreateData(void)
{
LastSphere=500;
InitialZoom=0.8;
Zoom=InitialZoom;
Boundary[0]=HalfWidth;
Boundary[1]=HalfWidth;
Boundary[2]=HalfHeight;

for(Sphere=0;Sphere<=LastSphere;Sphere++)
{
Radius[Sphere]=(10.0*rand())/(1.0*RAND_MAX)+1.0;
if((1.0*rand())/(1.0*RAND_MAX)<0.01)
{Colour[0][Sphere]=1.0;Colour[1][Sphere]=0.0;Colour[2][Sphere]=0.0;}
else
{Colour[0][Sphere]=1.0;Colour[1][Sphere]=1.0;Colour[2][Sphere]=1.0;}
for(xyz=0;xyz<=2;xyz++)
{
Coord[xyz][Sphere]=Boundary[xyz]*((2.0*rand())/(1.0*RAND_MAX)-1.0);
Vely[xyz][Sphere]=2.0*((2.0*rand())/(1.0*RAND_MAX)-1.0);
}
}
}
