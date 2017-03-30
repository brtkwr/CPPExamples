//Written by Chris J K Williams, University of Bath, UK

//The following lines determine from where to get glut.h for Macintosh, Windows etc.

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

int WindowHalfWidth;
int WindowHalfHeight;

void MyGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen);

static void Draw(void);
static void Key(unsigned char key, int x, int y);

void MouseClick(int button, int state, int x, int y);
void MouseMove(int x, int y);
void MouseDrag(int x, int y);

void MoveOrigin(void);
void Rotate(void);

int argc;
char **argv;

float rv[2][2];	
GLUquadricObj *quadObj;

float  InitialZoom,Zoom,ZoomIncrement,
xTrans=0.0,yTrans=0.0,zTrans=0.0,xTransInc,yTransInc,zTransInc,
zRotationIncrement,HorizontalRotationIncrement,zRotation=0.0,HorizontalRotation=0.0,cosz,cosH,sinz,sinH,PI,MaxHalfDimension,
InitialZoomButtonX,InitialZoomButtonY,ZoomButtonX,ZoomButtonY;

int clickx,clicky,ZoomButton;

bool LeftMouseDown=0,RightMouseDown=0;

void MyGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen)
{
    PI = 4.0 * atan(1.0);
	int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	
	cout<<"Your screen resolution is "<<screenWidth<<" x "<<screenHeight<<"\n";
	
	WindowHalfWidth =int(float(screenWidth)/2.0)-5;
	WindowHalfHeight=int(float(screenHeight)/2.0)-25;
	
	glutInitWindowSize(2*WindowHalfWidth, 2*WindowHalfHeight);
	
	InitialZoom=0.95*float(WindowHalfHeight)/MaxHalfDimension;
	Zoom=InitialZoom;
	glutInitWindowPosition(0,0);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Motion");
    glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
	glViewport(0, 0, 2*WindowHalfWidth, 2*WindowHalfHeight);
    glOrtho(-WindowHalfWidth,WindowHalfWidth,-WindowHalfHeight,WindowHalfHeight,-10.0*WindowHalfWidth,10.0*WindowHalfWidth);
	glClear(GL_COLOR_BUFFER_BIT);
	
    glLineWidth(0.5);
	
    glutKeyboardFunc(Key);
	glutIdleFunc(Draw);
	glutDisplayFunc(Draw);
	glutMouseFunc(MouseClick);
	glutMotionFunc(MouseDrag);
	glutPassiveMotionFunc(MouseMove);
	
	glEnable(GL_BLEND);
	glEnable(GL_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_TEXTURE);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glScaled(Zoom,Zoom,Zoom);
	zRotationIncrement=60.0;HorizontalRotationIncrement=-60.0;
	Rotate();
	
	InitialZoomButtonX = - 0.95;
	InitialZoomButtonY = - 0.0;
	
	ZoomButtonX = InitialZoomButtonX;
	ZoomButtonY = InitialZoomButtonY;
	
    glutMainLoop();
}

static void Key(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'q':gluDeleteQuadric(quadObj);exit(0);
		case 'Q':gluDeleteQuadric(quadObj);exit(0);
		case '\033':gluDeleteQuadric(quadObj);exit(0);
	}
}

void MouseMove(int x, int y) 
{
	if(
	   (float(x) / float(WindowHalfWidth) - 1.0 - ZoomButtonX) * (float(x) / float(WindowHalfWidth) - 1.0 - ZoomButtonX) +
	   (- (float(y) / float(WindowHalfHeight) - 1.0) + ZoomButtonY) * (- (float(y) / float(WindowHalfHeight) - 1.0) + ZoomButtonY) 
	   < 0.001)ZoomButton = 1;else
	   {
		   ZoomButton = 0;
		   ZoomButtonX = InitialZoomButtonX;
		   ZoomButtonY = InitialZoomButtonY;
	   }
}

void MouseClick(int button, int state, int x, int y) 
{
	if(state==GLUT_UP)glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
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
	
	if(LeftMouseDown||RightMouseDown)
	{
		if(ZoomButton == 1)
		{
			ZoomButtonX = float(x) / float(WindowHalfWidth) - 1.0;
			ZoomButtonY = - (float(y) / float(WindowHalfHeight) - 1.0);
			
			glutSetCursor(GLUT_CURSOR_CROSSHAIR);
			
			ZoomIncrement=(1.0+0.001*float(x-clickx))*(1.0-0.001*float(y-clicky));
			glScaled(ZoomIncrement,ZoomIncrement,ZoomIncrement);
			Zoom=Zoom*ZoomIncrement;
		}
		else
		{
			if(LeftMouseDown)
			{
				glutSetCursor(GLUT_CURSOR_INFO);
				
				TransInc=+float(x-clickx)/Zoom;
				xTransInc=TransInc*cosz;
				yTransInc=-TransInc*sinz;
				zTransInc=0.0;
				MoveOrigin();
				
				TransInc=-float(y-clicky)/Zoom;
				xTransInc=TransInc*sinz*cosH;
				yTransInc=TransInc*cosz*cosH;
				zTransInc=-TransInc*sinH;		
				MoveOrigin();
			}
			else
			{
				glutSetCursor(GLUT_CURSOR_CYCLE);
				
				HorizontalRotationIncrement = float(y-clicky) / 10.0;
				zRotationIncrement = float(x-clickx) / 10.0;
				Rotate();
			}
		}
		clickx=x;
		clicky=y;
	}
}

void Rotate(void)
{
	zRotation += zRotationIncrement;
	HorizontalRotation += HorizontalRotationIncrement;
	cosz = cos(zRotation*PI/180.0);cosH = cos(HorizontalRotation*PI/180.0);
	sinz = sin(zRotation*PI/180.0);sinH = sin(HorizontalRotation*PI/180.0);
	glRotated(zRotationIncrement,0.0,0.0,1.0);
	glRotated(HorizontalRotationIncrement,cosz,-sinz,0.0);
}

void MoveOrigin(void)
{
	xTrans+=xTransInc*Zoom/InitialZoom;
	yTrans+=yTransInc*Zoom/InitialZoom;
	zTrans+=zTransInc*Zoom/InitialZoom;
	glTranslated(xTransInc,yTransInc,zTransInc);
}

