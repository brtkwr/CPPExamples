//Written by Chris J K Williams, University of Bath, UK

//The following lines determine from where to get glut.h for Macintosh, Windows etc.

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define lastControlPoint 8

int WindowHalfWidth;
int WindowHalfHeight;

void MyGraphics(double BackGroundRed,double BackGroundBlue,double BackGroundGreen);

static void Draw(void);
static void Key(unsigned char key, int x, int y);

void MouseClick(int button, int state, int x, int y);
void MouseMove(int x, int y);
void MouseDrag(int x, int y);

void MoveOrigin(void);
void Rotate(void);

void writeDXF(void);

int argc;
char **argv;

double rv[2][2];	
GLUquadricObj *quadObj;

double  Zoom,a,ControlPointPos[lastControlPoint + 1][2];

int clickx,clicky,ControlPointOn[lastControlPoint + 1],StopStart;

bool LeftMouseDown=0,RightMouseDown=0;

void MyGraphics(double BackGroundRed,double BackGroundBlue,double BackGroundGreen)
{
    int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	
	cout<<"Your screen resolution is "<<screenWidth<<" x "<<screenHeight<<"\n";
	
	WindowHalfWidth =int(double(screenWidth)/2.0)-5;
	WindowHalfHeight=int(double(screenHeight)/2.0)-25;
	
	glutInitWindowSize(2*WindowHalfWidth, 2*WindowHalfHeight);
	
	Zoom=1.5*double(WindowHalfHeight)/a;
	glutInitWindowPosition(0,0);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Triangle");
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
	
    glutMainLoop();
}

static void Key(unsigned char key, int x, int y)
{
	switch (key)
	{
		case ' ':if(StopStart == 0)StopStart = 1;else StopStart = 0;break;
		case 'q'://writeDXF();
			gluDeleteQuadric(quadObj);exit(0);
		case 'Q'://writeDXF();
			gluDeleteQuadric(quadObj);exit(0);
		case '\033'://writeDXF();
			gluDeleteQuadric(quadObj);exit(0);			
	}
}

void MouseMove(int x, int y) 
{
	for(int ControlPoint = 0; ControlPoint <= lastControlPoint; ControlPoint ++)
	{
		if(
		   (double(x) - double(WindowHalfWidth) - ControlPointPos[ControlPoint][0] * Zoom) *
		   (double(x) - double(WindowHalfWidth) - ControlPointPos[ControlPoint][0] * Zoom) +
		   (double(y) - double(WindowHalfHeight) + ControlPointPos[ControlPoint][1] * Zoom) *
		   (double(y) - double(WindowHalfHeight) + ControlPointPos[ControlPoint][1] * Zoom)
		   < 10.0)ControlPointOn[ControlPoint] = 1;else
			   ControlPointOn[ControlPoint] = 0;
	}
}

void MouseClick(int button, int state, int x, int y) 
{
	if(state==GLUT_UP)glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
	if ((button==GLUT_LEFT_BUTTON)&&(state==GLUT_DOWN)) {
		LeftMouseDown=1;
	}
	if ((button==GLUT_LEFT_BUTTON)&&(state==GLUT_UP)) {
		LeftMouseDown=0;
	}
	
	if ((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_DOWN)) {
		RightMouseDown=1;
	}
	if ((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_UP)) {
		RightMouseDown=0;
	}
}

void MouseDrag(int x, int y) 
{
	if(LeftMouseDown||RightMouseDown)
	{
		for(int ControlPoint = 0;ControlPoint <= lastControlPoint; ControlPoint ++)
		{
			if(ControlPointOn[ControlPoint] == 1)
			{
				glutSetCursor(GLUT_CURSOR_CROSSHAIR);
				ControlPointPos[ControlPoint][0] = (double(x) - double(WindowHalfWidth)) / Zoom;
				ControlPointPos[ControlPoint][1] = - (double(y) - double(WindowHalfHeight)) / Zoom;
			}
		}
	}
	
	else glutSetCursor(GLUT_CURSOR_INFO);
}


