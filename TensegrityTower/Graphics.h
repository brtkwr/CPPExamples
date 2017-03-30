//Written by Chris J K Williams, University of Bath, UK

//The following lines determine from where to get glut.h for Macintosh, Windows etc.

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define LastButton 9

int WindowHalfWidth;
int WindowHalfHeight;

void MyGraphics(double,double,double);

static void Draw(void);
static void Key(unsigned char, int, int);

void MouseClick(int, int, int, int);
void MouseMove(int, int);
void MouseDrag(int, int);

void MoveOrigin(void);
void Rotate(void);

int argc;
char **argv;

GLUquadricObj *quadObj;

double  InitialZoom,Zoom,RelativeRingLength,RelativePrestress,ShrinkPerLevel,SizeMultiplier,
xTrans=0.0,yTrans=0.0,zTrans=0.0,xTransInc,yTransInc,zTransInc,
zRotationIncrement,HorizontalRotationIncrement,zRotation=0.0,HorizontalRotation=0.0,cosz,cosH,sinz,sinH,PI,MaxHalfDimension,
InitialButtonX,InitialButtonY[LastButton + 1],ButtonX,ButtonY[LastButton + 1];

int clickx,clicky,OverThisButton[LastButton + 1],RestartWhenReady,InitialStateWhenReady,FormfindUnstressLoad,WriteDXF_Now;

bool LeftMouseDown=0,RightMouseDown=0;

void MyGraphics(double BackGroundRed,double BackGroundBlue,double BackGroundGreen)
{
    PI = 4.0 * atan(1.0);
	int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	
	cout<<"Your screen resolution is "<<screenWidth<<" x "<<screenHeight<<".\n";
	
	WindowHalfWidth =int(double(screenWidth)/2.0)-5;
	WindowHalfHeight=int(double(screenHeight)/2.0)-25;
	
	glutInitWindowSize(2*WindowHalfWidth, 2*WindowHalfHeight);
	
	InitialZoom=0.8*double(WindowHalfHeight)/MaxHalfDimension;
	Zoom=InitialZoom;
	glutInitWindowPosition(0,0);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Tensegrity tower");
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
	glTranslated(0.0,- 0.8 * MaxHalfDimension,0.0);
	zRotationIncrement=90.0;HorizontalRotationIncrement=-90.0;Rotate();
	
	InitialButtonX = - 0.95;
	ButtonX = InitialButtonX;
	
	for(int Button = 0;Button <= LastButton; Button ++)
	{
		InitialButtonY[Button] = - 0.1 * double(Button);
		
		ButtonY[Button] = InitialButtonY[Button];
	}
	
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
	for(int Button = 0;Button <= LastButton; Button ++)
	{
		if(
		   (   double(x) / double(WindowHalfWidth)  - 1.0  - ButtonX       ) * (    double(x) / double(WindowHalfWidth)  - 1.0  - ButtonX) +
		   (- (double(y) / double(WindowHalfHeight) - 1.0) - ButtonY[Button]) * (- (double(y) / double(WindowHalfHeight) - 1.0) - ButtonY[Button]) 
		   < 0.001)
			OverThisButton[Button] = 1;
		else
		{
			OverThisButton[Button] = 0;
			ButtonX = InitialButtonX;
			ButtonY[Button] = InitialButtonY[Button];
		}
	}
}

void MouseClick(int ButtonDown, int state, int x, int y) 
{
	if(state==GLUT_UP)glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
	if ((ButtonDown==GLUT_LEFT_BUTTON)&&(state==GLUT_DOWN)) {
		clickx=x;
		clicky=y;
		LeftMouseDown=1;
	}
	if ((ButtonDown==GLUT_LEFT_BUTTON)&&(state==GLUT_UP)) {
		LeftMouseDown=0;
		if(OverThisButton[LastButton - 4] == 1)FormfindUnstressLoad = 1;
		if(OverThisButton[LastButton - 3] == 1)FormfindUnstressLoad = 2;
		if(OverThisButton[LastButton - 2] == 1)RestartWhenReady = 1;
		if(OverThisButton[LastButton - 1] == 1)InitialStateWhenReady = 1;
		if(OverThisButton[LastButton - 0] == 1)WriteDXF_Now = 1;
	}
	
	if ((ButtonDown==GLUT_RIGHT_BUTTON)&&(state==GLUT_DOWN)) {
		clickx=x;
		clicky=y;
		RightMouseDown=1;
	}
	if ((ButtonDown==GLUT_RIGHT_BUTTON)&&(state==GLUT_UP)) {
		RightMouseDown=0;
	}
}

void MouseDrag(int x, int y) 
{
	double TransInc;
	
	if(LeftMouseDown||RightMouseDown)
	{
		int OverAnyButton = 0;
		for(int Button = 0;Button <= LastButton; Button ++)
		{
			if(OverThisButton[Button] == 1)
			{
				OverAnyButton = 1;
				ButtonX = double(x) / double(WindowHalfWidth) - 1.0;
				ButtonY[Button] = - (double(y) / double(WindowHalfHeight) - 1.0);
				
				glutSetCursor(GLUT_CURSOR_CROSSHAIR);
				
				if(Button == 0)
				{
					double ZoomIncrement=(1.0+0.001*double(x-clickx))*(1.0-0.001*double(y-clicky));
					glScaled(ZoomIncrement,ZoomIncrement,ZoomIncrement);
					Zoom=Zoom*ZoomIncrement;
				}
				if(Button == 1)
				{
					double RelativeRingLengthIncrement=1.0+0.0001*double(y-clicky);
					RelativeRingLength *= RelativeRingLengthIncrement;
				}
				if(Button == 2)
				{
					double RelativePrestressIncrement=1.0+0.001*double(y-clicky);
					RelativePrestress *= RelativePrestressIncrement;
				}
				
				if(Button == 3)
				{
					double RelativeShrinkPerLevel=1.0-0.001*double(y-clicky);
					ShrinkPerLevel *= RelativeShrinkPerLevel;
				}
				
				if(Button == 4)
				{
					double RelativeSizeMultiplier=1.0-0.001*double(y-clicky);
					SizeMultiplier *= RelativeSizeMultiplier;
				}
			}
		}
		
		if(OverAnyButton == 0)
		{
			if(LeftMouseDown)
			{
				glutSetCursor(GLUT_CURSOR_INFO);
				
				TransInc=+double(x-clickx)/Zoom;
				xTransInc=TransInc*cosz;
				yTransInc=-TransInc*sinz;
				zTransInc=0.0;
				MoveOrigin();
				
				TransInc=-double(y-clicky)/Zoom;
				xTransInc=TransInc*sinz*cosH;
				yTransInc=TransInc*cosz*cosH;
				zTransInc=-TransInc*sinH;		
				MoveOrigin();
			}
			else
			{
				glutSetCursor(GLUT_CURSOR_CYCLE);
				
				HorizontalRotationIncrement = double(y-clicky) / 10.0;
				zRotationIncrement = double(x-clickx) / 10.0;
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

