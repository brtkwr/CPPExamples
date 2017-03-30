//The following lines determine from where to get glut.h for Macintosh, Windows etc.

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

int HalfWidth,HalfHeight;

void ChrisGraphics(double BackGroundRed,double BackGroundGreen,double BackGroundBlue);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
void MouseClick(int button, int state, int x, int y);
void MouseDrag(int x, int y);
void MoveOrigin(void);
void Rotate(void);

GLUquadricObj *quadObj;

#define   m 300
#define   n 200


double x[m+1][n+1],y[m+1][n+1],z[m+1][n+1];

double  Scale,InitialZoom,Zoom,ZoomIncrement,xTrans,xTransInc,yTrans,yTransInc;

double  zRotationIncrement,HorizontalRotationIncrement,zRotation,HorizontalRotation,
cosz,cosH,sinz,sinH,PI;

int argc,clickx,clicky,ZoomOrRotate=0;

char **argv;

bool LeftMouseDown=0,RightMouseDown=0;

void SetUpGraphics(double BackGroundRed,double BackGroundGreen,double BackGroundBlue)
{
	
    int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	
	cout<<"Your screen resolution is "<<screenWidth<<" x "<<screenHeight<<"\n";
	
	HalfWidth =int(double(screenWidth)/2.0)-5;
	HalfHeight=int(double(screenHeight)/2.0)-25;
	
	for(int i=0;i<=m;i++)
	{
		for(int j=0;j<=n;j++)
		{
			x[i][j]=HalfWidth*(2.0*double(i)/double(m)-1.0);
			y[i][j]=HalfWidth*(2.0*double(j)/double(m)-double(n)/double(m));
		}
	}
	
	InitialZoom=0.9;
	Zoom=InitialZoom;
	glutInitWindowSize(2*HalfWidth, 2*HalfHeight);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Water waves");
	glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
	glViewport(0, 0, 2*HalfWidth, 2*HalfHeight);
    glOrtho(-HalfWidth,HalfWidth,-HalfHeight,HalfHeight,-2.0*HalfWidth,2.0*HalfWidth);
	glClear(GL_COLOR_BUFFER_BIT);
	glutKeyboardFunc(Key);
	glutIdleFunc(Draw);
	glutDisplayFunc(Draw);
	glutMouseFunc(MouseClick);
	glutMotionFunc(MouseDrag);
	glEnable(GL_BLEND);
	glEnable(GL_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glScalef(Zoom,Zoom,Zoom);
	zRotationIncrement=10.0;HorizontalRotationIncrement=-60.0;Rotate();
	glutMainLoop();
}
static void Key(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'z':if(ZoomOrRotate==1)ZoomOrRotate=0;else ZoomOrRotate=1;break;
			
		case 'o':glScalef(InitialZoom/Zoom,InitialZoom/Zoom,InitialZoom/Zoom);Zoom=InitialZoom;break;
		case 'h':xTransInc=-xTrans;yTransInc=-yTrans;MoveOrigin();break;
		case 'r':zRotationIncrement=-zRotation;HorizontalRotationIncrement=-HorizontalRotation;Rotate();break;
		case 'q':gluDeleteQuadric(quadObj);exit(0);
		case 'Q':gluDeleteQuadric(quadObj);exit(0);
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
	double TransInc;
	
	if (LeftMouseDown)
    {
		TransInc=+static_cast<double>(x-clickx)/Zoom;
		xTransInc=TransInc;
		yTransInc=0.0;
		MoveOrigin();
		
		TransInc=-static_cast<double>(y-clicky)/Zoom;
		xTransInc=0.0;
		yTransInc=TransInc;		
		MoveOrigin();
    	clickx=x;
		clicky=y;
	}
	if (RightMouseDown)
	{
		if(ZoomOrRotate==0)
		{
			double X=x-double(HalfWidth);
			double Y=y-double(HalfHeight);
			
			zRotationIncrement=
				double(x-clickx)/10.0;
			HorizontalRotationIncrement=
				double(y-clicky)/10.0;
			Rotate();
		}
		
		else
		{
			ZoomIncrement=(1.0+0.001*static_cast<double>(x-clickx))*
			(1.0-0.001*static_cast<double>(y-clicky));
			glScaled(ZoomIncrement,ZoomIncrement,ZoomIncrement);
			Zoom=Zoom*ZoomIncrement;
		}
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
	glRotated(zRotationIncrement,0.0,0.0,1.0);
	glRotated(HorizontalRotationIncrement,cosz,-sinz,0.0);
}

void MoveOrigin(void)
{
	xTrans+=xTransInc*Zoom/InitialZoom;
	yTrans+=yTransInc*Zoom/InitialZoom;
	glTranslated(xTransInc,yTransInc,0.0);
}
