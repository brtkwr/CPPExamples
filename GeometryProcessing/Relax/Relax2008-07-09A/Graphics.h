//Written by Chris J K Williams, University of Bath, UK


//The following lines determine from where to get glut.h for Macintosh, Windows etc.

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

int WindowHalfWidth;
int WindowHalfHeight;

void SetUpGraphics(double,double,double);
static void Draw(void);
static void Key(unsigned char, int, int);
void MouseClick(int, int, int, int);
double SquareDist(int,int,int);
void MouseDrag(int x, int y);
void MouseMove(int x, int y);
void MoveOrigin(void);
void Controls(void);
void Rotate(void);
void PlanView(void);
void PrintTheText(void);
void SavePicture(void);

GLUquadricObj *quadObj;

double  InitialZoom,Zoom,ZoomIncrement,
xTrans,xTransInc,yTrans,yTransInc,zTrans,zTransInc,
zRotationIncrement,HorizontalRotationIncrement,zRotation,HorizontalRotation,
cosz,cosH,sinz,sinH,
PI,
Leftx,Lefty,
RadioLeftx=-0.95,RadioTopy=0.9,RadioSpacingx=0.025,RadioSpacingy=0.09,HalfLineHeight=0.02,TextRightShft=0.02,TextDownshift=0.008,
MenuBoundaryx=-0.77,
MaxSquareDist;

int clickx,clicky,FileNumber,ZoomPanRotate=2,ShowMemberLocalAxis=0,MemberTypeToShow=-1,ShowTensionOnlyElements=1,ShowLoads=0,ShowPeak=1,ShowInitial=1;

bool LeftMouseDown=0,RightMouseDown=0;

typedef struct {
	char id_len;                 // ID Field (Number of bytes - max 255)
	char map_type;               // Colormap Field (0 or 1)
	char img_type;               // Image Type (7 options - color vs. compression)
	int  map_first;              // Color Map stuff - first entry index
	int  map_len;                // Color Map stuff - total entries in file
	char map_entry_size;         // Color Map stuff - number of bits per entry
	int  x;                      // X-coordinate of origin 
	int  y;                      // Y-coordinate of origin
	int  width;                  // Width in Pixels
	int  height;                 // Height in Pixels
	char bpp;                    // Number of bits per pixel
	char misc;                   // Other stuff - scan origin and alpha bits
} targa_header;

void writeheader(targa_header h, FILE *tga);

void SetUpGraphics(double BackGroundRed,double BackGroundGreen,double BackGroundBlue)
{
	int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	
	cout<<"Your screen resolution is "<<screenWidth<<" x "<<screenHeight<<"\n";
	
	WindowHalfWidth =int(double(screenWidth)/2.0)-5;
	WindowHalfHeight=int(double(screenHeight)/2.0)-25;
	
	glutInitWindowSize(2*WindowHalfWidth, 2*WindowHalfHeight);
	
	InitialZoom=0.95*double(WindowHalfHeight)/MaxHalfDimension;
	Zoom=InitialZoom;
	glutInitWindowPosition(0,0);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Relax");
	glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
	glOrtho(-WindowHalfWidth,WindowHalfWidth,-WindowHalfHeight,WindowHalfHeight,-10.0*WindowHalfWidth,10.0*WindowHalfWidth);
	glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
	glClear(GL_COLOR_BUFFER_BIT);
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
	Rotate();
	glutMainLoop();
}

static void Key(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'q':WriteResults();gluDeleteQuadric(quadObj);exit(0);
		case 'Q':WriteResults();gluDeleteQuadric(quadObj);exit(0);
		case '\033':WriteResults();gluDeleteQuadric(quadObj);exit(0);
	}
}

void MouseClick(int button, int state, int x, int y) 
{
	if((button==GLUT_LEFT_BUTTON)&&(state==GLUT_DOWN))
	{
		MaxSquareDist=RadioSpacingy*RadioSpacingy/10.0;
		
		clickx=x;
		clicky=y;
		LeftMouseDown=1;
		
		Leftx=double(x)/double(WindowHalfWidth)-1.0;
		Lefty=1.0-double(y)/double(WindowHalfHeight);
		
		int DropButton=-1;
		
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist)ZoomPanRotate=0;
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist)ZoomPanRotate=1;
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist)ZoomPanRotate=2;
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist)PlanView();
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist)
		{
			PlanView();
			zRotationIncrement=0.0;
			HorizontalRotationIncrement=-90.0;
			Rotate();
		}
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist)
		{
			PlanView();
			zRotationIncrement=+90.0;
			HorizontalRotationIncrement=-90.0;
			Rotate();
		}
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist)
		{
			PlanView();
			zRotationIncrement=-90.0;
			HorizontalRotationIncrement=-90.0;
			Rotate();
		}
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist)
		{
			PlanView();
			zRotationIncrement=45.0;
			HorizontalRotationIncrement=-180.0*acos(1.0/sqrt(3.0))/PI;
			Rotate();
		}
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist)
		{
			xTransInc=-xTrans*InitialZoom/Zoom;
			yTransInc=-yTrans*InitialZoom/Zoom;
			zTransInc=-zTrans*InitialZoom/Zoom;
			MoveOrigin();
			ZoomIncrement=InitialZoom/Zoom;
			glScaled(ZoomIncrement,ZoomIncrement,ZoomIncrement);
			Zoom=Zoom*ZoomIncrement;
		}
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist){if(ShowInitial==0)ShowInitial=1;else ShowInitial=0;}
		DropButton++;if(SquareDist(0,DropButton,-1)<MaxSquareDist)MemberTypeToShow--;
		             if(SquareDist(0,DropButton,+1)<MaxSquareDist)MemberTypeToShow++;
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist){if(ShowTensionOnlyElements==0)ShowTensionOnlyElements=1;else ShowTensionOnlyElements=0;}
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist){if(ShowMemberLocalAxis==0)ShowMemberLocalAxis=1;else ShowMemberLocalAxis=0;}
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist){if(ShowLoads==0)ShowLoads=1;else ShowLoads=0;}
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist){if(NonLinear==0)NonLinear=1;else NonLinear=0;}
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist){if(ShowPeak==0)ShowPeak=1;else ShowPeak=0;}
		DropButton++;if(SquareDist(0,DropButton,-1)<MaxSquareDist)OwnWtSafetyFact-=0.1;
		             if(SquareDist(0,DropButton,+1)<MaxSquareDist)OwnWtSafetyFact+=0.1;
		DropButton++;if(SquareDist(0,DropButton,-1)<MaxSquareDist)AppliedSafetyFact-=0.1;
		             if(SquareDist(0,DropButton,+1)<MaxSquareDist)AppliedSafetyFact+=0.1;
		DropButton++;if(SquareDist(0,DropButton,0)<MaxSquareDist){WriteResults();gluDeleteQuadric(quadObj);exit(0);}		
	}
	if((button==GLUT_LEFT_BUTTON)&&(state==GLUT_UP)) {
		LeftMouseDown=0;
	}
	
	if((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_DOWN)) {
		clickx=x;
		clicky=y;
		RightMouseDown=1;
	}
	if((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_UP)) {
		RightMouseDown=0;
	}
}

double SquareDist(int RadioH,int RadioV,int RadioLine)
{
return
(Leftx-(RadioLeftx+double(RadioH)*RadioSpacingx))*
(Leftx-(RadioLeftx+double(RadioH)*RadioSpacingx))+
(Lefty-(RadioTopy-double(RadioV)*RadioSpacingy+double(RadioLine)*HalfLineHeight))*
(Lefty-(RadioTopy-double(RadioV)*RadioSpacingy+double(RadioLine)*HalfLineHeight));
}

void MouseDrag(int x, int y) 
{
	double ZoomIncrement,TransInc;
	
		if(LeftMouseDown||RightMouseDown)
	{
		if((LeftMouseDown&&ZoomPanRotate==0)||RightMouseDown)
		{
			glutSetCursor(GLUT_CURSOR_UP_DOWN);
			
			ZoomIncrement=(1.0+0.001*double(x-clickx))*(1.0-0.001*double(y-clicky));
			glScaled(ZoomIncrement,ZoomIncrement,ZoomIncrement);
			Zoom=Zoom*ZoomIncrement;
		}
		
		if(LeftMouseDown&&ZoomPanRotate==1)
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
		
		if(LeftMouseDown&&ZoomPanRotate==2)
		{
			glutSetCursor(GLUT_CURSOR_CYCLE);
			
			double X=x-double(WindowHalfWidth);
			double Y=y-double(WindowHalfHeight);
			
			zRotationIncrement=
			(double(x-clickx)*Y-double(y-clicky)*X)/1000.0;
			HorizontalRotationIncrement=
			(double(x-clickx)*X+double(y-clicky)*Y)*fabs(Y)/(1000.0*sqrt(X*X+Y*Y));
			Rotate();
		}
		clickx=x;
		clicky=y;
	}
}

void MouseMove(int x, int y) 
{
	if(double(x)/double(WindowHalfWidth)>1.0+MenuBoundaryx)glutSetCursor(GLUT_CURSOR_LEFT_ARROW);else glutSetCursor(GLUT_CURSOR_CROSSHAIR);
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

void PlanView(void)
{
	zRotationIncrement=-zRotation;
	HorizontalRotationIncrement=-HorizontalRotation;
	Rotate();
}

void MoveOrigin(void)
{
	xTrans+=xTransInc*Zoom/InitialZoom;
	yTrans+=yTransInc*Zoom/InitialZoom;
	zTrans+=zTransInc*Zoom/InitialZoom;
	glTranslated(xTransInc,yTransInc,zTransInc);
}

void Controls(void)
{
	double PointPosition[2];
	int NumberLength;
	
	glPushMatrix();
	glLoadIdentity();
	int LastRadioButton=18;
	
	glColor4f(1.0,1.0,1.0,0.7);
	glBegin(GL_QUADS);
	PointPosition[0]=MenuBoundaryx;PointPosition[1]=1.0;glVertex2dv(PointPosition);
	PointPosition[1]=-1.0;glVertex2dv(PointPosition);
	PointPosition[0]=-1.0;glVertex2dv(PointPosition);
	PointPosition[1]=+1.0;glVertex2dv(PointPosition);
	glEnd();
	
	glLineWidth(1.0);
	glColor4f(0.0,0.0,0.0,1.0);
	glBegin(GL_LINE_STRIP);
	PointPosition[0]=MenuBoundaryx;PointPosition[1]=1.0;glVertex2dv(PointPosition);
	PointPosition[1]=-1.0;glVertex2dv(PointPosition);
	glEnd();

    double TotalLoadMag=0.0;
	for(int i=0;i<=2;i++)TotalLoadMag+=TotalFactoredLoad[i]*TotalFactoredLoad[i];
	TotalLoadMag=sqrt(TotalLoadMag);
	if(TotalLoadMag!=0.0)
	{
	double IndicatorLength=0.9;
	for(int i=0;i<=2;i++)
	{
	glLineWidth(4.0);
	glColor4f(0.0,0.0,0.0,0.5);
	glBegin(GL_LINE_STRIP);
	PointPosition[0]=MenuBoundaryx*0.95+0.01*double(i);
	PointPosition[1]=-IndicatorLength;glVertex2dv(PointPosition);
	PointPosition[1]=0.0;glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.5,1.0,0.5);
	glBegin(GL_LINE_STRIP);
	glVertex2dv(PointPosition);
	PointPosition[1]=+IndicatorLength;glVertex2dv(PointPosition);
	glEnd();
	glLineWidth(2.0);
	glColor4f(0.9,1.0,0.0,1.0);
	glBegin(GL_LINE_STRIP);
	PointPosition[1]=0.0;glVertex2dv(PointPosition);
	PointPosition[1]=IndicatorLength*(TotalReaction[i]+TotalFactoredLoad[i])/TotalLoadMag;glVertex2dv(PointPosition);
	glEnd();
	glPointSize(3.0);
	glColor4f(1.0,0.0,0.0,1.0);
	glBegin(GL_POINTS);
	glVertex2dv(PointPosition);
	glEnd();
	}
	}
	
	glPointSize(7.0);
	glColor4f(0.0,0.0,0.0,1.0);
	glBegin(GL_POINTS);
	for(int RadioButton=0;RadioButton<=LastRadioButton;RadioButton++)
	{
	PointPosition[0]=RadioLeftx;
	if(RadioButton!=10&&RadioButton!=16&&RadioButton!=17)
	{
	PointPosition[1]=RadioTopy-double(RadioButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	}
	else
	{
	PointPosition[1]=RadioTopy-double(RadioButton)*RadioSpacingy+HalfLineHeight;
	glVertex2dv(PointPosition);
	PointPosition[1]=RadioTopy-double(RadioButton)*RadioSpacingy-HalfLineHeight;
	glVertex2dv(PointPosition);
	}
	}
	glEnd();
	
	glPointSize(5.0);
	glColor4f(1.0,1.0,1.0,1.0);
	glBegin(GL_POINTS);
	for(int RadioButton=0;RadioButton<=LastRadioButton;RadioButton++)
	{
	PointPosition[0]=RadioLeftx;
	if(RadioButton!=10&&RadioButton!=11&&RadioButton!=15)
	{
	PointPosition[1]=RadioTopy-double(RadioButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	}
	else
	{
	PointPosition[1]=RadioTopy-double(RadioButton)*RadioSpacingy+HalfLineHeight;
	glVertex2dv(PointPosition);
	PointPosition[1]=RadioTopy-double(RadioButton)*RadioSpacingy-HalfLineHeight;
	glVertex2dv(PointPosition);
	}
	}
	glEnd();
	
	int DropButton=-1;
	
	if(ZoomPanRotate==0)glColor4f(1.0,0.0,0.0,0.8);else glColor4f(0.0,0.0,1.0,0.8);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift);
	MyText="Zoom";PrintTheText();
	
	if(ZoomPanRotate==1)glColor4f(1.0,0.0,0.0,0.8);else glColor4f(0.0,0.0,1.0,0.8);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift);
	MyText="Pan";PrintTheText();
	
	if(ZoomPanRotate==2)glColor4f(1.0,0.0,0.0,0.8);else glColor4f(0.0,0.0,1.0,0.8);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift);
	MyText="Rotate";PrintTheText();
	
	glColor4f(0.0,0.0,0.0,0.5);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Plan";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="view";PrintTheText();
	
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Front";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="view";PrintTheText();
	
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Left";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="view";PrintTheText();
	
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Right";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="view";PrintTheText();
	
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Isometric";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="view";PrintTheText();

	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Zoom";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="extents";PrintTheText();
	
	if(ShowInitial==0)glColor4f(1.0,1.0,0.0,0.8);else glColor4f(0.0,1.0,1.0,0.8);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Show initial";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="shape";PrintTheText();
	
	if(NonLinear==1&&ShowPeak==1)glColor4f(1.0,1.0,1.0,1.0);else glColor4f(0.0,0.0,0.0,0.5);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;PointPosition[1]+=HalfLineHeight;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+0.5*TextRightShft,PointPosition[1]-TextDownshift);
	MyText="+ Show member";PrintTheText();
	if(NonLinear==1&&ShowPeak==1)glColor4f(1.0,1.0,1.0,1.0);else glColor4f(0.0,0.0,0.0,0.5);
	glBegin(GL_POINTS);
	PointPosition[1]-=2.0*HalfLineHeight;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+0.5*TextRightShft,PointPosition[1]-TextDownshift);
	MyText="- type ";PrintTheText();
	sprintf(NumericalValue,"%d",MemberTypeToShow);NumberLength=strlen(NumericalValue);
	for(int MyChar=0;MyChar<NumberLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyChar]);
	
	if(NonLinear==1&&ShowPeak==1)glColor4f(1.0,1.0,1.0,1.0);
	else{if(ShowTensionOnlyElements==0)glColor4f(1.0,1.0,0.0,0.8);else glColor4f(0.0,1.0,1.0,0.8);}
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Show tenson";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="only members";PrintTheText();
	
	if(ShowMemberLocalAxis==0)glColor4f(1.0,1.0,0.0,0.8);else glColor4f(0.0,1.0,1.0,0.8);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Show member";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="local axis";PrintTheText();
	
	if(ShowLoads==0)glColor4f(1.0,1.0,0.0,0.8);else glColor4f(0.0,1.0,1.0,0.8);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Show factored";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="loads";PrintTheText();
	
	if(NonLinear==0)glColor4f(1.0,1.0,0.0,0.8);else glColor4f(0.0,1.0,1.0,0.8);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Nonlinear";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="elastic";PrintTheText();
	
	if(NonLinear==1){if(ShowPeak==0)glColor4f(1.0,1.0,0.0,0.8);else glColor4f(0.0,1.0,1.0,0.8);}
	else glColor4f(1.0,1.0,1.0,1.0);
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift+HalfLineHeight);
	MyText="Show members";PrintTheText();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift-HalfLineHeight);
	MyText="near peak";PrintTheText();
	
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;PointPosition[1]+=HalfLineHeight;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+0.5*TextRightShft,PointPosition[1]-TextDownshift);
	MyText="+ Own weight";PrintTheText();
	glBegin(GL_POINTS);
	PointPosition[1]-=2.0*HalfLineHeight;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+0.5*TextRightShft,PointPosition[1]-TextDownshift);
	MyText="- load factor ";PrintTheText();
	sprintf(NumericalValue,"%.1f",OwnWtSafetyFact);NumberLength=strlen(NumericalValue);
	for(int MyChar=0;MyChar<NumberLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyChar]);
	
	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;PointPosition[1]+=HalfLineHeight;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+0.5*TextRightShft,PointPosition[1]-TextDownshift);
	MyText="+ Applied load";PrintTheText();
	glBegin(GL_POINTS);
	PointPosition[1]-=2.0*HalfLineHeight;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+0.5*TextRightShft,PointPosition[1]-TextDownshift);
	MyText="- load factor ";PrintTheText();
	sprintf(NumericalValue,"%.1f",AppliedSafetyFact);NumberLength=strlen(NumericalValue);
	for(int MyChar=0;MyChar<NumberLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyChar]);

	glBegin(GL_POINTS);
	PointPosition[0]=RadioLeftx;
	DropButton++;PointPosition[1]=RadioTopy-double(DropButton)*RadioSpacingy;
	glVertex2dv(PointPosition);
	glEnd();
	glRasterPos2d(PointPosition[0]+TextRightShft,PointPosition[1]-TextDownshift);
	MyText="Quit";PrintTheText();
	
	PointPosition[1]-=1.0*RadioSpacingy;
	glRasterPos2d(PointPosition[0],PointPosition[1]);
	MyText="Mean resultant force";PrintTheText();
	PointPosition[1]-=0.75*RadioSpacingy;
	glRasterPos2d(PointPosition[0],PointPosition[1]);
	sprintf(NumericalValue,"%7.1e",SumForce/double(NumberOfActiveNodes));NumberLength=strlen(NumericalValue);
	for(int MyChar=0;MyChar<NumberLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,NumericalValue[MyChar]);

	glPopMatrix();
}

void PrintTheText(void)
{
int StringLength=strlen(MyText);for(int MyChar=0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyChar]);
}

void SavePicture(void)
{
	char   ColourRed[2*WindowHalfWidth];
	char ColourGreen[2*WindowHalfWidth];
	char  ColourBlue[2*WindowHalfWidth];
	
	int pixelx,pixely;
	
	FILE *tga;               // Pointer to a FILE
	targa_header header;     // Variable of targa_header type
	
	/* First, set all the fields in the header to appropriate values */
	header.id_len = 0;          /* no ID field */
	header.map_type = 0;        /* no colormap */
	header.img_type = 2;        /* trust me */
	header.map_first = 0;       /* not used */
	header.map_len = 0;         /* not used */
	header.map_entry_size = 0;  /* not used */
	header.x = 0;               /* image starts at (0,0) */
	header.y = 0;
	header.width = 2*WindowHalfWidth;         
	header.height = 2*WindowHalfHeight;
	header.bpp = 24;            /* 24 bits per pixel */
	header.misc = 0x20;         /* scan from upper left corner */
	
	/* Open a file for writing targa data.  Call the file "test.tga" and
		write in binary mode (wb) so that nothing is lost as characters
		are written to the file */
	
	char FileName[40];
	char  zero[5]="0";
	char   one[5]="1";
	char   two[5]="2";
	char three[5]="3";
	char  four[5]="4";
	char  five[5]="5";
	char   six[5]="6";
	char seven[5]="7";
	char eight[5]="8";
	char  nine[5]="9";
	int NewFileNumber,ThisNumber,TenToThePower,Divisor,Power;
	FileNumber++;
	NewFileNumber=FileNumber;
	strcpy(FileName,"TestNo");
	for(TenToThePower=3;TenToThePower>=0;TenToThePower-=1)
	{
		Divisor=1;
		for(Power=1;Power<=TenToThePower;Power++)Divisor=Divisor*10;
		ThisNumber=NewFileNumber/Divisor;
		NewFileNumber-=ThisNumber*Divisor;
		if(ThisNumber==0)strcat(FileName,zero);
		if(ThisNumber==1)strcat(FileName,one);
		if(ThisNumber==2)strcat(FileName,two);
		if(ThisNumber==3)strcat(FileName,three);
		if(ThisNumber==4)strcat(FileName,four);
		if(ThisNumber==5)strcat(FileName,five);
		if(ThisNumber==6)strcat(FileName,six);
		if(ThisNumber==7)strcat(FileName,seven);
		if(ThisNumber==8)strcat(FileName,eight);
		if(ThisNumber==9)strcat(FileName,nine);
	}
	strcat (FileName,".tga");
	
	tga = fopen(FileName, "wb"); /* Write the header information  */
				
				writeheader(header, tga);  
				
				for(pixely=2*WindowHalfHeight-1;pixely>=0;pixely-=1)
				{
					glReadPixels(0,pixely,2*WindowHalfWidth,1,GL_BLUE,  GL_UNSIGNED_BYTE,&ColourBlue);
					glReadPixels(0,pixely,2*WindowHalfWidth,1,GL_GREEN, GL_UNSIGNED_BYTE,&ColourGreen);
					glReadPixels(0,pixely,2*WindowHalfWidth,1,GL_RED,   GL_UNSIGNED_BYTE,&ColourRed);
					
					for(pixelx=0;pixelx<=2*WindowHalfWidth-1;pixelx++)
					{
						fputc(ColourBlue[pixelx], tga);
						fputc(ColourGreen[pixelx],tga);
						fputc(ColourRed[pixelx],  tga);
					}
				}
				fclose(tga);
}

void writeheader(targa_header h, FILE *tga) 
{
	fputc(h.id_len, tga);          // Write chars for ID, map, and image type
	fputc(h.map_type, tga);
	fputc(h.img_type, tga);
	fputc(h.map_first % 256, tga); // Write integer, low order byte first
	fputc(h.map_first / 256, tga); // Write second byte of integer, high order
	fputc(h.map_len % 256, tga);   // Another integer 
	fputc(h.map_len / 256, tga);
	fputc(h.map_entry_size, tga);  // Write a char - only one byte
	fputc(h.x % 256, tga);         // More integers
	fputc(h.x / 256, tga);
	fputc(h.y % 256, tga);
	fputc(h.y / 256, tga);
	fputc(h.width % 256, tga);     // Even more integers
	fputc(h.width / 256, tga);
	fputc(h.height % 256, tga);
	fputc(h.height / 256, tga);
	fputc(h.bpp, tga);             // Write two chars
	fputc(h.misc, tga);
}
