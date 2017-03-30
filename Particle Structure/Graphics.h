//The following lines determine from where to get glut.h for Macintosh, Windows etc.

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define MaxHalfWidth 2000
int HalfWidth,HalfHeight;

void SetUpGraphics(float BackGroundRed,float BackGroundGreen,float BackGroundBlue);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
void MouseClick(int button, int state, int x, int y);
void MouseDrag(int x, int y);
void MoveOrigin(void);
void Rotate(void);
void SaveState(void);
void SavePicture(void);

GLUquadricObj *quadObj;

float AspectRatio,BorderHalfWidth,InitialZoom,Zoom,ZoomIncrement,xTrans,xTransInc,yTrans,yTransInc;

float zRotationIncrement,HorizontalRotationIncrement,zRotation,HorizontalRotation,cosz,cosH,sinz,sinH,PI;

int argc,clickx,clicky,FileNumber,VorticityColour,ShowArms,DrawGrid,DrawBoundary,ZoomOrRotate = 0;

char **argv;

bool LeftMouseDown = 0,RightMouseDown = 0;

typedef struct {
	char id_len;   // ID Field (Number of bytes - max 255)
	char map_type;  // Colormap Field (0 or 1)
	char img_type;  // Image Type (7 options - color vs. compression)
	int map_first;  // Color Map stuff - first entry index
	int map_len;  // Color Map stuff - total entries in file
	char map_entry_size;  // Color Map stuff - number of bits per entry
	int x;   // X-coordinate of origin 
	int y;   // Y-coordinate of origin
	int width;   // Width in Pixels
	int height;   // Height in Pixels
	char bpp;   // Number of bits per pixel
	char misc;   // Other stuff - scan origin and alpha bits
} targa_header;

void writeheader(targa_header h, FILE *tga);

void SetUpGraphics(float BackGroundRed,float BackGroundGreen,float BackGroundBlue)
{
	
 int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	
	cout<<"Your screen resolution is "<<screenWidth<<" x "<<screenHeight<<"\n";
	
	float Hori = float(screenWidth)-5.0;
	float Vert = float(screenHeight)-20.0;
	if(Hori<AspectRatio*Vert)HalfWidth = int(Hori/2.0);else HalfWidth = int(AspectRatio*Vert/2.0);
	if(HalfWidth>MaxHalfWidth)
	{
		cout<<"Screen half width = "<<HalfWidth<<" which is greater than the maximum allowed which is "<<MaxHalfWidth<<"\n";
		HalfWidth = MaxHalfWidth;
		cout<<"Screen half width set to this value\n";
	}
	HalfHeight = int(float(HalfWidth)/AspectRatio);
	
	VorticityColour = 1;ShowArms = 1;DrawGrid = 0;DrawBoundary = 0;
	
	InitialZoom = 0.95*float(HalfWidth)/BorderHalfWidth;
	//InitialZoomNoInfo = InitialZoom;
	//if(Information == 1)
	Zoom = InitialZoom;//else Zoom = InitialZoomNoInfo;
	
	glutInitWindowSize(2*HalfWidth, 2*HalfHeight);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Particle structure");
	glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
	glViewport(0, 0, 2*HalfWidth, 2*HalfHeight);
	glOrtho(-HalfWidth,HalfWidth,-HalfHeight,HalfHeight,-HalfHeight,HalfHeight);
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
	PI = 4.0 * atan(1.0);
	Rotate();
	glutMainLoop();
}
static void Key(unsigned char key, int x, int y)
{
	switch (key)
	{
		case '+':ZoomIncrement = sqrt(2.0)/1.0;glScalef(ZoomIncrement,ZoomIncrement,ZoomIncrement);Zoom = Zoom*ZoomIncrement;break;
		case '-':ZoomIncrement = 1.0/sqrt(2.0);glScalef(ZoomIncrement,ZoomIncrement,ZoomIncrement);Zoom = Zoom*ZoomIncrement;break;
		case 'r':xTransInc = +10;MoveOrigin();break;
		case 'l':xTransInc = -10;MoveOrigin();break;
		case 'z':if(ZoomOrRotate == 1)ZoomOrRotate = 0;else ZoomOrRotate = 1;break;
			case 'c':if(VorticityColour == 1)VorticityColour = 0;else VorticityColour = 1;break;
			case 'h':xTransInc = -xTrans;yTransInc = -yTrans;MoveOrigin();break;
			case 'a':ShowArms ++;if(ShowArms > 2) ShowArms = 0;break;
			case 'g':if(DrawGrid == 0)DrawGrid = 1;else DrawGrid = 0;break;
			case 'q':SaveState();gluDeleteQuadric(quadObj);exit(0);
			case 'Q':SaveState();gluDeleteQuadric(quadObj);exit(0);
			case '\033':SaveState();gluDeleteQuadric(quadObj);exit(0);
	}
}

void MouseClick(int button, int state, int x, int y) 
{
	if ((button == GLUT_LEFT_BUTTON)&&(state == GLUT_DOWN)) {
		clickx = x;
		clicky = y;
		LeftMouseDown = 1;
	}
	if ((button == GLUT_LEFT_BUTTON)&&(state == GLUT_UP)) {
		LeftMouseDown = 0;
	}
	
	if ((button == GLUT_RIGHT_BUTTON)&&(state == GLUT_DOWN)) {
		clickx = x;
		clicky = y;
		RightMouseDown = 1;
	}
	if ((button == GLUT_RIGHT_BUTTON)&&(state == GLUT_UP)) {
		RightMouseDown = 0;
	}
}

void MouseDrag(int x, int y) 
{
	float TransInc;
	
	if (LeftMouseDown)
 {
		TransInc = +static_cast<float>(x-clickx)/Zoom;
		xTransInc = TransInc;
		yTransInc = 0.0;
		MoveOrigin();
		
		TransInc = -static_cast<float>(y-clicky)/Zoom;
		xTransInc = 0.0;
		yTransInc = TransInc;		
		MoveOrigin();
 	clickx = x;
		clicky = y;
	}
	if (RightMouseDown)
	{
		if(ZoomOrRotate == 0)
		{
			//float X = x-float(HalfWidth);
			float Y = y-float(HalfHeight);
			
			zRotationIncrement = float(x-clickx) * Y / 5000.0;
			HorizontalRotationIncrement = float(y-clicky) / 2.0;
			Rotate();
			
			
		}
		
		else
		{
			ZoomIncrement = (1.0+0.001*static_cast<float>(x-clickx))*
			(1.0-0.001*static_cast<float>(y-clicky));
			glScaled(ZoomIncrement,ZoomIncrement,ZoomIncrement);
			Zoom = Zoom*ZoomIncrement;
		}
		clickx = x;
		clicky = y;
	}
}

void Rotate(void)
{
	zRotation += zRotationIncrement;
	HorizontalRotation += HorizontalRotationIncrement;
	cosz = cos(zRotation*PI/180.0);cosH = cos(HorizontalRotation*PI/180.0);
	sinz = sin(zRotation*PI/180.0);sinH = sin(HorizontalRotation*PI/180.0);
	glRotatef(zRotationIncrement,0.0,0.0,1.0);
	glRotatef(HorizontalRotationIncrement,cosz,-sinz,0.0);
}

void MoveOrigin(void)
{
	xTrans += xTransInc*Zoom/InitialZoom;
	yTrans += yTransInc*Zoom/InitialZoom;
	glTranslatef(xTransInc,yTransInc,0.0);
}

void SavePicture(void)
{
	char ColourRed[2*MaxHalfWidth];
	char ColourGreen[2*MaxHalfWidth];
	char ColourBlue[2*MaxHalfWidth];
	
	int pixelx,pixely;
	
	FILE *tga;  // Pointer to a FILE
	targa_header header; // Variable of targa_header type
	
	/* First, set all the fields in the header to appropriate values */
	header.id_len = 0;  /* no ID field */
	header.map_type = 0; /* no colormap */
	header.img_type = 2; /* trust me */
	header.map_first = 0; /* not used */
	header.map_len = 0;  /* not used */
	header.map_entry_size = 0; /* not used */
	header.x = 0;  /* image starts at (0,0) */
	header.y = 0;
	header.width = 2*HalfWidth;  
	header.height = 2*HalfHeight;
	header.bpp = 24;  /* 24 bits per pixel */
	header.misc = 0x20;  /* scan from upper left corner */
	
	/* Open a file for writing targa data. Call the file "test.tga" and
	 write in binary mode (wb) so that nothing is lost as characters
	 are written to the file */
	
	char FileName[40];
	char zero[5] = "0";
	char one[5] = "1";
	char two[5] = "2";
	char three[5] = "3";
	char four[5] = "4";
	char five[5] = "5";
	char six[5] = "6";
	char seven[5] = "7";
	char eight[5] = "8";
	char nine[5] = "9";
	int NewFileNumber,ThisNumber,TenToThePower,Divisor,Power;
	FileNumber++;
	NewFileNumber = FileNumber;
	strcpy(FileName,"Fibres");
	for(TenToThePower = 3;TenToThePower >= 0;TenToThePower -= 1)
	{
		Divisor = 1;
		for(Power = 1;Power <= TenToThePower;Power++)Divisor = Divisor*10;
		ThisNumber = NewFileNumber/Divisor;
		NewFileNumber -= ThisNumber*Divisor;
		if(ThisNumber == 0)strcat(FileName,zero);
		if(ThisNumber == 1)strcat(FileName,one);
		if(ThisNumber == 2)strcat(FileName,two);
		if(ThisNumber == 3)strcat(FileName,three);
		if(ThisNumber == 4)strcat(FileName,four);
		if(ThisNumber == 5)strcat(FileName,five);
		if(ThisNumber == 6)strcat(FileName,six);
		if(ThisNumber == 7)strcat(FileName,seven);
		if(ThisNumber == 8)strcat(FileName,eight);
		if(ThisNumber == 9)strcat(FileName,nine);
	}
	strcat (FileName,".tga");
	
	tga = fopen(FileName, "wb"); /* Write the header information */
	
	writeheader(header, tga); 
	
	for(pixely = 2*HalfHeight-1;pixely >= 0;pixely -= 1)
	{
		glReadPixels(0,pixely,2*HalfWidth,1,GL_BLUE, GL_UNSIGNED_BYTE,&ColourBlue);
		glReadPixels(0,pixely,2*HalfWidth,1,GL_GREEN, GL_UNSIGNED_BYTE,&ColourGreen);
		glReadPixels(0,pixely,2*HalfWidth,1,GL_RED, GL_UNSIGNED_BYTE,&ColourRed);
		
		for(pixelx = 0;pixelx <= 2*HalfWidth-1;pixelx++)
		{
			fputc(ColourBlue[pixelx], tga);
			fputc(ColourGreen[pixelx],tga);
			fputc(ColourRed[pixelx], tga);
		}
	}
	fclose(tga);
}

void writeheader(targa_header h, FILE *tga) 
{
	fputc(h.id_len, tga);  // Write chars for ID, map, and image type
	fputc(h.map_type, tga);
	fputc(h.img_type, tga);
	fputc(h.map_first % 256, tga); // Write integer, low order byte first
	fputc(h.map_first / 256, tga); // Write second byte of integer, high order
	fputc(h.map_len % 256, tga); // Another integer 
	fputc(h.map_len / 256, tga);
	fputc(h.map_entry_size, tga); // Write a char - only one byte
	fputc(h.x % 256, tga);  // More integers
	fputc(h.x / 256, tga);
	fputc(h.y % 256, tga);
	fputc(h.y / 256, tga);
	fputc(h.width % 256, tga); // Even more integers
	fputc(h.width / 256, tga);
	fputc(h.height % 256, tga);
	fputc(h.height / 256, tga);
	fputc(h.bpp, tga);  // Write two chars
	fputc(h.misc, tga);
}
