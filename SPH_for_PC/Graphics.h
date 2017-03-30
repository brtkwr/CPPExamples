#include "MacSunPC.h"

#define HalfWidth     640
#define HalfHeight    375

void ChrisGraphics(float BackGroundRed,float BackGroundGreen,float BackGroundBlue);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
void MouseClick(int button, int state, int x, int y);
void MouseDrag(int x, int y);
void MoveOrigin(void);
void SaveState(void);
void SavePicture(void);

GLUquadricObj *quadObj;

float  InitialZoom,InitialZoomNoInfo,Zoom,ZoomIncrement,xTrans,xTransInc,yTrans,yTransInc;

int argc,clickx,clicky,FileNumber,VorticityColour,Streak,Information;

char **argv;

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

void ChrisGraphics(float BackGroundRed,float BackGroundGreen,float BackGroundBlue)
{
	glutInitWindowSize(2*HalfWidth, 2*HalfHeight);
	if(Macintosh==1)    glutInit(&argc, argv);               // Remove this line for PC or Sun
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("2D Flow");
	glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
	glViewport(0, 0, 2*HalfWidth, 2*HalfHeight);
	gluOrtho2D(-HalfWidth,HalfWidth,-HalfHeight,HalfHeight);
	glClear(GL_COLOR_BUFFER_BIT);
	glLineWidth(1.0);
	glPointSize(2.0);
	glutKeyboardFunc(Key);
	glutIdleFunc(Draw);
	glutDisplayFunc(Draw);
	glutMouseFunc(MouseClick);
	glutMotionFunc(MouseDrag);
	glEnable(GL_POINT_SMOOTH);// glDisable(GL_POINT_SMOOTH) is the reverse. Enable means antialiasing 
	glEnable (GL_BLEND);glBlendFunc (GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glScalef(Zoom,Zoom,Zoom);
	glutMainLoop();
}
static void Key(unsigned char key, int x, int y)
{
switch (key)
{
case '+':ZoomIncrement=sqrt(2.0)/1.0;glScalef(ZoomIncrement,ZoomIncrement,ZoomIncrement);Zoom=Zoom*ZoomIncrement;break;
case '-':ZoomIncrement=1.0/sqrt(2.0);glScalef(ZoomIncrement,ZoomIncrement,ZoomIncrement);Zoom=Zoom*ZoomIncrement;break;
case 'r':xTransInc=+10;MoveOrigin();break;
case 'l':xTransInc=-10;MoveOrigin();break;
case 'c':if(VorticityColour==1)VorticityColour=0;else VorticityColour=1;break;
case 's':if(Streak==1)Streak=0;else Streak=1;break;
case 'i':if(Information==1){Information=0;glScalef(InitialZoomNoInfo/Zoom,InitialZoomNoInfo/Zoom,InitialZoomNoInfo/Zoom);Zoom=InitialZoomNoInfo;}
					  else {Information=1;glScalef(InitialZoom/Zoom,InitialZoom/Zoom,InitialZoom/Zoom);Zoom=InitialZoom;}
					  break;
case 'f':SavePicture();break;
case 'o':if(Information==1){glScalef(InitialZoom/Zoom,InitialZoom/Zoom,InitialZoom/Zoom);Zoom=InitialZoom;}
					  else {glScalef(InitialZoomNoInfo/Zoom,InitialZoomNoInfo/Zoom,InitialZoomNoInfo/Zoom);Zoom=InitialZoomNoInfo;}
                      break;
case 't':xTransInc=-xTrans;
		 yTransInc=-yTrans;MoveOrigin();break;
case 'q':SaveState();gluDeleteQuadric(quadObj);exit(0);
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
float TransInc;

	if (LeftMouseDown)
    {
		TransInc=+static_cast<float>(x-clickx)/Zoom;
		xTransInc=TransInc;
		yTransInc=0.0;
		MoveOrigin();
		
		TransInc=-static_cast<float>(y-clicky)/Zoom;
		xTransInc=0.0;
		yTransInc=TransInc;		
		MoveOrigin();
    	clickx=x;
		clicky=y;
   }
	if (RightMouseDown)
   {
		ZoomIncrement=(1.0+0.001*static_cast<float>(x-clickx))*
		              (1.0-0.001*static_cast<float>(y-clicky));
		glScalef(ZoomIncrement,ZoomIncrement,ZoomIncrement);
		Zoom=Zoom*ZoomIncrement;

		clickx=x;
		clicky=y;
	}
	}

void MoveOrigin(void)
{
		if(Information==1)
		{
		xTrans+=xTransInc*Zoom/InitialZoom;
		yTrans+=yTransInc*Zoom/InitialZoom;
		}
		else
		{
		xTrans+=xTransInc*Zoom/InitialZoomNoInfo;
		yTrans+=yTransInc*Zoom/InitialZoomNoInfo;
		}
		glTranslatef(xTransInc,yTransInc,0.0);
}

void SavePicture(void)
{
char   ColourRed[2*HalfWidth];
char ColourGreen[2*HalfWidth];
char  ColourBlue[2*HalfWidth];

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
   header.width = 2*HalfWidth;         
   header.height = 2*HalfHeight;
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

for(pixely=2*HalfHeight-1;pixely>=0;pixely-=1)
{
glReadPixels(0,pixely,2*HalfWidth,1,GL_BLUE,  GL_UNSIGNED_BYTE,&ColourBlue);
glReadPixels(0,pixely,2*HalfWidth,1,GL_GREEN, GL_UNSIGNED_BYTE,&ColourGreen);
glReadPixels(0,pixely,2*HalfWidth,1,GL_RED,   GL_UNSIGNED_BYTE,&ColourRed);

for(pixelx=0;pixelx<=2*HalfWidth-1;pixelx++)
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
