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

//#include "SavePicture.h" //Necessary only if you plan to save images, for example to make a movie. See SavePicture(..,..) below

#define   LastNode       1000
#define MaxLastZone0      27
#define MaxLastZone1      16
#define MaxLastInZone    500
#define LastPrevious      5
void Calculation(void);
void SortIntoZones(void);
void Contribution(int,int);
void Invert(int);
void MaudGraphics(float,float,float);
void MouseClick(int, int, int, int);
void MouseDrag(int x, int y);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
void WriteDXF(void);
double LineLengthSq(int,int);

int argc;
char **argv;
GLUquadricObj *quadObj;
bool LeftMouseDown=0,RightMouseDown=0;

int   LastInZone[MaxLastZone0+1][MaxLastZone1+1],
InZone[MaxLastZone0+1][MaxLastZone1+1][MaxLastInZone+1],
Zone[2],LastZone[2],WhichZone[2],HalfDimension[2],
TensorMass,ShowNodes,ShowLines,clickx,clicky,WhiteBackGround=1,ForceType=0;

float x[2][LastNode+1],Force[2][LastNode+1],Velocity[2][LastNode+1],xPrevious[2][LastNode+1][LastPrevious+1],
Stiffness[2][2][LastNode+1],Flex[2][2],deltax[2],HalfPicture[2],ZoneSize[2],
PI,a,aSq,ToverL,dTbydL,LSq,tempForce,tempStiffness,Slow,CarryOver,CentralAttraction,Bounce,L,MyExp,MinZoneSizeSq,
Zoom;

int main(int argc, char* argv[]){
	glutInit(&argc,argv);
	int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	
	cout<<"Your screen resolution is "<<screenWidth<<" x "<<screenHeight<<"\n";
	
	HalfDimension[0] =int(double(screenWidth)/2.0)-5;
	HalfDimension[1]=int(double(screenHeight)/2.0)-25;
	
	for(int i=0;i<=1;i++)HalfPicture[i]=1.0*HalfDimension[i]-25.0;
	cout<<"This program runs for ever until the window is closed or 'q' is pressed on the keyboard.\n";
	cout<<"Typing 'n' makes the nodes dissapear or reappear.\n";
	cout<<"Typing 'l' makes the nodes dissapear or reappear.\n";
	cout<<"Typing 't' makes the masses sclara or 'tensor'.\n";
	PI=4.0*atan(1.0);srand(time(NULL));
	Slow=0.0002;CarryOver=0.96;CentralAttraction=1.0/200000.0;Bounce=0.95;
	a=25.0;aSq=a*a;
	ShowNodes=1;ShowLines=0;
	TensorMass=1;
	for(int i=0;i<=1;i++)
	{
		LastZone[i]=int(2.0*HalfDimension[i]/(3.0*a));
		if((i==0&&LastZone[i]>MaxLastZone0)||(i==1&&LastZone[i]>MaxLastZone1))
		{
			cout<<LastZone[i]<<" is too many zones in the "<<i<<" direction\n";return 0;
		}
		ZoneSize[i]=2.0*HalfDimension[i]/(1.0*LastZone[i]);
		cout<<"ZoneSize["<<i<<"]/a ="<<ZoneSize[i]/a<<"\n";
	}
	cout<<"Number of zones is "<<LastZone[0]+1<<" x "<<LastZone[1]+1<<"\n";
	if(ZoneSize[0]<ZoneSize[1])MinZoneSizeSq=ZoneSize[0]*ZoneSize[0];
	else MinZoneSizeSq=ZoneSize[1]*ZoneSize[1];
	cout<<"Value of exponential at edge of zone = "<<exp(-MinZoneSizeSq/aSq)<<"\n";
	for(int Node=0;Node<=LastNode;Node++)
	{
		for(int i=0;i<=1;i++)
		{
			x[i][Node]=0.75*HalfPicture[i]*((2.0*rand())/(1.0*RAND_MAX)-1.0);
			Velocity[i][Node]=0.0;
			for(int Previous=0;Previous<=LastPrevious;Previous++)xPrevious[i][Node][Previous]=x[i][Node];
		}
	}
	if(WhiteBackGround==1)MaudGraphics(1.0,1.0,1.0);else MaudGraphics(0.0,0.0,0.0);
	return 0;
}
static void Draw(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	Calculation();
	glutSwapBuffers();
}
void Calculation(void)
{
	if(WhiteBackGround==1)glColor4f(0.0,0.0,0.0,0.5);else glColor4f(1.0,1.0,1.0,0.5);
	glRasterPos2f(-HalfDimension[0]+10,-HalfDimension[1]+20);
	char *MyText;
	int LengthOfString;
	MyText="Maud Algorithm - Chris Williams, Bath University";
	LengthOfString=strlen(MyText);
	for (int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
	
	glColor4f(1.0,0.0,0.0,0.5);
	glRasterPos2f(-HalfDimension[0]+10,-HalfDimension[1]+50);
	if(TensorMass==1)MyText="Tensor masses";else MyText="Scalar masses";
	LengthOfString=strlen(MyText);
	for (int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
	
	
	for(int Node=0;Node<=LastNode;Node++)
	{
		for(int i=0;i<=1;i++)
		{
			Force[i][Node]=0.0;
			for(int j=0;j<=1;j++)Stiffness[i][j][Node]=0.0;
		}
	}
	SortIntoZones();
	for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++)
	{
		for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++)
		{
			for(int NodeInZone=0;NodeInZone<=LastInZone[Zone[0]][Zone[1]];NodeInZone++)
			{
				int Node=InZone[Zone[0]][Zone[1]][NodeInZone];
				for(WhichZone[0]=Zone[0]-1;WhichZone[0]<=Zone[0]+1;WhichZone[0]++)
				{
					if(0<=WhichZone[0]&&WhichZone[0]<=LastZone[0])
					{
						for(WhichZone[1]=Zone[1]-1;WhichZone[1]<=Zone[1]+1;WhichZone[1]++)
						{
							if(0<=WhichZone[1]&&WhichZone[1]<=LastZone[1])
							{
								for(int OtherNodeInZone=0;OtherNodeInZone<=LastInZone[WhichZone[0]][WhichZone[1]];OtherNodeInZone++)
								{
									int OtherNode=InZone[WhichZone[0]][WhichZone[1]][OtherNodeInZone];
									if(OtherNode<Node)
									{
										LSq=0.0;
										for(int i=0;i<=1;i++)
										{
											deltax[i]=x[i][OtherNode]-x[i][Node];
											LSq+=deltax[i]*deltax[i];
										}
										if(LSq/aSq>1.0e-24&&LSq<MinZoneSizeSq)
										{
											Contribution(Node,OtherNode);
											if(ShowLines==1)
											{
												glColor4f(MyExp,0.0,1.0-MyExp,MyExp);
												glBegin(GL_LINES);
												glVertex2f(x[0][Node],x[1][Node]);
												glVertex2f(x[0][OtherNode],x[1][OtherNode]);
												glEnd();
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	for(int Node=0;Node<=LastNode;Node++)
	{
		glBegin(GL_LINES);
		for(int Previous=0;Previous<=LastPrevious;Previous++)
		{
			float factor=float(LastPrevious+1-Previous)/float(LastPrevious+1);
			factor=factor*factor;
			if(WhiteBackGround==1)glColor4f(0.0,0.0,0.0,factor);
			else glColor4f(1.0,1.0,1.0,factor);
			if(Previous==0)glVertex2f(x[0][Node],x[1][Node]);
			else glVertex2f(xPrevious[0][Node][Previous-1],xPrevious[1][Node][Previous-1]);
			glVertex2f(xPrevious[0][Node][Previous],xPrevious[1][Node][Previous]);
		}
		glEnd();
		if(ShowNodes==1)
		{
			if(WhiteBackGround==1)glColor4f(0.0,0.0,0.0,1.0);
			else glColor4f(1.0,1.0,1.0,1.0);
			glBegin(GL_POINTS);
			glVertex2f(x[0][Node],x[1][Node]);
			glEnd();
		}
		for(int Previous=LastPrevious;Previous>=0;Previous--)
		{
			for(int i=0;i<=1;i++)
			{
				if(Previous==0)xPrevious[i][Node][0]=x[i][Node];
				else xPrevious[i][Node][Previous]=xPrevious[i][Node][Previous-1];
			}
		}
	}
	//SavePicture(HalfDimension[0],HalfDimension[1]);//Necessary only if you plan to save images, for example to make a movie.
	for(int Node=0;Node<=LastNode;Node++)
	{
		Invert(Node);
		for(int i=0;i<=1;i++)
		{
			Velocity[i][Node]=CarryOver*Velocity[i][Node]-CentralAttraction*x[i][Node];
			for(int j=0;j<=1;j++)Velocity[i][Node]+=Slow*Flex[i][j]*Force[j][Node];
			x[i][Node]+=Velocity[i][Node];
		}
		for(int i=0;i<=1;i++)
		{
			if(x[i][Node]>HalfPicture[i])
			{
				x[i][Node]=HalfPicture[i];
				if(Velocity[i][Node]>0.0)Velocity[i][Node]=-Bounce*Velocity[i][Node];
			}
			if(x[i][Node]<-HalfPicture[i])
			{
				x[i][Node]=-HalfPicture[i];
				if(Velocity[i][Node]<0.0)Velocity[i][Node]=-Bounce*Velocity[i][Node];
			}
		}
	}
}
void MaudGraphics(float BackGroundRed,float BackGroundGreen,float BackGroundBlue)
{
    glutInitWindowSize(2*HalfDimension[0], 2*HalfDimension[1]);
	glutInitWindowPosition(0,0);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Maud");
	glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
	glViewport(0, 0, 2*HalfDimension[0], 2*HalfDimension[1]);
	gluOrtho2D(-HalfDimension[0],HalfDimension[0],-HalfDimension[1],HalfDimension[1]);
	glClear(GL_COLOR_BUFFER_BIT);
	glLineWidth(0.5);
	glPointSize(3.0);
	glutKeyboardFunc(Key);
	glutIdleFunc(Draw);
	glutDisplayFunc(Draw);
	glutMouseFunc(MouseClick);
	glutMotionFunc(MouseDrag);
	glEnable(GL_POINT_SMOOTH);
	glEnable (GL_BLEND);glBlendFunc (GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	Zoom=1.0;
	glutMainLoop();	
}
static void Key(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'n':if(ShowNodes==0)ShowNodes=1;else ShowNodes=0;break;
			case 'l':if(ShowLines==0)ShowLines=1;else ShowLines=0;break;
			case 't':if(TensorMass==0)TensorMass=1;else TensorMass=0;break;
			case 'q':gluDeleteQuadric(quadObj);WriteDXF();exit(0);
			case '\033':gluDeleteQuadric(quadObj);WriteDXF();exit(0);
	}
}
void SortIntoZones(void)
{
	for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++)
	{
		for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++)
			LastInZone[Zone[0]][Zone[1]]=-1;
	}
	for(int Node=0;Node<=LastNode;Node++)
	{
		for(int i=0;i<=1;i++)
		{
			Zone[i]=int((x[i][Node]+HalfDimension[i])/ZoneSize[i]);
			if(Zone[i]<0)Zone[i]=0;
			if(Zone[i]>LastZone[i])Zone[i]=LastZone[i];
		}
		LastInZone[Zone[0]][Zone[1]]++;
		if(LastInZone[Zone[0]][Zone[1]]<=MaxLastInZone)InZone[Zone[0]][Zone[1]][LastInZone[Zone[0]][Zone[1]]]=Node;
	}
	for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++)
	{
		for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++)
		{
			if(LastInZone[Zone[0]][Zone[1]]>MaxLastInZone)
			{
				cout<<"Number in Zone "<<Zone[0]<<" "<<Zone[1]<<" is "<<LastInZone[Zone[0]][Zone[1]]<<"\n";
				LastInZone[Zone[0]][Zone[1]]=MaxLastInZone;
			}
		}
	}
}
void Contribution(int Node,int OtherNode)
{
	L=sqrt(LSq);
	MyExp=exp(-LSq/aSq);
	if(ForceType==0)
	{
		ToverL=-(aSq/LSq)*MyExp;
		dTbydL=(aSq/LSq)*(1.0+2.0*LSq/aSq)*MyExp;
	}
	else
	{
		ToverL=-(a/L)*MyExp;
		dTbydL=(a/L)*(2.0*LSq/aSq)*MyExp;
	}
	
	for(int i=0;i<=1;i++)
	{
		tempForce=ToverL*deltax[i];
		Force[i][Node]+=tempForce;
		Force[i][OtherNode]-=tempForce;
		Stiffness[i][i][Node]+=ToverL;
		Stiffness[i][i][OtherNode]+=ToverL;
		for(int j=0;j<=1;j++)
		{
			tempStiffness=(dTbydL-ToverL)*deltax[i]*deltax[j]/LSq;
			Stiffness[i][j][Node]+=tempStiffness;
			Stiffness[i][j][OtherNode]+=tempStiffness;
		}
	}
}
void Invert(int Node)
{
	if(TensorMass==1)
	{
		float determinant;
		determinant=Stiffness[0][0][Node]*Stiffness[1][1][Node]
		-Stiffness[0][1][Node]*Stiffness[1][0][Node];
		if(fabs(determinant)>1.0e-48)
		{
			Flex[0][0]=Stiffness[1][1][Node]/determinant;
			Flex[1][1]=Stiffness[0][0][Node]/determinant;
			Flex[0][1]=-Stiffness[0][1][Node]/determinant;
			Flex[1][0]=-Stiffness[1][0][Node]/determinant;
		}
		else
		{
			Flex[0][0]=0.0;
			Flex[1][1]=0.0;
			Flex[0][1]=0.0;
			Flex[1][0]=0.0;
		}
	}
	else
	{
		Flex[0][0]=1.0/(Stiffness[0][0][Node]+Stiffness[1][1][Node]);
		Flex[1][1]=Flex[0][0];
		Flex[0][1]=0.0;
		Flex[1][0]=0.0;
	}
}

void MouseDrag(int x, int y)
{
	float ZoomIncrement;
	if(RightMouseDown)
	{
		ZoomIncrement=(1.0+0.001*float(x-clickx))*(1.0-0.001*float(y-clicky));
		glScalef(ZoomIncrement,ZoomIncrement,ZoomIncrement);
		Zoom=Zoom*ZoomIncrement;
		clickx=x;clicky=y;
	}
	if(LeftMouseDown)
	{
		float xTransInc=+float(x-clickx)/Zoom;
		float yTransInc=-float(y-clicky)/Zoom;
		glTranslatef(xTransInc,yTransInc,0.0);
		clickx=x;clicky=y;
	}
}

void MouseClick(int button, int state, int x, int y)
{
	if((button==GLUT_LEFT_BUTTON)&&(state==GLUT_DOWN)){LeftMouseDown=1;clickx=x;clicky=y;}
	if((button==GLUT_LEFT_BUTTON)&&(state==GLUT_UP))LeftMouseDown=0;
	if((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_DOWN)){RightMouseDown=1;clickx=x;clicky=y;}
	if((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_UP))RightMouseDown=0;
}

void WriteDXF(void)
{
	int dxfGreyMap[6];
	dxfGreyMap[0]=254;
	dxfGreyMap[1]=253;
	dxfGreyMap[2]=9;
	dxfGreyMap[3]=8;
	dxfGreyMap[4]=250;
	dxfGreyMap[5]=19;
	double MaxLineLengthSq=0.0;
	for(int Node=0;Node<=LastNode-1;Node++){
		for(int OtherNode=Node+1;OtherNode<=LastNode;OtherNode++){
		if(MaxLineLengthSq<LineLengthSq(Node,OtherNode))MaxLineLengthSq=LineLengthSq(Node,OtherNode);}}
	double MaxLineLength=sqrt(MaxLineLengthSq);
	double PointRadius=MaxLineLength/500.0;
	ofstream Julia;
	Julia.open("Maud.dxf");
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	for(int Node=0;Node<=LastNode-1;Node++){
		for(int OtherNode=Node+1;OtherNode<=LastNode;OtherNode++){
			int LineControl=5-int(20.0*sqrt(LineLengthSq(Node,OtherNode))/MaxLineLength);
			if(LineControl>5)LineControl=5;
			if(LineControl>=0){
				Julia<<"0\nLINE\n8\nLines\n";
				Julia<<"10\n"<<x[0][Node]     <<"\n20\n"<<x[1][Node]     <<"\n";
				Julia<<"11\n"<<x[0][OtherNode]<<"\n21\n"<<x[1][OtherNode]<<"\n";
			Julia<<"62\n"<<dxfGreyMap[LineControl]<<"\n";}}}
	for(int Node=0;Node<=LastNode-1;Node++){
		Julia<<"0\nCIRCLE\n8\nPoints\n";
		Julia<<"10\n"<<x[0][Node]     <<"\n20\n"<<x[1][Node]     <<"\n";
		Julia<<"40\n"<<PointRadius<<"\n";
	Julia<<"62\n"<<1<<"\n";}
	Julia<<"0\nENDSEC\n0\nEOF\n";
	Julia.close();
}

double LineLengthSq(int Node,int OtherNode)
{
	return
	(x[0][OtherNode]-x[0][Node])*(x[0][OtherNode]-x[0][Node])+
	(x[1][OtherNode]-x[1][Node])*(x[1][OtherNode]-x[1][Node]);
}
