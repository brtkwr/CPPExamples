#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#include <GLUT/glut.h>  // For Macintosh
//#include <GL/glut.h>    // For PC or Sun
#define Macintosh 1     // 1 for Macintosh, any other integer for PC

#define    halfW        630
#define    halfH        370
#define LastPerson    80000
#define MaxLastMaudZone   50
#define MaxLastInMaudZone 500

void MaudGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen);
void WriteDXFFile(void);

static void Draw(void);
static void Key(unsigned char key, int x, int y);
int argc;
char **argv;

float rv[2][2];	
GLUquadricObj *quadObj;

void StartingCoords(void);
void Collision(void);
void BounceOff(void);

double SizeDist(void);
double MaudRandom(void);

float Coord[2][LastPerson+1],Vely[2][LastPerson+1],
      m[LastPerson+1],r[LastPerson+1],
      Lower[2][MaxLastMaudZone+1],Upper[2][MaxLastMaudZone+1],
	  deltaCoord[2],deltaVely[2],
	  Overlap[2],
      MaudVertex[2],PreviousCoord[2],
	  PI,SeparationSq,Separation,OverlapFactor,MeanDiameter,
	  deltat,Walking,EWalking,NWalking,restitution,Drag,InOrOut,Stiffness,fineness,
	  ShuteAngle,StartAngle,CriticalStartAngle,Theta,ThetaOver2,MassSum,bOver2;

int    ColourOrHidden[LastPerson+1],InMaudZone[MaxLastInMaudZone+1],
       LastInMaudZone,MaudZone[2],LastMaudZone[2],
       Person,LastPersonUsed,otherPerson,PersonInMaudZone,otherPersonInMaudZone,
	   i,j,SwitchDirection,stop,NewPerson,NewPersonsPerCycle,FoundOne,
	   LengthOfString,MyCharacter,tracks;

char NumericalValue[101];
char *MyText;

ofstream Maud("People.dxf");

int main(void)
{
PI=4.0*atan(1.0);

/*cout<<"This program runs forever until the letter 'q' is pressed on the keyboard. ";
cout<<"A dxf file is then written and the program stops for good.\n\n";
cout<<"'s' will stop and restart the program.\n\n";
cout<<"Type the angle between the two flows in degrees. Must be in the range 0 to 180. Then press 'return'.\n\n";
for(;;)
{
cin>>Theta;
if(0.0<=Theta&&Theta<=180.0)break;
cout<<"Must be in the range 0 to 180\n";
}
Theta=Theta*PI/180.0;
cout<<"\nType a number to control how many people. The maximum value is 2.0 which gives a ";
cout<<"large number, but is slow. 0.5 is recommended. Then press 'return'.\n\n";
for(;;)
{
cin>>fineness;if(0.0<fineness&&fineness<=2.0)break;
else cout<<"Must be greater than 0 and less than or equal to 2.0\n";
}
cout<<"\nType the drag pulling people along. 0.5 recommended. Then press 'return'.\n\n";
cin>>Drag;
cout<<"\nType 1 if you wants tracks in the .dxf file. Otherwise any other number. Then press 'return'.\n\n";
cin>>tracks;*/
Theta=120.0*PI/180.0;fineness=0.5;Drag=1.0;tracks=0;

LastMaudZone[0]=int(24*fineness);if(LastMaudZone[0]>MaxLastMaudZone){cout<<"Too many zones\n";return 0;}
LastMaudZone[1]=int(18*fineness);if(LastMaudZone[1]>MaxLastMaudZone){cout<<"Too many zones\n";return 0;}

MeanDiameter=2.5/fineness;
deltat=1.0;
Walking=1.0/fineness;

ThetaOver2=Theta/2.0;

if(Theta>PI/2.0)
{
SwitchDirection=-1;
EWalking=Walking*sin(ThetaOver2);NWalking=-Walking*cos(ThetaOver2);
}
else
{
SwitchDirection=1;
EWalking=Walking*cos(ThetaOver2);NWalking=Walking*sin(ThetaOver2);
}

NewPersonsPerCycle=int(5.0*fineness);
restitution=0.7;

ShuteAngle=20.0*PI/180.0;
CriticalStartAngle=atan((1.0*halfH)/(1.0*halfW));

OverlapFactor=0.5;

Overlap[0]=(2.0*halfW)*OverlapFactor/(1.0*LastMaudZone[0]+1.0);
Overlap[1]=(2.0*halfH)*OverlapFactor/(1.0*LastMaudZone[1]+1.0);

for(MaudZone[0]=0;MaudZone[0]<=LastMaudZone[0];MaudZone[0]++)
{
Lower[0][MaudZone[0]]=-halfW+(2.0*halfW)*(1.0*MaudZone[0]+0.0)/(1.0*LastMaudZone[0]+1.0);
Upper[0][MaudZone[0]]=-halfW+(2.0*halfW)*(1.0*MaudZone[0]+1.0)/(1.0*LastMaudZone[0]+1.0);
}

for(MaudZone[1]=0;MaudZone[1]<=LastMaudZone[1];MaudZone[1]++)
{
Lower[1][MaudZone[1]]=-halfH+(2.0*halfH)*(1.0*MaudZone[1]+0.0)/(1.0*LastMaudZone[1]+1.0);
Upper[1][MaudZone[1]]=-halfH+(2.0*halfH)*(1.0*MaudZone[1]+1.0)/(1.0*LastMaudZone[1]+1.0);
}
for(Person=0;Person<=LastPerson;Person++)ColourOrHidden[Person]=2;
stop=0;
LastPersonUsed=-1;
Maud<<"0\nSECTION\n2\nENTITIES\n";
MaudGraphics(1.0,1.0,1.0);
return 0;
}
    
static void Draw(void)
{
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
if(stop!=1)
{
for(NewPerson=0;NewPerson<=NewPersonsPerCycle;NewPerson++)
{
FoundOne=0;
for(Person=0;Person<=LastPerson;Person++)
{
if(LastPersonUsed<Person)LastPersonUsed=Person;
if(ColourOrHidden[Person]==2)
{
r[Person]=SizeDist();
m[Person]=PI*r[Person]*r[Person];
if(MaudRandom()>0.5)
ColourOrHidden[Person]=0;
else
ColourOrHidden[Person]=1;

if(ColourOrHidden[Person]==0)
{
StartAngle=ThetaOver2+ShuteAngle*(1.0*MaudRandom()-0.5);
if(Theta>PI/2.0)
{
StartAngle+=PI/2.0;
Vely[0][Person]=+EWalking*SwitchDirection;
Vely[1][Person]=-NWalking;
}
else
{
Vely[0][Person]=+EWalking*SwitchDirection;
Vely[1][Person]=+NWalking;
}
StartingCoords();
}
else
{
StartAngle=-ThetaOver2+ShuteAngle*(1.0*MaudRandom()-0.5);
if(Theta>PI/2.0)
{
StartAngle+=PI/2.0;
Vely[0][Person]=+EWalking;
Vely[1][Person]=-NWalking;
}
else
{
Vely[0][Person]=+EWalking;
Vely[1][Person]=-NWalking;
}
StartingCoords();
}
FoundOne=1;
}
if(FoundOne==1)break;
}
}

for(Person=0;Person<=LastPersonUsed;Person++)
{
if(ColourOrHidden[Person]!=2)
{
if(ColourOrHidden[Person]==0)
{
Vely[0][Person]+=(+EWalking*SwitchDirection-Vely[0][Person])*Drag*deltat;
Vely[1][Person]+=(+NWalking*SwitchDirection-Vely[1][Person])*Drag*deltat;
}
else
{
Vely[0][Person]+=(+EWalking-Vely[0][Person])*Drag*deltat;
Vely[1][Person]+=(-NWalking-Vely[1][Person])*Drag*deltat;
}
if(Coord[0][Person]>+(halfW-r[Person])
 ||Coord[0][Person]<-(halfW-r[Person])
 ||Coord[1][Person]>+(halfH-r[Person])
 ||Coord[1][Person]<-(halfH-r[Person]))ColourOrHidden[Person]=2;
}
}

for(MaudZone[0]=0;MaudZone[0]<=LastMaudZone[0];MaudZone[0]++)
{
for(MaudZone[1]=0;MaudZone[1]<=LastMaudZone[1];MaudZone[1]++)
{
LastInMaudZone=-1;
for(Person=0;Person<=LastPersonUsed;Person++)
{
if(ColourOrHidden[Person]!=2)
{
if(Lower[0][MaudZone[0]]-Overlap[0]<Coord[0][Person]&&Coord[0][Person]<Upper[0][MaudZone[0]]+Overlap[0]
 &&Lower[1][MaudZone[1]]-Overlap[1]<Coord[1][Person]&&Coord[1][Person]<Upper[1][MaudZone[1]]+Overlap[1])
{
if(LastInMaudZone<MaxLastInMaudZone)
{
LastInMaudZone++;
InMaudZone[LastInMaudZone]=Person;
}
}
}
}

for(PersonInMaudZone=0;PersonInMaudZone<=LastInMaudZone-1;PersonInMaudZone++)
{
Person=InMaudZone[PersonInMaudZone];
if(ColourOrHidden[Person]!=2)
{
for(otherPersonInMaudZone=PersonInMaudZone+1;otherPersonInMaudZone<=LastInMaudZone;otherPersonInMaudZone++)
{
otherPerson=InMaudZone[otherPersonInMaudZone];
if(ColourOrHidden[otherPerson]!=2)Collision();
}
}
}
}
}

for(Person=0;Person<=LastPersonUsed;Person++)
{
if(ColourOrHidden[Person]!=2)
{
for(i=0;i<=1;i++)
{
PreviousCoord[i]=Coord[i][Person];
Coord[i][Person]+=Vely[i][Person]*deltat;
}
if(tracks==1)
{
Maud<<"0\nLINE\n8\nTrack "<<ColourOrHidden[Person]<<"\n";
Maud<<"10\n"<<PreviousCoord[0]<<"\n";
Maud<<"20\n"<<PreviousCoord[1]<<"\n";
Maud<<"11\n"<<Coord[0][Person]<<"\n";
Maud<<"21\n"<<Coord[1][Person]<<"\n";
Maud<<"62\n"<<ColourOrHidden[Person]<<"\n";
}
}
}
}

for(Person=0;Person<=LastPersonUsed;Person++)
{
if(ColourOrHidden[Person]!=2)
{
if(ColourOrHidden[Person]==0)glColor3f(0.25,0.0,1.0);
else glColor3f(0.3,0.5,0.0);
    glPushMatrix();
	glTranslatef(Coord[0][Person],Coord[1][Person],0.0);
	glutSolidSphere(r[Person],8,8);
	glPopMatrix();
}
}

glColor3f(0.0,0.0,0.0);

glRasterPos2f(-halfW+20.0,-halfH+20.0);
MyText="% used ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%i",int((100.0*LastPersonUsed)/LastPerson));
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);

glutSwapBuffers();
}

void Collision(void)
{
float MinSeparation,a,c,rootThing,tc;

SeparationSq=0.0;
for(i=0;i<=1;i++)
{
deltaCoord[i]=Coord[i][Person]-Coord[i][otherPerson];
deltaVely[i]=Vely[i][Person]-Vely[i][otherPerson];
SeparationSq+=deltaCoord[i]*deltaCoord[i];
}

MinSeparation=r[Person]+r[otherPerson];

a=deltaVely[0]*deltaVely[0]+deltaVely[1]*deltaVely[1];
bOver2=deltaCoord[0]*deltaVely[0]+deltaCoord[1]*deltaVely[1];
c=SeparationSq-MinSeparation*MinSeparation;
MassSum=m[Person]+m[otherPerson];
rootThing=bOver2*bOver2-a*c;

if(rootThing>=0.0)
{
tc=(-bOver2-sqrt(rootThing))/a;
if(0.0<=tc&&tc<=deltat&&c>0.0)
BounceOff();
}

if(c<0.0)
{
Separation=sqrt(SeparationSq);
for(i=0;i<=1;i++)
{
Coord[i][Person]     +=deltaCoord[i]*(MinSeparation/Separation-1.0)*m[otherPerson]/MassSum;
Coord[i][otherPerson]-=deltaCoord[i]*(MinSeparation/Separation-1.0)*m[Person]/MassSum;
BounceOff();
}
}
}

void BounceOff(void)
{
float I;

I=(1.0+restitution)*m[Person]*m[otherPerson]*bOver2/(MassSum*SeparationSq);
if(I<0.0)
{
for(i=0;i<=1;i++)
{
Vely[i][Person]-=I*deltaCoord[i]/m[Person];
Vely[i][otherPerson]+=I*deltaCoord[i]/m[otherPerson];
}
}
}

double MaudRandom(void)
{
return (1.0*rand())/(1.0*RAND_MAX);
}

double SizeDist(void)
{
double Thingy;
for(;;)
{
for(;;)
{
Thingy=1.0-MaudRandom();
if(Thingy>0.0)break;
}
Thingy=-MeanDiameter*MeanDiameter*log(Thingy)/PI;
if(Thingy>=0.0)break;
}
return 0.25*MeanDiameter+0.75*sqrt(Thingy);
}

void StartingCoords(void)
{
if(fabs(StartAngle)<CriticalStartAngle)
{
Coord[0][Person]=-(halfW-r[Person]);
Coord[1][Person]=Coord[0][Person]*tan(StartAngle);
}
else
{
if(fabs(StartAngle)<PI-CriticalStartAngle)
{
if(StartAngle>0.0)Coord[1][Person]=-(halfH-r[Person]);else Coord[1][Person]=(halfH-r[Person]);
Coord[0][Person]=Coord[1][Person]/tan(StartAngle);
}
else
{
Coord[0][Person]=(halfW-r[Person]);
Coord[1][Person]=Coord[0][Person]*tan(StartAngle);
}
}
}

void WriteDXFFile(void)
{
cout<<"Writing dxf file\n";

for(Person=0;Person<=LastPersonUsed;Person++)
{
if(ColourOrHidden[Person]!=2)
{
Maud<<"0\nCIRCLE\n8\nPerson "<<ColourOrHidden[Person]<<"\n";
Maud<<"10\n"<<Coord[0][Person]<<"\n";
Maud<<"20\n"<<Coord[1][Person]<<"\n";
Maud<<"40\n"<<r[Person]<<"\n";
Maud<<"62\n"<<ColourOrHidden[Person]<<"\n";
}
}

Maud<<"0\nENDSEC\n0\nEOF\n";
Maud.close();
cout<<"Finished\n";
}


void MaudGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen)
{
    glutInitWindowSize(2*halfW, 2*halfH);
	if(Macintosh==1)    glutInit(&argc, argv);               // Remove this line for PC or Sun
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("People flow");
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

	glViewport(0,0,2*halfW,2*halfH);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

    glOrtho(-halfW,halfW,-halfW,halfW,-halfW,halfW);
    glClear(GL_COLOR_BUFFER_BIT);
	
    glutKeyboardFunc(Key);
    glutIdleFunc(Draw);
    glutDisplayFunc(Draw);
	
	glEnable(GL_POINT_SMOOTH);// glDisable(GL_POINT_SMOOTH) is the reverse. Enable means antialiasing 
	
	glScalef(1.0,(1.0*halfW)/(1.0*halfH),1.0);
	
	glutMainLoop();	
}

static void Key(unsigned char key, int x, int y)
{
switch (key)
{
case 's':if(stop==0)stop=1;else stop=0;break;
case 'q':gluDeleteQuadric(quadObj);WriteDXFFile(); exit(0);
}
}
