#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
using namespace std;
#include <GLUT/glut.h>  // For Macintosh
						//#include <GL/glut.h>    // For PC or Sun
#define   Macintosh 1  //1 for Macintosh any other number for PC or Sun
#define   LastNode       500
#define MaxLastZone0      15
#define MaxLastZone1       5
#define MaxLastInZone    500
void Calculation(void);
void SortIntoZones(void);
void Contribution(void);
void Invert(void);
void ChrisGraphics(float BackGroundRed,float BackGroundBlue,float BackGroundGreen);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
int argc;
char **argv;
GLUquadricObj *quadObj;
int   LastInZone[MaxLastZone0+1][MaxLastZone1+1],
InZone[MaxLastZone0+1][MaxLastZone1+1][MaxLastInZone+1],
HalfDimension[2],Zone[2],LastZone[2],WhichZone[2],
Node,OtherNode,i,j,NodeInZone,OtherNodeInZone,ShowNodes,ShowLines;
float x[2][LastNode+1],Force[2][LastNode+1],Velocity[2][LastNode+1],
Stiffness[2][2][LastNode+1],Flex[2][2],deltax[2],HalfPicture[2],ZoneSize[2],
PI,a,aSq,ToverL,dTbydL,LSq,tempForce,tempStiffness,slow,CarryOver,Bounce,L,MyExp,MinZoneSizeSq;
int main(void){
	HalfDimension[0]=630;
	HalfDimension[1]=370;
	for(i=0;i<=1;i++)HalfPicture[i]=1.0*HalfDimension[i]-25.0;
	cout<<"This program runs for ever until the window is closed or 'q' is pressed on the keyboard.\n";
	cout<<"Typing 'n' makes the nodes dissapear or reappear.\n";
	cout<<"Typing 'l' makes the nodes dissapear or reappear.\n";
	PI=4.0*atan(1.0);srand(time(NULL));
	slow=0.001;CarryOver=0.95;Bounce=0.5;
	a=100.0;aSq=a*a;
	ShowNodes=1;ShowLines=1;
	for(i=0;i<=1;i++){
		LastZone[i]=int(2.0*HalfDimension[i]/(3.0*a));
		if(LastZone[i]>MaxLastZone0){cout<<LastZone[i]<<" is too many zones in the "<<i<<" direction\n";return 0;}
		ZoneSize[i]=2.0*HalfDimension[i]/(1.0*LastZone[i]);
		cout<<"ZoneSize["<<i<<"]/a ="<<ZoneSize[i]/a<<"\n";
	}
	cout<<"Number of zones is "<<LastZone[0]+1<<" x "<<LastZone[1]+1<<"\n";
	if(ZoneSize[0]<ZoneSize[1])MinZoneSizeSq=ZoneSize[0]*ZoneSize[0];
	else MinZoneSizeSq=ZoneSize[1]*ZoneSize[1];
	cout<<"Value of exponential at edge of zone = "<<exp(-MinZoneSizeSq/aSq)<<"\n";
	for(Node=0;Node<=LastNode;Node++){
		for(i=0;i<=1;i++)x[i][Node]=0.75*HalfPicture[i]*((2.0*rand())/(1.0*RAND_MAX)-1.0);
	}
	ChrisGraphics(1.0,1.0,1.0);
	return 0;
}
static void Draw(void){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	Calculation();
	glutSwapBuffers();
}
void Calculation(void){
	for(Node=0;Node<=LastNode;Node++){
		for(i=0;i<=1;i++){
			Force[i][Node]=0.0;
			for(j=0;j<=1;j++)Stiffness[i][j][Node]=0.0;
		}}
	SortIntoZones();
	for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++){
		for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++){
			for(NodeInZone=0;NodeInZone<=LastInZone[Zone[0]][Zone[1]];NodeInZone++){
				Node=InZone[Zone[0]][Zone[1]][NodeInZone];
				for(WhichZone[0]=Zone[0]-1;WhichZone[0]<=Zone[0]+1;WhichZone[0]++){
					if(0<=WhichZone[0]&&WhichZone[0]<=LastZone[0]){
						for(WhichZone[1]=Zone[1]-1;WhichZone[1]<=Zone[1]+1;WhichZone[1]++){
							if(0<=WhichZone[1]&&WhichZone[1]<=LastZone[1]){
								for(OtherNodeInZone=0;OtherNodeInZone<=LastInZone[WhichZone[0]][WhichZone[1]];OtherNodeInZone++){
									OtherNode=InZone[WhichZone[0]][WhichZone[1]][OtherNodeInZone];
									if(OtherNode<Node){
										LSq=0.0;
										for(i=0;i<=1;i++){
											deltax[i]=x[i][OtherNode]-x[i][Node];
											LSq+=deltax[i]*deltax[i];
										}
										if(LSq/aSq>1.0e-24&&LSq<MinZoneSizeSq){
											Contribution();
											if(ShowLines==1){
												glColor4f(0.0,0.0,1.0-MyExp,MyExp);
												glBegin(GL_LINES);
												glVertex2f(x[0][Node],x[1][Node]);
												glVertex2f(x[0][OtherNode],x[1][OtherNode]);
												glEnd();
											}}}}}}}}}}}
	if(ShowNodes==1){glColor4f(1.0,0.0,0.0,0.5);glBegin(GL_POINTS);}
	for(Node=0;Node<=LastNode;Node++){
		if(ShowNodes==1)glVertex2f(x[0][Node],x[1][Node]);
		Invert();
		for(i=0;i<=1;i++){
			Velocity[i][Node]=CarryOver*Velocity[i][Node];
			for(j=0;j<=1;j++)Velocity[i][Node]+=slow*Flex[i][j]*Force[j][Node];
			x[i][Node]+=Velocity[i][Node];
		}
		for(i=0;i<=1;i++){
			if(x[i][Node]>HalfPicture[i]){
				x[i][Node]=HalfPicture[i];
				if(Velocity[i][Node]>0.0)Velocity[i][Node]=-Bounce*Velocity[i][Node];
			}
			if(x[i][Node]<-HalfPicture[i]){
				x[i][Node]=-HalfPicture[i];
				if(Velocity[i][Node]<0.0)Velocity[i][Node]=-Bounce*Velocity[i][Node];
			}}}
	if(ShowNodes==1)glEnd();
}
void ChrisGraphics(float BackGroundRed,float BackGroundGreen,float BackGroundBlue){
	glutInitWindowSize(2*HalfDimension[0], 2*HalfDimension[1]);
	if(Macintosh==1)    glutInit(&argc, argv);               // Remove this line for PC or Sun
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Maud");
	glClearColor(BackGroundRed,BackGroundGreen,BackGroundBlue,0);
	glViewport(0, 0, 2*HalfDimension[0], 2*HalfDimension[1]);
	gluOrtho2D(-HalfDimension[0],HalfDimension[0],-HalfDimension[1],HalfDimension[1]);
	glClear(GL_COLOR_BUFFER_BIT);
	glLineWidth(0.5);
	glPointSize(5.0);
	glutKeyboardFunc(Key);
	glutIdleFunc(Draw);
	glutDisplayFunc(Draw);
	glEnable(GL_POINT_SMOOTH);
	glEnable (GL_BLEND);glBlendFunc (GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glutMainLoop();
}
static void Key(unsigned char key, int x, int y){
	switch (key){
		case 'n':if(ShowNodes==0)ShowNodes=1;else ShowNodes=0;break;
		case 'l':if(ShowLines==0)ShowLines=1;else ShowLines=0;break;
		case 'q':gluDeleteQuadric(quadObj);exit(0);
	}}
void SortIntoZones(void){
	for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++){
		for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++)
			LastInZone[Zone[0]][Zone[1]]=-1;
	}
	for(Node=0;Node<=LastNode;Node++){
		for(i=0;i<=1;i++){
			Zone[i]=int((x[i][Node]+HalfDimension[i])/ZoneSize[i]);
			if(Zone[i]<0)Zone[i]=0;
			if(Zone[i]>LastZone[i])Zone[i]=LastZone[i];
		}
		LastInZone[Zone[0]][Zone[1]]++;
		if(LastInZone[Zone[0]][Zone[1]]<=MaxLastInZone)InZone[Zone[0]][Zone[1]][LastInZone[Zone[0]][Zone[1]]]=Node;
	}
	for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++){
		for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++){
			if(LastInZone[Zone[0]][Zone[1]]>MaxLastInZone){
				cout<<"Number in Zone "<<Zone[0]<<" "<<Zone[1]<<" is "<<LastInZone[Zone[0]][Zone[1]]<<"\n";
				LastInZone[Zone[0]][Zone[1]]=MaxLastInZone;
			}}}}
void Contribution(void){
	L=sqrt(LSq);
	MyExp=exp(-LSq/aSq);
	ToverL=-(aSq/LSq)*MyExp;
	dTbydL=(aSq/LSq)*(1.0+2.0*LSq/aSq)*MyExp;
	for(i=0;i<=1;i++){
		tempForce=ToverL*deltax[i];
		Force[i][Node]+=tempForce;
		Force[i][OtherNode]-=tempForce;
		Stiffness[i][i][Node]+=ToverL;
		Stiffness[i][i][OtherNode]+=ToverL;
		for(j=0;j<=1;j++){
			tempStiffness=(dTbydL-ToverL)*deltax[i]*deltax[j]/LSq;
			Stiffness[i][j][Node]+=tempStiffness;
			Stiffness[i][j][OtherNode]+=tempStiffness;
		}}}
void Invert(void){
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
	else{
		Flex[0][0]=0.0;
		Flex[1][1]=0.0;
		Flex[0][1]=0.0;
		Flex[1][0]=0.0;
	}}
