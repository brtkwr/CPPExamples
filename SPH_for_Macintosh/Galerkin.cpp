#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#include <stdio.h>
#include <string.h>

#include "Graphics.h"

#define AbsoluteLastNode 1785000

#define MaxLastZone0 299
#define MaxLastZone1 169

#define MaxLastInZone 200

#define LastExp 10000

#define LastCylinderNode 100
#define LastTimeBeforeNow 10
#define MaxLastFabricNode 201

void DoThePicture(void);
void SortIntoZones(void);
void Boundaries(void);
void OtherNodePart(void);
void Interaction(void);
void FabricForcesAndMotion(void);
void Motion(void);
void StartNodeOnLeft(void);
void SplitNode(void);
int CreateData(void);
float My_exp(float MyArg);
void SetUpExp(void);

float  Coord[LastTimeBeforeNow+1][2][AbsoluteLastNode+1],v[2][AbsoluteLastNode+1],Force[2][AbsoluteLastNode+1],
	   Mass[AbsoluteLastNode+1],Omega[AbsoluteLastNode+1],
	   Stored_exp[LastExp+1],
	    FabricElementLength[MaxLastFabricNode],//Last fabric element = LastFabricNode - 1
	   RateOfLengthIncrease[MaxLastFabricNode],//Last fabric element = LastFabricNode - 1
	         FabricForce[2][MaxLastFabricNode+1],
	   
	   CylinderCoord[2][LastCylinderNode+1],

	   Boundary[2],delta_r[2],delta_v[2],NetForce[2],PreviousNetForce[2],
       
	   PI,ZoneSize,Deltat,FabricDeltat,ReynoldsNumber,mu_over_ro,
	   c,cSq,WindSpeed,
	   aiSqElastic,ajSqElastic,BetaSqElastic,
	   aiSqViscous,ajSqViscous,BetaSqViscous,
	   ZetaElastic,ZetaViscous,ThisForce,
	   ArgMultiplier,BoundaryZone,
	   roAverage,BasicMass,InitialTotalMass,TotalMass,
	   rij_Sq,MomentOfInertiaRatio,
	   DistanceFromCentreSq,DistanceFromInkCentreSq,MinDistanceSq,
	   OmegaArm,Moment_i,Moment_j,
	   DrawRadius,DistSqTimesBasicMassOverMassElastic,DistSqTimesBasicMassOverMassViscous,
	   SeparationAngle,SeparationDistance,
	   Scale,CylinderRadius,CylinderRadiusSq,ControlCylinderRadiusSq,CylinderPosition,
	   InkHalfDimension,InkHalfDimensionSq,InkPosition,TimeFactor,DragLiftFactor,
	   SecondMomentOfAreaBit,AngularBit,BaseOmega,tempforce,
	   ScalarProductOver_rij,
	   ForceSmooth,FabricPreStress,AverageFabricStrain,
	   FabricSlackLength,FabricMassPerNode,MassOfThisNode,
	   FabricEA,FabricPreStressFactor,FabricDamping,
	   FabricElementSlackLength;

int    NodeType[AbsoluteLastNode+1],
       InTank[AbsoluteLastNode+1],StreakOrNot[AbsoluteLastNode+1],
       LastInZone[MaxLastZone0+1][MaxLastZone1+1],
       InZone[MaxLastZone0+1][MaxLastZone1+1][MaxLastInZone+1],
	   Zone[2],LastZone[2],WhichZone[2],ThisWhichZone1,
	   InitialLastNodeOutOfTank,InitialNumberInTank,NumberInTank,
	   NewNode,FirstSplitNode,SecondSplitNode,Found,
	   Node,LastNode,LastNodeOutOfTank,otherNode,xyz,otherxyz,
	   NodeInZone,otherNodeInZone,
	   Cycle,PictureSaveInterval,SavePictureOrNot,
	   LengthOfString,MyCharacter,angle,
	   NumberToBeInjected,TimeBeforeNow,
	   Fabric,Sail,LastFabricNode,FabricCycle,LastFabricCycle;

char NumericalValue[101];
char *MyText;

ofstream Constance("DragLift.dxf");

int main(void)
{
PI=4.0*atan(1.0);
SetUpExp();

if(CreateData())return 0;

ChrisGraphics(1.0,1.0,1.0);

return 0;
}
    
static void Draw(void)
{
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
Cycle++;

NumberInTank=0;TotalMass=0.0;
for(Node=0;Node<=LastNode;Node++)
{
NumberInTank+=InTank[Node];
if(InTank[Node]!=1)
{
TotalMass+=Mass[Node];
LastNodeOutOfTank=Node;
}
}

NumberToBeInjected=int((InitialTotalMass-TotalMass)/BasicMass);

if(NumberToBeInjected>0)
{
for(NewNode=1;NewNode<=NumberToBeInjected;NewNode++)StartNodeOnLeft();
}

Boundaries();
SortIntoZones();
Interaction();
if(Fabric==1)FabricForcesAndMotion();
Motion();
DoThePicture();
if(Cycle==PictureSaveInterval*int(Cycle/PictureSaveInterval)&&SavePictureOrNot==1)SavePicture();
}

void DoThePicture(void)
{
float rv[2],Proportion,MyRed,MyGreen,MyBlue,ColourFactor;

glLineWidth(0.2);

for(Node=0;Node<=LastNodeOutOfTank;Node++)
{
if(Node>LastFabricNode)
{
if(InTank[Node]==0)
{
if(StreakOrNot[Node]==1||Streak==0)
{
if(VorticityColour==1)
{
if(NodeType[Node]==2)
{
  MyRed=0.5;
MyGreen=0.5;
 MyBlue=0.5;
}
else
{
ColourFactor=Omega[Node]/BaseOmega;
if(ColourFactor>+1.0)ColourFactor=+1.0;
if(ColourFactor<-1.0)ColourFactor=-1.0;
  MyRed=0.0;
MyGreen=0.0;
 MyBlue=0.0;
if(ColourFactor>0.0)MyRed=ColourFactor;else MyBlue=-ColourFactor;
}
}
else{MyRed=0.0;MyGreen=0.0;MyBlue=0.0;}

glBegin(GL_LINE_STRIP);                                  
for(TimeBeforeNow=0;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)
{
Proportion=(1.0*TimeBeforeNow)/(1.0*LastTimeBeforeNow);
Proportion=Proportion*Proportion;
//glColor4f(Proportion+MyRed*(1.0-Proportion),Proportion+MyGreen*(1.0-Proportion),Proportion+MyBlue*(1.0-Proportion),1.0-1.0*Proportion);                                                       
glColor4f(MyRed,MyGreen,MyBlue,1.0-1.0*Proportion);                                                       
glBegin(GL_LINE_STRIP);
rv[0]=Coord[TimeBeforeNow][0][Node];
rv[1]=Coord[TimeBeforeNow][1][Node];
glVertex2fv(rv);
}
glEnd();
}
}
}
}

glLineWidth(1.0);
glColor3f(0.0,0.0,0.0);

glBegin(GL_LINE_STRIP);
for(Node=0;Node<=LastCylinderNode;Node++)
{
rv[0]=CylinderCoord[0][Node];
rv[1]=CylinderCoord[1][Node];
glVertex2fv(rv);
}

rv[0]=CylinderCoord[0][0];
rv[1]=CylinderCoord[1][0];
glVertex2fv(rv);

glEnd();

if(Fabric==1)
{
glLineWidth(2.0);
glBegin(GL_LINE_STRIP);
for(Node=0;Node<=LastFabricNode;Node++)
{
rv[0]=Coord[0][0][Node];
rv[1]=Coord[0][1][Node];
glVertex2fv(rv);
}
glEnd();
}

glLineWidth(1.0);
glBegin(GL_LINE_STRIP);
rv[0]=-Boundary[0];rv[1]=-Boundary[1];glVertex2fv(rv);
rv[0]=+Boundary[0];rv[1]=-Boundary[1];glVertex2fv(rv);
rv[0]=+Boundary[0];rv[1]=+Boundary[1];glVertex2fv(rv);
rv[0]=-Boundary[0];rv[1]=+Boundary[1];glVertex2fv(rv);
rv[0]=-Boundary[0];rv[1]=-Boundary[1];glVertex2fv(rv);
glEnd();

if(Information==1)
{
glColor3f(0.0,0.0,0.0);

glRasterPos2f(-550.0*Scale,-450.0*Scale);
MyText="Reynolds number ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%i",int(ReynoldsNumber));
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);

MyText="   Cycle ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%i",Cycle);
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);

MyText="   In tank";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%6.2f %%",(100.0*NumberInTank)/(1.0*LastNode));
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);

if(Fabric==1)
{
MyText="   Average fabric strain ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%6.2f %%",100.0*AverageFabricStrain);
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);

MyText="   Fabric pretension ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%8.4f",
FabricPreStress/((1.0/2.0)*roAverage*WindSpeed*WindSpeed*FabricSlackLength));
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
}

/*MyText="   Total mass  ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%f",TotalMass);
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);

MyText="   Number of nodes ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%i",LastNode-NumberInTank);
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);

MyText="   Time ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%f",Cycle*TimeFactor);
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);

MyText="   Drag coefficient ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%f",NetForce[0]);
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);

MyText="   Lift coefficient ";
LengthOfString=strlen(MyText);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
sprintf(NumericalValue,"%f",NetForce[1]);
LengthOfString=strlen(NumericalValue);
for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);*/
}

glutSwapBuffers();
}

void Boundaries(void)
{
for(xyz=0;xyz<=1;xyz++)
{
PreviousNetForce[xyz]=NetForce[xyz];
NetForce[xyz]=0.0;
}

for(Node=0;Node<=LastNodeOutOfTank;Node++)
{
if(InTank[Node]==0)
{
DistanceFromInkCentreSq
=(Coord[0][0][Node]+InkPosition)*(Coord[0][0][Node]+InkPosition)+Coord[0][1][Node]*Coord[0][1][Node];

DistanceFromCentreSq
=(Coord[0][0][Node]+CylinderPosition)*(Coord[0][0][Node]+CylinderPosition)+Coord[0][1][Node]*Coord[0][1][Node];
if(DistanceFromInkCentreSq>DistanceFromCentreSq)MinDistanceSq=DistanceFromCentreSq;
										   else MinDistanceSq=DistanceFromInkCentreSq;
if(MinDistanceSq<9.0*ControlCylinderRadiusSq&&Mass[Node]>=0.9*BasicMass&&Node>LastFabricNode)SplitNode();
if(MinDistanceSq<6.0*ControlCylinderRadiusSq&&Mass[Node]>=0.3*BasicMass&&Node>LastFabricNode)SplitNode();
}
}

for(Node=0;Node<=LastNodeOutOfTank;Node++)
{
if(InTank[Node]==0)
{
for(xyz=0;xyz<=1;xyz++)Force[xyz][Node]=0.0;

if(Node>LastFabricNode)
{
DistanceFromInkCentreSq
=(Coord[0][0][Node]+InkPosition)*(Coord[0][0][Node]+InkPosition)+Coord[0][1][Node]*Coord[0][1][Node];

DistanceFromCentreSq
=(Coord[0][0][Node]+CylinderPosition)*(Coord[0][0][Node]+CylinderPosition)+Coord[0][1][Node]*Coord[0][1][Node];

if(DistanceFromInkCentreSq<InkHalfDimensionSq)StreakOrNot[Node]=1;

if(DistanceFromCentreSq<CylinderRadiusSq)
{
NodeType[Node]=1;
StreakOrNot[Node]=0;
for(xyz=0;xyz<=1;xyz++)
{
tempforce=0.9*v[xyz][Node]*Mass[Node]/Deltat;
Force[xyz][Node]-=tempforce;
NetForce[xyz]+=tempforce;
}
}
}
else NodeType[Node]=0;

if(Coord[0][0][Node]>+Boundary[0])InTank[Node]=1;
if(Coord[0][0][Node]<-Boundary[0])InTank[Node]=1;
if(Coord[0][1][Node]>+Boundary[1])
{
Coord[0][1][Node]-=2.0*Boundary[1];
for(TimeBeforeNow=1;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)  
{
for(xyz=0;xyz<=1;xyz++)Coord[TimeBeforeNow][xyz][Node]=Coord[0][xyz][Node];
}
}
if(Coord[0][1][Node]<-Boundary[1])
{
Coord[0][1][Node]+=2.0*Boundary[1];
for(TimeBeforeNow=1;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)  
{
for(xyz=0;xyz<=1;xyz++)Coord[TimeBeforeNow][xyz][Node]=Coord[0][xyz][Node];
}
}
if(NodeType[Node]!=1)
{
if((Boundary[0]+Coord[0][0][Node]<BoundaryZone)
 ||(Boundary[0]-Coord[0][0][Node]<BoundaryZone))
{
NodeType[Node]=2;
v[0][Node]=WindSpeed;
v[1][Node]=0.0;
}
else NodeType[Node]=0;
}
}
}

for(xyz=0;xyz<=1;xyz++)NetForce[xyz]=ForceSmooth*DragLiftFactor*NetForce[xyz]+(1.0-ForceSmooth)*PreviousNetForce[xyz];

Constance<<"0\nLINE\n8\nDrag Coefficient\n";
Constance<<"10\n"<<TimeFactor*(Cycle-1)<<"\n";
Constance<<"20\n"<<PreviousNetForce[0]<<"\n";
Constance<<"11\n"<<TimeFactor*Cycle<<"\n";
Constance<<"21\n"<<NetForce[0]<<"\n";

Constance<<"0\nLINE\n8\nLift Coefficient\n";
Constance<<"10\n"<<TimeFactor*(Cycle-1)<<"\n";
Constance<<"20\n"<<PreviousNetForce[1]<<"\n";
Constance<<"11\n"<<TimeFactor*Cycle<<"\n";
Constance<<"21\n"<<NetForce[1]<<"\n";
}

void SplitNode(void)
{
Found=0;
for(FirstSplitNode=0;FirstSplitNode<=LastNode;FirstSplitNode++)
{
if(InTank[FirstSplitNode]==1)Found=1;
if(Found==1)break;
}
if(Found==1)InTank[FirstSplitNode]=0;

Found=0;
for(SecondSplitNode=0;SecondSplitNode<=LastNode;SecondSplitNode++)
{
if(InTank[SecondSplitNode]==1)Found=1;
if(Found==1)break;
}
if(Found!=1)
InTank[FirstSplitNode]=1;
else
{
if(LastNodeOutOfTank<FirstSplitNode)LastNodeOutOfTank=FirstSplitNode;
if(LastNodeOutOfTank<SecondSplitNode)LastNodeOutOfTank=SecondSplitNode;

InTank[SecondSplitNode]=0;

SeparationAngle=(2.0*PI*rand())/(1.0*RAND_MAX);
SeparationDistance=2.0*sqrt(DistSqTimesBasicMassOverMassViscous*Mass[Node]/BasicMass);

MomentOfInertiaRatio=3.0/27.0+
SeparationDistance*SeparationDistance/(Mass[Node]*DistSqTimesBasicMassOverMassViscous);
//cout<<MomentOfInertiaRatio<<"\n";
Omega[Node]=Omega[Node]/MomentOfInertiaRatio;

Coord[0][0][FirstSplitNode]=Coord[0][0][Node]+SeparationDistance*cos(SeparationAngle);
Coord[0][1][FirstSplitNode]=Coord[0][1][Node]+SeparationDistance*sin(SeparationAngle);
for(TimeBeforeNow=1;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)  
{
for(xyz=0;xyz<=1;xyz++)Coord[TimeBeforeNow][xyz][FirstSplitNode]=Coord[0][xyz][FirstSplitNode];
}
v[0][FirstSplitNode]=v[0][Node]-Omega[Node]*SeparationDistance*sin(SeparationAngle);
v[1][FirstSplitNode]=v[1][Node]+Omega[Node]*SeparationDistance*cos(SeparationAngle);

SeparationAngle+=2.0*PI/3.0;
Coord[0][0][SecondSplitNode]=Coord[0][0][Node]+SeparationDistance*cos(SeparationAngle);
Coord[0][1][SecondSplitNode]=Coord[0][1][Node]+SeparationDistance*sin(SeparationAngle);
for(TimeBeforeNow=1;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)  
{
for(xyz=0;xyz<=1;xyz++)Coord[TimeBeforeNow][xyz][SecondSplitNode]=Coord[0][xyz][SecondSplitNode];
}
v[0][SecondSplitNode]=v[0][Node]-Omega[Node]*SeparationDistance*sin(SeparationAngle);
v[1][SecondSplitNode]=v[1][Node]+Omega[Node]*SeparationDistance*cos(SeparationAngle);

SeparationAngle+=2.0*PI/3.0;
Coord[0][0][Node]=Coord[0][0][Node]+SeparationDistance*cos(SeparationAngle);
Coord[0][1][Node]=Coord[0][1][Node]+SeparationDistance*sin(SeparationAngle);
for(TimeBeforeNow=1;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)  
{
for(xyz=0;xyz<=1;xyz++)Coord[TimeBeforeNow][xyz][Node]=Coord[0][xyz][Node];
}
v[0][Node]=v[0][Node]-Omega[Node]*SeparationDistance*sin(SeparationAngle);
v[1][Node]=v[1][Node]+Omega[Node]*SeparationDistance*cos(SeparationAngle);

Mass[Node]=Mass[Node]/3.0;
Mass[FirstSplitNode]=Mass[Node];
Mass[SecondSplitNode]=Mass[Node];

StreakOrNot[Node]=0;
StreakOrNot[FirstSplitNode]=0;
StreakOrNot[SecondSplitNode]=0;
}
}

void StartNodeOnLeft(void)
{
Found=0;
for(Node=0;Node<=LastNode;Node++)
{
if(InTank[Node]==1)
{
if(LastNodeOutOfTank<Node)LastNodeOutOfTank=Node;

Coord[0][0][Node]=-Boundary[0]+WindSpeed*Deltat*(1.0*rand())/(1.0*RAND_MAX);
Coord[0][1][Node]=Boundary[1]*((2.0*rand())/(1.0*RAND_MAX)-1.0);
for(TimeBeforeNow=1;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)  
{
for(xyz=0;xyz<=1;xyz++)Coord[TimeBeforeNow][xyz][Node]=Coord[0][xyz][Node];
}
v[0][Node]=WindSpeed;v[1][Node]=0.0;
Mass[Node]=BasicMass;
Omega[Node]=0.0;
InTank[Node]=0;
StreakOrNot[Node]=0;
Found=1;
}
if(Found==1)break;
}
}

void SortIntoZones(void)
{
for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++)
{
for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++)
LastInZone[Zone[0]][Zone[1]]=-1;
}

for(Node=0;Node<=LastNodeOutOfTank;Node++)
{
if(InTank[Node]==0)
{
for(xyz=0;xyz<=1;xyz++)
Zone[xyz]=int((Coord[0][xyz][Node]+Boundary[xyz])/ZoneSize);
if(Zone[0]<0||Zone[0]>LastZone[0]||Zone[1]<0||Zone[1]>LastZone[1])InTank[Node]=1;
else
{
LastInZone[Zone[0]][Zone[1]]++;
if(LastInZone[Zone[0]][Zone[1]]<=MaxLastInZone)InZone[Zone[0]][Zone[1]][LastInZone[Zone[0]][Zone[1]]]=Node;
}
}
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

void Interaction(void)
{
for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++)
{
for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++)
{
for(NodeInZone=0;NodeInZone<=LastInZone[Zone[0]][Zone[1]];NodeInZone++)
{
Node=InZone[Zone[0]][Zone[1]][NodeInZone];
aiSqElastic=DistSqTimesBasicMassOverMassElastic*Mass[Node];
aiSqViscous=DistSqTimesBasicMassOverMassViscous*Mass[Node];
TotalMass+=Mass[Node];

OtherNodePart();
}
}
}
}

void OtherNodePart(void)
{
for(WhichZone[0]=Zone[0]-1;WhichZone[0]<=Zone[0]+1;WhichZone[0]++)
{
for(WhichZone[1]=Zone[1]-1;WhichZone[1]<=Zone[1]+1;WhichZone[1]++)
{
if(0<=WhichZone[0]&&WhichZone[0]<=LastZone[0])
{
ThisWhichZone1=WhichZone[1];
if(ThisWhichZone1>LastZone[1])ThisWhichZone1=0;
if(ThisWhichZone1<0)ThisWhichZone1=LastZone[1];
for(otherNodeInZone=0;otherNodeInZone<=LastInZone[WhichZone[0]][ThisWhichZone1];otherNodeInZone++)
{
otherNode=InZone[WhichZone[0]][ThisWhichZone1][otherNodeInZone];
if(otherNode<Node)
{
ajSqElastic=DistSqTimesBasicMassOverMassElastic*Mass[otherNode];
ajSqViscous=DistSqTimesBasicMassOverMassViscous*Mass[otherNode];

for(xyz=0;xyz<=1;xyz++)
{
delta_r[xyz]=Coord[0][xyz][otherNode]-Coord[0][xyz][Node];
delta_v[xyz]=v[xyz][otherNode]-v[xyz][Node];
}
if(delta_r[1]>+Boundary[1])delta_r[1]-=2.0*Boundary[1];
if(delta_r[1]<-Boundary[1])delta_r[1]+=2.0*Boundary[1];

rij_Sq=delta_r[0]*delta_r[0]+delta_r[1]*delta_r[1];
BetaSqElastic=rij_Sq/(aiSqElastic+ajSqElastic);
BetaSqViscous=rij_Sq/(aiSqViscous+ajSqViscous);
if(My_exp(-BetaSqElastic)!=0.0)
{
ZetaElastic=Mass[Node]*Mass[otherNode]/(PI*(aiSqElastic+ajSqElastic)*(aiSqElastic+ajSqElastic))*My_exp(-BetaSqElastic);
ZetaViscous=Mass[Node]*Mass[otherNode]/(PI*(aiSqViscous+ajSqViscous)*(aiSqViscous+ajSqViscous))*My_exp(-BetaSqViscous);

ScalarProductOver_rij=(delta_r[0]*delta_v[0]+delta_r[1]*delta_v[1])/sqrt(rij_Sq);

OmegaArm=(delta_r[0]*delta_v[1]-delta_r[1]*delta_v[0])/rij_Sq;

Moment_i=-8.0*(mu_over_ro/roAverage)*ZetaViscous*rij_Sq
*(Omega[Node     ]-OmegaArm)*aiSqViscous/(aiSqViscous+ajSqViscous);

Moment_j=-8.0*(mu_over_ro/roAverage)*ZetaViscous*rij_Sq
*(Omega[otherNode]-OmegaArm)*ajSqViscous/(aiSqViscous+ajSqViscous);

Omega[Node     ]+=Moment_i/(Mass[Node     ]*aiSqViscous);
Omega[otherNode]+=Moment_j/(Mass[otherNode]*ajSqViscous);

for(xyz=0;xyz<=1;xyz++)
{
ThisForce=-cSq*ZetaElastic*delta_r[xyz]/roAverage;

ThisForce+=(mu_over_ro/roAverage)*ZetaViscous*ScalarProductOver_rij*delta_r[xyz];

if(xyz==0)ThisForce-=delta_r[1]*(Moment_i+Moment_j)/rij_Sq;
	 else ThisForce+=delta_r[0]*(Moment_i+Moment_j)/rij_Sq;

Force[xyz][Node]     +=ThisForce;
Force[xyz][otherNode]-=ThisForce;
}
}
}
}
}
}
}
}

void FabricForcesAndMotion(void)
{
float FabricTensionCoeficient,TotalFabricLength;
for(FabricCycle=1;FabricCycle<=LastFabricCycle;FabricCycle++)
{
TotalFabricLength=0.0;
for(Node=0;Node<=LastFabricNode-1;Node++)
{
FabricElementLength[Node]=sqrt((Coord[0][0][Node+1]-Coord[0][0][Node])
                              *(Coord[0][0][Node+1]-Coord[0][0][Node])
				              +(Coord[0][1][Node+1]-Coord[0][1][Node])
				              *(Coord[0][1][Node+1]-Coord[0][1][Node]));
TotalFabricLength+=FabricElementLength[Node];
RateOfLengthIncrease[Node]
=((v[0][Node+1]-v[0][Node])*(Coord[0][0][Node+1]-Coord[0][0][Node])
 +(v[1][Node+1]-v[1][Node])*(Coord[0][1][Node+1]-Coord[0][1][Node]))/FabricElementLength[Node];
}
AverageFabricStrain=(TotalFabricLength-FabricSlackLength)/FabricSlackLength;
FabricPreStress+=FabricPreStressFactor*FabricEA*LastFabricNode*AverageFabricStrain;

for(Node=0;Node<=LastFabricNode;Node++)
{
for(xyz=0;xyz<=1;xyz++)FabricForce[xyz][Node]=Force[xyz][Node];
}
for(Node=0;Node<=LastFabricNode-1;Node++)
{
FabricTensionCoeficient=
(FabricPreStress
+FabricEA*(FabricElementLength[Node]-FabricElementSlackLength)/FabricElementSlackLength
+FabricDamping*RateOfLengthIncrease[Node]
)/FabricElementLength[Node];

for(xyz=0;xyz<=1;xyz++)
{
ThisForce=FabricTensionCoeficient*(Coord[0][xyz][Node+1]-Coord[0][xyz][Node]);
FabricForce[xyz][Node  ]+=ThisForce;
FabricForce[xyz][Node+1]-=ThisForce;
}
}

for(Node=1;Node<=LastFabricNode;Node++)
{
if(Node!=LastFabricNode||Sail==0)
{
MassOfThisNode=Mass[Node]+FabricMassPerNode;
for(xyz=0;xyz<=1;xyz++)v[xyz][Node]+=FabricForce[xyz][Node]*FabricDeltat/MassOfThisNode;

for(xyz=0;xyz<=1;xyz++)Coord[0][xyz][Node]+=v[xyz][Node]*FabricDeltat;
}
}
}
}

void Motion(void)
{
for(Node=LastFabricNode+1;Node<=LastNodeOutOfTank;Node++)
{
if(InTank[Node]!=1)
{
if(NodeType[Node]!=2)
{
for(xyz=0;xyz<=1;xyz++)v[xyz][Node]+=Force[xyz][Node]*Deltat/Mass[Node];
}

for(xyz=0;xyz<=1;xyz++)Coord[0][xyz][Node]+=v[xyz][Node]*Deltat;
for(TimeBeforeNow=LastTimeBeforeNow;TimeBeforeNow>=1;TimeBeforeNow--)  
{
for(xyz=0;xyz<=1;xyz++)Coord[TimeBeforeNow][xyz][Node]=Coord[TimeBeforeNow-1][xyz][Node];
}
}
}
}

void SaveState(void)
{
int Grey,dxfGreyMap[6];

ofstream Julia("State.txt");
Julia<<InitialTotalMass<<"  "<<Cycle<<"\n";
if(Fabric==1)Julia<<FabricPreStress<<"\n";
for(Node=0;Node<=LastNode;Node++)
{
Julia<<InTank[Node]<<"\n";
if(InTank[Node]!=1)
{
Julia<<NodeType[Node]<<" "<<Mass[Node]<<" "<<StreakOrNot[Node]<<" "<<Omega[Node]<<"\n";
for(xyz=0;xyz<=1;xyz++)Julia<<Coord[0][xyz][Node]<<"  "<<v[xyz][Node]<<"\n";
}
}
Julia.close();
cout<<"Text file written\n";

dxfGreyMap[0]=254;
dxfGreyMap[1]=253;
dxfGreyMap[2]=9;
dxfGreyMap[3]=8;
dxfGreyMap[4]=250;
dxfGreyMap[5]=19;

ofstream Madeleine("Flow.dxf");
Madeleine<<"0\nSECTION\n2\nENTITIES\n";
for(Node=0;Node<=LastNodeOutOfTank;Node++)
{
if(Node>LastFabricNode)
{
if(InTank[Node]==0)
{
if(StreakOrNot[Node]==1)
{
DrawRadius=0.2*sqrt(DistSqTimesBasicMassOverMassElastic*Mass[Node]);
Madeleine<<"0\nCIRCLE\n8\nStreak\n";
Madeleine<<"10\n"<<Coord[0][0][Node]<<"\n";
Madeleine<<"20\n"<<Coord[0][1][Node]<<"\n";
Madeleine<<"40\n"<<DrawRadius<<"\n";
}
for(TimeBeforeNow=1;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)  
{
Madeleine<<"0\nLINE\n8\nMini paths\n";
Madeleine<<"10\n"<<Coord[TimeBeforeNow-1][0][Node]<<"\n";
Madeleine<<"20\n"<<Coord[TimeBeforeNow-1][1][Node]<<"\n";
Madeleine<<"11\n"<<Coord[TimeBeforeNow  ][0][Node]<<"\n";
Madeleine<<"21\n"<<Coord[TimeBeforeNow  ][1][Node]<<"\n";
Grey=int((5*(LastTimeBeforeNow-TimeBeforeNow))/LastTimeBeforeNow);
if(Grey<0)Grey=0;
if(Grey>5)Grey=5;
Madeleine<<"62\n"<<dxfGreyMap[Grey]<<"\n";
}
}
}
}
Madeleine<<"0\nCIRCLE\n8\nCylinder\n";
Madeleine<<"10\n"<<-CylinderPosition<<"\n";
Madeleine<<"20\n"<<0.0<<"\n";
Madeleine<<"40\n"<<CylinderRadius<<"\n";
Madeleine<<"62\n1\n";

for(Node=1;Node<=LastFabricNode;Node++)
{
Madeleine<<"0\nLINE\n8\nFabric\n";
Madeleine<<"10\n"<<Coord[0][0][Node-1]<<"\n";
Madeleine<<"20\n"<<Coord[0][1][Node-1]<<"\n";
Madeleine<<"11\n"<<Coord[0][0][Node ]<<"\n";
Madeleine<<"21\n"<<Coord[0][1][Node ]<<"\n";
Madeleine<<"62\n1\n";
}

Madeleine<<"0\nENDSEC\n0\nEOF\n";
Madeleine.close();
cout<<"DXF file written\n";

Constance<<"0\nENDSEC\n0\nEOF\n";
Constance.close();
cout<<"Lift and drag DXF file written, end of program\n";
}

int CreateData(void)
{
int OldState;
float AngleOfAttack,FabricInitialRadius,TotalAngle,xFabricRoot,yFabricRoot;

Fabric=0;Sail=0;SavePictureOrNot=0;

Scale=0.1;//1.0 max
c=200.0;
cSq=c*c;

ZoneSize=5.0;
DistSqTimesBasicMassOverMassElastic=5.0;
DistSqTimesBasicMassOverMassViscous=3.0;

BasicMass=1.0;
cout<<"Value at edge of zone = "<<exp(-ZoneSize*ZoneSize/(DistSqTimesBasicMassOverMassElastic*BasicMass))<<"\n";

TimeBeforeNow=0;
Deltat=0.005;

LastFabricCycle=50;
FabricDeltat=Deltat/(1.0*LastFabricCycle);

if(Fabric==0)ReynoldsNumber=1000.0*(Scale/0.5);
        else ReynoldsNumber=5000.0*(Scale/0.5);
WindSpeed=0.2*c;

VorticityColour=1;Streak=1;Information=1;

InitialZoom=0.8/Scale;
InitialZoomNoInfo=0.9/Scale;
if(Information==1)Zoom=InitialZoom;else Zoom=InitialZoomNoInfo;

LastZone[0]=int(300*Scale)-1;
if(LastZone[0]>MaxLastZone0){cout<<LastZone[0]<<" is too many zones in the 0 direction\n";return 1;}
LastZone[1]=int(170*Scale)-1;
if(LastZone[1]>MaxLastZone1){cout<<LastZone[1]<<" is too many zones in the 1 direction\n";return 1;}

for(xyz=0;xyz<=1;xyz++)
{
if(LastZone[xyz]<0)LastZone[xyz]=0;
Boundary[xyz]=0.5*(LastZone[xyz]+1)*ZoneSize;
}
cout<<"Number of zones is "<<LastZone[0]+1<<" x "<<LastZone[1]+1<<"\n";

CylinderRadius=0.075*Boundary[0];
CylinderRadiusSq=CylinderRadius*CylinderRadius;
ControlCylinderRadiusSq=CylinderRadiusSq;
CylinderPosition=0.5*Boundary[0];
BaseOmega=2.0*WindSpeed/CylinderRadius;

InkPosition=CylinderPosition+0.0*CylinderRadius;
InkHalfDimension=1.05*CylinderRadius;

if(Fabric==1)
{
CylinderRadius=0.005*Boundary[0];
CylinderRadiusSq=CylinderRadius*CylinderRadius;
InkHalfDimension=3.0*CylinderRadius;
}

InkHalfDimensionSq=InkHalfDimension*InkHalfDimension;

if(Fabric==1)
{
FabricSlackLength=0.5*Boundary[0];
LastFabricNode=int(2.0*FabricSlackLength/sqrt(DistSqTimesBasicMassOverMassElastic));
if(LastFabricNode>MaxLastFabricNode){cout<<"Too many fabric nodes, last fabric node = "<<LastFabricNode<<"\n";return 1;}
cout<<"Last fabric node = "<<LastFabricNode<<"\n";
FabricElementSlackLength=FabricSlackLength/(1.0*LastFabricNode);
FabricPreStress=0.0;
}
else LastFabricNode=-1;

if(Fabric==0)mu_over_ro=WindSpeed*2.0*CylinderRadius/ReynoldsNumber;
        else mu_over_ro=WindSpeed*FabricSlackLength/ReynoldsNumber;

BoundaryZone=0.5*ZoneSize;

PictureSaveInterval=10;ForceSmooth=0.1;

InitialLastNodeOutOfTank=int(5.0*(LastZone[0]+1)*(LastZone[1]+1));
InitialNumberInTank=6*InitialLastNodeOutOfTank;
LastNode=InitialLastNodeOutOfTank+InitialNumberInTank;
cout<<"Last particle is "<<LastNode<<"\n";
if(LastNode>AbsoluteLastNode){cout<<"which is too many particles\n";return 1;}
cout<<"Initial average number of particles per zone is "<<InitialLastNodeOutOfTank/(1.0*(LastZone[0]+1)*(LastZone[1]+1))<<"\n";

cout<<"Total mass out of tank = "<<(1.0+1.0*InitialLastNodeOutOfTank)*BasicMass<<"\n";

OldState=0;//cout<<"To read from existing file type 1, otherwise 0\n";cin>>OldState;

for(Node=0;Node<=LastCylinderNode;Node++)
{
CylinderCoord[0][Node]=-CylinderPosition+CylinderRadius*cos(2.0*PI*Node/(LastCylinderNode+1));
CylinderCoord[1][Node]=					 CylinderRadius*sin(2.0*PI*Node/(LastCylinderNode+1));
}

if(OldState==1)
{
cout<<"Reading from file\n";
ifstream Maud("State.txt");
Maud>>InitialTotalMass>>Cycle;
if(Fabric==1)Maud>>FabricPreStress;
for(Node=0;Node<=LastNode;Node++)
{
Maud>>InTank[Node];
if(InTank[Node]!=1)
{
Maud>>NodeType[Node]>>Mass[Node]>>StreakOrNot[Node]>>Omega[Node];
for(xyz=0;xyz<=1;xyz++)Maud>>Coord[0][xyz][Node]>>v[xyz][Node];
}
}
Maud.close();
cout<<"Finished reading from file\n";
}
else
{
Cycle=-1;
InitialTotalMass=0.0;

if(Fabric==1&&Sail==1)
{
AngleOfAttack=10.0*PI/180.0;
TotalAngle=30.0*PI/180.0;
FabricInitialRadius=FabricSlackLength/TotalAngle;
xFabricRoot=0.0;
yFabricRoot=0.0;
}

for(Node=0;Node<=LastNode;Node++)
{
if(Node<=LastFabricNode)
{
InTank[Node]=0;
NodeType[Node]=0;
StreakOrNot[Node]=0;

if(Sail==0)
{
Coord[0][0][Node]=-CylinderPosition+CylinderRadius+(FabricSlackLength*Node)/(LastFabricNode);
Coord[0][1][Node]=0.0;
}
else
{
Coord[0][0][Node]=-xFabricRoot+
FabricInitialRadius*sin(AngleOfAttack+TotalAngle*(2.0*Node-1.0*LastFabricNode)/(2.0*LastFabricNode));
Coord[0][1][Node]=-yFabricRoot+
FabricInitialRadius*cos(AngleOfAttack+TotalAngle*(2.0*Node-1.0*LastFabricNode)/(2.0*LastFabricNode));
if(Node==0)
{
xFabricRoot=-(-CylinderPosition+CylinderRadius)+Coord[0][0][Node];
yFabricRoot=Coord[0][1][Node];
Coord[0][0][Node]-=xFabricRoot;
Coord[0][1][Node]-=yFabricRoot;
}
}

for(xyz=0;xyz<=1;xyz++)v[xyz][Node]=0.0;
}
else
{
if(Node>InitialLastNodeOutOfTank)InTank[Node]=1;
else
{
InTank[Node]=0;
NodeType[Node]=0;
Mass[Node]=BasicMass;
StreakOrNot[Node]=0;
for(xyz=0;xyz<=1;xyz++)
Coord[0][xyz][Node]=Boundary[xyz]*((2.0*rand())/(1.0*RAND_MAX)-1.0);

InitialTotalMass+=Mass[Node];
for(xyz=0;xyz<=1;xyz++)
{
if(xyz==0&&NodeType[Node]!=1)v[xyz][Node]=WindSpeed;
else v[xyz][Node]=0.0;
}
}
}
}
}

for(Node=0;Node<=LastNode;Node++)
{
for(TimeBeforeNow=1;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)  
{
for(xyz=0;xyz<=1;xyz++)Coord[TimeBeforeNow][xyz][Node]=Coord[0][xyz][Node];
}
}

if(Cycle==-1)FileNumber=-1;else FileNumber=int(Cycle/PictureSaveInterval);
roAverage=InitialTotalMass/(4.0*Boundary[0]*Boundary[1]);

if(Fabric==1)
{
FabricMassPerNode=2.0*roAverage*FabricSlackLength/(1.0*LastFabricNode);
for(Node=0;Node<=LastFabricNode;Node++)Mass[Node]=5.0*roAverage*FabricSlackLength*FabricSlackLength/(1.0*LastFabricNode*LastFabricNode);
if(Sail==1)
{
FabricEA=10.0*(FabricMassPerNode/(Deltat*Deltat))*FabricElementSlackLength;
FabricPreStressFactor=0.001;
FabricDamping=10.0*sqrt(FabricMassPerNode*FabricEA/FabricElementSlackLength);
}
else
{
//FabricEA=0.05*(FabricMassPerNode/(Deltat*Deltat))*FabricElementSlackLength;
//FabricPreStressFactor=0.01;
FabricEA=1.0*(FabricMassPerNode/(Deltat*Deltat))*FabricElementSlackLength;
FabricPreStressFactor=0.0001;
FabricDamping=10.0*sqrt(FabricMassPerNode*FabricEA/FabricElementSlackLength);
}
}

TimeFactor=Deltat*WindSpeed/(2.0*CylinderRadius);
DragLiftFactor=1.0/(0.5*roAverage*WindSpeed*WindSpeed*2.0*CylinderRadius);
for(xyz=0;xyz<=1;xyz++)NetForce[xyz]=0.0;

Constance<<"0\nSECTION\n2\nENTITIES\n";
return 0;
}

float My_exp(float MyArg)
{
int MyIntArg;
MyIntArg=int(-ArgMultiplier*MyArg);
if(MyIntArg<0)return 1.0;
if(MyIntArg<=LastExp)return Stored_exp[MyIntArg];
else return 0.0;
}

void SetUpExp(void)
{
int MyIntArg;
ArgMultiplier=0.05*LastExp;
for(MyIntArg=0;MyIntArg<=LastExp;MyIntArg++)Stored_exp[MyIntArg]=exp(-(0.5+MyIntArg)/ArgMultiplier);
}
