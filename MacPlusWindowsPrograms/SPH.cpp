//Written by Chris J K Williams, University of Bath, UK

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include <stdio.h>
#include <string.h>
#include "Graphics.h"
#define AbsoluteLastParticle 127500
#define MaxLastZone0 299
#define MaxLastZone1   0
#define MaxLastZone2 169
#define MaxLastInZone 200
#define LastExp 10000
#define LastTimeBeforeNow 1
void DoTheDrawing(void);
void SortIntoZones(void);
void Boundaries(void);
//void Wall(void);
void TheCalculation(int,int);
void Interaction(void);
void Motion(void);
//void StartParticleOnLeft(int);
int CreateData(void);
float My_exp(float);
float PressureOverDensitySq(float,float);
void SetUpExp(void);
int xyzPlus1(int);
int xyzPlus2(int);
float  PastCoords[LastTimeBeforeNow+1][3][AbsoluteLastParticle+1],Coord[3][AbsoluteLastParticle+1],Velocity[3][AbsoluteLastParticle+1],
Force[3][AbsoluteLastParticle+1],Moment[3][AbsoluteLastParticle+1],Omega[3][AbsoluteLastParticle+1],
Mass[AbsoluteLastParticle+1],Density[AbsoluteLastParticle+1],hSq[AbsoluteLastParticle+1],h[AbsoluteLastParticle+1],
Stored_exp[LastExp+1],
Boundary[3],
PI3over2,ZoneSize,BoundaryOffset,Deltat,
ReynoldsNumber,MachNumber,TotalMass,
FundamentalSpeed,Gravity,
MaxExpArg,ArgMultiplier,minRadiusSquared,hSq_ViscosityOverPressure,hSq_ViscosityOverPressureRaised,
c0Sq,ro0,ro_ZeroPressure,mu_over_roSq,
//DistanceFromCentreSq,DistanceFromInkCentreSq,
FundamentalLength,//CylinderRadiusSq,ControlCylinderRadiusSq,CylinderPosition,
//InkHalfDimension,InkHalfDimensionSq,InkPosition,
ZoneSizeSq,TwelveZoneSizeSq;//,WallPosition;

int    ParticleType[AbsoluteLastParticle+1],
StreakOrNot[AbsoluteLastParticle+1],
LastInZone[MaxLastZone0+1][MaxLastZone1+1][MaxLastZone2+1],
InZone[MaxLastZone0+1][MaxLastZone1+1][MaxLastZone2+1][MaxLastInZone+1],
Zone[3],LastZone[3],
LastParticle,Cycle,
DrawingSaveInterval,
LengthOfString,MyCharacter,PathUpdate,BlackBackground;
char NumericalValue[101];
char *MyText;
int main(int argc, char* argv[])
{
	float PI=4.0*atan(1.0);
	PI3over2=pow(double(PI),1.5);
	if(CreateData())return 0;
	glutInit(&argc,argv);
	BlackBackground=0;
	if(BlackBackground==1)SetUpGraphics(0.0,0.0,0.0);else SetUpGraphics(1.0,1.0,1.0);
	return 0;
}

static void Draw(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    Cycle++;
	Boundaries();
	//Wall();
	DoTheDrawing();
	if(DrawingSaveInterval>=0&&Cycle==DrawingSaveInterval*int(Cycle/DrawingSaveInterval))SavePicture();
	SortIntoZones();
	Interaction();
	Motion();
	//DoTheDrawing();
	//if(DrawingSaveInterval>=0&&Cycle==DrawingSaveInterval*int(Cycle/DrawingSaveInterval))SavePicture();
}

void DoTheDrawing(void)
{
	float rv[3],Proportion,ColourFactor,MyRed,MyGreen,MyBlue;
	glLineWidth(1.0);
	glBegin(GL_LINE_STRIP);
	if(BlackBackground==1)glColor4f(1.0,1.0,1.0,1.0);else glColor4f(0.0,0.0,0.0,1.0);
	
	rv[0]=-Boundary[0];rv[1]=-Boundary[1];rv[2]=-Boundary[2];glVertex3fv(rv);
	
	rv[0]=-rv[0];glVertex3fv(rv);
	rv[1]=-rv[1];glVertex3fv(rv);
	rv[0]=-rv[0];glVertex3fv(rv);
	rv[1]=-rv[1];glVertex3fv(rv);
	
	rv[2]=-rv[2];glVertex3fv(rv);
	
	rv[0]=-rv[0];glVertex3fv(rv);
	rv[1]=-rv[1];glVertex3fv(rv);
	rv[0]=-rv[0];glVertex3fv(rv);
	rv[1]=-rv[1];glVertex3fv(rv);
	
	glEnd();
	
	glBegin(GL_LINE_STRIP);rv[0]=-rv[0];glVertex3fv(rv);rv[2]=-rv[2];glVertex3fv(rv);glEnd();
	glBegin(GL_LINE_STRIP);rv[1]=-rv[1];glVertex3fv(rv);rv[2]=-rv[2];glVertex3fv(rv);glEnd();
	glBegin(GL_LINE_STRIP);rv[0]=-rv[0];glVertex3fv(rv);rv[2]=-rv[2];glVertex3fv(rv);glEnd();
	
	if(PathTuftStreak==0)
	{
		glLineWidth(1.0);
		for(int Particle=0;Particle<=LastParticle;Particle++)
		{
			if(ParticleType[Particle]!=1)
			{
				glBegin(GL_LINE_STRIP);
				for(int TimeBeforeNow=0;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)
				{
					Proportion=(1.0*TimeBeforeNow)/(1.0*LastTimeBeforeNow);
					if(BlackBackground==1)glColor4f(1.0,1.0,1.0,1.0-1.0*Proportion);
					else glColor4f(0.0,0.0,0.0,1.0-1.0*Proportion);
					for(int xyz=0;xyz<=2;xyz++)rv[xyz]=PastCoords[TimeBeforeNow][xyz][Particle];
					glVertex3fv(rv);
				}
				glEnd();
			}
		}
		/*for(int Particle=0;Particle<=LastParticle;Particle++)
		 {
		 if(ParticleType[Particle]!=1)
		 {
		 float ThisPointSize=0.5*h[Particle]/Scale;
		 if(ThisPointSize<1.0)ThisPointSize=1.0;
		 glPointSize(ThisPointSize);
		 glColor4f(1.0,0.0,0.0,0.5);
		 glBegin(GL_POINTS);
		 for(int xyz=0;xyz<=2;xyz++)rv[xyz]=Coord[xyz][Particle];
		 glVertex3fv(rv);
		 glEnd();
		 }
		 }*/
	}
	if(PathTuftStreak==1)
	{
		glLineWidth(1.0);
		for(int Particle=0;Particle<=LastParticle;Particle++)
		{
			if(ParticleType[Particle]!=1)
			{
				float SpeedSq=0.0;
				for(int xyz=0;xyz<=2;xyz++)SpeedSq+=Velocity[xyz][Particle]*Velocity[xyz][Particle];
				float TuftValue=0.5/sqrt(SpeedSq);
				glBegin(GL_LINE_STRIP);
				for(int End=-1;End<=1;End+=2)
				{
					if(BlackBackground==1)glColor4f(1.0,1.0,1.0,1.0);
					else glColor4f(0.0,0.0,0.0,1.0);
					for(int xyz=0;xyz<=2;xyz++)rv[xyz]=Coord[xyz][Particle]+TuftValue*End*Velocity[xyz][Particle];
					glVertex3fv(rv);
				}
				glEnd();
			}
		}
	}
	if(PathTuftStreak==2)
	{
		for(int Particle=0;Particle<=LastParticle;Particle++)
		{
			//if(StreakOrNot[Particle]==1&&ParticleType[Particle]==0)
			{
				//float ThisPointSize=1.0*h[Particle]*Zoom;
				float ThisPointSize=1.0*h[Particle]*Zoom;
				if(ThisPointSize<1.0)ThisPointSize=1.0;
				glPointSize(ThisPointSize);
				if(VorticityColour==1)
				{
					ColourFactor=Omega[2][Particle]*FundamentalLength/FundamentalSpeed;
					if(ColourFactor>+1.0)ColourFactor=+1.0;
					if(ColourFactor<-1.0)ColourFactor=-1.0;
					if(BlackBackground==1)
					{
						if(ColourFactor>0.0)
						{
							MyRed=1.0;MyBlue=1.0-ColourFactor;MyGreen=MyBlue;
						}
						else
						{
							MyBlue=1.0;MyRed=1.0+ColourFactor;MyGreen=MyRed;
						}
					}
					else
					{
						if(ColourFactor>0.0)
						{
							MyRed=ColourFactor;MyGreen=0.0;MyBlue=0.0;
						}
						else
						{
							MyRed=0.0;MyGreen=0.0;MyBlue=-ColourFactor;
						}
					}
				}
				else
				{
					if(BlackBackground==1)
					{
						MyRed=1.0;MyGreen=1.0;MyBlue=1.0;
					}
					else
					{
						MyRed=0.0;MyGreen=0.0;MyBlue=1.0;
						if(Density[Particle]<ro_ZeroPressure){MyRed=0.5;MyGreen=0.5;MyBlue=0.5;}
					}
				}
				glBegin(GL_POINTS);
				glColor4f(MyRed,MyGreen,MyBlue,1.0);
				for(int xyz=0;xyz<=2;xyz++)rv[xyz]=Coord[xyz][Particle];
				glVertex3fv(rv);
				glEnd();
			}
		}
	}
	if(Information==1)
	{
		glPushMatrix();
		glLoadIdentity();
		if(BlackBackground==1)glColor4f(1.0,1.0,1.0,0.5);else glColor4f(0.0,0.0,0.0,0.5);
		glRasterPos2f(-0.9,-0.9);
		/*MyText="Mach number ";
		 LengthOfString=strlen(MyText);
		 for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
		 sprintf(NumericalValue,"%8.2f",MachNumber);
		 LengthOfString=strlen(NumericalValue);
		 for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
		 MyText="         Reynolds number  ";
		 LengthOfString=strlen(MyText);
		 for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
		 sprintf(NumericalValue,"%i",int(ReynoldsNumber));
		 LengthOfString=strlen(NumericalValue);
		 for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);*/
		MyText="Cycle  ";
		LengthOfString=strlen(MyText);
		for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
		sprintf(NumericalValue,"%i",Cycle);
		LengthOfString=strlen(NumericalValue);
		for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
		/*MyText="         Non-dimensional time ";
		 LengthOfString=strlen(MyText);
		 for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
		 sprintf(NumericalValue,"%8.2f",Deltat*Cycle*WindSpeed/(2.0*CylinderRadius));
		 LengthOfString=strlen(NumericalValue);
		 for (MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);*/
		
		glPopMatrix();
	}
	glutSwapBuffers();
}

void Boundaries(void)
{
	/*for(int Particle=0;Particle<=LastParticle;Particle++)
	 {
	 if(Coord[0][Particle]>+Boundary[0]||Coord[0][Particle]<-Boundary[0])StartParticleOnLeft(Particle);
	 }*/
	for(int Particle=0;Particle<=LastParticle;Particle++)
	{
		ParticleType[Particle]=0;
		/*Density[Particle]=Mass[Particle]/(PI3over2*pow(double(2.0*hSq[Particle]),1.5));
		 for(int xyz=0;xyz<=2;xyz++)
		 {
		 Force[xyz][Particle]=0.0;
		 Moment[xyz][Particle]=0.0;
		 }*/
		/*DistanceFromInkCentreSq
		 =(Coord[0][Particle]+InkPosition)*(Coord[0][Particle]+InkPosition)+Coord[1][Particle]*Coord[1][Particle];
		 DistanceFromCentreSq
		 =(Coord[0][Particle]+CylinderPosition)*(Coord[0][Particle]+CylinderPosition)+Coord[1][Particle]*Coord[1][Particle];
		 if(DistanceFromInkCentreSq<InkHalfDimensionSq)StreakOrNot[Particle]=1;
		 if(DistanceFromCentreSq<CylinderRadiusSq)ParticleType[Particle]=1;
		 if((Boundary[0]+Coord[0][Particle]<ZoneSize)
		 ||(Boundary[0]-Coord[0][Particle]<ZoneSize))
		 {
		 ParticleType[Particle]=2;
		 for(int xyz=0;xyz<=2;xyz++)
		 {
		 Velocity[xyz][Particle]=0.0;Omega[xyz][Particle]=0.0;
		 }
		 Velocity[0][Particle]=WindSpeed;
		 }*/
		for(int xyz=0;xyz<=2;xyz++)
		{
			/*while(Coord[xyz][Particle]>+Boundary[xyz])
			 {
			 Coord[xyz][Particle]-=2.0*Boundary[xyz];
			 for(int TimeBeforeNow=0;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)
			 PastCoords[TimeBeforeNow][xyz][Particle]-=2.0*Boundary[xyz];
			 }
			 while(Coord[xyz][Particle]<-Boundary[xyz])
			 {
			 Coord[xyz][Particle]+=2.0*Boundary[xyz];
			 for(int TimeBeforeNow=0;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)
			 PastCoords[TimeBeforeNow][xyz][Particle]+=2.0*Boundary[xyz];
			 }*/
			if(Coord[xyz][Particle]>+Boundary[xyz])
			{
				Coord[xyz][Particle]=+(2.0*Boundary[xyz]-BoundaryOffset)-Coord[xyz][Particle];
				if(Coord[xyz][Particle]<-Boundary[xyz])Coord[xyz][Particle]=+(Boundary[xyz]-BoundaryOffset);
				if(Velocity[xyz][Particle]>0.0)Velocity[xyz][Particle]=-Velocity[xyz][Particle];
			}
			if(Coord[xyz][Particle]<-Boundary[xyz])
			{
				Coord[xyz][Particle]=-(2.0*Boundary[xyz]-BoundaryOffset)-Coord[xyz][Particle];
				if(Coord[xyz][Particle]>+Boundary[xyz])Coord[xyz][Particle]=-(Boundary[xyz]-BoundaryOffset);
				if(Velocity[xyz][Particle]<0.0)Velocity[xyz][Particle]=-Velocity[xyz][Particle];
			}
		}
	}
}

/*void Wall(void)
 {
 for(int Particle=0;Particle<=LastParticle;Particle++)
 {
 if(Coord[0][Particle]>WallPosition)
 {
 Coord[0][Particle]=2.0*WallPosition-BoundaryOffset-Coord[0][Particle];
 if(Coord[0][Particle]<-Boundary[0])Coord[0][Particle]=+WallPosition-BoundaryOffset;
 if(Velocity[0][Particle]>0.0)Velocity[0][Particle]=-Velocity[0][Particle];
 }
 }
 }*/

/*void StartParticleOnLeft(int Particle)
 {
 Coord[0][Particle]=-Boundary[0]+WindSpeed*Deltat*(1.0*rand())/(1.0*RAND_MAX);
 Coord[1][Particle]=Boundary[1]*((2.0*rand())/(1.0*RAND_MAX)-1.0);
 for(int TimeBeforeNow=0;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)
 {
 for(int xyz=0;xyz<=2;xyz++)PastCoords[TimeBeforeNow][xyz][Particle]=Coord[xyz][Particle];
 }
 Velocity[0][Particle]=WindSpeed;Velocity[1][Particle]=0.0;
 for(int xyz=0;xyz<=2;xyz++)Omega[xyz][Particle]=0.0;
 StreakOrNot[Particle]=0;
 }*/

void SortIntoZones(void)
{
	for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++)
	{
		for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++)
		{
			for(Zone[2]=0;Zone[2]<=LastZone[2];Zone[2]++)LastInZone[Zone[0]][Zone[1]][Zone[2]]=-1;
		}
	}
	for(int Particle=0;Particle<=LastParticle;Particle++)
	{
		for(int xyz=0;xyz<=2;xyz++)
		{
			float floatZone=(Coord[xyz][Particle]+Boundary[xyz])/ZoneSize;
			if(floatZone<0.0)
			{
				Zone[xyz]=0;cout<<"Negative floatZone\n";
			}
			else
			{
				Zone[xyz]=int(floatZone);
				if(Zone[xyz]<0)Zone[xyz]=0;
				if(Zone[xyz]>LastZone[xyz])Zone[xyz]=LastZone[xyz];
			}
		}
		LastInZone[Zone[0]][Zone[1]][Zone[2]]++;
		if(LastInZone[Zone[0]][Zone[1]][Zone[2]]<=MaxLastInZone)
			InZone[Zone[0]][Zone[1]][Zone[2]][LastInZone[Zone[0]][Zone[1]][Zone[2]]]=Particle;
	}
	for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++)
	{
		for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++)
		{
			for(Zone[2]=0;Zone[2]<=LastZone[2];Zone[2]++)
			{
				if(LastInZone[Zone[0]][Zone[1]][Zone[2]]>MaxLastInZone)
				{
					cout<<"Number in Zone "<<Zone[0]<<" "<<Zone[1]<<" "<<Zone[2]<<" is "<<LastInZone[Zone[0]][Zone[1]][Zone[2]]<<"\n";
					LastInZone[Zone[0]][Zone[1]][Zone[2]]=MaxLastInZone;
				}
			}
		}
	}
}

void Interaction(void)
{
	for(int Particle=0;Particle<=LastParticle;Particle++)
	{
		ParticleType[Particle]=0;
		Density[Particle]=Mass[Particle]/(PI3over2*pow(double(2.0*hSq[Particle]),1.5));
		for(int xyz=0;xyz<=2;xyz++)
		{
			Force[xyz][Particle]=0.0;
			Moment[xyz][Particle]=0.0;
		}
	}
	for(int DensityCalcIs0=0;DensityCalcIs0<=1;DensityCalcIs0++)
	{
		for(Zone[0]=0;Zone[0]<=LastZone[0];Zone[0]++)
		{
			for(Zone[1]=0;Zone[1]<=LastZone[1];Zone[1]++)
			{
				for(Zone[2]=0;Zone[2]<=LastZone[2];Zone[2]++)
				{
					if(LastInZone[Zone[0]][Zone[1]][Zone[2]]>=0)
					{
						for(int ParticleInZone=0;ParticleInZone<=LastInZone[Zone[0]][Zone[1]][Zone[2]];ParticleInZone++)
						{
							int Particle_a=InZone[Zone[0]][Zone[1]][Zone[2]][ParticleInZone];
							TheCalculation(DensityCalcIs0,Particle_a);
						}
					}
				}
			}
		}
	}
}

void TheCalculation(int DensityCalcIs0,int Particle_a)
{
	float delta_r_ab[3],delta_Velocity_ab[3],delta_Velocity_abNET[3],Moment_a[3],Moment_b[3],Pab[3],Qab[3],rv[3];
	int WhichZone[3],ThisWhichZone[3],MirrorBoundary[3];
	
	float ha2=hSq[Particle_a];
	
	for(WhichZone[0]=Zone[0]-1;WhichZone[0]<=Zone[0]+1;WhichZone[0]++)
	{
		for(WhichZone[1]=Zone[1]-1;WhichZone[1]<=Zone[1]+1;WhichZone[1]++)
		{
			for(WhichZone[2]=Zone[2]-1;WhichZone[2]<=Zone[2]+1;WhichZone[2]++)
			{
				int AnyMirror=0;
				for(int xyz=0;xyz<=2;xyz++)
				{
					ThisWhichZone[xyz]=WhichZone[xyz];MirrorBoundary[xyz]=0;
					if(ThisWhichZone[xyz]<0)
					{
						//ThisWhichZone[xyz]=LastZone[xyz];MirrorBoundary[xyz]=-1;
						ThisWhichZone[xyz]=0;MirrorBoundary[xyz]=-1;AnyMirror=1;
					}
					if(ThisWhichZone[xyz]>LastZone[xyz])
					{
						//ThisWhichZone[xyz]=0;MirrorBoundary[xyz]=+1;
						ThisWhichZone[xyz]=LastZone[xyz];MirrorBoundary[xyz]=+1;AnyMirror=1;
					}
				}
				if(LastInZone[ThisWhichZone[0]][ThisWhichZone[1]][ThisWhichZone[2]]>=0)
				{
					for(int otherParticleInZone=0;otherParticleInZone<=LastInZone[ThisWhichZone[0]][ThisWhichZone[1]][ThisWhichZone[2]];otherParticleInZone++)
					{
						int Particle_b=InZone[ThisWhichZone[0]][ThisWhichZone[1]][ThisWhichZone[2]][otherParticleInZone];
						if(Particle_a>Particle_b||AnyMirror!=0)
						{
							float hb2=hSq[Particle_b];
							float rab2=0.0;
							for(int xyz=0;xyz<=2;xyz++)
							{
								if(MirrorBoundary[xyz]==0)delta_r_ab[xyz]=Coord[xyz][Particle_a]-Coord[xyz][Particle_b];
								else delta_r_ab[xyz]=Coord[xyz][Particle_a]+Coord[xyz][Particle_b]-2.0*MirrorBoundary[xyz]*Boundary[xyz];
								rab2+=delta_r_ab[xyz]*delta_r_ab[xyz];
							}
							if(rab2>TwelveZoneSizeSq)
							{
								cout<<"\nSeparation bigger than it should be.";
								cout<<"\nDeltax/ZoneSize = "<<delta_r_ab[0]/ZoneSize;
								cout<<"\nDeltay/ZoneSize = "<<delta_r_ab[1]/ZoneSize;
								cout<<"\nDeltaz/ZoneSize = "<<delta_r_ab[2]/ZoneSize<<"\n";
							}
							if(rab2<ZoneSizeSq&&rab2>minRadiusSquared)
							{
								float ratio_ab_Pressure=rab2/(ha2+hb2);
								if(DensityCalcIs0==0)
									Density[Particle_a]+=Mass[Particle_b]*My_exp(-ratio_ab_Pressure)/(PI3over2*pow(double(ha2+hb2),1.5));
								else
								{
									float ScalarProduct=0.0;
									for(int xyz=0;xyz<=2;xyz++)
									{
										if(MirrorBoundary[xyz]==0)delta_Velocity_ab[xyz]=Velocity[xyz][Particle_a]-Velocity[xyz][Particle_b];
										else delta_Velocity_ab[xyz]=Velocity[xyz][Particle_a]+Velocity[xyz][Particle_b];
										ScalarProduct+=delta_r_ab[xyz]*delta_Velocity_ab[xyz];
									}
									float ha2Plus_hb2Raised=pow(double(ha2+hb2),2.5);
									float eta_ab_Pressure=Mass[Particle_a]*Mass[Particle_b]*My_exp(-ratio_ab_Pressure)/(PI3over2*ha2Plus_hb2Raised);
									float ratio_ab_Viscosity=rab2/(hSq_ViscosityOverPressure*(ha2+hb2));
									float eta_ab_Viscosity=Mass[Particle_a]*Mass[Particle_b]
									*My_exp(-ratio_ab_Viscosity)/(PI3over2*hSq_ViscosityOverPressureRaised*ha2Plus_hb2Raised);
									
									if(ShowArms==1)
									{
										if(MirrorBoundary[0]==0&&MirrorBoundary[1]==0)
										{
											glColor4f(0.5,0.5,0.5,My_exp(-ratio_ab_Viscosity));
											glBegin(GL_LINE_STRIP);
											for(int xyz=0;xyz<=2;xyz++)rv[xyz]=Coord[xyz][Particle_a];glVertex3fv(rv);
											for(int xyz=0;xyz<=2;xyz++)rv[xyz]=Coord[xyz][Particle_b];glVertex3fv(rv);
											glEnd();
										}
									}
									
									for(int xyz=0;xyz<=2;xyz++)
									{
										int xyzp1=xyzPlus1(xyz);
										int xyzp2=xyzPlus2(xyz);
										delta_Velocity_abNET[xyz]=delta_Velocity_ab[xyz]-
										(+(Omega[xyzp1][Particle_a]+Omega[xyzp1][Particle_b])*delta_r_ab[xyzp2]
										 -(Omega[xyzp2][Particle_a]+Omega[xyzp2][Particle_b])*delta_r_ab[xyzp1])/2.0;
									}
									
									float Damping=2.0*mu_over_roSq*eta_ab_Viscosity;
									
									for(int xyz=0;xyz<=2;xyz++)Qab[xyz]=-Damping*2.0*
										(delta_Velocity_abNET[xyz]*ratio_ab_Viscosity-ScalarProduct*delta_r_ab[xyz]/(hSq_ViscosityOverPressure*(ha2+hb2)));
									
									for(int xyz=0;xyz<=2;xyz++)
									{
										int xyzp1=xyzPlus1(xyz);
										int xyzp2=xyzPlus2(xyz);
										float ArmEndMomentSumOver2=(Qab[xyzp1]*delta_r_ab[xyzp2]-Qab[xyzp2]*delta_r_ab[xyzp1])/2.0;
										float ArmEndMomentDifOver2=Damping*rab2*ratio_ab_Viscosity
										*(Omega[xyz][Particle_a]-Omega[xyz][Particle_b])/2.0;
										Moment_a[xyz]=ArmEndMomentSumOver2-ArmEndMomentDifOver2;
										Moment_b[xyz]=ArmEndMomentSumOver2+ArmEndMomentDifOver2;
									}
									for(int xyz=0;xyz<=2;xyz++)
									{
										Pab[xyz]=4.0*PressureOverDensitySq(Density[Particle_a],Density[Particle_b])*eta_ab_Pressure*delta_r_ab[xyz];
										Force[xyz][Particle_a]+=Pab[xyz]+Qab[xyz];
										Moment[xyz][Particle_a]+=Moment_a[xyz];
										if(AnyMirror==0)
										{
											Force[xyz][Particle_b]-=Pab[xyz]+Qab[xyz];
											Moment[xyz][Particle_b]+=Moment_b[xyz];
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

void Motion(void)
{
	for(int Particle=0;Particle<=LastParticle;Particle++)
	{
		if(ParticleType[Particle]!=2)
		{
			Velocity[2][Particle]-=Gravity*Deltat;
			for(int xyz=0;xyz<=2;xyz++)
			{
				Velocity[xyz][Particle]+=Force[xyz][Particle]*Deltat/Mass[Particle];
				Omega[xyz][Particle]+=Moment[xyz][Particle]*Deltat/(Mass[Particle]*hSq[Particle]);
				//Omega[xyz][Particle]+=Moment[xyz][Particle]*Deltat/(Mass[Particle]*hSq_ViscosityOverPressure*hSq[Particle]);//Moment of inertia is based on viscous size
			}
		}
		if(ParticleType[Particle]==1)
		{
			for(int xyz=0;xyz<=2;xyz++)Velocity[xyz][Particle]=0.01*Velocity[xyz][Particle];
		}
		for(int xyz=0;xyz<=2;xyz++)Coord[xyz][Particle]+=Velocity[xyz][Particle]*Deltat;
		//if(PathUpdate>=12)
		{
			PathUpdate=0;
			for(int TimeBeforeNow=LastTimeBeforeNow;TimeBeforeNow>=1;TimeBeforeNow--)
			{
				for(int xyz=0;xyz<=2;xyz++)PastCoords[TimeBeforeNow][xyz][Particle]=PastCoords[TimeBeforeNow-1][xyz][Particle];
			}
			for(int xyz=0;xyz<=2;xyz++)PastCoords[0][xyz][Particle]=Coord[xyz][Particle];
		}
		PathUpdate++;
	}
}

void SaveState(void)
{
	float DrawRadius;
	int Grey,dxfGreyMap[6];
	ofstream Julia("State.txt");
	cout<<"Saving files\n";
	Julia<<TotalMass<<"  "<<Cycle<<"\n";
	for(int Particle=0;Particle<=LastParticle;Particle++)
	{
		Julia<<"\n"<<ParticleType[Particle]<<" "<<Mass[Particle]<<" "<<hSq[Particle]<<" "<<StreakOrNot[Particle]<<"\n";
		for(int xyz=0;xyz<=2;xyz++)Julia<<Coord[xyz][Particle]<<"  "<<Velocity[xyz][Particle]<<" "<<Omega[xyz][Particle]<<"\n";
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
	for(int Particle=0;Particle<=LastParticle;Particle++)
	{
		if(StreakOrNot[Particle]==1)
		{
			for(Grey=4;Grey>=0;Grey--)
			{
				DrawRadius=0.1*sqrt(hSq[Particle])*(1.0+float(Grey));
				Madeleine<<"0\nCIRCLE\n8\nStreak "<<Grey<<"\n";
				Madeleine<<"10\n"<<Coord[0][Particle]<<"\n";
				Madeleine<<"20\n"<<Coord[1][Particle]<<"\n";
				Madeleine<<"40\n"<<DrawRadius<<"\n";
				Madeleine<<"62\n"<<dxfGreyMap[Grey]<<"\n";
			}
		}
		for(int TimeBeforeNow=0;TimeBeforeNow<=LastTimeBeforeNow-1;TimeBeforeNow++)
		{
			Madeleine<<"0\nLINE\n8\nMini paths\n";
			Madeleine<<"10\n"<<PastCoords[TimeBeforeNow][0][Particle]<<"\n";
			Madeleine<<"20\n"<<PastCoords[TimeBeforeNow][1][Particle]<<"\n";
			Madeleine<<"11\n"<<PastCoords[TimeBeforeNow+1][0][Particle]<<"\n";
			Madeleine<<"21\n"<<PastCoords[TimeBeforeNow+1][1][Particle]<<"\n";
			Grey=int((5*(LastTimeBeforeNow-TimeBeforeNow))/LastTimeBeforeNow);
			if(Grey<0)Grey=0;
			if(Grey>5)Grey=5;
			Madeleine<<"62\n"<<dxfGreyMap[Grey]<<"\n";
		}
	}
	/*Madeleine<<"0\nCIRCLE\n8\nCylinder\n";
	 Madeleine<<"10\n"<<-CylinderPosition<<"\n";
	 Madeleine<<"20\n"<<0.0<<"\n";
	 Madeleine<<"40\n"<<CylinderRadius<<"\n";
	 Madeleine<<"62\n0\n";*/
	Madeleine<<"0\nLINE\n8\nBorder\n";
	Madeleine<<"10\n"<<-Boundary[0]<<"\n";Madeleine<<"20\n"<<-Boundary[1]<<"\n";
	Madeleine<<"11\n"<<+Boundary[0]<<"\n";Madeleine<<"21\n"<<-Boundary[1]<<"\n";
	Madeleine<<"62\n0\n";
	Madeleine<<"0\nLINE\n8\nBorder\n";
	Madeleine<<"10\n"<<+Boundary[0]<<"\n";Madeleine<<"20\n"<<-Boundary[1]<<"\n";
	Madeleine<<"11\n"<<+Boundary[0]<<"\n";Madeleine<<"21\n"<<+Boundary[1]<<"\n";
	Madeleine<<"62\n0\n";
	Madeleine<<"0\nLINE\n8\nBorder\n";
	Madeleine<<"10\n"<<+Boundary[0]<<"\n";Madeleine<<"20\n"<<+Boundary[1]<<"\n";
	Madeleine<<"11\n"<<-Boundary[0]<<"\n";Madeleine<<"21\n"<<+Boundary[1]<<"\n";
	Madeleine<<"62\n0\n";
	Madeleine<<"0\nLINE\n8\nBorder\n";
	Madeleine<<"10\n"<<-Boundary[0]<<"\n";Madeleine<<"20\n"<<+Boundary[1]<<"\n";
	Madeleine<<"11\n"<<-Boundary[0]<<"\n";Madeleine<<"21\n"<<-Boundary[1]<<"\n";
	Madeleine<<"62\n0\n";
	Madeleine<<"0\nENDSEC\n0\nEOF\n";
	Madeleine.close();
	cout<<"DXF file written\n";
	//cout<<"Non-dimensional time = "<<Deltat*Cycle*WindSpeed/(2.0*CylinderRadius)
	cout<<"\nEnd of program\n";
}

int CreateData(void)
{
	float BasicDistSq,BasicMass,c0,mu_over_ro,factor;
	int OldState;
	Scale=0.2;//1.5 for Re = 300
	DrawingSaveInterval=-10;//Never saves if negative
	c0=1000.0;
	c0Sq=c0*c0;
	ZoneSize=5.0;
	BoundaryOffset=0.01*ZoneSize;
	//int NumberInZone=20;
	int NumberInZone=8;
	float DensityRatio=4.0;
	ZoneSizeSq=ZoneSize*ZoneSize;
	TwelveZoneSizeSq=12.0*ZoneSizeSq;
	BasicDistSq=2.0;
	BasicMass=1.0;
	hSq_ViscosityOverPressure=4.0;
	hSq_ViscosityOverPressureRaised=pow(double(hSq_ViscosityOverPressure),2.5);
	MaxExpArg=ZoneSize*ZoneSize/(BasicDistSq*pow(double(BasicMass),2.0/3.0));
	if(hSq_ViscosityOverPressure<1.0)MaxExpArg=MaxExpArg/hSq_ViscosityOverPressure;
	cout<<"Nominal value at edge of zone for pressure = "<<exp(-MaxExpArg)<<"\n";
	cout<<"Nominal value at edge of zone for viscosity = "<<exp(-MaxExpArg/hSq_ViscosityOverPressure)<<"\n";
	SetUpExp();
	Deltat=0.002;
	MachNumber=0.5;
	FundamentalSpeed=MachNumber*c0;
	//ReynoldsNumber=1.0*FundamentalSpeed*Scale;
	ReynoldsNumber=1.0*FundamentalSpeed*Scale;
	LastZone[0]=int(300*Scale)-1;
	if(LastZone[0]>MaxLastZone0)
	{
		cout<<LastZone[0]<<" is too large a last zone in the 0 direction\n";return 1;
	}
	//LastZone[1]=int(30*Scale)-1;
	LastZone[1]=0;
	if(LastZone[1]<0)
	{
		cout<<LastZone[1]<<" is negative\n";return 1;
	}
	if(LastZone[1]>MaxLastZone1)
	{
		cout<<LastZone[1]<<" is too large a last zone in the 1 direction\n";return 1;
	}
	LastZone[2]=int(170*Scale)-1;
	if(LastZone[2]>MaxLastZone2)
	{
		cout<<LastZone[2]<<" is too large a last zone in the 2 direction\n";return 1;
	}
	for(int xyz=0;xyz<=2;xyz++)
	{
		if(LastZone[xyz]<0)LastZone[xyz]=0;
		Boundary[xyz]=0.5*(LastZone[xyz]+1)*ZoneSize;
	}
	cout<<"Number of zones is "<<LastZone[0]+1<<" x "<<LastZone[1]+1<<" x "<<LastZone[2]+1<<"\n";
	FundamentalLength=2.0*Boundary[2];
	Gravity=FundamentalSpeed*FundamentalSpeed/(2.0*FundamentalLength);
	//CylinderRadiusSq=CylinderRadius*CylinderRadius;
	//ControlCylinderRadiusSq=CylinderRadiusSq;
	//CylinderPosition=0.5*Boundary[0];
	//InkPosition=CylinderPosition+0.0*CylinderRadius;
	//InkHalfDimension=1.05*CylinderRadius;
	//InkHalfDimensionSq=InkHalfDimension*InkHalfDimension;
	mu_over_ro=FundamentalSpeed*FundamentalLength/ReynoldsNumber;
	/*cout<<"Fundamentalspeed = "<<FundamentalSpeed<<"\n";
	cout<<"Fundamental length = "<<FundamentalLength<<"\n";
	cout<<"ReynoldsNumber = "<<ReynoldsNumber<<"\n";
	cout<<"Viscosity over density = "<<mu_over_ro<<"\n";*/
	LastParticle=int(float(NumberInZone*(LastZone[0]+1)*(LastZone[1]+1)*(LastZone[2]+1))/DensityRatio);
	cout<<"Last particle is "<<LastParticle<<"\n";
	if(LastParticle>AbsoluteLastParticle)
	{
		cout<<"which is too many particles\n";return 1;
	}
	//cout<<"\nType 0 to generate new data, or any other number to read from file.\n";cin>>OldState;
	OldState=0;
	if(OldState!=0)
	{
		cout<<"Reading from file\n";
		ifstream Maud("State.txt");
		Maud>>TotalMass>>Cycle;cout<<"Cycle = "<<Cycle<<"\n";
		for(int Particle=0;Particle<=LastParticle;Particle++)
		{
			Maud>>ParticleType[Particle]>>Mass[Particle]>>hSq[Particle]>>StreakOrNot[Particle];
			//StreakOrNot[Particle]=0;//To reset streaks
			for(int xyz=0;xyz<=2;xyz++)Maud>>Coord[xyz][Particle]>>Velocity[xyz][Particle]>>Omega[xyz][Particle];
		}
		Maud.close();
		cout<<"Finished reading from file\n";
	}
	else
	{
		Cycle=-1;
		TotalMass=0.0;
		float LowestMass=BasicMass;
		float HighestMass=0.0;
		//WallPosition=-0.5*Boundary[0];
		for(int Particle=0;Particle<=LastParticle;Particle++)
		{
			ParticleType[Particle]=0;
			for(;;)
			{
				factor=((1.0*rand())/(2.0*RAND_MAX)-1.0);
				if((1.0*rand())/(1.0*RAND_MAX)<exp(-factor*factor))break;
			}
			Mass[Particle]=BasicMass*(1.0-0.5*(factor+1.0));
			if(LowestMass>Mass[Particle])LowestMass=Mass[Particle];
			if(HighestMass<Mass[Particle])HighestMass=Mass[Particle];
			hSq[Particle]=BasicDistSq*pow(double(Mass[Particle]),2.0/3.0);
			StreakOrNot[Particle]=0;
			Coord[0][Particle]=-Boundary[0]+0.5*Boundary[0]*(1.0*rand())/(1.0*RAND_MAX);
			Coord[1][Particle]=Boundary[1]*((2.0*rand())/(1.0*RAND_MAX)-1.0);
			Coord[2][Particle]=-Boundary[2]+1.9*Boundary[2]*(1.0*rand())/(1.0*RAND_MAX);
			TotalMass+=Mass[Particle];
			for(int xyz=0;xyz<=2;xyz++)
				//{
				//if(xyz==0)Velocity[xyz][Particle]=WindSpeed;
				//else
				Velocity[xyz][Particle]=0.0;
			//}
		}
		cout<<"Lowest mass = "<<LowestMass<<"\nHighest mass = "<<HighestMass<<"\n";
	}
	for(int Particle=0;Particle<=LastParticle;Particle++)
	{
		h[Particle]=sqrt(hSq[Particle]);
		for(int TimeBeforeNow=0;TimeBeforeNow<=LastTimeBeforeNow;TimeBeforeNow++)
		{
			for(int xyz=0;xyz<=2;xyz++)PastCoords[TimeBeforeNow][xyz][Particle]=Coord[xyz][Particle];
		}
	}
	if(Cycle==-1)FileNumber=-1;else FileNumber=int(Cycle/DrawingSaveInterval);
	ro0=TotalMass/(8.0*Boundary[0]*Boundary[1]*Boundary[2]);
	ro_ZeroPressure=1.0*DensityRatio*ro0;
	mu_over_roSq=mu_over_ro/ro_ZeroPressure;
	minRadiusSquared=pow(double(ZoneSize/100.0),2);
	PathUpdate=1;
	return 0;
}

float My_exp(float MyArg)
{
	int MyIntArg;
	MyIntArg=int(-ArgMultiplier*MyArg);
	if(MyIntArg<0){cout<<"Positive exponential\n";return 0.0;}
	if(MyIntArg<=LastExp)return Stored_exp[MyIntArg];
	else {cout<<"Exponential argument too big and negative\n";return 0.0;}
}

void SetUpExp(void)
{
	int MyIntArg;
	ArgMultiplier=LastExp/MaxExpArg;
	for(MyIntArg=0;MyIntArg<=LastExp;MyIntArg++)Stored_exp[MyIntArg]=exp(-(0.5+MyIntArg)/ArgMultiplier);
	cout<<"Smallest exponential = "<<Stored_exp[LastExp]<<"\n\n";
}

int xyzPlus1(int this_xyz)
{
	if(this_xyz<2)return this_xyz+1;else return this_xyz-2;
}

int xyzPlus2(int this_xyz)
{
	if(this_xyz<1)return this_xyz+2;else return this_xyz-1;
}

float PressureOverDensitySq(float Density_a,float Density_b)
{
	/*float ro_ab=(ha2*Density_b+hb2*Density_a)/(ha2+hb2);
	if(ro_ab>ro_ZeroPressure)return c0Sq*(ro_ab-ro_ZeroPressure)/(ro_ab*ro_ab);
	else return 0.0;*/
	/*float x=Density_a/ro_ZeroPressure-1.0;
	float y=Density_b/ro_ZeroPressure-1.0;
	float a=0.01;
	float b=0.1*a;
	float X=sqrt(a*a+x*x)+x;
	float Y=sqrt(a*a+y*y)+y;
	return 0.707*c0Sq*(X+Y-sqrt((X-Y)*(X-Y)+4.0*b*b))/ro_ZeroPressure;*/
	float x=Density_a/ro_ZeroPressure-1.0;
	float y=Density_b/ro_ZeroPressure-1.0;
	if(x<0.0||y<0.0)return 0.0;
	else
	{
	if(x<y)return c0Sq*x/ro_ZeroPressure;
	else return c0Sq*y/ro_ZeroPressure;
	}
}

