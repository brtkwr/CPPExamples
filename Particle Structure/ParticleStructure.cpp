//Written by Chris J K Williams, University of Bath, UK

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include <stdio.h>
#include <string.h>
#include "Graphics.h"

#define MaxLastZone0             59
#define MaxLastZone1             4
#define MaxLastZone2               19
#define AbsoluteLastParticle 102000

#define MaxLastMember        3000000
#define MaxLastInZone             50
#define LastCylinder               3
#define Scale                    0.12
#define DrawingSaveInterval      -10
void DoTheDrawing(void);
void SortIntoZones(void);
int DefineMembers(int);
int LookForMembers(void);
void MemberForces(void);
void Motion(void);
int CreateData(void);
float Coord[3][AbsoluteLastParticle + 1],Velocity[3][AbsoluteLastParticle + 1],
PreviousVelocity[3][AbsoluteLastParticle + 1],
Force[3][AbsoluteLastParticle + 1],
ParticleMass[AbsoluteLastParticle + 1],h_Multiplier_Squared[AbsoluteLastParticle + 1],
BasicMemberStiffness[MaxLastMember + 1],UnstressedLength[MaxLastMember + 1],
Boundary[3],
Basic_h,Basic_h_Squared,
mTension,mCompression,strainJoin,mJoin,eta_Squared,
TensionConstant,CompressionConstant,ColourBasis,
ZoneSize,Maximum_r_over_h_Squared,ZoneSizeSquared,TwelveZoneSizeSquared,
Deltat,DistanceFromCentreSq,
CylinderRadius,CylinderRadiusSq,CylinderHorizontalPosition[LastCylinder + 1],CylinderVerticalPosition[LastCylinder + 1];

int  LastInZone[MaxLastZone0 + 1][MaxLastZone1 + 1][MaxLastZone2 + 1],
InZone[MaxLastZone0 + 1][MaxLastZone1 + 1][MaxLastZone2 + 1][MaxLastInZone + 1],
MemberEnd[2][MaxLastMember + 1],
Zone[3],LastZone[3],dxfGreyMap[6],dxf,
LastParticle,Cycle,LastMember,OldState,BlackOrWhite,ScreenUpdate;
char NumericalValue[101];
char * MyText;
ofstream Madeleine("Fibres.dxf");
int main(int argc, char * argv[])
{
	if(CreateData())return 0;
	if(OldState == 0)
	{
		SortIntoZones();
		if(LookForMembers() == 1)return 1;
	}
	cout<<"Last member = "<<LastMember<<"\n";
	glutInit(&argc,argv);
	zRotation=0.0;
	HorizontalRotation=0.0;
	zRotationIncrement=0.0;HorizontalRotationIncrement=-90.0;
	ScreenUpdate = 0;
	if(BlackOrWhite == 1)SetUpGraphics(1.0,1.0,1.0);
	else SetUpGraphics(0.0,0.0,0.0);
	return 0;
}

static void Draw(void)
{
	Cycle ++;
	if(ScreenUpdate == 0 || dxf ==1)glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	DoTheDrawing();
	ScreenUpdate ++;
	if(ScreenUpdate > 20)ScreenUpdate = 0;
	Motion();
	if(DrawingSaveInterval >= 0&&Cycle == DrawingSaveInterval * int(Cycle / DrawingSaveInterval))SavePicture();
}

void DoTheDrawing(void)
{
	float rv[3];
	
	for(int Particle = 0; Particle <= LastParticle; Particle ++)
	{
		for(int Cylinder = 0;Cylinder <= LastCylinder;Cylinder ++)
		{
			float delta_h = Coord[0][Particle] - CylinderHorizontalPosition[Cylinder];
			float delta_v = Coord[2][Particle] - CylinderVerticalPosition[Cylinder];
			float DistanceFromAxisSquared = delta_h * delta_h + delta_v * delta_v;
			if(DistanceFromAxisSquared < CylinderRadiusSq)
			{
				float ratio = sqrt(CylinderRadiusSq / DistanceFromAxisSquared);
				Coord[0][Particle] = CylinderHorizontalPosition[Cylinder] + ratio * delta_h;
				Coord[2][Particle] =   CylinderVerticalPosition[Cylinder] + ratio * delta_v;
				
				float radial_velocity_over_radius = (delta_h * Velocity[0][Particle] + delta_v * Velocity[2][Particle]) / DistanceFromAxisSquared;
				Velocity[0][Particle] -= radial_velocity_over_radius * delta_h;
				Velocity[2][Particle] -= radial_velocity_over_radius * delta_v;
			}
		}
	}
	
	MemberForces();
	
	/*for(int Particle = 0; Particle <= LastParticle; Particle ++)
	{
		float ThisPointSize = 0.1 * Basic_h * sqrt(h_Multiplier_Squared[Particle]) * Zoom;
		if(ThisPointSize<1.0)ThisPointSize = 1.0;
		glPointSize(ThisPointSize);
		glColor4f(0.0,0.0,0.0,0.1);
		glBegin(GL_POINTS);
		
		for(int xyz = 0; xyz <= 2; xyz ++)rv[xyz] = Coord[xyz][Particle];
		glVertex3fv(rv);
		glEnd();
	}*/
	
	if(ScreenUpdate == 0 || dxf ==1)
	{
		if(DrawGrid == 1)
		{
			glLineWidth(1.0); glColor4f(1.0,0.0,0.0,1.0);
			glBegin(GL_LINES);
			for(int gridDirection = 0;gridDirection <=2; gridDirection++)
			{
				int      next_gridDirection = gridDirection + 1;if(     next_gridDirection > 2)     next_gridDirection -= 3;
				int next_next_gridDirection = gridDirection + 2;if(next_next_gridDirection > 2)next_next_gridDirection -= 3;
				for(int ZoneLine = 0; ZoneLine <= LastZone[next_gridDirection] + 1; ZoneLine ++)
				{
					for(int otherZoneLine = 0; otherZoneLine <= LastZone[next_next_gridDirection] + 1; otherZoneLine ++)
					{
						rv[     next_gridDirection] = float(     ZoneLine) * ZoneSize - Boundary[     next_gridDirection];
						rv[next_next_gridDirection] = float(otherZoneLine) * ZoneSize - Boundary[next_next_gridDirection];
						rv[gridDirection] = - Boundary[gridDirection]; glVertex3fv(rv);
						rv[gridDirection] = + Boundary[gridDirection]; glVertex3fv(rv);
					}
				}
			}
			glEnd();
		}
		if(DrawBoundary == 1)
		{
			glLineWidth(1.0); glColor4f(0.0,0.0,0.0,0.5);
			glBegin(GL_LINE_STRIP);
			rv[0] = - Boundary[0]; rv[1] = - Boundary[1]; rv[2] = Boundary[2]; glVertex3fv(rv);
			rv[0] = + Boundary[0]; glVertex3fv(rv);
			rv[1] = + Boundary[1]; glVertex3fv(rv);
			rv[0] = - Boundary[0]; glVertex3fv(rv);
			rv[1] = - Boundary[1]; glVertex3fv(rv);
			glEnd();
		}
		glutSwapBuffers();
	}
}

void SortIntoZones(void)
{
	for(Zone[0] = 0; Zone[0] <= LastZone[0]; Zone[0] ++)
	{
		for(Zone[1] = 0; Zone[1] <= LastZone[1]; Zone[1] ++)
		{
			for(Zone[2] = 0; Zone[2] <= LastZone[2]; Zone[2] ++)LastInZone[Zone[0]][Zone[1]][Zone[2]] = - 1;
		}
	}
	for(int Particle = 0; Particle <= LastParticle; Particle ++)
	{
		for(int xyz = 0; xyz <= 2; xyz ++)
		{
			float floatZone = (Coord[xyz][Particle] + Boundary[xyz]) / ZoneSize;
			if(floatZone<0.0)
			{
				Zone[xyz] = 0;
				cout<<"Negative floatZone, xyz = "<<xyz<<" floatZone = "<<floatZone<<"\n";
			}
			else 
			{
				Zone[xyz] = int(floatZone);
				if(Zone[xyz]<0)
				{
					Zone[xyz] = 0;
					cout<<"Negative integer zone\n";
				}
				if(Zone[xyz]>LastZone[xyz])
				{
					Zone[xyz] = LastZone[xyz];
					cout<<"Integer zone too large, xyz = "<<xyz<<" floatZone = "<<floatZone<<"\n";
				}
			}
		}
		LastInZone[Zone[0]][Zone[1]][Zone[2]] ++;
		if(LastInZone[Zone[0]][Zone[1]][Zone[2]] <= MaxLastInZone)
			InZone[Zone[0]][Zone[1]][Zone[2]][LastInZone[Zone[0]][Zone[1]][Zone[2]]] = Particle;
	}
	for(Zone[0] = 0; Zone[0] <= LastZone[0]; Zone[0] ++)
	{
		for(Zone[1] = 0; Zone[1] <= LastZone[1]; Zone[1] ++)
		{
			for(Zone[2] = 0; Zone[2] <= LastZone[2]; Zone[2] ++)
			{
				if(LastInZone[Zone[0]][Zone[1]][Zone[2]]>MaxLastInZone)
				{
					cout<<"Number in Zone "<<Zone[0]<<" "<<Zone[1]<<" "<<Zone[2]<<" is "<<LastInZone[Zone[0]][Zone[1]][Zone[2]]<<"\n";
					LastInZone[Zone[0]][Zone[1]][Zone[2]] = MaxLastInZone;
				}
			}
		}
	}
}

int LookForMembers(void)
{
	for(Zone[0] = 0; Zone[0] <= LastZone[0]; Zone[0] ++)
	{
		for(Zone[1] = 0; Zone[1] <= LastZone[1]; Zone[1] ++)
		{
			for(Zone[2] = 0; Zone[2] <= LastZone[2]; Zone[2] ++)
			{
				if(LastInZone[Zone[0]][Zone[1]][Zone[2]] >= 0)
				{
					for(int ParticleInZone = 0; ParticleInZone <= LastInZone[Zone[0]][Zone[1]][Zone[2]]; ParticleInZone ++)
					{
						int Particle_A = InZone[Zone[0]][Zone[1]][Zone[2]][ParticleInZone];
						if(DefineMembers(Particle_A) == 1)return 1;
					}
				}
			}
		}
	}
	return 0;
}

int DefineMembers(int Particle_A)
{
	float delta_r[3];
	int WhichZone[3],UpperZone[3],LowerZone[3];
	
	for(int xyz = 0; xyz <= 2; xyz ++)
	{
		UpperZone[xyz] = Zone[xyz] + 1;if(UpperZone[xyz] > LastZone[xyz])UpperZone[xyz] = LastZone[xyz];
		LowerZone[xyz] = Zone[xyz] - 1;if(LowerZone[xyz] < 0)LowerZone[xyz] = 0;
	}
	
	for(WhichZone[0] = LowerZone[0]; WhichZone[0] <= UpperZone[0]; WhichZone[0] ++)
	{
		for(WhichZone[1] = LowerZone[1]; WhichZone[1] <= UpperZone[1]; WhichZone[1] ++)
		{
			for(WhichZone[2] = LowerZone[2]; WhichZone[2] <= UpperZone[2]; WhichZone[2] ++)
			{
				if(LastInZone[WhichZone[0]][WhichZone[1]][WhichZone[2]] >= 0)
				{
					for(int otherParticleInZone = 0; otherParticleInZone <= LastInZone[WhichZone[0]][WhichZone[1]][WhichZone[2]]; otherParticleInZone ++)
					{
						int Particle_B = InZone[WhichZone[0]][WhichZone[1]][WhichZone[2]][otherParticleInZone];
						if(Particle_B<Particle_A)
						{
							for(int xyz = 0; xyz <= 2; xyz ++)delta_r[xyz] = Coord[xyz][Particle_B] - Coord[xyz][Particle_A];
							
							float r_AB_Squared = delta_r[0] * delta_r[0] + delta_r[1] * delta_r[1] + delta_r[2] * delta_r[2];
							if(r_AB_Squared>TwelveZoneSizeSquared)
							{
								cout<<"\nSeparation bigger than it should be.";
								cout<<"\nDeltax / L = "<<delta_r[0] / ZoneSize;
								cout<<"\nDeltay / L = "<<delta_r[1] / ZoneSize;
								cout<<"\nDeltaz / L = "<<delta_r[2] / ZoneSize<<"\n";
							}
							
							float h_Multiplier_AB_Squared = (h_Multiplier_Squared[Particle_A] + h_Multiplier_Squared[Particle_B]) / 2.0;
							
							float r_AB_Squared_over_h_AB_Squared = r_AB_Squared / (Basic_h_Squared * h_Multiplier_AB_Squared);
							
							if(r_AB_Squared_over_h_AB_Squared < Maximum_r_over_h_Squared && r_AB_Squared < ZoneSizeSquared)
							{
								LastMember++;if(LastMember > MaxLastMember){cout<<"Too many members\n";return 1;}
								MemberEnd[0][LastMember] = Particle_A;
								MemberEnd[1][LastMember] = Particle_B;
								UnstressedLength[LastMember] = sqrt(r_AB_Squared);
								BasicMemberStiffness[LastMember] = r_AB_Squared_over_h_AB_Squared * exp(- r_AB_Squared_over_h_AB_Squared);
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

void MemberForces(void)
{
	float rv[3],delta_r[3],delta_Velocity[3],Velocity_A[3],Velocity_B[3];
	float Proportion = 1.0 / 3.0;
	float OneMinusProportion = 1.0 - Proportion;
	
	for(int Particle = 0; Particle <= LastParticle; Particle ++)
	{
		for(int xyz = 0; xyz <= 2; xyz ++)Force[xyz][Particle] = 0.0;
	}
	for(int Member = 0; Member <= LastMember; Member ++)
	{
		int Particle_A = MemberEnd[0][Member];
		int Particle_B = MemberEnd[1][Member];
		float delat_r_dot_delta_Velocity = 0.0;
		for(int xyz = 0; xyz <= 2; xyz ++)
		{
			delta_r[xyz] = Coord[xyz][Particle_B] - Coord[xyz][Particle_A];
			Velocity_A[xyz] = (Velocity[xyz][Particle_A] - Proportion * PreviousVelocity[xyz][Particle_A]) / OneMinusProportion;
			Velocity_B[xyz] = (Velocity[xyz][Particle_B] - Proportion * PreviousVelocity[xyz][Particle_B]) / OneMinusProportion;
			delta_Velocity[xyz] = Velocity_B[xyz] - Velocity_A[xyz];
			delat_r_dot_delta_Velocity += delta_r[xyz] * delta_Velocity[xyz];
		}
		float r_AB_Squared = delta_r[0] * delta_r[0] + delta_r[1] * delta_r[1] + delta_r[2] * delta_r[2];
		
		float r_AB = sqrt(r_AB_Squared);
		
		float VelocityComponent = 0.0;
		if(r_AB > 0.1 * ZoneSize)VelocityComponent = delat_r_dot_delta_Velocity / r_AB;
		
		float elongation = r_AB - UnstressedLength[Member];
		
		float strain = 0.0;
		
		if(UnstressedLength[Member] > 0.001 * Basic_h)strain = elongation / UnstressedLength[Member];
		
		float temp0 = strain - strainJoin;
		float temp1 = TensionConstant * strain - strainJoin;
		float strainFunction = mTension * strain / 2.0
			+ (mJoin / 2.0) * (sqrt(temp0 * temp0 + eta_Squared) - sqrt(temp1 * temp1 + eta_Squared));
		
		temp0 = strain - strainJoin;
		temp1 = CompressionConstant * strain - strainJoin;
		strainFunction += mCompression * strain / 2.0
			- (mJoin / 2.0) * (sqrt(temp0 * temp0 + eta_Squared) - sqrt(temp1 * temp1 + eta_Squared));
		
		float TensionCoefficient = 0.05 * (strainFunction * UnstressedLength[Member] + 0.0 * VelocityComponent) * BasicMemberStiffness[Member];
		//float TensionCoefficient = 0.02 * (elongation + 0.1 * VelocityComponent) * BasicMemberStiffness[Member];
		
		for(int xyz = 0; xyz <= 2; xyz ++)
		{
			Force[xyz][Particle_A] += TensionCoefficient * delta_r[xyz];
			Force[xyz][Particle_B] -= TensionCoefficient * delta_r[xyz];
		}
		
		if(ScreenUpdate == 0 || dxf ==1)
		{
			float ColourDensity;
			if(strainFunction > 0.0)
			{
				float ColourFactor = ColourBasis / (strainJoin * mTension);
				ColourDensity = tanh(ColourFactor * strainFunction);
				float ColourControl = tanh(ColourFactor * strainFunction);
				glColor4f(0.0,0.0,ColourControl,ColourDensity);
			}
			else
			{
				float ColourFactor = ColourBasis / (strainJoin * mCompression);
				ColourDensity = tanh(- ColourFactor * strainFunction);
				float ColourControl = tanh(- ColourFactor * strainFunction);
				glColor4f(ColourControl,0.0,0.0,ColourDensity);
			}
			
			glBegin(GL_LINE_STRIP);
			for(int xyz = 0; xyz <= 2; xyz ++)rv[xyz] = Coord[xyz][Particle_A]; glVertex3fv(rv);
			for(int xyz = 0; xyz <= 2; xyz ++)rv[xyz] = Coord[xyz][Particle_B]; glVertex3fv(rv);
			glEnd();
			
			if(dxf == 1)
			{
				int Layer = 0;
				if(ColourDensity > 0.1)
				{
					int Grey = int(6.0 * ColourDensity);
					if(Grey < 0)Grey = 0;
					if(Grey > 5)Grey = 5;
					Madeleine<<"0\nLINE\n8\n"<<Layer<<"\n";
					Madeleine<<"10\n"<<   Coord[0][Particle_A]<<"\n";
					Madeleine<<"20\n"<<   Coord[2][Particle_A]<<"\n";
					Madeleine<<"30\n"<< - Coord[1][Particle_A]<<"\n";
					Madeleine<<"11\n"<<   Coord[0][Particle_B]<<"\n";
					Madeleine<<"21\n"<<   Coord[2][Particle_B]<<"\n";
					Madeleine<<"31\n"<< - Coord[1][Particle_B]<<"\n";
					Madeleine<<"62\n"<<dxfGreyMap[Grey]<<"\n";
				}
			}
		}
	}
}

void Motion(void)
{
	for(int Particle = 0; Particle <= LastParticle; Particle ++)
	{
		for(int xyz = 0; xyz <= 2; xyz ++)PreviousVelocity[xyz][Particle] = Velocity[xyz][Particle];
		
		for(int xyz = 0; xyz <= 2; xyz ++)Velocity[xyz][Particle] = 0.99 * Velocity[xyz][Particle] + Force[xyz][Particle] * Deltat / ParticleMass[Particle];
		for(int xyz = 0; xyz <= 2; xyz ++)Coord[xyz][Particle] += Velocity[xyz][Particle] * Deltat;
	}
}

void SaveState(void)
{
	ofstream Julia("State.txt");
	cout<<"Saving files\n";
	Julia<<Cycle<<"\n";
	for(int Particle = 0; Particle <= LastParticle; Particle ++)
	{
		for(int xyz = 0; xyz <= 2; xyz ++)Julia<<h_Multiplier_Squared[Particle]<<" "<<Coord[xyz][Particle]<<" "<<Velocity[xyz][Particle]<<" "<<PreviousVelocity[xyz][Particle]<<"\n";
	}
	Julia<<LastMember<<"\n";
	for(int Member = 0; Member <= LastMember; Member ++)
	{
		Julia<<MemberEnd[0][Member]<<" "<<MemberEnd[1][Member]<<" "<<UnstressedLength[Member]<<" "<<BasicMemberStiffness[Member]<<"\n";
	}
	Julia.close();
	cout<<"Text file written\n";
	dxfGreyMap[0] = 254;
	dxfGreyMap[1] = 253;
	dxfGreyMap[2] = 9;
	dxfGreyMap[3] = 8;
	dxfGreyMap[4] = 250;
	dxfGreyMap[5] = 19;
	Madeleine<<"0\nSECTION\n2\nENTITIES\n";
	dxf = 1;MemberForces();
	Madeleine<<"0\nENDSEC\n0\nEOF\n";
	Madeleine.close();
	cout<<"DXF file written\n";
}

int CreateData(void)
{
	BlackOrWhite = 0;
	dxf = 0;
	
	ZoneSize = 5.0;
	
	float Maximum_r_over_h = 1.0 / 0.35;
	Basic_h = ZoneSize / Maximum_r_over_h;
	Maximum_r_over_h_Squared = Maximum_r_over_h * Maximum_r_over_h;
	cout << "Maximum_r_over_h x Basic_h / ZoneSize = " << Maximum_r_over_h * Basic_h / ZoneSize << "\n";
	float h_Multiplier_Maximum = 1.0;
	float h_Multiplier_Minimum = 0.8;
	Basic_h_Squared = Basic_h * Basic_h;
	float AverageNumberParticlesPerZone = (MaxLastInZone + 1) / 3;
	
	mTension = 1.0;
	mCompression = mTension;
	mJoin = 10.0;
	strainJoin = 0.01;
	float eta = 0.1;
	eta_Squared = eta * eta;
	ColourBasis = 1.0;
	
	TensionConstant     = 1.0 - mTension     / mJoin;
	CompressionConstant = 1.0 - mCompression / mJoin;
	
	cout<<"Average number of particles in zone = "<<AverageNumberParticlesPerZone<<"\n";
	cout<<"Maximum number of particles in zone = "<<MaxLastInZone + 1<<"\n";
	ZoneSizeSquared = ZoneSize * ZoneSize;
	TwelveZoneSizeSquared = 12.0 * ZoneSizeSquared;
	float Basic_ParticleMass = 1.0;
	
	Deltat = 1.0;
	LastZone[0] = int(300 * Scale) - 1;
	LastZone[1] = int(25 * Scale) - 1;
	LastZone[2] = int(100 * Scale) - 1;
	for(int xyz = 0; xyz <= 2; xyz ++)
	{
		if(LastZone[xyz]<0)LastZone[xyz] = 0;
	}
	cout << "\nCorrect array sizes:\n\n";
	cout << "#define MaxLastZone0             " << LastZone[0] << "\n";
	cout << "#define MaxLastZone1             " << LastZone[1] << "\n";
	cout << "#define MaxLastZone2               " << LastZone[2] << "\n";
	if(LastZone[0]>MaxLastZone0||LastZone[1]>MaxLastZone1||LastZone[2]>MaxLastZone2){cout<<"Too many zones in at least one direction.\n"; return 1; }
	
	for(int xyz = 0; xyz <= 2; xyz ++)Boundary[xyz] = 0.5 * float(LastZone[xyz] + 1) * ZoneSize;
	AspectRatio = 2.0;
	BorderHalfWidth = Boundary[0];
	
	CylinderRadius = 0.2 * Boundary[0];
	CylinderRadiusSq = CylinderRadius * CylinderRadius;
	
	CylinderHorizontalPosition[0] = - 0.9 * Boundary[0];
	CylinderHorizontalPosition[1] = - CylinderHorizontalPosition[0];
	CylinderHorizontalPosition[2] = 0.3 * Boundary[0];
	CylinderHorizontalPosition[3] = -CylinderHorizontalPosition[2];
	
	CylinderVerticalPosition[0] = - (Boundary[2] + 0.9 * CylinderRadius);
	CylinderVerticalPosition[1] = CylinderVerticalPosition[0];
	CylinderVerticalPosition[2] = - CylinderVerticalPosition[0];
	CylinderVerticalPosition[3] = CylinderVerticalPosition[2];
	
	LastParticle = int(AverageNumberParticlesPerZone * float((LastZone[0] + 1) * (LastZone[1] + 1) * (LastZone[2] + 1)));
	cout << "#define AbsoluteLastParticle " << LastParticle << "\n\n";
	if(LastParticle>AbsoluteLastParticle)
	{
		cout<<"which is too many particles\n"; return 1;
	}
	cout<<"Type 0 to generate new data, or any other number to read from file.\n";cin>>OldState;
	//OldState = 0;
	if(OldState != 0)
	{
		cout<<"Reading from file\n";
		ifstream Maud("State.txt");
		Maud>>Cycle; cout<<"Cycle = "<<Cycle<<"\n";
		for(int Particle = 0; Particle <= LastParticle; Particle ++)
		{
			for(int xyz = 0; xyz <= 2; xyz ++)Maud>>h_Multiplier_Squared[Particle]>>Coord[xyz][Particle]>>Velocity[xyz][Particle]>>PreviousVelocity[xyz][Particle];
		}
		Maud>>LastMember;
		for(int Member = 0; Member <= LastMember; Member ++)
		{
			Maud>>MemberEnd[0][Member]>>MemberEnd[1][Member]>>UnstressedLength[Member]>>BasicMemberStiffness[Member];
		}
		Maud.close();
		cout<<"Finished reading from file\n";
	}
	else
	{
		Cycle = - 1;
		float PIby2 = 2.0 * atan(1.0);
		float Max_Actual_Multiplier = 0.0;
		float Min_Actual_Multiplier;
		double Sum_Actual_Multiplier = 0.0;
		double Sum_Square_Actual_Multiplier = 0.0;
		for(int Particle = 0; Particle <= LastParticle; Particle ++)
		{
			float h_Multiplier = (h_Multiplier_Maximum + h_Multiplier_Minimum) / 2.0
			+ (asin(0.99 * (1.0 - 2.0 * float(rand()) / float(RAND_MAX))) / PIby2) * (h_Multiplier_Maximum - h_Multiplier_Minimum) / 2.0;
			
			h_Multiplier_Squared[Particle] = h_Multiplier * h_Multiplier;
			
			if(Min_Actual_Multiplier > h_Multiplier || Particle == 0)Min_Actual_Multiplier = h_Multiplier;
			if(Max_Actual_Multiplier < h_Multiplier)Max_Actual_Multiplier = h_Multiplier;
			Sum_Actual_Multiplier += double(h_Multiplier);
			Sum_Square_Actual_Multiplier += double(h_Multiplier_Squared[Particle]);
			
			for(int xyz = 0; xyz <= 2; xyz ++)
				Coord[xyz][Particle] = Boundary[xyz] * ((2.0 * rand()) / (1.0 * RAND_MAX) - 1.0);
			for(int xyz = 0; xyz <= 2; xyz ++)
			{
				Velocity[xyz][Particle] = 0.0;
				PreviousVelocity[xyz][Particle] = 0.0;
			}
		}
		cout<<"Minimum actual h multiplier = "<<Min_Actual_Multiplier<<"\n";
		cout<<"Maximum actual h multiplier = "<<Max_Actual_Multiplier<<"\n";
		float Mean_h_Multiplier = Sum_Actual_Multiplier / float(LastParticle + 1);
		cout<<"Mean actual h multiplier = "<<Mean_h_Multiplier<<"\n";
		cout<<"Standard deviation of h multiplier = "<<sqrt(Sum_Square_Actual_Multiplier  / float(LastParticle + 1) - Mean_h_Multiplier * Mean_h_Multiplier)<<"\n";
	}
	if(Cycle == - 1)FileNumber = - 1; else FileNumber = int(Cycle / DrawingSaveInterval);
	
	float TotalMass = 0.0;
	for(int Particle = 0; Particle <= LastParticle; Particle ++)
	{
		ParticleMass[Particle] = Basic_ParticleMass * h_Multiplier_Squared[Particle] * sqrt(h_Multiplier_Squared[Particle]);
		TotalMass += ParticleMass[Particle];
	}	
	return 0;
}
