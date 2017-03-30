//Written by Chris J K Williams, University of Bath, UK

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include "Graphics.h"

#define MaxLastZone0             59
#define MaxLastZone1             0
#define MaxLastZone2             9

#define AbsoluteLastParticle     (11 * 1000)

#define MaxLastFibre             (250 * 1000)
#define MaxLastInZone            50
#define LastCylinder             3
#define Scale                    0.2
#define DrawingSaveInterval      2
#define ScreenUpdateInterval     50
#define LAST_THREAD              1

pthread_t CalcThread[LAST_THREAD + 1];
pthread_mutex_t OverlapMutex[LAST_THREAD + 1];
int FirstZone0InThread[LAST_THREAD + 1];
int LastZone0InThread[LAST_THREAD + 1];
int FirstParticleInThread[LAST_THREAD + 1];
int LastParticleInThread[LAST_THREAD + 1];

void DoTheDrawing(void);
void SortIntoZones(void);
int DefineFibres(int);
int LookForFibres(void);
void FibreForces(void);
void MoveParticles(void);
void FibreForcesOnThread(int, int);
void Motion(int, int);
int CreateData(void);

void ThreadsInit(void);

void* ForceThreadFunction(void *threadid);
void ForceThreadsRun(void);

void* MotionThreadFunction(void *threadid);
void MotionThreadsRun(void);

float Coord[3][AbsoluteLastParticle + 1],Velocity[3][AbsoluteLastParticle + 1],
PreviousVelocity[3][AbsoluteLastParticle + 1],
Force[3][AbsoluteLastParticle + 1],
ParticleMass[AbsoluteLastParticle + 1],h_Multiplier_Squared[AbsoluteLastParticle + 1],
BasicFibreStiffness[MaxLastFibre + 1],UnstressedLength[MaxLastFibre + 1],
CylinderHorizontalPosition[LastCylinder + 1],CylinderVerticalPosition[LastCylinder + 1],
Boundary[3],
Basic_h,Basic_h_Squared,
ZoneSize,Maximum_r_over_h_Squared,ZoneSizeSquared,TwelveZoneSizeSquared,
Deltat,DistanceFromCentreSq,
CylinderRadius,CylinderRadiusSq,CylinderMotionperCycle,FractureStrain;

float Strain[MaxLastFibre];

int  LastInZone[MaxLastZone0 + 1][MaxLastZone1 + 1][MaxLastZone2 + 1],
InZone[MaxLastZone0 + 1][MaxLastZone1 + 1][MaxLastZone2 + 1][MaxLastInZone + 1],
FibreEnd[2][MaxLastFibre + 1],LastFibreInZone0[MaxLastZone0+1],
FibreBroken[MaxLastFibre + 1],
Zone[3],LastZone[3],dxfGreyMap[6],dxf,
LastParticle,
LastFibre,OldState,BlackOrWhite,ScreenUpdate,Cycle,CyclesSinceDrawingSave,TwoD;
char NumericalValue[101];
char * MyText;
//ofstream Madeleine("Fibres.dxf");
int main(int argc, char * argv[])
{
    if(CreateData())return 0;
	if(OldState == 0)
	{
		SortIntoZones();
		if(LookForFibres() == 1)return 1;
	}
	cout<<"Last fibre = "<<LastFibre<<"\n";
	for(int Particle = 0; Particle <= LastParticle; Particle ++)
	{
		for(int xyz = 0; xyz <= 2; xyz ++)Force[xyz][Particle] = 0.0;
	}
	ThreadsInit();
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
	CylinderVerticalPosition[0] = - (Boundary[2] + 1.0 * CylinderRadius) + float(Cycle) * CylinderMotionperCycle;
	CylinderVerticalPosition[1] =   CylinderVerticalPosition[0];
	CylinderVerticalPosition[2] = - CylinderVerticalPosition[0];
	CylinderVerticalPosition[3] = - CylinderVerticalPosition[0];
	
	FibreForces();
	MoveParticles();	
	
	Cycle ++;
	if(ScreenUpdate == 0 || dxf ==1)
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		DoTheDrawing();
		if(DrawingSaveInterval >= 0 && CyclesSinceDrawingSave == 0)SavePicture();
		CyclesSinceDrawingSave ++;
		if(CyclesSinceDrawingSave >= DrawingSaveInterval)CyclesSinceDrawingSave = 0;
		glutSwapBuffers();
	}	
	
	ScreenUpdate ++;
	if(ScreenUpdate >= ScreenUpdateInterval) ScreenUpdate = 0;
}

void DoTheDrawing(void)
{
	float rv[3];
	
	if(ShowParticles == 1)
	{
		for(int Particle = 0; Particle <= LastParticle; Particle ++)
		{
			float ThisPointSize = 1.0 * Basic_h * sqrt(h_Multiplier_Squared[Particle]) * Zoom;
			if(ThisPointSize<1.0)ThisPointSize = 1.0;
			glPointSize(ThisPointSize);
			if(BlackOrWhite == 0)glColor4f(1.0,1.0,1.0,0.2);else glColor4f(0.0,0.0,0.0,0.2);
			glBegin(GL_POINTS);
			for(int xyz = 0; xyz <= 2; xyz ++)rv[xyz] = Coord[xyz][Particle];
			glVertex3fv(rv);
			glEnd();
		}
	}
	
	for(int Fibre = 0; Fibre <= LastFibre; Fibre ++)
	{
		if(FibreBroken[Fibre] == 0)
		{
			int Particle_A = FibreEnd[0][Fibre];
			int Particle_B = FibreEnd[1][Fibre];
			
			float StrainFactor = 10.0 * Strain[Fibre] / FractureStrain;
			float ColourFactor = (1.0 - tanh(StrainFactor)) / 2.0;
			
			if(ShowFibres == 1)
			{
				float ColourArgument = - StrainFactor * StrainFactor;
				
				glColor4f(ColourFactor,0.0,1.0 - ColourFactor,0.1 * (1.0 - exp(ColourArgument)));
				
				glBegin(GL_LINE_STRIP);
				for(int xyz = 0; xyz <= 2; xyz ++)rv[xyz] = Coord[xyz][Particle_A]; glVertex3fv(rv);
				for(int xyz = 0; xyz <= 2; xyz ++)rv[xyz] = Coord[xyz][Particle_B]; glVertex3fv(rv);
				glEnd();
			}
			/*if(dxf == 1)
			 {
			 int Layer = 0;
			 if(ColourFactor > 0.1)
			 {
			 int Grey = int(6.0 * ColourFactor);
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
			 }*/
		}
	}
	
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

int LookForFibres(void)
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
						if(DefineFibres(Particle_A) == 1)return 1;
					}
				}
			}
		}
		// make it easier to find which fibres are in each zone.
		LastFibreInZone0[Zone[0]] = LastFibre;
	}
	return 0;
}

int DefineFibres(int Particle_A)
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
								LastFibre++;if(LastFibre > MaxLastFibre){cout<<"Too many fibres, "<<LastFibre<<" so far\n";return 1;}
								FibreEnd[0][LastFibre] = Particle_A;
								FibreEnd[1][LastFibre] = Particle_B;
								UnstressedLength[LastFibre] = sqrt(r_AB_Squared);
								BasicFibreStiffness[LastFibre] = r_AB_Squared_over_h_AB_Squared * exp(- r_AB_Squared_over_h_AB_Squared);
								FibreBroken[LastFibre] = 0;
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

void FibreForces(void)
{
#if LAST_THREAD > 0
	ForceThreadsRun();
#else
	FibreForcesOnThread(0, LastFibre);
#endif
}

void MoveParticles(void)
{
#if LAST_THREAD > 0
	MotionThreadsRun();
#else
	Motion(0, LastParticle);
#endif
}

void FibreForcesOnThread(int FirstFibreInThread, int LastFibreInThread)
{
   	float delta_r[3],delta_Velocity[3],Velocity_A[3],Velocity_B[3];
	float Proportion = 1.0 / 3.0;
	float OneMinusProportion = 1.0 - Proportion;
	for(int Fibre = FirstFibreInThread; Fibre <= LastFibreInThread; Fibre ++)
	{
		if(FibreBroken[Fibre] == 0)
		{
			int Particle_A = FibreEnd[0][Fibre];
			int Particle_B = FibreEnd[1][Fibre];
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
			if(r_AB > 0.01 * ZoneSize)VelocityComponent = delat_r_dot_delta_Velocity / r_AB;
			
			float elongation = r_AB - UnstressedLength[Fibre];
			
			Strain[Fibre] = 0.0;
			
			if(UnstressedLength[Fibre] > 0.01 * Basic_h)Strain[Fibre] = elongation / UnstressedLength[Fibre];
			
			if(Strain[Fibre] > FractureStrain)FibreBroken[Fibre] = 1;
			
			float TensionCoefficient = 0.05 * (Strain[Fibre] * UnstressedLength[Fibre] + 0.1 * VelocityComponent) * BasicFibreStiffness[Fibre];
			
			for(int xyz = 0; xyz <= 2; xyz ++)
			{
				Force[xyz][Particle_A] += TensionCoefficient * delta_r[xyz];
				Force[xyz][Particle_B] -= TensionCoefficient * delta_r[xyz];
			}
		}
	}
}

void Motion(int FirstParticleInThread,int LastParticleInThread)
{
	for(int Particle = FirstParticleInThread; Particle <= LastParticleInThread; Particle ++)
	{
		for(int xyz = 0; xyz <= 2; xyz ++)PreviousVelocity[xyz][Particle] = Velocity[xyz][Particle];
		
		for(int xyz = 0; xyz <= 2; xyz ++)Velocity[xyz][Particle] = 0.99 * Velocity[xyz][Particle] + Force[xyz][Particle] * Deltat / ParticleMass[Particle];
		for(int xyz = 0; xyz <= 2; xyz ++)Coord[xyz][Particle] += Velocity[xyz][Particle] * Deltat;
		if(TwoD == 1)Coord[1][Particle] = 0.0;
		
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
		for(int xyz = 0; xyz <= 2; xyz ++)Force[xyz][Particle] = 0.0;
	}
}

int CreateData(void)
{
	BlackOrWhite = 1;
	dxf = 0;
	CyclesSinceDrawingSave = 0;
	CylinderMotionperCycle = 0.005;
	FractureStrain = 0.01;
	TwoD = 1;
	
	ZoneSize = 5.0;
	
	float Maximum_r_over_h = 1.0 / 0.35;
	Basic_h = ZoneSize / Maximum_r_over_h;
	Maximum_r_over_h_Squared = Maximum_r_over_h * Maximum_r_over_h;
	cout << "Maximum_r_over_h x Basic_h / ZoneSize = " << Maximum_r_over_h * Basic_h / ZoneSize << "\n";
	float h_Multiplier_Maximum = 1.0;
	float h_Multiplier_Minimum = 0.8;
	Basic_h_Squared = Basic_h * Basic_h;
	float AverageNumberParticlesPerZone = (MaxLastInZone + 1) / 3;
	
	cout<<"Average number of particles in zone = "<<AverageNumberParticlesPerZone<<"\n";
	cout<<"Maximum number of particles in zone = "<<MaxLastInZone + 1<<"\n";
	ZoneSizeSquared = ZoneSize * ZoneSize;
	TwelveZoneSizeSquared = 12.0 * ZoneSizeSquared;
	float Basic_ParticleMass = 1.0;
	
	if(TwoD == 1)Deltat = 1.0; else Deltat = 2.0;
	LastZone[0] = int(300 * Scale) - 1;
	if(TwoD == 1)LastZone[1] = 0;else LastZone[1] = int(25 * Scale) - 1;
	LastZone[2] = int(50 * Scale) - 1;
	for(int xyz = 0; xyz <= 2; xyz ++)
	{
		if(LastZone[xyz]<0)LastZone[xyz] = 0;
	}
	cout << "\nCorrect array sizes:\n\n";
	cout << "#define MaxLastZone0             " << LastZone[0] << "\n";
	cout << "#define MaxLastZone1             " << LastZone[1] << "\n";
	cout << "#define MaxLastZone2             " << LastZone[2] << "\n\n";
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
	
	LastParticle = int(AverageNumberParticlesPerZone * float((LastZone[0] + 1) * (LastZone[1] + 1) * (LastZone[2] + 1)));
	cout << "#define AbsoluteLastParticle     " << LastParticle << "\n\n";
	if(LastParticle>AbsoluteLastParticle)
	{
		cout<<"which is too many particles\n"; return 1;
	}
	
	OldState = 0;
	//cout<<"Type 0 to generate new data, or any other number to read from file.\n";cin>>OldState;
	
	if(OldState != 0)
	{
		cout<<"Reading from file\n";
		ifstream Maud("State.txt");
		Maud>>Cycle>>CyclesSinceDrawingSave>>FileNumber;
		cout<<"Cycle = "<<Cycle<<"\n";
		for(int Particle = 0; Particle <= LastParticle; Particle ++)
		{
			Maud>>h_Multiplier_Squared[Particle];
			for(int xyz = 0; xyz <= 2; xyz ++)Maud>>Coord[xyz][Particle]>>Velocity[xyz][Particle]>>PreviousVelocity[xyz][Particle];
		}
		Maud>>LastFibre;
		for(int Fibre = 0; Fibre <= LastFibre; Fibre ++)
			Maud>>FibreEnd[0][Fibre]>>FibreEnd[1][Fibre]>>UnstressedLength[Fibre]>>BasicFibreStiffness[Fibre]>>FibreBroken[Fibre];
		
		for(Zone[0] = 0; Zone[0] <= LastZone[0]; Zone[0] ++)Maud>>LastFibreInZone0[Zone[0]];
		
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
			
			if(Particle == 0)Min_Actual_Multiplier = h_Multiplier;
			if(Min_Actual_Multiplier > h_Multiplier || Particle == 0)Min_Actual_Multiplier = h_Multiplier;
			if(Max_Actual_Multiplier < h_Multiplier)Max_Actual_Multiplier = h_Multiplier;
			Sum_Actual_Multiplier += double(h_Multiplier);
			Sum_Square_Actual_Multiplier += double(h_Multiplier_Squared[Particle]);
			
			for(int xyz = 0; xyz <= 2; xyz ++)
				Coord[xyz][Particle] = Boundary[xyz] * ((2.0 * rand()) / (1.0 * RAND_MAX) - 1.0);
			if(TwoD == 1)Coord[1][Particle] = 0.0;
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
	if(Cycle == - 1)FileNumber = - 1;
	
	for(int Particle = 0; Particle <= LastParticle; Particle ++)
		ParticleMass[Particle] = Basic_ParticleMass * h_Multiplier_Squared[Particle] * sqrt(h_Multiplier_Squared[Particle]);
	
	return 0;
}

void ThreadsInit(void)
{
	//allocate a range of Zones to each partition and initialise the mutexes.
	for(int thread = 0; thread <= LAST_THREAD; thread++)
	{
		pthread_mutex_init(&OverlapMutex[thread], NULL);
		
		LastZone0InThread[thread] = ((thread + 1) * LastZone[0]) / (LAST_THREAD + 1);
		if(thread !=0)FirstZone0InThread[thread] = LastZone0InThread[thread - 1] + 1;
		else FirstZone0InThread[thread] = 0;
		
		LastParticleInThread[thread] = ((thread + 1) * LastParticle) / (LAST_THREAD + 1);
		if(thread != 0)FirstParticleInThread[thread] = LastParticleInThread[thread - 1] + 1;
		else FirstParticleInThread[thread] = 0;
	}
	//make sure last thread really runs to the end of the domain.
	LastZone0InThread[LAST_THREAD] = LastZone[0];
	LastParticleInThread[LAST_THREAD] = LastParticle;
	
	for(int thread = 0; thread <= LAST_THREAD; thread++)
		cout<<"Thread = "<<thread<<" First zone 0 = "<<FirstZone0InThread[thread]<<" Last zone 0 = "<<LastZone0InThread[thread]<<"\n";
	for(int thread = 0; thread <= LAST_THREAD; thread++)
		cout<<"Thread = "<<thread<<" First particle = "<<FirstParticleInThread[thread]<<" Last Particle = "<<LastParticleInThread[thread]<<"\n";
}

void ForceThreadsRun(void)
{
	int rc;
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//start all the threads
	for(int thread = 0; thread <= LAST_THREAD; thread++)
	{
		rc = pthread_create(&CalcThread[thread], &attr, ForceThreadFunction, (void*)thread);
		if (rc)
		{
            cerr << "ERROR: Couldn't launch a thread,\n";
            exit(-1);        
		}  
	}
	
	//wait for all the threads to finish
	for(int thread = 0; thread <= LAST_THREAD; thread++)
	{
		rc = pthread_join(CalcThread[thread], NULL);
		if (rc)
		{
            cerr << "ERROR: Couldn't join a thread,\n";
            exit(-1);        
		}
	}     
}

void *ForceThreadFunction(void *threadid)
{
    int id = (int)threadid;
	
	if(id > 0) pthread_mutex_lock(&OverlapMutex[id] - 1);
	
    for(int zone = FirstZone0InThread[id];  zone <= LastZone0InThread[id]; zone++)
	{
        if(id > 0           && zone == FirstZone0InThread[id] + 1) pthread_mutex_unlock(&OverlapMutex[id - 1]);
		if(id < LAST_THREAD && zone ==  LastZone0InThread[id] - 1)   pthread_mutex_lock(&OverlapMutex[id + 1]);
        
        int FirstFibreInThisZone, LastFibreInThisZone;
		
        if(zone==0) FirstFibreInThisZone = 0;else FirstFibreInThisZone = LastFibreInZone0[zone - 1] + 1;
		
        LastFibreInThisZone = LastFibreInZone0[zone];
		
		FibreForcesOnThread(FirstFibreInThisZone, LastFibreInThisZone);
    }
	if(id < LAST_THREAD) pthread_mutex_unlock(&OverlapMutex[id + 1]);
	
    pthread_exit(NULL);
	return 0;
}

void MotionThreadsRun(void)
{
	int rc;
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//start all the threads
	for(int thread = 0; thread <= LAST_THREAD; thread++)
	{
		rc = pthread_create(&CalcThread[thread], &attr, MotionThreadFunction, (void*)thread);
		if (rc)
		{
            cerr << "ERROR: Couldn't launch a thread,\n";
            exit(-1);        
		}  
	}
	
	//wait for all the threads to finish
	for(int thread = 0; thread <= LAST_THREAD; thread++)
	{
		rc = pthread_join(CalcThread[thread], NULL);
		if (rc)
		{
            cerr << "ERROR: Couldn't join a thread,\n";
            exit(-1);        
		}
	}     
}

void *MotionThreadFunction(void *threadid)
{
    int id = (int)threadid;
	Motion(FirstParticleInThread[id], LastParticleInThread[id]);
    pthread_exit(NULL);
	return 0;
}
