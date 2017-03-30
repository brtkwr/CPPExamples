//Written by Chris J K Williams, University of Bath, UK
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include <pthread.h>

#include "Graphics.h"

void Calculation(void);
void CalculateForces(int);
void Axial(int,int,int,int);
void Initialise(int);
void CalculateMotion(int);
void FinaliseForces(void);
void CalclateDiagonalLengths(int,int);
void FindNormal(int,int);

void ThreadsInit(void);

void* ThreadFunction(void *threadid);
void ThreadsRun(void);

void* ThreadMotionFunction(void *threadid);
void ThreadsMotionRun(void);

void JoinUp(void);
void MakeLine(double,double,double,double,double,double,double);

void CableColourControl(int);

#define  fine 6
#define  cableSpacing 2

#define SHAPE_CONTROL 1
#if SHAPE_CONTROL == 0
#define  q ( 1 * fine * cableSpacing)
#define  r ( 4 * fine * cableSpacing)
#define  s ( 5 * fine * cableSpacing)
#define  t ( 4 * fine * cableSpacing)
#define  u ( 6 * fine * cableSpacing)
#else
#define  q ( 1 * fine * cableSpacing)
#define  r ( 3 * fine * cableSpacing)
#define  s (10 * fine * cableSpacing)
#define  t ( 8 * fine * cableSpacing)
#define  u ( 7 * fine * cableSpacing)
#endif

#define  n (q + r + s + t + u)

#define FLAT_CONTROL 0
#if FLAT_CONTROL == 0
#define  SupportHeightMultiplier 0.25
#define  SupportCanMoveHorizontally 0
#else
#define  SupportHeightMultiplier 0.0
#define  SupportCanMoveHorizontally 1
#endif

#define  CarryOver 0.99
#define  MovementFactor 0.2

#define LAST_THREAD  1

pthread_t CalcThread[LAST_THREAD + 1];
pthread_mutex_t OverlapMutex[LAST_THREAD + 1];
int Partition_iStart[LAST_THREAD + 1];
int Partition_iStop[LAST_THREAD + 1];

int NodeType[n + 1][n + 1],SquareExists[n + 1][n + 1],LastNodeTypeExisting,NumberOfExistingSquares;

double x[n + 1][n + 1][3],Force[n + 1][n + 1][3],
velocity[n + 1][n + 1][3],unitNormal[n + 1][n + 1][3],
Stiffness[n + 1][n + 1],DiagonalLength[n + 1][n + 1][2],NetCableLength[n + 1][n + 1][2],
SupportPosition[3],SupportVelocity[3],SupportStiffness,SupportForce[3],
PointPosition[2],totalDiagonalLength[2],
EyeCableTension,EAEyeCable,EyeCableSlackLength,EyeCableLength,CableNetTensionCoefficient,
netCableRadius,eyeCableRadius,a,b,sumOfLengthErrorsSq,sumOfDiagonalLengthErrorsSq;

char NumericalValue[101];
string MyText;
//ofstream Madeleine("MinimalSurface.dxf");

int main(int argc, char * argv[])
{
	cout<<"The program will now run forever until 'q' or ESCAPE is pressed on the keyboard\n";
	cout<<"Left mouse button for pan\n";
	cout<<"Right mouse button (or CONTROL left button) for rotate\n";
	cout<<"Move spot to zoom\n";
	
	supportRadius = 100.0;
	SupportHeight = SupportHeightMultiplier * supportRadius;
	a = supportRadius;
	b = a / sqrt(2.0);
	
	netCableRadius = supportRadius / 10000.0;
	eyeCableRadius = supportRadius / 5000.0;
	
	CableNetTensionCoefficient = 0.1;
	EyeCableSlackLength = 0.5 * supportRadius;
	
	EAEyeCable = 0.1 * (double)fine * CableNetTensionCoefficient * supportRadius;
	
	double setupmultiplier = 8.0 * atan(1.0) / double(2 * (n - q));
	for(int i = 0;i <= n;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			double theta = 4.0 * atan(1.0) + double(i - j) * setupmultiplier;
			double phi = double(i + j - (n + q + r)) * setupmultiplier;
			double radius = supportRadius * exp(phi);
			x[i][j][0] = radius * cos(theta) - 0.4 * supportRadius;
			x[i][j][1] = radius * sin(theta);
			x[i][j][2] = SupportHeight * double((n + q + r) - (i + j)) / double(n + q + r);
			
			LastNodeTypeExisting = 8;
			NodeType[i][j] = 0;
			if(i - j > n - q || j - i > n - q
			   ||(i > n - s && j > q + r) || (j > n - s && i > q + r)
			   || (i > q + r + t && j > q + r + t)
			   )NodeType[i][j] = LastNodeTypeExisting + 1;//That is removed
			else
			{
				if(i == n || j == n)NodeType[i][j] = 1;
				if(j >= n - s && i >= q + r)        {if(NodeType[i][j] != 1)NodeType[i][j] = 2;else NodeType[i][j] = 6;}
				if(i >= q + r + t && j >= q + r + t){if(NodeType[i][j] != 2)NodeType[i][j] = 3;else NodeType[i][j] = 7;}
				if(i >= n - s && j >= q + r)        {if(NodeType[i][j] != 3)NodeType[i][j] = 4;else NodeType[i][j] = 8;}
				if((i == n || j == n) && NodeType[i][j] == 4)NodeType[i][j] = 5;
			}
			
			for(int xyz = 0; xyz <= 2;xyz ++)velocity[i][j][xyz] = 0.0;
		}
		Initialise(i);
		SupportPosition[0] = x[0][0][0];
		SupportPosition[1] = x[0][0][1];
		SupportPosition[2] = x[0][0][2];
		
		SupportHorizontalPosition = SupportPosition[0];
		
		for(int xyz = 0; xyz <= 2;xyz ++)SupportVelocity[xyz] = 0.0;
	}
	
	NumberOfExistingSquares = 0;
	for(int i = 0;i <= n;i ++)
	{
		for(int j = 0;j <= n;j ++)
	{
		SquareExists[i][j] = 0;
		if(i <= n - 1 && j <= n - 1)
		{
		if(NodeType[i    ][j    ] <= LastNodeTypeExisting &&
		   NodeType[i + 1][j + 1] <= LastNodeTypeExisting &&
		   NodeType[i    ][j + 1] <= LastNodeTypeExisting)
		{
			if(NodeType[i + 1][j] > LastNodeTypeExisting && i - j != n - q){cout<<"\nUnexpected square condition\n";return 0;}
			
			SquareExists[i][j] = 1;
			NumberOfExistingSquares ++;
		}
	}
	}
}

#if LAST_THREAD > 0
	ThreadsInit();
#endif
	glutInit(&argc,argv);
	MyGraphics(1.0,1.0,1.0);
	return 0;
}

static void Draw(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	glLoadIdentity();
	
	glColor4f(0.0,0.0,0.0,0.2);
	if(totalDiagonalLength[1] != 0.0)
	{
		glRasterPos2d(ButtonX[0],0.65);
		MyText="Ratio of total cable lengths = ";
		int LengthOfString = MyText.length();
		for(int MyCharacter = 0; MyCharacter < LengthOfString; MyCharacter++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyCharacter]);
		sprintf(NumericalValue,"%8.6f",totalDiagonalLength[0] / totalDiagonalLength[1]);
		LengthOfString=strlen(NumericalValue);
		for(int MyCharacter = 0; MyCharacter < LengthOfString; MyCharacter++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	}
	
	{
		glRasterPos2d(ButtonX[0],0.6);
	MyText="Sum of length errors squared = ";
	int LengthOfString = MyText.length();
	for(int MyCharacter = 0; MyCharacter < LengthOfString; MyCharacter++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyCharacter]);
	sprintf(NumericalValue,"%8.6f",sumOfLengthErrorsSq / double(NumberOfExistingSquares));
	LengthOfString=strlen(NumericalValue);
	for(int MyCharacter = 0; MyCharacter < LengthOfString; MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	}
	
	{
	glRasterPos2d(ButtonX[0],0.55);
	MyText="Sum of angle errors squared = ";
	int LengthOfString = MyText.length();
	for(int MyCharacter = 0; MyCharacter < LengthOfString; MyCharacter++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyCharacter]);
	sprintf(NumericalValue,"%8.6f",sumOfDiagonalLengthErrorsSq / double(NumberOfExistingSquares));
	LengthOfString=strlen(NumericalValue);
	for(int MyCharacter = 0; MyCharacter < LengthOfString; MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	}
	
	glPointSize(5.0);
	for(int Button = 0; Button <= lastButton - SupportCanMoveHorizontally; Button ++)
	{
		glBegin(GL_POINTS);
		
		if(ButtonOn[Button] == 1)glColor4f(1.0,0.0,0.0,1.0);else glColor4f(0.0,0.0,0.0,1.0);
		PointPosition[0] = ButtonX[Button];
		PointPosition[1] = ButtonY[Button];
		glVertex2dv(PointPosition);
		glEnd();
		
		glColor4f(0.0,0.0,0.0,0.2);
		if(Button == 0)MyText = "Zoom";
		if(Button == 1)MyText = "Support height";
		if(Button == 2)MyText = "Support horizontal position";
		
		glRasterPos2d(ButtonX[Button] + 0.02,ButtonY[Button] - 0.01);
		
		int StringLength = MyText.length();
		for(int MyChar = 0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
	}
	
	glPopMatrix();
	
	for(int j = 0;j <= n;j += cableSpacing)
	{
		if(j == 0)CableColourControl(0);else CableColourControl(1);
		glBegin(GL_LINE_STRIP);
		for(int i = 0;i <= n;i ++)
		{
			if(NodeType[i][j] <= LastNodeTypeExisting)glVertex3dv(x[i][j]);
		}
		glEnd();
	}
	
	for(int i = 0;i <= n;i += cableSpacing)
	{
		if(i == 0)CableColourControl(0);else CableColourControl(2);
		glBegin(GL_LINE_STRIP);
		for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] <= LastNodeTypeExisting)glVertex3dv(x[i][j]);
		}
		glEnd();
	}
	
	int mini = 1;
	int maxi = 1;
	int minj = 1;
	int maxj = 1;
	for(int i = 1;i <= n - 1;i ++)
	{
		for(int j = 1;j <= n - 1;j ++)
		{
			if(NodeType[i][j] == 0 && (NodeType[i + 1][j] == 2 || NodeType[i][j + 1] == 2 || NodeType[i + 1][j] == 4 || NodeType[i][j + 1] == 4))
			{
				FindNormal(i,j);
				if(unitNormal[i][j][1] < unitNormal[mini][minj][1]){mini = i;minj = j;}
				if(unitNormal[i][j][1] > unitNormal[maxi][maxj][1]){maxi = i;maxj = j;}
			}
		}
	}
	
	glPointSize(5.0);
	glBegin(GL_POINTS);
	for(int i = 1;i <= n - 1;i ++)
	{
		glBegin(GL_LINE_STRIP);
		for(int j = 1;j <= n - 1;j ++)
		{
			if(NodeType[i][j] == 0 && (NodeType[i + 1][j] == 2 || NodeType[i][j + 1] == 2 || NodeType[i + 1][j] == 4 || NodeType[i][j + 1] == 4))
			{
				if((i == mini && j == minj) || (i == maxi && j == maxj))glColor4f(1.0,0.0,0.0,1.0);else glColor4f(0.0,0.0,0.0,0.5);
				glVertex3dv(x[i][j]);
			}
		}
	}
	
	glEnd();
	
	glutSwapBuffers();
	
	for(int CalculationLoop = 0;CalculationLoop <= 0;CalculationLoop++)
	{
		sumOfLengthErrorsSq = 0.0;
		
		for(int i = 0;i <= n - 1;i ++)
		{
			for(int j = 0;j <= n - 1;j ++)
			{
				if(SquareExists[i][j] == 1)
				{
				double Xlength = NetCableLength[i][j + 1][0];
				double Ylength = NetCableLength[i][j    ][1];
					
					if(i - j == n - q)
					{
						Xlength += NetCableLength[j][i][1];
					Ylength += NetCableLength[j][i + 1][0];
					}
					else //This is the standard case:
					{
						Xlength += NetCableLength[i][j][0];
						Ylength += NetCableLength[i + 1][j][1];
					}
					
						sumOfLengthErrorsSq += 4.0 * ((Xlength - Ylength) * (Xlength - Ylength)) / ((Xlength + Ylength) * (Xlength + Ylength));
				}
				}
		}
		
		for(int DiagonalDirection = 0; DiagonalDirection <= 1;DiagonalDirection ++)
		{
			totalDiagonalLength[DiagonalDirection] = 0.0;
			for(int i = 0;i <= n - 1;i ++)
			{
				for(int j = 0;j <= n - 1;j ++)
					if(SquareExists[i][j] == 1)totalDiagonalLength[DiagonalDirection] += DiagonalLength[i][j][DiagonalDirection];
			}
		}
			sumOfDiagonalLengthErrorsSq = 0.0;
			for(int i = 0;i <= n - 1;i ++)
			{
				for(int j = 0;j <= n - 1;j ++)
				{
					if(SquareExists[i][j] == 1)
					{
						sumOfDiagonalLengthErrorsSq += 4.0 *
						((DiagonalLength[i][j][1] - DiagonalLength[i][j][0]) * (DiagonalLength[i][j][1] - DiagonalLength[i][j][0])) /
						((DiagonalLength[i][j][1] + DiagonalLength[i][j][0]) * (DiagonalLength[i][j][1] + DiagonalLength[i][j][0]));
					}
			}
			}
		
		EyeCableSlackLength -= 0.00001 * (totalDiagonalLength[0] - totalDiagonalLength[1]);
		
		EyeCableLength = 0.0;
		for(int i = 0;i <= n - 1;i ++)
		{
			double lengthsq = 0.0;
			if(NodeType[i + 1][0] == 0)
			{
				for(int xyz = 0; xyz <= 2;xyz ++)lengthsq += 
					(x[i + 1][0][xyz] - x[i][0][xyz]) *
					(x[i + 1][0][xyz] - x[i][0][xyz]);
				EyeCableLength += sqrt(lengthsq);
			}
		}
		for(int j = 0;j <= n - 1;j ++)
		{
			double lengthsq = 0.0;
			if(NodeType[0][j + 1] == 0)
			{
				for(int xyz = 0; xyz <= 2;xyz ++)lengthsq += 
					(x[0][j + 1][xyz] - x[0][j][xyz]) *
					(x[0][j + 1][xyz] - x[0][j][xyz]);
				EyeCableLength += sqrt(lengthsq);
			}
		}
		
		EyeCableTension = EAEyeCable * (EyeCableLength - EyeCableSlackLength) / EyeCableSlackLength;
		
		Calculation();
		
		for(int i = n - q; i <= n; i ++)
		{
			for(int xyz = 0;xyz <= 2;xyz ++)
			{
				x[i][i - (n - q)][xyz] += x[i - (n - q)][i][xyz];
				x[i][i - (n - q)][xyz] /= 2.0;
				x[i - (n - q)][i][xyz] = x[i][i - (n - q)][xyz];
			}	
		}
		
		if(SupportCanMoveHorizontally == 1)
		{
			for(int xy = 0;xy <= 1; xy ++)
			{
				SupportVelocity[xy] = 0.5 * CarryOver * SupportVelocity[xy]
				+ MovementFactor * SupportForce[xy] * EyeCableSlackLength / max(SupportStiffness,1000.0 * EAEyeCable);
				SupportPosition[xy] += SupportVelocity[xy];
			}
		}
		else
		{
			SupportPosition[0] = SupportHorizontalPosition;
			SupportPosition[1] = 0.0;
		}
		
		SupportPosition[2] = SupportHeight;
		
		for(int xyz = 0; xyz <= 2;xyz ++)x[0][0][xyz] = SupportPosition[xyz];
	}
}

void Calculation(void)
{
#if LAST_THREAD > 0
	ThreadsRun();
	FinaliseForces();
	ThreadsMotionRun();
#else
	for(int i = 0;i <= n;i ++)CalculateForces(i);
	FinaliseForces();
	for(int i = 0;i <= n;i ++)CalculateMotion(i);
#endif
}

void Initialise(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		Stiffness[i][j] = 0.0;
		for(int xyz = 0; xyz <= 2;xyz ++)Force[i][j][xyz] = 0.0;
	}
}

void CalculateForces(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		Axial(i,i + 1,j,j);
		Axial(i,i,j,j + 1);
		CalclateDiagonalLengths(i,j);
	}
}

void CalculateMotion(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		if(NodeType[i][j] <= LastNodeTypeExisting)
		{
			for(int xyz = 0;xyz <= 2;xyz ++)
			{
				if(NodeType[i][j] == 0 || 
				   ((NodeType[i][j] == 2 || NodeType[i][j] == 4) && xyz == 0) ||
				   ((NodeType[i][j] == 1 || NodeType[i][j] == 3) && xyz == 1))
				{
					velocity[i][j][xyz] = CarryOver * velocity[i][j][xyz] + MovementFactor * Force[i][j][xyz] / Stiffness[i][j];
					x[i][j][xyz] += velocity[i][j][xyz];
				}
			}
			double InitialMotionFactor = 0.0001;
			if(NodeType[i][j] == 1 || NodeType[i][j] == 5 || NodeType[i][j] == 6)
			{
				velocity[i][j][0] = CarryOver * velocity[i][j][0] + InitialMotionFactor * (+ a - x[i][j][0]);
				x[i][j][0] += velocity[i][j][0];
				velocity[i][j][2] = CarryOver * velocity[i][j][2] - InitialMotionFactor * x[i][j][2];x[i][j][2] += velocity[i][j][2];
			}
			if(NodeType[i][j] == 3 || NodeType[i][j] == 7 || NodeType[i][j] == 8)
			{
				velocity[i][j][0] = CarryOver * velocity[i][j][0] + InitialMotionFactor * (- a - x[i][j][0]);
				x[i][j][0] += velocity[i][j][0];
				velocity[i][j][2] = CarryOver * velocity[i][j][2] - InitialMotionFactor * x[i][j][2];x[i][j][2] += velocity[i][j][2];
			}
			if(NodeType[i][j] == 2 || NodeType[i][j] == 6 || NodeType[i][j] == 7)
			{
				velocity[i][j][1] = CarryOver * velocity[i][j][1] + InitialMotionFactor * (+ b - x[i][j][1]);
				x[i][j][1] += velocity[i][j][1];
				velocity[i][j][2] = CarryOver * velocity[i][j][2] - InitialMotionFactor * x[i][j][2];x[i][j][2] += velocity[i][j][2];
			}
			if(NodeType[i][j] == 4 || NodeType[i][j] == 5 || NodeType[i][j] == 8)
			{
				velocity[i][j][1] = CarryOver * velocity[i][j][1] + InitialMotionFactor * (- b - x[i][j][1]);
				x[i][j][1] += velocity[i][j][1];
				velocity[i][j][2] = CarryOver * velocity[i][j][2] - InitialMotionFactor * x[i][j][2];x[i][j][2] += velocity[i][j][2];
			}
		}
	}
	Initialise(i);
}

void Axial(int i,int inext,int j,int jnext)
{
	double deltax[3];
	if(inext <= n && jnext <= n)
	{
		if(NodeType[i][j] <= LastNodeTypeExisting && NodeType[inext][jnext] <= LastNodeTypeExisting)
		{
			for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[inext][jnext][xyz] - x[i][j][xyz];
			
			double tensioncoefficient = CableNetTensionCoefficient;
			
			double thisLengthSq = 0.0;
			for(int xyz = 0;xyz <= 2;xyz ++)thisLengthSq += deltax[xyz] * deltax[xyz];
			double thisLength = sqrt(thisLengthSq);
			
			if((i == 0 && inext == 0) || (j == 0 && jnext ==0))tensioncoefficient = EyeCableTension / thisLength;
			
			if(i == inext)NetCableLength[i][j][0] = thisLength;else NetCableLength[i][j][1] = thisLength;
			
			Stiffness[i][j] += tensioncoefficient;
			Stiffness[inext][jnext] += tensioncoefficient;
			
			for(int xyz = 0;xyz <= 2;xyz ++)
			{
				double thisForce = tensioncoefficient * deltax[xyz];
				Force[i][j][xyz] += thisForce;
				Force[inext][jnext][xyz] -= thisForce;
			}
		}
	}
}

void CalclateDiagonalLengths(int i,int j)
{
	double deltax[3];
	if(i <= n - 1 && j <= n - 1)
	{
		if(SquareExists[i][j] == 1)
			{
	for(int DiagonalDirection = 0; DiagonalDirection <= 1;DiagonalDirection ++)
	{
		if(NodeType[i + 1][j] > LastNodeTypeExisting && DiagonalDirection == 0)
		{
		for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[i + 1][j + DiagonalDirection][xyz] - x[j + DiagonalDirection][i + 1][xyz];
		}
		else //this is the standard case:
			{
				for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[i + 1][j + DiagonalDirection][xyz] - x[i][j + 1 - DiagonalDirection][xyz];
			}
		
		double lengthsq = 0.0;
				for(int xyz = 0;xyz <= 2;xyz ++)lengthsq += deltax[xyz] * deltax[xyz];
				DiagonalLength[i][j][DiagonalDirection] = sqrt(lengthsq);
			}
		}
	}
}

void FinaliseForces(void)
{
	SupportStiffness = Stiffness[0][0];
	for(int xyz = 0;xyz <= 2;xyz ++)
	{
		SupportForce[xyz] = Force[0][0][xyz];
	}
	
	for(int i = n - q; i <= n; i ++)
	{
		for(int xyz = 0;xyz <= 2;xyz ++)
		{
			Force[i][i - (n - q)][xyz] += Force[i - (n - q)][i][xyz];
			Force[i - (n - q)][i][xyz] = Force[i][i - (n - q)][xyz];
		}
	}
}

void ThreadsInit(void)
{
	Partition_iStart[0] = 0;
	
	for(int myThread = 0; myThread <= LAST_THREAD - 1; myThread++)
	{
		Partition_iStop[myThread] = ((myThread + 1) * n) / (LAST_THREAD + 1);
		Partition_iStart[myThread + 1] = Partition_iStop[myThread] + 1;
	}
	
	Partition_iStop[LAST_THREAD] = n;
	
	for(int myThread = 0; myThread <= LAST_THREAD; myThread++)pthread_mutex_init(&OverlapMutex[myThread], NULL);
}

void ThreadsRun(void)
{
	int rc;
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//start all the threads
	for(int myThread = 0; myThread <= LAST_THREAD; myThread++)
	{
		rc = pthread_create(&CalcThread[myThread], &attr, ThreadFunction, (void*)myThread);
		if(rc)
		{
			cerr << "ERROR: Couldn't launch a thread,\n";
			exit(-1);    
		}
	}
	
	//wait for all the threads to finish
	for(int myThread = 0; myThread <= LAST_THREAD; myThread++)
	{
		rc = pthread_join(CalcThread[myThread], NULL);
		if(rc)
		{
			cerr << "ERROR: Couldn't join a thread,\n";
			exit(-1);    
		}
	}
}

void *ThreadFunction(void *threadid)
{
	long id = (long)threadid;
	
	for(int i = Partition_iStart[id]; i <= Partition_iStop[id]; i++)
	{
		if(i == 0                     )pthread_mutex_lock(&OverlapMutex[LAST_THREAD]);
		if(i == Partition_iStart[id]  && id > 0     )pthread_mutex_lock(&OverlapMutex[id-1]);
		if(i == Partition_iStop[id] - 2 && id < LAST_THREAD)pthread_mutex_lock(&OverlapMutex[id+1]);
		if(i == n - 2                   )pthread_mutex_lock(&OverlapMutex[0]);
		
		CalculateForces(i);
		
		if(i == 2                      )pthread_mutex_unlock(&OverlapMutex[LAST_THREAD]);		
		if(i == Partition_iStart[id] + 2 && id > 0     )pthread_mutex_unlock(&OverlapMutex[id-1]);
		if(i == Partition_iStop[id]   && id < LAST_THREAD)pthread_mutex_unlock(&OverlapMutex[id+1]);
		if(i == n                      )pthread_mutex_unlock(&OverlapMutex[0]);
	} 
	pthread_exit(NULL);  
}

void ThreadsMotionRun(void)
{
	int rc;
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//start all the threads
	for(int myThread = 0; myThread <= LAST_THREAD; myThread++)
	{
		rc = pthread_create(&CalcThread[myThread], &attr, ThreadMotionFunction, (void*)myThread);
		if(rc)
		{
			cerr << "ERROR: Couldn't launch a thread,\n";
			exit(-1);    
		}
	}
	
	//wait for all the threads to finish
	for(int myThread = 0; myThread <= LAST_THREAD; myThread++)
	{
		rc = pthread_join(CalcThread[myThread], NULL);
		if(rc)
		{
			cerr << "ERROR: Couldn't join a thread,\n";
			exit(-1);    
		}
	}
}

void *ThreadMotionFunction(void *threadid)
{
	long id = (long)threadid;
	
	for(int i = Partition_iStart[id]; i <= Partition_iStop[id]; i++)CalculateMotion(i);
	
	pthread_exit(NULL);  
}

void writeDXF(void)
{
	ofstream Julia("MinimalSurfaceNet.dxf");
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	
	for(int j = 0;j <= n;j += cableSpacing)
	{
		Julia<<"0\nPOLYLINE\n";
		if(j == 0)Julia<<"100\nAcDbEntity\n8\nEyeCable one way\n";else Julia<<"100\nAcDbEntity\n8\nNet one way\n";
		Julia<<"100\nAcDb3dPolyline\n66\n1\n";
		Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
		if(j == 0)Julia<<"62\n1\n";else Julia<<"62\n0\n";
		
		for(int i = 0;i <= n;i ++)
		{
			if(NodeType[i][j] <= LastNodeTypeExisting)
			{
				Julia<<"0\nVERTEX\n";
				Julia<<"10\n"<<x[i][j][0]<<"\n";
				Julia<<"20\n"<<x[i][j][1]<<"\n";
				Julia<<"30\n"<<x[i][j][2]<<"\n";
			}
		}
		Julia<<"0\nSEQEND\n";
	}
	
	for(int i = 0;i <= n;i += cableSpacing)
	{
		Julia<<"0\nPOLYLINE\n";
		if(i == 0)Julia<<"100\nAcDbEntity\n8\nEyeCable other way\n";else Julia<<"100\nAcDbEntity\n8\nNet other way\n";
		Julia<<"100\nAcDb3dPolyline\n66\n1\n";
		Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
		if(i == 0)Julia<<"62\n1\n";else Julia<<"62\n0\n";
		
		for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] <= LastNodeTypeExisting)
			{
				Julia<<"0\nVERTEX\n";
				Julia<<"10\n"<<x[i][j][0]<<"\n";
				Julia<<"20\n"<<x[i][j][1]<<"\n";
				Julia<<"30\n"<<x[i][j][2]<<"\n";
			}
		}
		Julia<<"0\nSEQEND\n";
	}
	
	Julia<<"0\nPOLYLINE\n";
	Julia<<"100\nAcDbEntity\n8\nRectangularBoundary\n";
	Julia<<"100\nAcDb3dPolyline\n66\n1\n";
	Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
	Julia<<"62\n1\n";
	
	Julia<<"0\nVERTEX\n10\n"<<- a<<"\n20\n"<<- b<<"\n30\n"<<0.0<<"\n";
	Julia<<"0\nVERTEX\n10\n"<<+ a<<"\n20\n"<<- b<<"\n30\n"<<0.0<<"\n";
	Julia<<"0\nVERTEX\n10\n"<<+ a<<"\n20\n"<<+ b<<"\n30\n"<<0.0<<"\n";
	Julia<<"0\nVERTEX\n10\n"<<- a<<"\n20\n"<<+ b<<"\n30\n"<<0.0<<"\n";
	Julia<<"0\nVERTEX\n10\n"<<- a<<"\n20\n"<<- b<<"\n30\n"<<0.0<<"\n";
	Julia<<"0\nSEQEND\n";
	
	Julia<<"0\nENDSEC\n0\nEOF\n\n";
	Julia.close();
	
	/*Madeleine<<"0\nSECTION\n2\nENTITIES\n";
	
	int SurfaceToRender = 1;
	if(SurfaceToRender == 1)
	{
		for(int j = 0;j <= n - 1;j ++)
		{
			for(int i = 0;i <= n - 1;i ++)
			{
				if(NodeType[i    ][j    ] <= LastNodeTypeExisting &&
				   NodeType[i + 1][j    ] <= LastNodeTypeExisting &&
				   NodeType[i + 1][j + 1] <= LastNodeTypeExisting &&
				   NodeType[i    ][j + 1] <= LastNodeTypeExisting
				   )
				{
					Madeleine<<"0\n3DFACE\n8\nSurface\n";
					Madeleine<<"10\n"<<x[i    ][j    ][0]<<"\n";
					Madeleine<<"20\n"<<x[i    ][j    ][1]<<"\n";
					Madeleine<<"30\n"<<x[i    ][j    ][2]<<"\n";
					Madeleine<<"11\n"<<x[i + 1][j    ][0]<<"\n";
					Madeleine<<"21\n"<<x[i + 1][j    ][1]<<"\n";
					Madeleine<<"31\n"<<x[i + 1][j    ][2]<<"\n";
					Madeleine<<"12\n"<<x[i + 1][j + 1][0]<<"\n";
					Madeleine<<"22\n"<<x[i + 1][j + 1][1]<<"\n";
					Madeleine<<"32\n"<<x[i + 1][j + 1][2]<<"\n";
					Madeleine<<"13\n"<<x[i    ][j + 1][0]<<"\n";
					Madeleine<<"23\n"<<x[i    ][j + 1][1]<<"\n";
					Madeleine<<"33\n"<<x[i    ][j + 1][2]<<"\n";
				}
			}
		}
	}
	
	for(int i = n - q;i <= n - 1;i ++)
	{
		int j = i - (n - q);
		if(NodeType[i    ][j    ] <= LastNodeTypeExisting &&
		   NodeType[i + 1 - (n - q)][j + (n - q)] <= LastNodeTypeExisting &&
		   NodeType[i + 1][j + 1] <= LastNodeTypeExisting &&
		   NodeType[i    ][j + 1] <= LastNodeTypeExisting
		   )
		{
			Madeleine<<"0\n3DFACE\n8\nSurface\n";
			Madeleine<<"10\n"<<x[i    ][j    ][0]<<"\n";
			Madeleine<<"20\n"<<x[i    ][j    ][1]<<"\n";
			Madeleine<<"30\n"<<x[i    ][j    ][2]<<"\n";
			Madeleine<<"11\n"<<x[i + 1 - (n - q)][j + (n - q)][0]<<"\n";
			Madeleine<<"21\n"<<x[i + 1 - (n - q)][j + (n - q)][1]<<"\n";
			Madeleine<<"31\n"<<x[i + 1 - (n - q)][j + (n - q)][2]<<"\n";
			Madeleine<<"12\n"<<x[i + 1][j + 1][0]<<"\n";
			Madeleine<<"22\n"<<x[i + 1][j + 1][1]<<"\n";
			Madeleine<<"32\n"<<x[i + 1][j + 1][2]<<"\n";
			Madeleine<<"13\n"<<x[i    ][j + 1][0]<<"\n";
			Madeleine<<"23\n"<<x[i    ][j + 1][1]<<"\n";
			Madeleine<<"33\n"<<x[i    ][j + 1][2]<<"\n";
		}
	}
	
	for(int j = 0;j <= n;j += cableSpacing)
	{
		for(int i = 0;i <= n - 1;i ++)
		{
			if(NodeType[i][j] <= LastNodeTypeExisting && NodeType[i + 1][j] <= LastNodeTypeExisting)
			{
				if(j >= 1)
					MakeLine(netCableRadius,x[i][j][0],x[i][j][1],x[i][j][2],x[i + 1][j][0],x[i + 1][j][1],x[i + 1][j][2]);
				else
					MakeLine(eyeCableRadius,x[i][j][0],x[i][j][1],x[i][j][2],x[i + 1][j][0],x[i + 1][j][1],x[i + 1][j][2]);
			}
		}
	}
	
	for(int i = 0;i <= n;i += cableSpacing)
	{
		for(int j = 0;j <= n - 1;j ++)
		{
			if(NodeType[i][j] <= LastNodeTypeExisting && NodeType[i][j + 1] <= LastNodeTypeExisting)
			{
				if(i >= 1)
					MakeLine(netCableRadius,x[i][j][0],x[i][j][1],x[i][j][2],x[i][j + 1][0],x[i][j + 1][1],x[i][j + 1][2]);
				else
					MakeLine(eyeCableRadius,x[i][j][0],x[i][j][1],x[i][j][2],x[i][j + 1][0],x[i][j + 1][1],x[i][j + 1][2]);
			}
		}
	}
	
	Madeleine<<"0\nENDSEC\n0\nEOF\n\n";
	Madeleine.close();*/
	
	cout<<"DXF files written, end of program\n\n";
}

void MakeLine(double radius,double x1,double y1,double z1,double x2,double y2,double z2)
{
	double memberlength,Xx,Xy,Xz,Yx,Yy,Yz,Zx,Zy,Zz,thislength;
	
	memberlength=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	
	Zx=(x1-x2)/memberlength;
	Zy=(y1-y2)/memberlength;
	Zz=(z1-z2)/memberlength;
	if(fabs(Zx)<1.0/64.0&&fabs(Zy)<1.0/64.0)
	{
		Xx=Zz;Xy=0.0;Xz=-Zx;
	}
	else
	{
		Xx=-Zy;Xy=+Zx;Xz=0.0;
	}
	thislength=sqrt(Xx*Xx+Xy*Xy+Xz*Xz);
	Xx=Xx/thislength;Xy=Xy/thislength;Xz=Xz/thislength;
	Yx=Zy*Xz-Zz*Xy;
	Yy=Zz*Xx-Zx*Xz;
	Yz=Zx*Xy-Zy*Xx;
	/*Madeleine<<"0\nCIRCLE\n8\nLines\n";
	Madeleine<<"39\n"<<memberlength<<"\n";
	Madeleine<<"10\n"<<x2*Xx+y2*Xy+z2*Xz<<"\n";
	Madeleine<<"20\n"<<x2*Yx+y2*Yy+z2*Yz<<"\n";
	Madeleine<<"30\n"<<x2*Zx+y2*Zy+z2*Zz<<"\n";
	Madeleine<<"40\n"<<radius<<"\n";
	Madeleine<<"210\n"<<Zx<<"\n";
	Madeleine<<"220\n"<<Zy<<"\n";
	Madeleine<<"230\n"<<Zz<<"\n";*/
}

void CableColourControl(int CableColour)
{
	if(CableColour == 0)
	{
		glLineWidth(1.0);
		glColor4f(0.0,0.0,0.0,0.5);
	}
	else
	{
		glLineWidth(1.0);
		glColor4f(0.0,0.0,1.0,0.2);
	}		
}

void FindNormal(int i,int j)
{
	if(i >= 1 && i <= n - 1 && j >= 1 && j <= n - 1)
	{
		double normalLengthSq = 0.0;
		for(int xyz = 0;xyz <= 2;xyz ++)
		{
			int xyzp1 = xyz + 1;if(xyzp1 > 2)xyzp1 -= 3;
			int xyzp2 = xyz + 2;if(xyzp2 > 2)xyzp2 -= 3;
			unitNormal[i][j][xyzp2] = 
			(x[i + 1][j][xyz  ] -  x[i - 1][j][xyz  ]) * (x[i][j + 1][xyzp1] -  x[i][j - 1][xyzp1]) - 
			(x[i + 1][j][xyzp1] -  x[i - 1][j][xyzp1]) * (x[i][j + 1][xyz  ] -  x[i][j - 1][xyz  ]);
			normalLengthSq += unitNormal[i][j][xyzp2] * unitNormal[i][j][xyzp2];
		}
		if(normalLengthSq != 0.0)
		{
			double normalLength = sqrt(normalLengthSq);
			for(int xyz = 0;xyz <= 2;xyz ++)unitNormal[i][j][xyz] /= normalLength;
		}
	}
}