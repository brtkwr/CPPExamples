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
void Axial(int,int,int,int,int);
void Initialise(int);
void CalculateMotion(int);

void ThreadsInit(void);

void* ThreadFunction(void *threadid);
void ThreadsRun(void);

void* ThreadMotionFunction(void *threadid);
void ThreadsMotionRun(void);

int jOnCable(int);

void JoinUp(void);
void MakeLine(double,double,double,double,double,double,double);

void CableColourControl(int);

#define  fine 5
#define  cableSpacingcontrol 2
#define  cableSpacing (2 * cableSpacingcontrol + 1)
#define  mOver2 (4 * fine * cableSpacing + cableSpacingcontrol)
#define  m (mOver2 * 2)
#define  n (mOver2 + 2 * cableSpacing * fine)

#define  CarryOver 0.99
#define  MovementFactor 0.2
#define  SupportCanMoveHorizontally 0

#define LAST_THREAD  1

pthread_t CalcThread[LAST_THREAD + 1];
pthread_mutex_t OverlapMutex[LAST_THREAD + 1];
int Partition_iStart[LAST_THREAD + 1];
int Partition_iStop[LAST_THREAD + 1];

int NodeType[m + 1][n + 1];

double x[m + 1][n + 1][3],CableNetForce[m + 1][n + 1][3],EyeCableForce[m + 1][n + 1][3],
velocity[m + 1][n + 1][3],
Stiffness[m + 1][n + 1],CableLength[m + 1][n + 1][2],
SupportPosition[3],SupportVelocity[3],SupportStiffness,SupportForce[3],
PointPosition[2],EyeCableTension,EAEyeCable,EyeCableSlackLength,EyeCableLength,CableNetTensionCoefficient,
totalLength_i,totalLength_j,netCableRadius,eyeCableRadius;

char NumericalValue[101];
string MyText;
//ofstream Madeleine("MinimalSurface.dxf");

int main(int argc, char * argv[])
{
	cout<<"m = "<<m<<"\n";
	cout<<"n = "<<n<<"\n";
	cout<<"The program will now run forever until 'q' or ESCAPE is pressed on the keyboard\n";
	cout<<"Left mouse button for pan\n";
	cout<<"Right mouse button (or CONTROL left button) for rotate\n";
	cout<<"Move spot to zoom\n";
	
	supportRadius = 100.0;
	SupportHeight = 0.11 * supportRadius;
	
	netCableRadius = supportRadius / 10000.0;
	eyeCableRadius = supportRadius / 5000.0;
	
	CableNetTensionCoefficient = 0.1;
	EyeCableSlackLength = 0.5 * supportRadius;
	
	EAEyeCable = 0.5 * (double)fine * CableNetTensionCoefficient * EyeCableSlackLength;
	
	double setupmultiplier = 8.0 * atan(1.0) / double(m + 1);
	for(int i = 0;i <= m;i ++)
	{
		double theta = double(i) * setupmultiplier;
		for(int j = 0;j <= n;j ++)
		{
			double phi = - double(j) * setupmultiplier;
			double radius = supportRadius * exp(phi);
			x[i][j][0] = radius * cos(theta);
			x[i][j][1] = radius * sin(theta);
			x[i][j][2] = SupportHeight * double(j) / double(n + 1);
			
			NodeType[i][j] = 0;
			if(j == jOnCable(i))NodeType[i][j] = 1;
			if(j > jOnCable(i))NodeType[i][j] = 2;
			
			for(int xyz = 0; xyz <= 2;xyz ++)velocity[i][j][xyz] = 0.0;
			for(int CableDirection = 0; CableDirection <= 1;CableDirection ++)CableLength[i][j][CableDirection] = 0.0;
		}
		Initialise(i);
		SupportPosition[0] = - supportRadius * exp(- double(n + 1) * setupmultiplier);
		SupportPosition[1] =  0.0;
		SupportPosition[2] = SupportHeight;
		
		SupportHorizontalPosition = SupportPosition[0];
		
		for(int xyz = 0; xyz <= 2;xyz ++)SupportVelocity[xyz] = 0.0;
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
	{
		glRasterPos2d(ButtonX[0],0.65);
		MyText="Ratio of total cable lengths = ";
		int LengthOfString = MyText.length();
		for(int MyCharacter = 0; MyCharacter < LengthOfString; MyCharacter++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyCharacter]);
		sprintf(NumericalValue,"%8.6f",totalLength_i / totalLength_j);
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
	
	int CableColour = jOnCable(0);
	do
	{
		CableColour -= cableSpacing;
	}
	while(CableColour >= cableSpacing);
	
	for(int j = 0;j <= n;j ++)
	{
		if(CableColour >= cableSpacing)CableColour -= cableSpacing;
		CableColourControl(CableColour);
		CableColour ++;
		
		if(j == 0)glColor4f(0.0,0.0,1.0,0.25);
		glBegin(GL_LINE_STRIP);
		for(int i = 0;i <= m;i ++)
		{
			if(NodeType[i][j] <= 1)glVertex3dv(x[i][j]);
		}
		if(NodeType[0][j] <= 1)glVertex3dv(x[0][j]);
		glEnd();
	}
	
	CableColour = 0;
	for(int i = 0;i <= m;i ++)
	{
		if(CableColour >= cableSpacing)CableColour -= cableSpacing;
		CableColourControl(CableColour);
		CableColour ++;
		glBegin(GL_LINE_STRIP);
		for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] <= 1)glVertex3dv(x[i][j]);
		}
		glEnd();
	}
	glColor4f(1.0,0.0,0.0,0.5);
	glBegin(GL_LINE_STRIP);
	for(int i = 0;i <= m;i ++)
	{
		glVertex3dv(x[i][jOnCable(i)]);
		if(i == mOver2)glVertex3dv(SupportPosition);
	}
	glVertex3dv(x[0][jOnCable(0)]);
	glEnd();
	
	glutSwapBuffers();
	for(int CalculationLoop = 0;CalculationLoop <= 50;CalculationLoop++)
	{
		totalLength_i = 0.0;
		totalLength_j = 0.0;
		
		for(int i = 0;i <= m;i ++)
		{
			for(int j = 0;j <= n;j ++)
			{
				if(NodeType[i][j] == 0)
				{
					totalLength_i += CableLength[i][j][0];
					totalLength_j += CableLength[i][j][1];
				}
				for(int CableDirection = 0; CableDirection <= 1;CableDirection ++)CableLength[i][j][CableDirection] = 0.0;
			}
		}
		
		EyeCableSlackLength += 0.00001 * (totalLength_i - totalLength_j);
		
		EyeCableLength = 0.0;
		for(int i = 0;i <= m;i ++)
		{
			double lengthsq = 0.0;
			int inext = i + 1;
			if(inext > m)inext = 0;
			if(i != mOver2)
			{
				for(int xyz = 0; xyz <= 2;xyz ++)lengthsq += 
					(x[i + 1][jOnCable(inext)][xyz] - x[i][jOnCable(i)][xyz]) *
					(x[i + 1][jOnCable(inext)][xyz] - x[i][jOnCable(i)][xyz]);
				EyeCableLength += sqrt(lengthsq);
			}
			else
			{
				for(int whichSide = 0;whichSide <= 1; whichSide++)
				{
					for(int xyz = 0; xyz <= 2;xyz ++)lengthsq += 
						(SupportPosition[xyz] - x[mOver2 + whichSide][jOnCable(inext)][xyz]) *
						(SupportPosition[xyz] - x[mOver2 + whichSide][jOnCable(inext)][xyz]);
					EyeCableLength += sqrt(lengthsq);
				}
			}
		}
		EyeCableTension = EAEyeCable * (EyeCableLength - EyeCableSlackLength) / EyeCableSlackLength;
		
		Calculation();
		
		if(SupportCanMoveHorizontally == 1)
		{
			for(int xy = 0;xy <= 1; xy ++)
			{
				SupportVelocity[xy] = CarryOver * SupportVelocity[xy] + MovementFactor * SupportForce[xy] / SupportStiffness;
				SupportPosition[xy] += SupportVelocity[xy];
			}
		}
		else SupportPosition[0] = SupportHorizontalPosition;
		
		SupportPosition[2] = SupportHeight;
	}
}

void Calculation(void)
{
#if LAST_THREAD > 0
	ThreadsRun();
	ThreadsMotionRun();
#else
	for(int i = 0;i <= m;i ++)CalculateForces(i);
	for(int i = 0;i <= m;i ++)CalculateMotion(i);
#endif
}

void Initialise(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		Stiffness[i][j] = 0.0;
		for(int xyz = 0; xyz <= 2;xyz ++)
		{
			CableNetForce[i][j][xyz] = 0.0;
			EyeCableForce[i][j][xyz] = 0.0;
		}
	}
}

void CalculateForces(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		Axial(i,i + 1,j,j,0);
		if(j + 1 <= n)Axial(i,i,j,j + 1,0);
	}
	Axial(i,i + 1,jOnCable(i),jOnCable(i + 1),1);
}

void CalculateMotion(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		if(j == 0)
		{
			float planDistanceFromOriginSq = 0.0;
			float dotProduct = 0.0;
			CableNetForce[i][j][2] = 0.0;
			//THIS ONLY REMOVES PLAN COMPONENT OF FORCE
			for(int xy = 0;xy <= 1;xy ++)
			{
				planDistanceFromOriginSq += x[i][j][xy] * x[i][j][xy];
				dotProduct += x[i][j][xy] * CableNetForce[i][j][xy];
			}
			
			for(int xy = 0;xy <= 1;xy ++)
				CableNetForce[i][j][xy] -= x[i][j][xy] * dotProduct / planDistanceFromOriginSq;
			
			float planDistanceFromOrigin = sqrt(planDistanceFromOriginSq);
			
			x[i][j][2] = 0.0;
			for(int xy = 0;xy <= 1;xy ++)
				x[i][j][xy] -= (planDistanceFromOrigin - supportRadius) * x[i][j][xy] / planDistanceFromOrigin;
			
		}
		if(NodeType[i][j] <= 1)
		{
			for(int xyz = 0;xyz <= 2;xyz ++)
			{
				velocity[i][j][xyz] = CarryOver * velocity[i][j][xyz] + MovementFactor * 
				(CableNetForce[i][j][xyz] + EyeCableForce[i][j][xyz]) / Stiffness[i][j];
				x[i][j][xyz] += velocity[i][j][xyz];
			}
		}
	}
	Initialise(i);
}

void Axial(int i,int inext_trial,int j,int jnext,int CableType)
{
	double deltax[3];
	int inext = inext_trial;
	if(inext > m)inext = 0;
	if(NodeType[i][j] <= 1 && NodeType[inext][jnext] <= 1)
	{
		if(CableType == 0)
		{
			for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[inext][jnext][xyz] - x[i][j][xyz];
			
			if(inext == i || jnext == j)
			{
				double lengthSq = 0.0;
				
				for(int xyz = 0;xyz <= 2;xyz ++)lengthSq += deltax[xyz] * deltax[xyz];
				
				double length = sqrt(lengthSq);
				
				int CableDirection = 0;
				
				if(jnext == j)CableDirection = 1;
				CableLength[i][j][CableDirection] += length;
				CableLength[inext][jnext][CableDirection] += length;
			}
			
			double tensioncoefficient = CableNetTensionCoefficient;
			
			if(j == 0 && jnext == 0)tensioncoefficient = tensioncoefficient / 2.0;//This is for the edge ring
			
			Stiffness[i][j] += tensioncoefficient;
			Stiffness[inext][jnext] += tensioncoefficient;
			
			for(int xyz = 0;xyz <= 2;xyz ++)
			{
				double thisForce = tensioncoefficient * deltax[xyz];
				CableNetForce[i][j][xyz] += thisForce;
				CableNetForce[inext][jnext][xyz] -= thisForce;
			}
		}
		else
		{
			if(i != mOver2)
			{
				for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[inext][jnext][xyz] - x[i][j][xyz];
				double lengthsq = 0.0;
				for(int xyz = 0;xyz <= 2;xyz ++)lengthsq += deltax[xyz] * deltax[xyz];
				double tensioncoefficient = EyeCableTension / sqrt(lengthsq);
				
				Stiffness[i][j] += tensioncoefficient;
				Stiffness[inext][jnext] += tensioncoefficient;
				for(int xyz = 0;xyz <= 2;xyz ++)
				{
					double thisForce = tensioncoefficient * deltax[xyz];
					EyeCableForce[i][j][xyz] += thisForce;
					EyeCableForce[inext][jnext][xyz] -= thisForce;
				}
			}
			else
			{
				for(int xyz = 0;xyz <= 2;xyz ++)SupportForce[xyz] = 0.0;
				SupportStiffness = 0.0;
				for(int whichSide = 0;whichSide <= 1; whichSide++)
				{
					for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = SupportPosition[xyz] - x[mOver2 + whichSide][jOnCable(inext)][xyz];
					double lengthsq = 0.0;
					for(int xyz = 0;xyz <= 2;xyz ++)lengthsq += deltax[xyz] * deltax[xyz];
					double tensioncoefficient = EyeCableTension / sqrt(lengthsq);
					
					Stiffness[i + whichSide][j] += tensioncoefficient;
					SupportStiffness += tensioncoefficient;
					for(int xyz = 0;xyz <= 2;xyz ++)
					{
						double thisForce = tensioncoefficient * deltax[xyz];
						EyeCableForce[i + whichSide][j][xyz] += thisForce;
						SupportForce[xyz] -= thisForce;
					}
				}
			}
		}
	}
}

void ThreadsInit(void)
{
	Partition_iStart[0] = 0;
	
	for(int myThread = 0; myThread <= LAST_THREAD - 1; myThread++)
	{
		Partition_iStop[myThread] = ((myThread + 1) * m) / (LAST_THREAD + 1);
		Partition_iStart[myThread + 1] = Partition_iStop[myThread] + 1;
	}
	
	Partition_iStop[LAST_THREAD] = m;
	
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
		if(i == m - 2                   )pthread_mutex_lock(&OverlapMutex[0]);
		
		CalculateForces(i);
		
		if(i == 2                      )pthread_mutex_unlock(&OverlapMutex[LAST_THREAD]);		
		if(i == Partition_iStart[id] + 2 && id > 0     )pthread_mutex_unlock(&OverlapMutex[id-1]);
		if(i == Partition_iStop[id]   && id < LAST_THREAD)pthread_mutex_unlock(&OverlapMutex[id+1]);
		if(i == m                      )pthread_mutex_unlock(&OverlapMutex[0]);
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

int jOnCable(int i)
{
	if(i <= mOver2)return n + (i - mOver2);
	else return n - (i - (mOver2 + 1));
}

void writeDXF(void)
{
	ofstream Julia("MinimalSurfaceNet.dxf");
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	
	for(int j = 0;j <= n;j ++)
	{
		Julia<<"0\nPOLYLINE\n";
		if(j == 0)Julia<<"100\nAcDbEntity\n8\nEdge\n";else Julia<<"100\nAcDbEntity\n8\nNet\n";
		Julia<<"100\nAcDb3dPolyline\n66\n1\n";
		Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
		if(j == 0)Julia<<"62\n2\n";else Julia<<"62\n0\n";
		
		for(int trial_i = 0;trial_i <= m + 1;trial_i ++)
		{
			int i = trial_i;
			if(i > m)i = 0;
			if(NodeType[i][j] <= 1)
			{
				
				Julia<<"0\nVERTEX\n";
				if(j == 0)Julia<<"100\nAcDbEntity\n8\nEdge\n";else Julia<<"100\nAcDbEntity\n8\nNet\n";
				Julia<<"100\nAcDb3dPolylineVertex\n";
				Julia<<"10\n"<<x[i][j][0]<<"\n";
				Julia<<"20\n"<<x[i][j][1]<<"\n";
				Julia<<"30\n"<<x[i][j][2]<<"\n";
				if(j == 0)Julia<<"62\n2\n";else Julia<<"62\n0\n";
			}
		}
		Julia<<"0\nSEQEND\n";
	}
	
	for(int i = 0;i <= m;i ++)
	{
		Julia<<"0\nPOLYLINE\n";
		Julia<<"100\nAcDbEntity\n8\nNet\n";
		Julia<<"100\nAcDb3dPolyline\n66\n1\n";
		Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
		Julia<<"62\n0\n";
		
		for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] <= 1)
			{
				
				Julia<<"0\nVERTEX\n";
				Julia<<"100\nAcDbEntity\n8\nNet\n";
				Julia<<"100\nAcDb3dPolylineVertex\n";
				Julia<<"10\n"<<x[i][j][0]<<"\n";
				Julia<<"20\n"<<x[i][j][1]<<"\n";
				Julia<<"30\n"<<x[i][j][2]<<"\n";
				Julia<<"62\n0\n";
			}
		}
		Julia<<"0\nSEQEND\n";
	}
	
	Julia<<"0\nPOLYLINE\n";
	Julia<<"100\nAcDbEntity\n8\nEyeCable\n";
	Julia<<"100\nAcDb3dPolyline\n66\n1\n";
	Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
	Julia<<"62\n1\n";
	
	Julia<<"0\nVERTEX\n";
	Julia<<"100\nAcDbEntity\n8\nEyeCable\n";
	Julia<<"100\nAcDb3dPolylineVertex\n";
	Julia<<"10\n"<<SupportPosition[0]<<"\n";
	Julia<<"20\n"<<SupportPosition[1]<<"\n";
	Julia<<"30\n"<<SupportPosition[2]<<"\n";
	Julia<<"62\n1\n";
	
	for(int this_i = mOver2;this_i >= - mOver2;this_i --)
	{
		int i = this_i;
		if(i < 0) i = m + i + 1;
		Julia<<"0\nVERTEX\n";
		Julia<<"100\nAcDbEntity\n8\nEyeCable\n";
		Julia<<"100\nAcDb3dPolylineVertex\n";
		Julia<<"10\n"<<x[i][jOnCable(i)][0]<<"\n";
		Julia<<"20\n"<<x[i][jOnCable(i)][1]<<"\n";
		Julia<<"30\n"<<x[i][jOnCable(i)][2]<<"\n";
		Julia<<"62\n1\n";
	}
	
	Julia<<"0\nVERTEX\n";
	Julia<<"100\nAcDbEntity\n8\nEyeCable\n";
	Julia<<"100\nAcDb3dPolylineVertex\n";
	Julia<<"10\n"<<SupportPosition[0]<<"\n";
	Julia<<"20\n"<<SupportPosition[1]<<"\n";
	Julia<<"30\n"<<SupportPosition[2]<<"\n";
	Julia<<"62\n1\n";
	
	Julia<<"0\nSEQEND\n";
	
	Julia<<"0\nENDSEC\n0\nEOF\n\n";
	Julia.close();
	
	/*Madeleine<<"0\nSECTION\n2\nENTITIES\n";
	
	int SurfaceToRender = 0;
	if(SurfaceToRender == 1)
	{
	for(int j = 0;j <= n - 1;j ++)
	{
		for(int i = 0;i <= m;i ++)
		{
			int nexti = i + 1;
			if(nexti > m)nexti = 0;
			if(NodeType[i    ][j    ] +
			   NodeType[nexti][j    ] +
			   NodeType[nexti][j + 1] +
			   NodeType[i    ][j    ] <= 4
			   )
			{
				Madeleine<<"0\n3DFACE\n8\nSurface\n";
				Madeleine<<"10\n"<<x[i    ][j    ][0]<<"\n";
				Madeleine<<"20\n"<<x[i    ][j    ][1]<<"\n";
				Madeleine<<"30\n"<<x[i    ][j    ][2]<<"\n";
				Madeleine<<"11\n"<<x[nexti][j    ][0]<<"\n";
				Madeleine<<"21\n"<<x[nexti][j    ][1]<<"\n";
				Madeleine<<"31\n"<<x[nexti][j    ][2]<<"\n";
				
				if(NodeType[nexti][j + 1] <= 1)
				{
					Madeleine<<"12\n"<<x[nexti][j + 1][0]<<"\n";
					Madeleine<<"22\n"<<x[nexti][j + 1][1]<<"\n";
					Madeleine<<"32\n"<<x[nexti][j + 1][2]<<"\n";
				}
				else
				{
					Madeleine<<"12\n"<<x[i    ][j + 1][0]<<"\n";
					Madeleine<<"22\n"<<x[i    ][j + 1][1]<<"\n";
					Madeleine<<"32\n"<<x[i    ][j + 1][2]<<"\n";
				}
				
				if(NodeType[i    ][j + 1] <= 1)
				{
					Madeleine<<"13\n"<<x[i    ][j + 1][0]<<"\n";
					Madeleine<<"23\n"<<x[i    ][j + 1][1]<<"\n";
					Madeleine<<"33\n"<<x[i    ][j + 1][2]<<"\n";
				}
				else
				{
					Madeleine<<"13\n"<<x[nexti][j + 1][0]<<"\n";
					Madeleine<<"23\n"<<x[nexti][j + 1][1]<<"\n";
					Madeleine<<"33\n"<<x[nexti][j + 1][2]<<"\n";
				}
			}
		}
	}
	
	Madeleine<<"0\n3DFACE\n8\nSurface\n";
	Madeleine<<"10\n"<<x[mOver2][n][0]<<"\n";
	Madeleine<<"20\n"<<x[mOver2][n][1]<<"\n";
	Madeleine<<"30\n"<<x[mOver2][n][2]<<"\n";
	Madeleine<<"11\n"<<SupportPosition[0]<<"\n";
	Madeleine<<"21\n"<<SupportPosition[1]<<"\n";
	Madeleine<<"31\n"<<SupportPosition[2]<<"\n";
	Madeleine<<"12\n"<<x[mOver2 + 1][n][0]<<"\n";
	Madeleine<<"22\n"<<x[mOver2 + 1][n][1]<<"\n";
	Madeleine<<"32\n"<<x[mOver2 + 1][n][2]<<"\n";
	Madeleine<<"13\n"<<x[mOver2][n][0]<<"\n";
	Madeleine<<"23\n"<<x[mOver2][n][1]<<"\n";
	Madeleine<<"33\n"<<x[mOver2][n][2]<<"\n";
	}
	
	int CableColour = jOnCable(0);
	do
	{
		CableColour -= cableSpacing;
	}
	while(CableColour >= cableSpacing);
	
	for(int j = 0;j <= n;j ++)
	{
		if(CableColour >= cableSpacing)CableColour -= cableSpacing;
		if(CableColour == 0)
		{
			for(int i = 0;i <= m;i ++)
			{
				int nexti = i + 1;
				if(nexti > m)nexti -= m + 1;
				if(NodeType[i][j] <= 1 && NodeType[nexti][j] <= 1)
				{
					if(j > 0)
						MakeLine(netCableRadius,x[i][j][0],x[i][j][1],x[i][j][2],x[nexti][j][0],x[nexti][j][1],x[nexti][j][2]);
					else
						MakeLine(eyeCableRadius,x[i][j][0],x[i][j][1],x[i][j][2],x[nexti][j][0],x[nexti][j][1],x[nexti][j][2]);
				}
			}
		}
		CableColour ++;
	}
	
	CableColour = 0;
	for(int i = 0;i <= m;i ++)
	{
		if(CableColour >= cableSpacing)CableColour -= cableSpacing;
		if(CableColour == 0)
		{
			for(int j = 0;j <= n - 1;j ++)
			{
				if(NodeType[i][j] <= 1 && NodeType[i][j + 1] <= 1)
					MakeLine(netCableRadius,x[i][j][0],x[i][j][1],x[i][j][2],x[i][j + 1][0],x[i][j + 1][1],x[i][j + 1][2]);
			}
		}
		CableColour ++;
	}
	
	for(int i = 0;i <= m;i ++)
	{
		int nexti = i + 1;
		if(nexti > m)nexti -= m + 1;
		if(i != mOver2)
			MakeLine(eyeCableRadius,
					 x[i][jOnCable(i)][0],
					 x[i][jOnCable(i)][1],
					 x[i][jOnCable(i)][2],
					 x[nexti][jOnCable(nexti)][0],
					 x[nexti][jOnCable(nexti)][1],
					 x[nexti][jOnCable(nexti)][2]);
		else
		{
			MakeLine(eyeCableRadius,
					 x[i][jOnCable(i)][0],
					 x[i][jOnCable(i)][1],
					 x[i][jOnCable(i)][2],
					 SupportPosition[0],
					 SupportPosition[1],
					 SupportPosition[2]);
			MakeLine(eyeCableRadius,
					 SupportPosition[0],
					 SupportPosition[1],
					 SupportPosition[2],
					 x[nexti][jOnCable(nexti)][0],
					 x[nexti][jOnCable(nexti)][1],
					 x[nexti][jOnCable(nexti)][2]);
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
	glLineWidth(0.5);
	glColor4f(0.0,0.0,0.0,0.5);
}		
}