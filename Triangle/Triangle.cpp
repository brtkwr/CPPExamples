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
void CalculateMotion(int);

void ThreadsInit(void);

void* ThreadFunction(void *threadid);
void ThreadsRun(void);

void* ThreadMotionFunction(void *threadid);
void ThreadsMotionRun(void);

void faint(void);
void bold(void);
double FindEdgePoint(int,int,int,int,int,double);

#define  Spacing 3
#define  n        (81 * Spacing)
#define  CarryOver 0.995
#define  MovementFactor 0.1

#define LAST_THREAD  1

pthread_t CalcThread[LAST_THREAD + 1];
pthread_mutex_t OverlapMutex[LAST_THREAD + 1];
int Partition_iStart[LAST_THREAD + 1];
int Partition_iStop[LAST_THREAD + 1];

int  DoNotDrawFaint;

double x[n + 1][n + 1][2],CableNetForce[n + 1][n + 1][2],
velocity[n + 1][n + 1][2],
u[n + 1],v[n + 1],w[n + 1],uDot[n + 1],vDot[n + 1],wDot[n + 1],
Tangent[3],
CableNetTensionCoefficient;

int main(int argc, char * argv[])
{
	cout<<"The program will now run forever until 'q' or ESCAPE is pressed on the keyboard\n";
	cout<<"The space bar starts and stops the analysis\n";
	
	CableNetTensionCoefficient = 0.1;
	
	a = 5000.0;
	
	ControlPointPos[0][0] = - a / 2.0;
	ControlPointPos[0][1] = - sqrt(3.0) * a / 4.0;
	
	ControlPointPos[3][0] = - ControlPointPos[0][0];
	ControlPointPos[3][1] = ControlPointPos[0][1];
	
	ControlPointPos[6][0] = 0.0;
	ControlPointPos[6][1] = - ControlPointPos[0][1];
	
	for(int xy = 0; xy <= 1;xy ++)
	{
		for(int ControlPoint = 1; ControlPoint <= 7; ControlPoint += 3)
		{
			int NextNextControlPoint = ControlPoint + 2;
			if(NextNextControlPoint > lastControlPoint)NextNextControlPoint -= lastControlPoint + 1;
			ControlPointPos[ControlPoint][xy] = (2.0 * ControlPointPos[ControlPoint - 1][xy] + ControlPointPos[NextNextControlPoint][xy]) / 3.0;
		}
		for(int ControlPoint = 2; ControlPoint <= 8; ControlPoint += 3)
		{
			int NextControlPoint = ControlPoint + 1;
			if(NextControlPoint > lastControlPoint)NextControlPoint -= lastControlPoint + 1;
			ControlPointPos[ControlPoint][xy] = (2.0 * ControlPointPos[NextControlPoint][xy] + ControlPointPos[ControlPoint - 2][xy]) / 3.0;
		}
	}
	
	for(int i = 0;i <= n;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			for(int xy = 0; xy <= 1;xy ++)x[i][j][xy]
				= (1.0 - (double) i / (double) n - (double) j / (double) n) * ControlPointPos[0][xy] + 
				((double) i / (double) n) * ControlPointPos[3][xy] + 
				((double) j / (double) n) * ControlPointPos[6][xy];
			
			for(int xy = 0; xy <= 1;xy ++)velocity[i][j][xy] = 0.0;
		}
	}
	
	for(int EdgePoint = 0;EdgePoint <=n; EdgePoint ++)
	{
		u[EdgePoint] = (double) EdgePoint / (double) n;
		v[EdgePoint] = (double) EdgePoint / (double) n;
		w[EdgePoint] = (double) EdgePoint / (double) n;
		
		uDot[EdgePoint] = 0.0;
		vDot[EdgePoint] = 0.0;
		wDot[EdgePoint] = 0.0;
	}
	
	StopStart = 0;
	
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
	glPointSize(5.0);
	for(int ControlPoint = 0; ControlPoint <= lastControlPoint; ControlPoint ++)
	{
		glBegin(GL_POINTS);
		
		if(ControlPointOn[ControlPoint] == 1)glColor4f(1.0,0.0,0.0,1.0);else faint();
		glVertex2dv(ControlPointPos[ControlPoint]);
		glEnd();
	}
	
	glLineWidth(1.0);
	
	for(int j = 0;j <= n - 1;j += Spacing)
	{
		glBegin(GL_LINE_STRIP);
		if(j == 0)bold();else faint();
		for(int i = 0;i <= n - j;i ++)glVertex2dv(x[i][j]);
		glEnd();
	}
	
	for(int i = 0;i <= n - 1;i += Spacing)
	{
		glBegin(GL_LINE_STRIP);
		if(i == 0)bold();else faint();
		for(int j = 0;j <= n - i;j ++)glVertex2dv(x[i][j]);
		glEnd();
	}
	
	for(int i = Spacing;i <= n;i += Spacing)
	{
		glBegin(GL_LINE_STRIP);
		if(i == n)bold();else faint();
		for(int j = 0;j <= i;j ++)glVertex2dv(x[i - j][j]);
		glEnd();
	}
	
	glutSwapBuffers();
	
	for(int CalculationLoop = 0;CalculationLoop <= 0;CalculationLoop++)
	{
		//if(StopStart == 1)
		{
			for(int i = 0;i <= n;i ++)
			{
				for(int j = 0;j <= n;j ++)
				{
					for(int xy = 0; xy <= 1;xy ++)CableNetForce[i][j][xy] = 0.0;
				}
			}
			
			Calculation();
			
			for(int EdgePoint = 0;EdgePoint <= n; EdgePoint ++)
			{
				double EdgeMovementFactor = 0.01;
				double Top = 0.0;
				double Bot = 0.0;
				int i = EdgePoint;
				int j = 0;
				for(int xy = 0; xy <= 1;xy ++)x[i][j][xy] = FindEdgePoint(0,1,2,3,xy,u[EdgePoint]);
				for(int xy = 0; xy <= 1;xy ++)
				{
					Top += CableNetForce[i][j][xy] * Tangent[xy];
					Bot += Tangent[xy] * Tangent[xy];
				}
				if(EdgePoint != 0 && EdgePoint != n)
				{
					uDot[EdgePoint] = CarryOver * uDot[EdgePoint] + EdgeMovementFactor * Top / Bot;
					u[EdgePoint] += uDot[EdgePoint];
				}
				
				Top = 0.0;
				Bot = 0.0;
				i = n - EdgePoint;
				j = EdgePoint;
				for(int xy = 0; xy <= 1;xy ++)x[i][j][xy] = FindEdgePoint(3,4,5,6,xy,v[EdgePoint]);
				for(int xy = 0; xy <= 1;xy ++)
				{
					Top += CableNetForce[i][j][xy] * Tangent[xy];
					Bot += Tangent[xy] * Tangent[xy];
				}
				if(EdgePoint != 0 && EdgePoint != n)
				{
					vDot[EdgePoint] = CarryOver * vDot[EdgePoint] + EdgeMovementFactor * Top / Bot;
					v[EdgePoint] += vDot[EdgePoint];
				}
				
				Top = 0.0;
				Bot = 0.0;
				i = 0;
				j = n - EdgePoint;
				for(int xy = 0; xy <= 1;xy ++)x[0][n - EdgePoint][xy] = FindEdgePoint(6,7,8,0,xy,w[EdgePoint]);
				for(int xy = 0; xy <= 1;xy ++)
				{
					Top += CableNetForce[i][j][xy] * Tangent[xy];
					Bot += Tangent[xy] * Tangent[xy];
				}
				if(EdgePoint != 0 && EdgePoint != n)
				{
					wDot[EdgePoint] = CarryOver * wDot[EdgePoint] + EdgeMovementFactor * Top / Bot;
					w[EdgePoint] += wDot[EdgePoint];
				}				
			}
		}
	}
}

void Calculation(void)
{
#if LAST_THREAD > 0
	ThreadsRun();
	ThreadsMotionRun();
#else
	for(int i = 0;i <= n;i ++)CalculateForces(i);
	for(int i = 0;i <= n;i ++)CalculateMotion(i);
#endif
}

void CalculateForces(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		Axial(i,i + 1,j,j);
		Axial(i,i,j,j + 1);
		Axial(i,i - 1,j,j + 1);
	}
}

void CalculateMotion(int i)
{
	for(int j = 0;j <= n;j ++)
	{			
		if(i > 0 && j > 0 && n - (i + j) > 0)
		{
			for(int xy = 0;xy <= 1;xy ++)
			{
				velocity[i][j][xy] = CarryOver * velocity[i][j][xy] + MovementFactor * CableNetForce[i][j][xy];
				x[i][j][xy] += velocity[i][j][xy];
			}
		}
	}
}

void Axial(int i,int inext,int j,int jnext)
{
	double deltax[3];
	
	
	if((
		i     >= 0 && j     >= 0 && n - (i     + j    ) >= 0 &&
		inext >= 0 && jnext >= 0 && n - (inext + jnext) >= 0) &&
	   (i != 0 || inext != 0) && (j != 0 || jnext != 0) && (i + j != n || inext + jnext != n)
	   )
	{
		for(int xy = 0;xy <= 1;xy ++)deltax[xy] = x[inext][jnext][xy] - x[i][j][xy];
		
		double tensioncoefficient = CableNetTensionCoefficient;
		
		for(int xy = 0;xy <= 1;xy ++)
		{
			double thisForce = tensioncoefficient * deltax[xy];
			CableNetForce[i][j][xy] += thisForce;
			CableNetForce[inext][jnext][xy] -= thisForce;
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
	ofstream Julia("Triangle.dxf");
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	
	int NowBold = 1;
	
	for(int j = 0;j <= n - 1;j += Spacing)
	{
		if(j == 0)NowBold = 1;else NowBold = 0;
		
		Julia<<"0\nPOLYLINE\n";
		if(NowBold == 1)Julia<<"100\nAcDbEntity\n8\nBold\n";else Julia<<"100\nAcDbEntity\n8\nFaint\n";
		Julia<<"100\nAcDb3dPolyline\n66\n1\n";
		Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
		if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";
		
		for(int i = 0;i <= n - j;i ++)
		{
			Julia<<"0\nVERTEX\n";
			Julia<<"10\n"<<x[i][j][0]<<"\n";
			Julia<<"20\n"<<x[i][j][1]<<"\n";
			Julia<<"30\n0.0\n";
		}
		Julia<<"0\nSEQEND\n";
	}
	
	for(int i = 0;i <= n - 1;i += Spacing)
	{
		if(i == 0)NowBold = 1;else NowBold = 0;
		
		Julia<<"0\nPOLYLINE\n";
		if(NowBold == 1)Julia<<"100\nAcDbEntity\n8\nBold\n";else Julia<<"100\nAcDbEntity\n8\nFaint\n";
		Julia<<"100\nAcDb3dPolyline\n66\n1\n";
		Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
		if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";
		
		for(int j = 0;j <= n - i;j ++)
		{
			Julia<<"0\nVERTEX\n";
			Julia<<"10\n"<<x[i][j][0]<<"\n";
			Julia<<"20\n"<<x[i][j][1]<<"\n";
			Julia<<"30\n0.0\n";
		}
		Julia<<"0\nSEQEND\n";
	}
	
	for(int i = Spacing;i <= n;i += Spacing)
	{
		if(i == n)NowBold = 1;else NowBold = 0;
		
		Julia<<"0\nPOLYLINE\n";
		if(NowBold == 1)Julia<<"100\nAcDbEntity\n8\nBold\n";else Julia<<"100\nAcDbEntity\n8\nFaint\n";
		Julia<<"100\nAcDb3dPolyline\n66\n1\n";
		Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
		if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";
		
		for(int j = 0;j <= i;j ++)
		{
			Julia<<"0\nVERTEX\n";
			Julia<<"10\n"<<x[i - j][j][0]<<"\n";
			Julia<<"20\n"<<x[i - j][j][1]<<"\n";
			Julia<<"30\n0.0\n";
		}		
		Julia<<"0\nSEQEND\n";
	}
	
	Julia<<"0\nENDSEC\n0\nEOF\n\n";
	Julia.close();
	cout<<"DXF file written, end of program\n\n";
}

void faint(void)
{
	glColor4f(0.0,0.0,1.0,0.25);
	DoNotDrawFaint = 0;
}

void bold(void)
{
	glColor4f(0.0,0.0,0.0,0.5);
	DoNotDrawFaint = 0;
}

double FindEdgePoint(int Point0,int Point1,int Point2,int Point3,int xy,double this_u)
{
	double p = 0.0 - 3.0 * this_u;
	double q = 1.0 - 3.0 * this_u;
	double r = 2.0 - 3.0 * this_u;
	double s = 3.0 - 3.0 * this_u;
	
	Tangent[xy] = - 3.0 * (
						   + (r * s + q * s + q * r) * ControlPointPos[Point0][xy] / 6.0 +
						   - (r * s + p * s + p * r) * ControlPointPos[Point1][xy] / 2.0
						   + (q * s + p * s + p * q) * ControlPointPos[Point2][xy] / 2.0 +
						   - (q * r + p * r + p * q) * ControlPointPos[Point3][xy] / 6.0);
	
	return
	+ q * r * s * ControlPointPos[Point0][xy] / 6.0 +
	- p * r * s * ControlPointPos[Point1][xy] / 2.0
	+ p * q * s * ControlPointPos[Point2][xy] / 2.0 +
	- p * q * r * ControlPointPos[Point3][xy] / 6.0;
}