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

void ThreadsInit(void);

void* ThreadFunction(void *threadid);
void ThreadsRun(void);

void* ThreadMotionFunction(void *threadid);
void ThreadsMotionRun(void);

void faint(void);
void bold(void);

#define  HalfSpacing 5
#define  Spacing (2 * HalfSpacing)
#define  Fineness 3
#define  n_over_6 (2 * Fineness * Spacing)
#define  n_over_3 (2 * n_over_6)
#define  n_over_2 (3 * n_over_6)
#define  n        (6 * n_over_6)
#define  CarryOver 0.99
#define  MovementFactor 1.0

#define  Hexagon 1
#define  Sphere 1

#define LAST_THREAD  1

pthread_t CalcThread[LAST_THREAD + 1];
pthread_mutex_t OverlapMutex[LAST_THREAD + 1];
int Partition_iStart[LAST_THREAD + 1];
int Partition_iStop[LAST_THREAD + 1];

int  NodeType[n + 1][n + 1],DoNotDrawFaint;

double x[n + 1][n + 1][3],CableNetForce[n + 1][n + 1][3],
velocity[n + 1][n + 1][3],
PointPosition[2],
CableNetTensionCoefficient;

char NumericalValue[101];
string MyText;

int main(int argc, char * argv[])
{
	cout<<"The program will now run forever until 'q' or ESCAPE is pressed on the keyboard\n";
	cout<<"Left mouse button for pan\n";
	cout<<"Right mouse button (or CONTROL left button) for rotate\n";
	cout<<"Move spot to zoom\n";
	cout<<"The space bar starts and stops the analysis\n";
	
	CableNetTensionCoefficient = 0.1;
	
	a = 5000.0;
	double b;
	if(Hexagon == 1)b = 3.0 * a / sqrt(10.0);else b = a / sqrt(2.0);
	
	for(int i = 0;i <= n;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			NodeType[i][j] = 0;
			
			if(i == 0 || i == n || j == 0 || j == n)NodeType[i][j] = 1;
			if(Hexagon == 1)
			{
				if(i + j  < n_over_3)NodeType[i][j] = - 1;
				if(i + j == n_over_3)NodeType[i][j] = + 1;
				
				if(i + j  > 2 * n - n_over_3)NodeType[i][j] = - 1;
				if(i + j == 2 * n - n_over_3)NodeType[i][j] = + 1;
				
				if(i - j  > n - n_over_3)NodeType[i][j] = - 1;
				if(i - j == n - n_over_3)NodeType[i][j] = + 1;
				
				if(j - i  > n - n_over_3)NodeType[i][j] = - 1;
				if(j - i == n - n_over_3)NodeType[i][j] = + 1;
			}
								}
	}
	
	for(int i = 0;i <= n;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			x[i][j][0] = b * (double) (2 * i - n) / (double) n;
			x[i][j][1] = b * (double) (2 * j - n) / (double) n;
			if(Sphere == 1)
			{
				double zSquared = a * a - x[i][j][0] * x[i][j][0] -  x[i][j][1] * x[i][j][1];
				if(zSquared > 0.0)x[i][j][2] = sqrt(zSquared);else x[i][j][2] = 0.0;
			}
			else x[i][j][2] = 0.0;
			
			for(int xyz = 0; xyz <= 2;xyz ++)velocity[i][j][xyz] = 0.0;
		}
	}
	
	StopStart = 0;
	
	for(int i = 0;i <= n;i ++)Initialise(i);
	
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
		
	glPointSize(5.0);
	for(int Button = 0; Button <= lastButton; Button ++)
	{
		glBegin(GL_POINTS);
		
		if(ButtonOn[Button] == 1)glColor4f(1.0,0.0,0.0,1.0);else faint();
		PointPosition[0] = ButtonX[Button];
		PointPosition[1] = ButtonY[Button];
		glVertex2dv(PointPosition);
		glEnd();
		
		glColor4f(0.0,0.0,0.0,0.2);
		if(Button == 0)MyText = "Zoom";
		if(Button == 1)MyText = "Circle size";
		if(Button == 2)MyText = "Difference in height";
		
		glRasterPos2d(ButtonX[Button] + 0.02,ButtonY[Button] - 0.01);
		
		int StringLength = MyText.length();
		for(int MyChar = 0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
	}
	
	glPopMatrix();
	
	glLineWidth(1.0);
	
	int mirrorContol;
	if(Sphere == 1) mirrorContol = 1;else mirrorContol = 0;
	for(int mirror = - mirrorContol; mirror <= mirrorContol; mirror += 2)
	{
	int LineDrawStep = HalfSpacing - 1;
	
		int i_stop;
		//if(Sphere == 1)i_stop = n_over_2;else 
			i_stop = n;
	for(int j = 0;j <= n;j ++)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			bold();
		}
		else faint();
		
		int NowDrawing = 0;
		for(int i = 0;i <= i_stop;i ++)
		{
			if(NodeType[i][j] >= 0)
			{
				if(NowDrawing == 0){NowDrawing = 1;glBegin(GL_LINE_STRIP);}
				glVertex3d(x[i][j][0],x[i][j][1],(double) mirror * x[i][j][2]);
			}
			else
			{
				if(NowDrawing == 1){NowDrawing = 0;glEnd();}
			}
		}
		glEnd();
	}

	LineDrawStep = HalfSpacing - 1;
	for(int i = 0;i <= i_stop;i ++)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			bold();
		}
		else faint();
		
		int NowDrawing = 0;
		for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] >= 0)
			{
				if(NowDrawing == 0){NowDrawing = 1;glBegin(GL_LINE_STRIP);}
				glVertex3d(x[i][j][0],x[i][j][1],(double) mirror * x[i][j][2]);
			}
			else
			{
				if(NowDrawing == 1){NowDrawing = 0;glEnd();}
			}
		}
		glEnd();
	}
}

	glutSwapBuffers();
	if(StopStart == 1)
	{
	for(int CalculationLoop = 0;CalculationLoop <= 0;CalculationLoop++)Calculation();
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

void Initialise(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		for(int xyz = 0; xyz <= 2;xyz ++)CableNetForce[i][j][xyz] = 0.0;
	}
}

void CalculateForces(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		if(i + 1 <= n)Axial(i,i + 1,j,j);
			if(j + 1 <= n)Axial(i,i,j,j + 1);
		}
	}

void CalculateMotion(int i)
{
	for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] >= 0)
			{
				if(NodeType[i][j] == 1 || Sphere == 1)
				{
					double RadiusSq = 0.0;
					double forcedotProduct = 0.0;
					double velocitydotProduct = 0.0;
					for(int xyz = 0;xyz <= 2;xyz ++)
					{
						RadiusSq += x[i][j][xyz] * x[i][j][xyz];
						forcedotProduct += x[i][j][xyz] * CableNetForce[i][j][xyz];
						velocitydotProduct += x[i][j][xyz] * velocity[i][j][xyz];
					}
					double Radius = sqrt(RadiusSq);
					for(int xyz = 0;xyz <= 2;xyz ++)
					{
						CableNetForce[i][j][xyz] -= forcedotProduct * x[i][j][xyz] / RadiusSq;
						velocity[i][j][xyz] -= velocitydotProduct * x[i][j][xyz] / RadiusSq;
						x[i][j][xyz] *= a / Radius;
					}
				}
					if(NodeType[i][j] == 1)
					{
						CableNetForce[i][j][2] = 0.0;
					velocity[i][j][2] = 0.0;
					x[i][j][2] = 0.0;
				}
				
				for(int xyz = 0;xyz <= 2;xyz ++)
				{
					if(NodeType[i][j] >= 0)
					{
						velocity[i][j][xyz] = CarryOver * velocity[i][j][xyz] + MovementFactor * CableNetForce[i][j][xyz];
						x[i][j][xyz] += velocity[i][j][xyz];
					}
				}
			}
		}
	Initialise(i);
}

void Axial(int i,int inext,int j,int jnext)
{
	double deltax[3];
	
	if(NodeType[i][j] >= 0 && NodeType[inext][jnext] >= 0)
	{
		for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[inext][jnext][xyz] - x[i][j][xyz];
				
		double tensioncoefficient = CableNetTensionCoefficient;
		
		if(NodeType[i][j] > 0 && NodeType[inext][jnext] > 0)tensioncoefficient /= 2.0;
		
		for(int xyz = 0;xyz <= 2;xyz ++)
		{
			double thisForce = tensioncoefficient * deltax[xyz];
			CableNetForce[i][j][xyz] += thisForce;
			CableNetForce[inext][jnext][xyz] -= thisForce;
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
	
	int mirrorContol;
	if(Sphere == 1) mirrorContol = 1;else mirrorContol = 0;
	for(int mirror = - mirrorContol; mirror <= mirrorContol; mirror += 2)
	{
		int NowBold = 1;
	int LineDrawStep = HalfSpacing - 1;
		int i_stop;
		//if(Sphere == 1)i_stop = n_over_2;else 
			i_stop = n;
	for(int j = 0;j <= n;j ++)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			NowBold = 1;
		}
		else NowBold = 0;
		
		int NowDrawing = 0;
		for(int i = 0;i <= i_stop;i ++)
		{
			if(NodeType[i][j] >= 0)
			{
				if(NowDrawing == 0)
				{
					NowDrawing = 1;
					Julia<<"0\nPOLYLINE\n";
					if(NowBold == 1)Julia<<"100\nAcDbEntity\n8\nBold\n";else Julia<<"100\nAcDbEntity\n8\nFaint\n";
					Julia<<"100\nAcDb3dPolyline\n66\n1\n";
					Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
					if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";
				}
				Julia<<"0\nVERTEX\n";
				Julia<<"10\n"<<x[i][j][0]<<"\n";
				Julia<<"20\n"<<x[i][j][1]<<"\n";
				Julia<<"30\n"<<(double) mirror * x[i][j][2]<<"\n";
			}
			else
			{
				if(NowDrawing == 1){NowDrawing = 0;Julia<<"0\nSEQEND\n";}
			}
		}
		Julia<<"0\nSEQEND\n";
	}
	
	LineDrawStep = HalfSpacing - 1;
	
	for(int i = 0;i <= i_stop;i ++)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			NowBold = 1;
		}
		else NowBold = 0;
		
		int NowDrawing = 0;
		for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] >= 0)
			{
				if(NowDrawing == 0)
				{
					NowDrawing = 1;
					Julia<<"0\nPOLYLINE\n";
					if(NowBold == 1)Julia<<"100\nAcDbEntity\n8\nBold\n";else Julia<<"100\nAcDbEntity\n8\nFaint\n";
					Julia<<"100\nAcDb3dPolyline\n66\n1\n";
					Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
					if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";
				}
				Julia<<"0\nVERTEX\n";
				Julia<<"10\n"<<x[i][j][0]<<"\n";
				Julia<<"20\n"<<x[i][j][1]<<"\n";
				Julia<<"30\n"<<(double) mirror * x[i][j][2]<<"\n";
			}
			else
			{
				if(NowDrawing == 1){NowDrawing = 0;Julia<<"0\nSEQEND\n";}
			}
		}
		Julia<<"0\nSEQEND\n";
	}
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