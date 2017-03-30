//Written by Chris J K Williams, University of Bath, UK
//Modified for multi-threading Rob Hart

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
void Bending(int,int,int,int,int,int);
void Initialise(int);
void CalculateMotion(int);

void ThreadsInit(void);

void* ThreadFunction(void *threadid);
void ThreadsRun(void);

void* ThreadMotionFunction(void *threadid);
void ThreadsMotionRun(void);

#define   m 50
#define   n 50

#define LAST_THREAD    4


pthread_t CalcThread[LAST_THREAD + 1];
pthread_mutex_t OverlapMutex[LAST_THREAD + 1];
int Partition_iStart[LAST_THREAD + 1];
int Partition_iStop[LAST_THREAD + 1];

int NodeType[m + 1][n + 1];

float x[m + 1][n + 1][3],force[m + 1][n + 1][3],velocity[m + 1][n + 1][3],
Stiffness[m + 1][n + 1],
SphereCentre[3],PointPosition[2],
EAOverTwoSlackCubed,EAOverSlack,
slackSq,
eta,weightpernode,factor,carryover,
SphereRadius,SphereRadiusSq,eighteta,root2,DragConstant,floatNumberOfNodes;

char *MyText;

int main(int argc, char *  argv[])
{
	cout<<"The program will now run forever until 'q' or ESCAPE is pressed on the keyboard\n";
	cout<<"Left mouse button for pan\n";
	cout<<"Right mouse button (or CONTROL left button) for rotate\n";
	cout<<"Move spot to zoom\n";
	int mcentre = int(0.45 * m);
	int ncentre = int(0.48 * n);
	float slack = 800.0/(1.0 * m);
	slackSq = slack * slack;
	EAOverTwoSlackCubed = 1.0;
	EAOverSlack = EAOverTwoSlackCubed * 2.0 * slack * slack;
	eta = 0.5 * EAOverSlack;
	weightpernode = 0.002 * EAOverSlack * slack;
	DragConstant = 0.0001 * weightpernode;
	factor = 0.01;
	carryover = 0.999;
	eighteta = 8.0 * eta;
	SphereRadius = 100.0;SphereCentre[0] = 0.0;SphereCentre[1] = 0.0;SphereCentre[2] = 200.0;
	SphereRadiusSq = SphereRadius * SphereRadius;
	root2 = sqrt(2.0);
	MaxHalfDimension = 0.3 * slack * sqrt(float(m * m + n * n));
	floatNumberOfNodes = float((m + 1) * (n + 1));
	for(int i = 0;i <= m;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			x[i][j][0] = slack * float(i - mcentre);
			x[i][j][1] = slack * float(j - ncentre);
			x[i][j][2] = SphereCentre[2] + SphereRadius;
			NodeType[i][j] = 0;
		}
		Initialise(i);
	}
#if LAST_THREAD > 0
	ThreadsInit();
#endif
	glutInit(&argc,argv);
	MyGraphics(0.0,0.0,0.0);
	return 0;
}

static void Draw(void)
{
for(int xyz = 0; xyz <= 2;xyz ++)
{
float xSum = 0.0;
	for(int i = 0;i <= m;i ++)
	{
	for(int j = 0;j <= n;j ++)
	{
xSum += x[i][j][xyz];
	}
	}
	float xAverage = xSum /floatNumberOfNodes;
	for(int i = 0;i <= m;i ++)
	{
	for(int j = 0;j <= n;j ++)x[i][j][xyz] -= xAverage;
	}
	SphereCentre[xyz] -= xAverage;
	}
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	glLoadIdentity();
	
	glPointSize(10.0);
	if(ZoomButton == 1)glColor4f(1.0,0.0,0.0,0.5);else glColor4f(1.0,1.0,1.0,0.2);
	glBegin(GL_POINTS);
	PointPosition[0] = ZoomButtonX;
	PointPosition[1] = ZoomButtonY;
	glVertex2fv(PointPosition);
	glEnd();
	
	MyText="Zoom";
	glColor4f(1.0,1.0,1.0,0.2);
	glRasterPos2d(ZoomButtonX + 0.02,ZoomButtonY - 0.01);
	
	int StringLength=strlen(MyText);
	for(int MyChar=0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
	
	glPopMatrix();
	glColor4f(1.0,1.0,1.0,0.25);
	glLineWidth(2.0);
	for(int j = 0;j <= n;j ++)
	{
		glBegin(GL_LINE_STRIP);
		for(int i = 0;i <= m;i ++)glVertex3fv(x[i][j]);
		glEnd();
	}
	for(int i = 0;i <= m;i ++)
	{
		glBegin(GL_LINE_STRIP);
		for(int j = 0;j <= n;j ++)
			glVertex3fv(x[i][j]);
		glEnd();
	}
	glutSwapBuffers();
	for(int CalculationLoop = 0;CalculationLoop <= 30;CalculationLoop++)Calculation();
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
		force[i][j][0] = 0.0;force[i][j][1] = 0.0;force[i][j][2] =  - weightpernode;
	}
}

void CalculateForces(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		if(i < m)Axial(i,i + 1,j,j);
		if(j < n)Axial(i,i,j,j + 1);
		if(i < m - 1)Bending(i,i + 1,i + 2,j,j,j);
		if(j < n - 1)Bending(i,i,i,j,j + 1,j + 2);
	}
}

void CalculateMotion(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		if(NodeType[i][j] == 0)
		{
			float Normal[3];
			if(i > 0 && i < m && j > 0 && j < n)
			{
				float NormalDotVelocity = 0.0;
				for(int xyz = 0; xyz <= 2;xyz ++)
				{
					int xyzp1 = xyz + 1;if(xyzp1 > 2)xyzp1 -=3;
					int xyzp2 = xyz + 2;if(xyzp2 > 2)xyzp2 -=3;
					Normal[xyz]=(x[i + 1][j][xyzp1] - x[i - 1][j][xyzp1]) * (x[i][j + 1][xyzp2] - x[i][j - 1][xyzp2]);
					NormalDotVelocity += Normal[xyz] * velocity[i][j][xyz];
				}
				for(int xyz = 0; xyz <= 2;xyz ++)
				{
					force[i][j][xyz] -= DragConstant * Normal[xyz] * NormalDotVelocity * fabs(NormalDotVelocity);
				}
			}
			float Radius[3];
			for(int xyz = 0;xyz <= 2;xyz ++)
			{
				velocity[i][j][xyz] = carryover * velocity[i][j][xyz] + factor * force[i][j][xyz]/(Stiffness[i][j] + eighteta);
				x[i][j][xyz] = x[i][j][xyz] + velocity[i][j][xyz];
				Radius[xyz] = x[i][j][xyz] - SphereCentre[xyz];
			}
			
			float RadSq = Radius[0] * Radius[0] + Radius[1] * Radius[1] + Radius[2] * Radius[2];
			if(RadSq<SphereRadiusSq)
			{
				float dot = velocity[i][j][0] * Radius[0] + velocity[i][j][1] * Radius[1] + velocity[i][j][2] * Radius[2];
				float Rad = sqrt(RadSq);
				for(int xyz = 0;xyz <= 2;xyz ++)
				{
					velocity[i][j][xyz] -= Radius[xyz] * dot/RadSq;
					x[i][j][xyz] += (SphereRadius - Rad) * Radius[xyz]/Rad;
				}
			}
		}
	}
	Initialise(i);
}

void Axial(int i,int inext,int j,int jnext)
{
	float deltax[3];
	for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[inext][jnext][xyz] - x[i][j][xyz];
	float lengthSq = deltax[0] * deltax[0] + deltax[1] * deltax[1] + deltax[2] * deltax[2];
	float tensioncoefficient = EAOverTwoSlackCubed * (lengthSq - slackSq);
	float thisStiffness = EAOverTwoSlackCubed * (3.0 * lengthSq - 2.0 * slackSq);
	Stiffness[i][j] += thisStiffness;
	Stiffness[inext][jnext] += thisStiffness;
	for(int xyz = 0;xyz <= 2;xyz ++)
	{
		float thisforce = tensioncoefficient * deltax[xyz];
		force[i][j][xyz] += thisforce;
		force[inext][jnext][xyz] -= thisforce;
	}
}

void Bending(int i,int inext,int inextnext,int j,int jnext,int jnextnext)
{
	for(int xyz = 0;xyz <= 2;xyz ++)
	{
		float deltaA = x[inext][jnext][xyz] - x[i][j][xyz];
		float deltaB = x[inextnext][jnextnext][xyz] - x[inext][jnext][xyz];
		float thisforce = eta * (deltaA - deltaB);
		force[i][j][xyz] += thisforce;
		force[inext][jnext][xyz] -= 2.0 * thisforce;
		force[inextnext][jnextnext][xyz] += thisforce;
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
    int id = (int)threadid;
	
	for(int i = Partition_iStart[id];  i <= Partition_iStop[id]; i++)
	{
		if(i == Partition_iStart[id]    && id > 0          )pthread_mutex_lock(&OverlapMutex[id-1]);
		if(i == Partition_iStop[id] - 2 && id < LAST_THREAD)pthread_mutex_lock(&OverlapMutex[id+1]);
		
        CalculateForces(i);
		
		if(i == Partition_iStart[id] + 2 && id > 0          )pthread_mutex_unlock(&OverlapMutex[id-1]);
		if(i == Partition_iStop[id]      && id < LAST_THREAD)pthread_mutex_unlock(&OverlapMutex[id+1]);
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
    int id = (int)threadid;
	
	for(int i = Partition_iStart[id];  i <= Partition_iStop[id]; i++)CalculateMotion(i);
	
	pthread_exit(NULL);   
}