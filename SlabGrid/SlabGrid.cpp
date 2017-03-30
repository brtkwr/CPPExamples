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
int LineExists(int,int,int,int);
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
#define  Fineness 1
#define  m (13 * Fineness * Spacing)
#define  nOver2 (6 * Fineness * Spacing)
#define  n (2 * nOver2)
#define  p (2 * Fineness * Spacing + HalfSpacing)

#define  CarryOver 0.99
#define  MovementFactor 0.2
#define LAST_THREAD  1

pthread_t CalcThread[LAST_THREAD + 1];
pthread_mutex_t OverlapMutex[LAST_THREAD + 1];
int Partition_iStart[LAST_THREAD + 1];
int Partition_iStop[LAST_THREAD + 1];

int  NodeType[m + 1][n + 1];

double x[m + 1][n + 1][3],CableNetForce[m + 1][n + 1][3],
velocity[m + 1][n + 1][3],
PointPosition[2];

char NumericalValue[101];
string MyText;

int main(int argc, char * argv[])
{
	cout<<"m = "<<m<<"\n";
	cout<<"n = "<<n<<"\n";
	cout<<"The program will now run forever until 'q' or ESCAPE is pressed on the keyboard\n";
	cout<<"Left mouse button for pan\n";
	cout<<"Right mouse button (or CONTROL left button) for rotate\n";
	cout<<"Move spot to zoom\n";
	cout<<"The space bar starts and stops the analysis\n";
	
	for(int i = 0;i <= m;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			NodeType[i][j] = - 1;
			
			if(j == nOver2 - (i - m))NodeType[i][j] = 1;
			if(j == nOver2 + (i - m))NodeType[i][j] = 2;
			if( i == m && j == nOver2)NodeType[i][j] = 3;
			
			if(i == m - p && j > nOver2)NodeType[i][j] = 1;
			if(i == m - p && j < nOver2)NodeType[i][j] = 2;
			
			if((j < nOver2 - (i - m)) && (j > nOver2 + (i - m)))NodeType[i][j] = 0;
			
			if(i < m - p)NodeType[i][j] = 0;
			
			if((j == 0) && (i < m - p))NodeType[i][j] = 1;
			if((j == n) && (i < m - p))NodeType[i][j] = 2;
			
			if((j == 0) && (i == m - p))NodeType[i][j] = 3;
			if((j == n) && (i == m - p))NodeType[i][j] = 3;
			
			if(i == 0)NodeType[i][j] = 4;
		}
	}
	
	a = 10.0;
	
	for(int i = 0;i <= m;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			x[i][j][2] = 0.0;
			
			x[i][j][0] = - a / 2.0 + ((double) (i - j + nOver2) / (double) m) * a / 4.0;
			x[i][j][1] = - a / 2.0 + ((double) (i + j - nOver2) / (double) m) * a / 4.0;
			
			if(i == 0 || j == n)x[i][j][0] = - a;
			if(i == 0 || j == 0)x[i][j][1] = - a;
			if(NodeType[i][j] == 2 && j != n)x[i][j][0] = 0.0;
			if(NodeType[i][j] == 1 && j != 0)x[i][j][1] = 0.0;
			
			if((j == 0) && (i == m - p))x[i][j][0] = 0.0;
			if((j == n) && (i == m - p))x[i][j][1] = 0.0;
			
			if( i == m && j == nOver2){x[i][j][0] = 0.0;x[i][j][1] = 0.0;}
		}
	}
	
	StopStart = 0;
	
	for(int i = 0;i <= m;i ++)Initialise(i);
	
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
		
		glRasterPos2d(ButtonX[Button] + 0.02,ButtonY[Button] - 0.01);
		
		int StringLength = MyText.length();
		for(int MyChar = 0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
	}
	
	glPopMatrix();
	
	
	for(int shiftX = -1;shiftX <= 1; shiftX += 2)
	{
		for(int shiftY = -1;shiftY <= 1; shiftY += 2)
		{
			for(int MirrorX = -1;MirrorX <= 1; MirrorX += 2)
	{
		for(int MirrorY = -1;MirrorY <= 1; MirrorY += 2)
		{
			glLineWidth(1.0);
			
			int LineDrawStep = Spacing - 1;
			
			for(int j = 0;j <= n;j ++)
			{
				LineDrawStep ++;
				
				if(LineDrawStep >= Spacing)
				{
					LineDrawStep = 0;
					bold();
				}
				else faint();
				glBegin(GL_LINES);
				for(int i = 0;i <= m - 1;i ++)
				{
					if(LineExists(i,i + 1,j,j) >= 1 && (j != n || MirrorX == - 1 || shiftX == - 1) && (j != 0 || MirrorY == - 1 || shiftY == - 1))
					{
						glVertex3d((double)shiftX * a + (double)MirrorX * x[i + 0][j + 0][0],(double)shiftY * a + (double)MirrorY * x[i + 0][j + 0][1],x[i + 0][j + 0][2]);
						glVertex3d((double)shiftX * a + (double)MirrorX * x[i + 1][j + 0][0],(double)shiftY * a + (double)MirrorY * x[i + 1][j + 0][1],x[i + 1][j + 0][2]);
					}
				}
				glEnd();
			}
			
			LineDrawStep = Spacing - 1;
			for(int i = 0;i <= m;i ++)
			{
				LineDrawStep ++;
				
				if(LineDrawStep >= Spacing)
				{
					LineDrawStep = 0;
					bold();
				}
				else faint();
				glBegin(GL_LINES);
				for(int j = 0;j <= n - 1;j ++)
				{
					if(LineExists(i,i,j,j + 1) >= 1)
					{
					if(LineExists(i,i,j,j + 1) == 1 || (MirrorX == 1 && MirrorY == 1) || (MirrorX == - 1 && MirrorY == - 1))
					{
						glVertex3d((double)shiftX * a + (double)MirrorX * x[i + 0][j + 0][0],(double)shiftY * a + (double)MirrorY * x[i + 0][j + 0][1],x[i + 0][j + 0][2]);
						glVertex3d((double)shiftX * a + (double)MirrorX * x[i + 0][j + 1][0],(double)shiftY * a + (double)MirrorY * x[i + 0][j + 1][1],x[i + 0][j + 1][2]);
					}
					}
				}
				glEnd();
			}
			
			int ShowNodes = 0;
			if(ShowNodes == 1)
			{
			glPointSize(2.0);
			glBegin(GL_POINTS);
			for(int i = 0;i <= m;i ++)
			{
				for(int j = 0;j <= n;j ++)
				{
					if(NodeType[i][j] >= 0)
					{
						if(NodeType[i][j] == 0)glColor4f(0.0,0.0,0.0,0.5);
						if(NodeType[i][j] == 1)glColor4f(1.0,0.0,0.0,1.0);
						if(NodeType[i][j] == 2)glColor4f(0.0,1.0,0.0,1.0);
						if(NodeType[i][j] == 3)glColor4f(0.0,0.0,1.0,1.0);
						if(NodeType[i][j] == 4)glColor4f(0.0,1.0,1.0,1.0);
						
						glVertex3d((double)shiftX * a + (double)MirrorX * x[i][j][0],(double)shiftY * a + (double)MirrorY * x[i][j][1],x[i][j][2]);
					}
				}
			}
			glEnd();
		}
	}
	}
		}
	}
	
	glutSwapBuffers();
	if(StopStart == 1)
	{
		for(int CalculationLoop = 0;CalculationLoop <= 0;CalculationLoop++)
		{
			Calculation();
		}
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
		for(int xyz = 0; xyz <= 2;xyz ++)CableNetForce[i][j][xyz] = 0.0;
	}
}

void CalculateForces(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		if(i < m)Axial(i,i + 1,j,j);
		if(j < n)Axial(i,i,j,j + 1);
	}
}

void CalculateMotion(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		if(NodeType[i][j] >= 0)
		{
			for(int xyz = 0;xyz <= 2;xyz ++)
			{
				if(NodeType[i][j] == 0 || (NodeType[i][j] == 1 && xyz != 1) || (NodeType[i][j] == 2 && xyz != 0) || (NodeType[i][j] == 3 && xyz == 2))
				{
					velocity[i][j][xyz] = CarryOver * velocity[i][j][xyz] + MovementFactor * CableNetForce[i][j][xyz];
					x[i][j][xyz] += velocity[i][j][xyz];
				}
			}
		}
	}
	Initialise(i);
}

int LineExists(int i,int inext,int j,int jnext)
{
	if(NodeType[i][j] == - 1 || NodeType[inext][jnext] == - 1 || (NodeType[i][j] == 3 && NodeType[inext][jnext] == 3))return 0;
	else
	{
if(NodeType[i][j] > 0 && NodeType[inext][jnext] > 0)return 2;
	   else return 1;
	}
}

void Axial(int i_trial,int inext_trial,int j,int jnext)
{
	double deltax[3];
	
	int i = i_trial;
	if(i < 0)i= m;
	if(i > m)i = 0;
	
	int inext = inext_trial;
	if(inext < 0)inext = m;
	if(inext > m)inext = 0;
	
	if(LineExists(i,inext,j,jnext) >= 1)
	{
		for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[inext][jnext][xyz] - x[i][j][xyz];
		
		double tensioncoefficient = 1.0;
		
		if(LineExists(i,inext,j,jnext) == 2)tensioncoefficient /= 2.0;
		
		for(int xyz = 0;xyz <= 2;xyz ++)
		{
			double thisForce = tensioncoefficient * deltax[xyz];
			CableNetForce[i][j][xyz] += thisForce;
			CableNetForce[inext][jnext][xyz] -= thisForce;
		}
		double ThisLoad = 0.0;
		for(int xy = 0;xy <= 1;xy ++)ThisLoad += 0.01 * deltax[xy] * deltax[xy];
		if(LineExists(i,inext,j,jnext) == 2)ThisLoad /= 2.0;
		CableNetForce[i][j][2] -= ThisLoad;
		CableNetForce[inext][jnext][2] -= ThisLoad;
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

void writeDXF(void)
{
	ofstream Julia("SlabRibs.dxf");
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	
	int BoldLines;
	for(int MirrorX = -1;MirrorX <= 1; MirrorX += 2)
	{
		for(int MirrorY = -1;MirrorY <= 1; MirrorY += 2)
		{
			int LineDrawStep = Spacing - 1;
			
			for(int j = 0;j <= n;j ++)
			{
				LineDrawStep ++;
				
				if(LineDrawStep >= Spacing)
				{
					LineDrawStep = 0;
					BoldLines = 1;
				}
				else BoldLines = 0;
				for(int i = 0;i <= m - 1;i ++)
				{
					if(LineExists(i,i + 1,j,j) >= 1)
					{
						if(BoldLines == 1)
							Julia<<"0\nLINE\n8\nMainRibs\n";
						else
							Julia<<"0\nLINE\n8\nSubsiduaryRibs\n";
						Julia<<"10\n"<<(double)MirrorX * x[i + 0][j + 0][0]<<"\n";
						Julia<<"20\n"<<(double)MirrorY * x[i + 0][j + 0][1]<<"\n";
						Julia<<"30\n"<<                  x[i + 0][j + 0][2]<<"\n";
						Julia<<"11\n"<<(double)MirrorX * x[i + 1][j + 0][0]<<"\n";
						Julia<<"21\n"<<(double)MirrorY * x[i + 1][j + 0][1]<<"\n";
						Julia<<"31\n"<<                  x[i + 1][j + 0][2]<<"\n";	
						if(BoldLines == 1)
							Julia<<"62\n0\n";
						else
							Julia<<"62\n1\n";
					}
				}
			}
			
			LineDrawStep = Spacing - 1;
			for(int i = 0;i <= m;i ++)
			{
				LineDrawStep ++;
				
				if(LineDrawStep >= Spacing)
				{
					LineDrawStep = 0;
					BoldLines = 1;
				}
				else BoldLines = 0;
				for(int j = 0;j <= n - 1;j ++)
				{
					if(LineExists(i,i,j,j + 1) >= 1)
					{
						if(LineExists(i,i,j,j + 1) == 1 || (MirrorX == 1 && MirrorY == 1) || (MirrorX == - 1 && MirrorY == - 1))
						{
							if(BoldLines == 1)
								Julia<<"0\nLINE\n8\nMainRibs\n";
							else
								Julia<<"0\nLINE\n8\nSubsiduaryRibs\n";							
							Julia<<"10\n"<<(double)MirrorX * x[i + 0][j + 0][0]<<"\n";
							Julia<<"20\n"<<(double)MirrorY * x[i + 0][j + 0][1]<<"\n";
							Julia<<"30\n"<<                  x[i + 0][j + 0][2]<<"\n";
							Julia<<"11\n"<<(double)MirrorX * x[i + 0][j + 1][0]<<"\n";
							Julia<<"21\n"<<(double)MirrorY * x[i + 0][j + 1][1]<<"\n";
							Julia<<"31\n"<<                  x[i + 0][j + 1][2]<<"\n";									}
						if(BoldLines == 1)
							Julia<<"62\n0\n";
						else
							Julia<<"62\n1\n";
					}
				}
			}
		}
	}
			
	Julia<<"0\nENDSEC\n0\nEOF\n\n";
	Julia.close();
	cout<<"DXF file written, end of program\n\n";
}

void faint(void)
{
	glColor4f(0.0,0.0,1.0,0.25);
}

void bold(void)
{
	glColor4f(0.0,0.0,0.0,0.5);
	}