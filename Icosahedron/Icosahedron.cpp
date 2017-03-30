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
int LineExists(int,int,int,int);

#define  Spacing 20
#define  InnerSpacing 2
#define  n_over_6 (Spacing * InnerSpacing)
#define  n        (6 * n_over_6)
#define  CarryOver 0.995
#define  MovementFactor 0.1

#define LAST_THREAD  1

pthread_t CalcThread[LAST_THREAD + 1];
pthread_mutex_t OverlapMutex[LAST_THREAD + 1];
int Partition_iStart[LAST_THREAD + 1];
int Partition_iStop[LAST_THREAD + 1];

int  NodeType[n + 1][n + 1],Sister_i[n + 1][n + 1],Sister_j[n + 1][n + 1],DoNotDrawFaint;

double x[n + 1][n + 1][3],CableNetForce[n + 1][n + 1][3],
velocity[n + 1][n + 1][3],
PointPosition[2],
CableNetTensionCoefficient;

char NumericalValue[101];
string MyText;

int main(int argc, char * argv[])
{
	PI = 4.0 * atan(1.0);
	
	cout<<"The program will now run forever until 'q' or ESCAPE is pressed on the keyboard\n";
	cout<<"Left mouse button for pan\n";
	cout<<"Right mouse button (or CONTROL left button) for rotate\n";
	cout<<"Move spot to zoom\n";
	cout<<"The space bar starts and stops the analysis\n";
	
	CableNetTensionCoefficient = 0.1;
	
	a = 5000.0;
	double b = 2.0 * a;
	
	for(int i = 0;i <= n;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			NodeType[i][j] = 0;
			Sister_i[i][j] = i;
			Sister_j[i][j] = j;
			if(
			   (i < 1 * n_over_6 && j < 4 * n_over_6) ||
			   (i < 2 * n_over_6 && j < 3 * n_over_6) ||
			   (i < 3 * n_over_6 && j < 2 * n_over_6) ||
			   (i < 4 * n_over_6 && j < 1 * n_over_6) ||
			   (i < 5 * n_over_6 && j < 0 * n_over_6) ||
			   (i > 1 * n_over_6 && j > 5 * n_over_6) ||
			   (i > 2 * n_over_6 && j > 4 * n_over_6) ||
			   (i > 3 * n_over_6 && j > 3 * n_over_6) ||
			   (i > 4 * n_over_6 && j > 2 * n_over_6) ||
			   (i > 5 * n_over_6))
				NodeType[i][j] = - 1;
		}
	}
	NodeType[0][0] = 0;
	NodeType[5 * n_over_6][n] = 0;
	
	for(int Notch = 0;Notch <= 4;Notch ++)
	{
		for(int Notch_i = 1;Notch_i <= n_over_6;Notch_i ++)
		{
			int i = Notch * n_over_6 + Notch_i;
			int j = (6 - Notch) * n_over_6;
			NodeType[i][j] = 1;
			NodeType[5 * n_over_6 - i][n - j] = 1;
			if(Notch > 0)
			{
				int Sis_i = Notch * n_over_6;
				int Sis_j = (6 - Notch) * n_over_6 + Notch_i;
				
				Sister_i[i][j] = Sis_i;
				Sister_j[i][j] = Sis_j;
				
				Sister_i[5 * n_over_6 - i][n - j] = 5 * n_over_6 - Sis_i;
				Sister_j[5 * n_over_6 - i][n - j] = n - Sis_j;
			}
			else
			{
				Sister_i[i][j] = 5 * n_over_6;
				Sister_j[i][j] = n_over_6 + Notch_i;
				
				Sister_i[5 * n_over_6 - i][n - j] = 0;
				Sister_j[5 * n_over_6 - i][n - j] = 5 * n_over_6 - Notch_i;
			}
		}
	}	
	
	for(int j = 5 * n_over_6;j <= n; j ++)
	{
		NodeType[0][j] = 1;
		Sister_i[0][j] = 5 * n_over_6;
		Sister_j[0][j] = j - 5 * n_over_6;
	}
	
	for(int i = n_over_6;i <= n;i ++)
	{
		int j = 7 * n_over_6 - i;
		if(NodeType[i][j] >= 0)
		{
			NodeType[i][j] = 1;
			Sister_i[i][j] = n;
			Sister_j[i][j] = n;
		}
	}
	
	for(int i = 0;i <= 5 * n_over_6;i ++)
	{
		int j = 4 * n_over_6 - i;
		if(NodeType[i][j] >= 0)
		{
			NodeType[i][j] = 1;
			Sister_i[i][j] = 0;
			Sister_j[i][j] = 0;
		}
	}
	
	for(int i = 0;i <= n;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			double theta = PI * (double) (i - j) / (double) (5 * n_over_6);
			x[i][j][0] = (b / (2.0 * PI)) * cos(theta);
			x[i][j][1] = (b / (2.0 * PI)) * sin(theta);
			x[i][j][2] = b * (2.0 / sqrt(3.0)) * (double) (2 * (i + j) - 11 * n_over_6) / (double) (2 * n);
			
			for(int xyz = 0; xyz <= 2;xyz ++)velocity[i][j][xyz] = 0.0;
		}
	}
	
	NodeType[0][0] = 0;
	x[0][0][0] = 0.0;
	x[0][0][1] = 0.0;
	x[0][0][2] = - b / 3.0;
	
	NodeType[n][n] = 0;
	x[n][n][0] = 0.0;
	x[n][n][1] = 0.0;
	x[n][n][2] = - x[0][0][2];
	
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
	
	int LineDrawStep = Spacing - 1;
	for(int j_value = 4 * n_over_6;j_value <= 7 * n_over_6;j_value += InnerSpacing)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			bold();
		}
		else faint();
		
		for(int i = 0;i <= n - 1;i ++)
		{
			int j = j_value - i;
			if(j >= 1 && j <= n)
			{
				if(LineExists(i,j,i + 1,j - 1))
				{
					glBegin(GL_LINES);
					glVertex3dv(x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]]);
					glVertex3dv(x[Sister_i[i + 1][j - 1]][Sister_j[i + 1][j - 1]]);
					glEnd();
				}
			}
		}
	}
	
	LineDrawStep = - 1;
	for(int j = 0;j <= n;j += InnerSpacing)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			bold();
		}
		else faint();
		
		for(int i = 0;i <= n - 1;i ++)
		{
			if(LineExists(i,j,i + 1,j + 0))
			{
				glBegin(GL_LINES);
				glVertex3dv(x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]]);
				glVertex3dv(x[Sister_i[i + 1][j + 0]][Sister_j[i + 1][j + 0]]);
				glEnd();
			}
		}
	}
	
	LineDrawStep = Spacing - 1;
	for(int i = 0;i <= n;i += InnerSpacing)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			bold();
		}
		else faint();
		
		for(int j = 0;j <= n - 1;j ++)
		{
			if(LineExists(i,j,i + 0,j + 1))
			{
				glBegin(GL_LINES);
				glVertex3dv(x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]]);
				glVertex3dv(x[Sister_i[i + 0][j + 1]][Sister_j[i + 0][j + 1]]);
				glEnd();
			}
		}
	}
	
	glutSwapBuffers();
	if(StopStart == 1)
	{
		for(int CalculationLoop = 0;CalculationLoop <= 0;CalculationLoop++)
		{
			Calculation();
			double lambda = 2.0;
			
			velocity[1 * n_over_6][4 * n_over_6][0] -= lambda * x[1 * n_over_6][4 * n_over_6][0];
			velocity[2 * n_over_6][4 * n_over_6][0] -= lambda * x[2 * n_over_6][4 * n_over_6][0];
			velocity[4 * n_over_6][1 * n_over_6][0] -= lambda * x[4 * n_over_6][1 * n_over_6][0];
			velocity[4 * n_over_6][2 * n_over_6][0] -= lambda * x[4 * n_over_6][2 * n_over_6][0];
			
			double average = (
							  x[1 * n_over_6][4 * n_over_6][1] + 
							  x[2 * n_over_6][4 * n_over_6][1]) / 2.0;
			velocity[1 * n_over_6][4 * n_over_6][1] += lambda * (average - x[1 * n_over_6][4 * n_over_6][1]);
			velocity[2 * n_over_6][4 * n_over_6][1] += lambda * (average - x[2 * n_over_6][4 * n_over_6][1]);
			
			average = (
					   x[4 * n_over_6][1 * n_over_6][1] + 
					   x[4 * n_over_6][2 * n_over_6][1]) / 2.0;
			velocity[4 * n_over_6][1 * n_over_6][1] += lambda * (average - x[4 * n_over_6][1 * n_over_6][1]);
			velocity[4 * n_over_6][2 * n_over_6][1] += lambda * (average - x[4 * n_over_6][2 * n_over_6][1]);	
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
		if(i + 1 <= n && j - 1 >= 0)Axial(i,i + 1,j,j - 1);
	}
}

void CalculateMotion(int i)
{
	for(int j = 0;j <= n;j ++)
	{
		if(NodeType[i][j] == 0)
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
			
			for(int xyz = 0;xyz <= 2;xyz ++)
			{
				velocity[i][j][xyz] = CarryOver * velocity[i][j][xyz] + MovementFactor * CableNetForce[i][j][xyz];
				x[i][j][xyz] += velocity[i][j][xyz];
			}
		}
	}
	
	Initialise(i);
}

void Axial(int original_i,int original_inext,int original_j,int original_jnext)
{
	double deltax[3];
	
	if(LineExists(original_i,original_j,original_inext,original_jnext))
	{
		int i = Sister_i[original_i][original_j];
		int j = Sister_j[original_i][original_j];
		int inext = Sister_i[original_inext][original_jnext];
		int jnext = Sister_j[original_inext][original_jnext];
		
		for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[inext][jnext][xyz] - x[i][j][xyz];
		
		double tensioncoefficient = CableNetTensionCoefficient;
		
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
	ofstream Julia("Icosahedron.dxf");
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	
	int NowBold = 1;
	int LineDrawStep = Spacing - 1;
	for(int j_value = 4 * n_over_6;j_value <= 7 * n_over_6;j_value += InnerSpacing)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			NowBold = 1;
		}
		else NowBold = 0;
		
		for(int i = 0;i <= n - 1;i ++)
		{
			int j = j_value - i;
			if(j >= 1 && j <= n)
			{
				if(LineExists(i,j,i + 1,j - 1))
				{
					if(NowBold == 1)Julia<<"0\nLINE\n8\nBold\n";else Julia<<"0\nLINE\n8\nFaint\n";
					Julia<<"10\n"<<x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]][0]<<"\n";
					Julia<<"20\n"<<x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]][1]<<"\n";
					Julia<<"30\n"<<x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]][2]<<"\n";
					Julia<<"11\n"<<x[Sister_i[i + 1][j - 1]][Sister_j[i + 1][j - 1]][0]<<"\n";
					Julia<<"21\n"<<x[Sister_i[i + 1][j - 1]][Sister_j[i + 1][j - 1]][1]<<"\n";
					Julia<<"31\n"<<x[Sister_i[i + 1][j - 1]][Sister_j[i + 1][j - 1]][2]<<"\n";
					if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";
				}
			}
		}
	}
	
	LineDrawStep = - 1;
	for(int j = 0;j <= n;j += InnerSpacing)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			NowBold = 1;
		}
		else NowBold = 0;
		
		for(int i = 0;i <= n - 1;i ++)
		{
			if(LineExists(i,j,i + 1,j + 0))
			{
				if(NowBold == 1)Julia<<"0\nLINE\n8\nBold\n";else Julia<<"0\nLINE\n8\nFaint\n";
				Julia<<"10\n"<<x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]][0]<<"\n";
				Julia<<"20\n"<<x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]][1]<<"\n";
				Julia<<"30\n"<<x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]][2]<<"\n";
				Julia<<"11\n"<<x[Sister_i[i + 1][j + 0]][Sister_j[i + 1][j + 0]][0]<<"\n";
				Julia<<"21\n"<<x[Sister_i[i + 1][j + 0]][Sister_j[i + 1][j + 0]][1]<<"\n";
				Julia<<"31\n"<<x[Sister_i[i + 1][j + 0]][Sister_j[i + 1][j + 0]][2]<<"\n";
				if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";
			}
		}
	}
	
	LineDrawStep = Spacing - 1;
	for(int i = 0;i <= n;i += InnerSpacing)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			NowBold = 1;
		}
		else NowBold = 0;
		
		for(int j = 0;j <= n - 1;j ++)
		{
			if(LineExists(i,j,i + 0,j + 1))
			{
				if(NowBold == 1)Julia<<"0\nLINE\n8\nBold\n";else Julia<<"0\nLINE\n8\nFaint\n";
				Julia<<"10\n"<<x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]][0]<<"\n";
				Julia<<"20\n"<<x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]][1]<<"\n";
				Julia<<"30\n"<<x[Sister_i[i + 0][j + 0]][Sister_j[i + 0][j + 0]][2]<<"\n";
				Julia<<"11\n"<<x[Sister_i[i + 0][j + 1]][Sister_j[i + 0][j + 1]][0]<<"\n";
				Julia<<"21\n"<<x[Sister_i[i + 0][j + 1]][Sister_j[i + 0][j + 1]][1]<<"\n";
				Julia<<"31\n"<<x[Sister_i[i + 0][j + 1]][Sister_j[i + 0][j + 1]][2]<<"\n";
				if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";	
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
	DoNotDrawFaint = 0;
}

void bold(void)
{
	glColor4f(0.0,0.0,0.0,0.5);
	DoNotDrawFaint = 0;
}

int LineExists(int i,int j,int next_i,int next_j)
{
	if((i == next_i && ((NodeType[i][j] == 0 && NodeType[next_i][next_j] >= 0) || (NodeType[i][j] >= 0 && NodeType[next_i][next_j] == 0)))
	   || (j == next_j && ((i ==0 || NodeType[i][j] == 0) && (next_i == 0 || NodeType[next_i][next_j] == 0))
		   || (i != next_i && j != next_j) && 
		   (NodeType[i][j] >= 0 && NodeType[next_i][next_j] >= 0 && NodeType[next_i][j] >= 0 && NodeType[i][next_j] >= 0)
		   )) return 1;
	else return 0;
}