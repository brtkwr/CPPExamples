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
void Axial(int,int,int,int,int);
void Initialise(int);
void CalculateMotion(int);

void ThreadsInit(void);

void* ThreadFunction(void *threadid);
void ThreadsRun(void);

void* ThreadMotionFunction(void *threadid);
void ThreadsMotionRun(void);

void faint(void);
void bold(void);

#define  ThirdCable 0

#define  Spacing 9
#define  Fineness 6
#define  eighth    (3 * Fineness * Spacing)
#define  n ((5 - 2 * ThirdCable) * Fineness * Spacing)
#define  nAll (eighth - Fineness * Spacing)

#define  m (8 * eighth - 1)

#define  CarryOver 0.99
#define  MovementFactor 1.0
#define  CircleRadiusFactor 0.005
#define  CheckGeometry 0

#define LAST_THREAD  1

pthread_t CalcThread[LAST_THREAD + 1];
pthread_mutex_t OverlapMutex[LAST_THREAD + 1];
int Partition_iStart[LAST_THREAD + 1];
int Partition_iStop[LAST_THREAD + 1];

int  NodeType[m + 1][n + 1],DoNotDrawFaint;

double x[m + 1][n + 1][3],CableNetForce[m + 1][n + 1][3],
velocity[m + 1][n + 1][3],CableLength[m + 1][n + 1][3],
PointPosition[2],totalLength[3],
CableNetTensionCoefficient,NetCableRadius;

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
	
	CableNetTensionCoefficient = 0.1;
	
	if(ThirdCable == 1)b = 1750.0;else b = 1950.0;
	a = 5000.0;
	NetCableRadius = 2.0;
	HalfHeightDifference = 0.0;
	
	double setupmultiplier = 8.0 * atan(1.0) / (double)(m + 1);
	
	int NodeRemove = 1;
	for(int i = 0;i <= m;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			NodeRemove = - NodeRemove;
			if(NodeRemove == 1)NodeType[i][j] = - 1;
			else
			{
				NodeType[i][j] = 0;
				if(j == n)NodeType[i][j] = 1;
				else
				{
					for(int WhichQuarter = 0;WhichQuarter <= 3; WhichQuarter ++)
					{
						int  first_i_control = 2 * (WhichQuarter + 0) * eighth;
						int second_i_control = 2 * (WhichQuarter + 1) * eighth;
						
						if(i >= first_i_control && i <= second_i_control)
						{
							if(j <= i - first_i_control && j <= second_i_control - i && j <= nAll)
							{
								if(j == i - first_i_control || j == second_i_control - i || j == nAll)
								{
									if(NodeType[i][j] != 0)NodeType[i][j] += 1;
									if(i == 0)NodeType[i][j] += 10;
									NodeType[i][j] += WhichQuarter + 2;
								}
								else NodeType[i][j] = - 1;
							}
						}
					}
				}
			}
		}
	}
	
	for(int i = 0;i <= m;i ++)
	{
		double theta = atan(1.0) + (double)(2 * i + 1) * setupmultiplier / 2.0;
		for(int j = 0;j <= n;j ++)
		{
			double phi = (double) (n - j) * setupmultiplier;
			if(ThirdCable == 1)phi *= sqrt(3.0);
			
			double radius = b * exp(phi);
			x[i][j][0] = radius * cos(theta);
			x[i][j][1] = radius * sin(theta);
			x[i][j][2] = SupportHeight * (double) j / (double) (n + 1);
			
			for(int xyz = 0; xyz <= 2;xyz ++)velocity[i][j][xyz] = 0.0;
			for(int CableDirection = 0; CableDirection <= 2;CableDirection ++)CableLength[i][j][CableDirection] = 0.0;
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
	if(CheckGeometry == 0 && StopStart == 1)
	{
		glRasterPos2d(ButtonX[0],0.65);
		MyText="Ratio of total cable lengths = ";
		int LengthOfString = MyText.length();
		for(int MyCharacter = 0; MyCharacter < LengthOfString; MyCharacter++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyCharacter]);
		if(ThirdCable == 1)
			sprintf(NumericalValue,"%8.6f",2.0 * totalLength[0] / (totalLength[1] + totalLength[2]));
		else
			sprintf(NumericalValue,"%8.6f",sqrt(2.0) * totalLength[0] / (totalLength[1] + totalLength[2]));
		LengthOfString=strlen(NumericalValue);
		for(int MyCharacter = 0; MyCharacter < LengthOfString; MyCharacter++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	}
	
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
	int i_LineBegin = 1;
	
	if(ThirdCable == 1)
	{
	for(int j = 0;j <= n;j ++)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			bold();
		}
		else faint();
		
		if(i_LineBegin == 0)i_LineBegin = 1;else i_LineBegin = 0;
		int NowDrawing = 0;
		int PreviousExisted = 0;
		{
			int nexti;
			for(int i = i_LineBegin;i <= m;i += 2)
			{
				nexti = i + 2;
				if(nexti > m)nexti -= m + 1;
				if(LineExists(i,nexti,j,j) && (DoNotDrawFaint == 0 || CheckGeometry == 1))
				{
					if(NowDrawing == 0){NowDrawing = 1;glBegin(GL_LINE_STRIP);}
					glVertex3dv(x[i][j]);
					PreviousExisted = 1;
				}
				else
				{
					if(NowDrawing == 1)
					{
						if(PreviousExisted = 1)glVertex3dv(x[i][j]);
						NowDrawing = 0;
						glEnd();
					}
					PreviousExisted = 0;
				}
			}
			if(NowDrawing == 1)
			{
				if(PreviousExisted = 1)glVertex3dv(x[nexti][j]);
				glEnd();
			}
		}
		}
	}
	
	LineDrawStep = Spacing - 1;
	for(int i = 0;i <= m;i += 2)
	{
		LineDrawStep ++;
		
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			bold();
		}
		else faint();
		glBegin(GL_LINE_STRIP);
		for(int j = 0;j <= n - 1;j ++)
		{
			int thisi = i + j;while(thisi > m)thisi -= m + 1;
			int nexti = thisi + 1;while(nexti > m)nexti -= m + 1;
			if(LineExists(thisi,nexti,j,j + 1) && (DoNotDrawFaint == 0 || CheckGeometry == 1))glVertex3dv(x[thisi][j]);
		}
		{
			int thisi = i + n;while(thisi > m)thisi -= m + 1;
			glVertex3dv(x[thisi][n]);
		}
		glEnd();
	}
	
	LineDrawStep = Spacing - 1;
	for(int i = 0;i <= m;i += 2)
	{
		LineDrawStep ++;
		
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			bold();
		}
		else faint();
		glBegin(GL_LINE_STRIP);
		for(int j = 0;j <= n - 1;j ++)
		{
			int thisi = i - j;while(thisi < 0)thisi += m + 1;
			int nexti = thisi - 1;while(nexti < 0)nexti += m + 1;
			if(LineExists(thisi,nexti,j,j + 1)  && (DoNotDrawFaint == 0 || CheckGeometry == 1))glVertex3dv(x[thisi][j]);
		}
		{
			int thisi = i - n;while(thisi < 0)thisi += m + 1;
			glVertex3dv(x[thisi][n]);
		}
		glEnd();
	}
	
	if(CheckGeometry == 1)
	{
		glPointSize(2.0);
		glBegin(GL_POINTS);
		for(int i = 0;i <= m;i ++)
		{
			for(int j = 0;j <= n;j ++)
			{
				if(NodeType[i][j] == 0)glColor4f(0.0,0.0,0.0,0.5);
				else
				{
					if(NodeType[i][j] > 0)
					{
						if(NodeType[i][j] == 1 || NodeType[i][j] == 2 || NodeType[i][j] == 3 || NodeType[i][j] == 4 || NodeType[i][j] == 5)glColor4f(1.0,0.0,0.0,1.0);
						else
						{
							if(NodeType[i][j] == 6 || NodeType[i][j] == 8 || NodeType[i][j] == 10 || NodeType[i][j] == 12)glColor4f(0.0,0.0,1.0,1.0);
							else glColor4f(0.0,0.0,0.0,1.0);
						}
					}
					else glColor4f(0.0,1.0,0.0,1.0);
				}
				glVertex3dv(x[i][j]);
			}
		}
	}
	glEnd();
	
	glutSwapBuffers();
	if(StopStart == 1)
	{
		for(int CalculationLoop = 0;CalculationLoop <= 0;CalculationLoop++)
		{
			for(int CableDirection = 0; CableDirection <= 2;CableDirection ++)totalLength[CableDirection] = 0.0;
			
			Calculation();
			
			if(CheckGeometry != 1)
			{
				for(int i = 0;i <= m;i ++)
				{
					for(int j = 0;j <= n;j ++)
					{
						if(NodeType[i][j] >= 0)
						{
							for(int CableDirection = 0; CableDirection <= 2;CableDirection ++)totalLength[CableDirection] += CableLength[i][j][CableDirection];
						}
						for(int CableDirection = 0; CableDirection <= 2;CableDirection ++)CableLength[i][j][CableDirection] = 0.0;
					}
				}
				if(ThirdCable == 1)
					b -= CircleRadiusFactor * a * (2.0 * totalLength[0] - totalLength[1] - totalLength[2]) / (totalLength[0] + totalLength[1] + totalLength[2]);
				else
					b -= CircleRadiusFactor * a * (sqrt(2.0) * totalLength[0] - totalLength[1] - totalLength[2]) / (totalLength[0] + totalLength[1] + totalLength[2]);
				
				double lambda = 0.2;
				for(int i = 0;i <= m;i ++)
				{
					for(int j = 0;j <= n;j ++)
					{
						if(NodeType[i][j] == 2 || NodeType[i][j] ==  6 || NodeType[i][j] == 12){x[i][j][1] += lambda * (+ a - x[i][j][1]);x[i][j][2] = - HalfHeightDifference;}
						if(NodeType[i][j] == 3 || NodeType[i][j] ==  6 || NodeType[i][j] ==  8){x[i][j][0] += lambda * (- a - x[i][j][0]);x[i][j][2] = - HalfHeightDifference;}
						if(NodeType[i][j] == 4 || NodeType[i][j] ==  8 || NodeType[i][j] == 10){x[i][j][1] += lambda * (- a - x[i][j][1]);x[i][j][2] = - HalfHeightDifference;}
						if(NodeType[i][j] == 5 || NodeType[i][j] == 10 || NodeType[i][j] == 12){x[i][j][0] += lambda * (+ a - x[i][j][0]);x[i][j][2] = - HalfHeightDifference;}
					}
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
		Axial(i - 1,i + 1,j,j,0);
		if(j + 1 <= n)
		{
			Axial(i,i + 1,j,j + 1,1);
			Axial(i,i - 1,j,j + 1,2);
		}
	}
}

void CalculateMotion(int i)
{
	if(CheckGeometry != 1)
	{
		for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] >= 0)
			{
				if(NodeType[i][j] == 1)
				{
					double RadiusSq = 0.0;
					double forcedotProduct = 0.0;
					double velocitydotProduct = 0.0;
					for(int xy = 0;xy <= 1;xy ++)
					{
						RadiusSq += x[i][j][xy] * x[i][j][xy];
						forcedotProduct += x[i][j][xy] * CableNetForce[i][j][xy];
						velocitydotProduct += x[i][j][xy] * velocity[i][j][xy];
					}
					double Radius = sqrt(RadiusSq);
					for(int xy = 0;xy <= 1;xy ++)
					{
						CableNetForce[i][j][xy] -= forcedotProduct * x[i][j][xy] / RadiusSq;
						velocity[i][j][xy] -= velocitydotProduct * x[i][j][xy] / RadiusSq;
						x[i][j][xy] *= b / Radius;
					}
					CableNetForce[i][j][2] = 0.0;
					velocity[i][j][2] = 0.0;
					x[i][j][2] = HalfHeightDifference;
				}
				
				for(int xyz = 0;xyz <= 2;xyz ++)
				{
					if((NodeType[i][j] == 0 || NodeType[i][j] == 1) ||
					   ((NodeType[i][j] == 2 || NodeType[i][j] == 4) && xyz == 0) ||
					   ((NodeType[i][j] == 3 || NodeType[i][j] == 5) && xyz == 1))
					{
						velocity[i][j][xyz] = CarryOver * velocity[i][j][xyz] + MovementFactor * CableNetForce[i][j][xyz];
						x[i][j][xyz] += velocity[i][j][xyz];
					}
				}
			}
		}
	}
	Initialise(i);
}

int LineExists(int i,int inext,int j,int jnext)
{
	if((NodeType[i][j] >= 0 && NodeType[inext][jnext] == 0) ||
	   (NodeType[i][j] == 0 && NodeType[inext][jnext] >= 0) ||
	   (NodeType[i][j] >= 0 && NodeType[inext][jnext] >= 0 && (j != jnext || j == 1 || j == nAll || j == n)))return 1;
	else return 0;
}
void Axial(int i_trial,int inext_trial,int j,int jnext,int CableDirection)
{
	double deltax[3];
	
	int i = i_trial;
	if(i < 0)i= m;
	if(i > m)i = 0;
	
	int inext = inext_trial;
	if(inext < 0)inext = m;
	if(inext > m)inext = 0;
	
	if(LineExists(i,inext,j,jnext))
	{
		double lengthSq = 0.0;
		
		for(int xyz = 0;xyz <= 2;xyz ++)
		{
			deltax[xyz] = x[inext][jnext][xyz] - x[i][j][xyz];
			lengthSq += deltax[xyz] * deltax[xyz];
		}
		double length = sqrt(lengthSq);
		
		CableLength[i][j][CableDirection] += length;
		CableLength[inext][jnext][CableDirection] += length;
		
		if(ThirdCable == 1 || j != jnext)
		{
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
	ofstream Julia("MinimalSurfaceNet.dxf");
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	
	int NowBold = 1;
	int LineDrawStep = Spacing - 1;
	int i_LineBegin = 1;
	if(ThirdCable == 1)
	{
	for(int j = 0;j <= n;j ++)
	{
		LineDrawStep ++;
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			NowBold = 1;
		}
		else NowBold = 0;
		
		if(i_LineBegin == 0)i_LineBegin = 1;else i_LineBegin = 0;
		int NowDrawing = 0;
		int PreviousExisted = 0;
		{
			int nexti;
			for(int i = i_LineBegin;i <= m;i += 2)
			{
				nexti = i + 2;
				if(nexti > m)nexti -= m + 1;
				if(LineExists(i,nexti,j,j) && (DoNotDrawFaint == 0 || CheckGeometry == 1))
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
					Julia<<"30\n"<<x[i][j][2]<<"\n";
					
					PreviousExisted = 1;
				}
				else
				{
					if(NowDrawing == 1)
					{
						if(PreviousExisted = 1)
						{
							Julia<<"0\nVERTEX\n";
							Julia<<"10\n"<<x[i][j][0]<<"\n";
							Julia<<"20\n"<<x[i][j][1]<<"\n";
							Julia<<"30\n"<<x[i][j][2]<<"\n";	
						}
						NowDrawing = 0;
						Julia<<"0\nSEQEND\n";
					}
					PreviousExisted = 0;
				}
			}
			if(NowDrawing == 1)
			{
				if(PreviousExisted = 1)
				{
					Julia<<"0\nVERTEX\n";
					Julia<<"10\n"<<x[nexti][j][0]<<"\n";
					Julia<<"20\n"<<x[nexti][j][1]<<"\n";
					Julia<<"30\n"<<x[nexti][j][2]<<"\n";
				}
				Julia<<"0\nSEQEND\n";
			}
		}
	}
	}
	
	LineDrawStep = Spacing - 1;
	for(int i = 0;i <= m;i += 2)
	{
		LineDrawStep ++;
		
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			NowBold = 1;
		}
		else NowBold = 0;
		Julia<<"0\nPOLYLINE\n";
		if(NowBold == 1)Julia<<"100\nAcDbEntity\n8\nBold\n";else Julia<<"100\nAcDbEntity\n8\nFaint\n";
		Julia<<"100\nAcDb3dPolyline\n66\n1\n";
		Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
		if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";
		for(int j = 0;j <= n - 1;j ++)
		{
			int thisi = i + j;while(thisi > m)thisi -= m + 1;
			int nexti = thisi + 1;while(nexti > m)nexti -= m + 1;
			if(LineExists(thisi,nexti,j,j + 1) && (DoNotDrawFaint == 0 || CheckGeometry == 1))
			{
				Julia<<"0\nVERTEX\n";
				Julia<<"10\n"<<x[thisi][j][0]<<"\n";
				Julia<<"20\n"<<x[thisi][j][1]<<"\n";
				Julia<<"30\n"<<x[thisi][j][2]<<"\n";	
			}
		}
		{
			int thisi = i + n;while(thisi > m)thisi -= m + 1;
			Julia<<"0\nVERTEX\n";
			Julia<<"10\n"<<x[thisi][n][0]<<"\n";
			Julia<<"20\n"<<x[thisi][n][1]<<"\n";
			Julia<<"30\n"<<x[thisi][n][2]<<"\n";	
		}
		Julia<<"0\nSEQEND\n";
	}
	
	LineDrawStep = Spacing - 1;
	for(int i = 0;i <= m;i += 2)
	{
		LineDrawStep ++;
		
		if(LineDrawStep >= Spacing)
		{
			LineDrawStep = 0;
			NowBold = 1;
		}
		else NowBold = 0;
		Julia<<"0\nPOLYLINE\n";
		if(NowBold == 1)Julia<<"100\nAcDbEntity\n8\nBold\n";else Julia<<"100\nAcDbEntity\n8\nFaint\n";
		Julia<<"100\nAcDb3dPolyline\n66\n1\n";
		Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
		if(NowBold == 1)Julia<<"62\n0\n";else Julia<<"62\n5\n";
		for(int j = 0;j <= n - 1;j ++)
		{
			int thisi = i - j;while(thisi < 0)thisi += m + 1;
			int nexti = thisi - 1;while(nexti < 0)nexti += m + 1;
			if(LineExists(thisi,nexti,j,j + 1)  && (DoNotDrawFaint == 0 || CheckGeometry == 1))
			{
				Julia<<"0\nVERTEX\n";
				Julia<<"10\n"<<x[thisi][j][0]<<"\n";
				Julia<<"20\n"<<x[thisi][j][1]<<"\n";
				Julia<<"30\n"<<x[thisi][j][2]<<"\n";	
			}
		}
		{
			int thisi = i - n;while(thisi < 0)thisi += m + 1;
			Julia<<"0\nVERTEX\n";
			Julia<<"10\n"<<x[thisi][n][0]<<"\n";
			Julia<<"20\n"<<x[thisi][n][1]<<"\n";
			Julia<<"30\n"<<x[thisi][n][2]<<"\n";	
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