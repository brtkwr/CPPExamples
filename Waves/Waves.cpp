//Written by Chris J K Williams, University of Bath, UK

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include <stdio.h>
#include <string.h>
#include "Graphics.h"

#define LastComponent 100

double omega[LastComponent+1],relativeMagnitude[LastComponent+1],Lambda[LastComponent+1],phi[LastComponent+1],
cosine[LastComponent+1],sine[LastComponent+1];
double t,rootdeltaomega;

int main(int argc, char* argv[])
{
	PI=4.0*atan(1.0);
	double g=9.81;
	double omegaMax=10.0;
	double deltaomega=omegaMax/double(LastComponent);
	rootdeltaomega=sqrt(deltaomega);
	for(int Component=0;Component<=LastComponent;Component++)
	{
		omega[Component]=double(Component)*deltaomega;
		Lambda[Component]=2.0*PI*g/(omega[Component]*omega[Component]);
		double Argument=10.0*omega[Component]/omegaMax;
		relativeMagnitude[Component]=exp(-Argument*Argument);
		phi[Component]=2.0*PI*double(rand())/double(RAND_MAX);
		for(;;)
		{
			double theta=2.0*PI*double(rand())/double(RAND_MAX)-PI;
			cosine[Component]=cos(theta);sine[Component]=sin(theta);
			if(double(rand())/double(RAND_MAX)<exp(-5.0*(theta+PI/2.0)*(theta+PI/2.0))>0.0)break;
		}
	}
	
	glutInit(&argc,argv);
	SetUpGraphics(1.0,1.0,1.0);
	
	cout<<"Finished.\n";
	return 0;
}

static void Draw(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	t+=0.2;
	for(int j=0;j<=n;j++)
	{
		for(int i=0;i<=m;i++)
		{
			z[i][j]=0.0;
			for(int Component=0;Component<=LastComponent;Component++)
				z[i][j]+=10.0*rootdeltaomega*relativeMagnitude[Component]
					*cos(2.0*PI*(x[i][j]*cosine[Component]+y[i][j]*sine[Component])/Lambda[Component]
						 -omega[Component]*t+phi[Component]);
		}
	}
	glLineWidth(2.0);
	glColor4d(0.0,0.0,1.0,0.5);
	double rv[3];
	for(int j=0;j<=n;j++)
	{
		glBegin(GL_LINE_STRIP);
		for(int i=0;i<=m;i++)
		{
			rv[0]=x[i][j];rv[1]=y[i][j];rv[2]=z[i][j];glVertex3dv(rv);glVertex3dv(rv);
		}
		glEnd();
	}
	for(int i=0;i<=m;i++)
	{
		glBegin(GL_LINE_STRIP);
		for(int j=0;j<=n;j++)
		{
			rv[0]=x[i][j];rv[1]=y[i][j];rv[2]=z[i][j];glVertex3dv(rv);glVertex3dv(rv);
		}
		glEnd();
	}
	glutSwapBuffers();	
}
