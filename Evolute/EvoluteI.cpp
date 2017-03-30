#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   PointsPerCircle 5
#define   LastCircle 600

double xCoordinate[LastCircle*PointsPerCircle+1],yCoordinate[LastCircle*PointsPerCircle+1],
xCentreOfCurvature[LastCircle*PointsPerCircle+1],yCentreOfCurvature[LastCircle*PointsPerCircle+1],
Curvature[LastCircle*PointsPerCircle+1],PI;

double x(double,int);
double y(double,int);
double xDash(double,int);
double yDash(double,int);
double xDashDash(double,int);
double yDashDash(double,int);

ofstream Chris("Evolute.dxf");

int main(void)
{
	PI=4.0*atan(1.0);
	int   n=4;
	
	int LastPoint=LastCircle*PointsPerCircle;
	double MinCurvature=1.0/5.0;
	
	Chris<<"0\nSECTION\n2\nENTITIES\n";
	
	for(int Point=0;Point<=LastPoint;Point++)
	{
		double theta=2.0*PI*double(Point)/double(LastPoint);
		xCoordinate[Point]=x(theta,n);
		yCoordinate[Point]=y(theta,n);
	}
	
	for(int Point=0;Point<=LastPoint-1;Point++)
	{
		Chris<<"0\nLINE\n8\nCurve\n";
		Chris<<"10\n"<<xCoordinate[Point]<<"\n";
		Chris<<"20\n"<<yCoordinate[Point]<<"\n";
		Chris<<"30\n"<<0.0<<"\n";
		Chris<<"11\n"<<xCoordinate[Point+1]<<"\n";
		Chris<<"21\n"<<yCoordinate[Point+1]<<"\n";
		Chris<<"31\n"<<0.0<<"\n";
		Chris<<"62\n2\n";
	}
	
	for(int Point=0;Point<=LastPoint;Point+=PointsPerCircle)
	{
        Chris<<"0\nCIRCLE\n8\nCircleOnCurve\n";
		Chris<<"10\n"<<xCoordinate[Point]<<"\n";
		Chris<<"20\n"<<yCoordinate[Point]<<"\n";
		Chris<<"30\n"<<0.0<<"\n";
		Chris<<"40\n"<<0.01<<"\n";
		Chris<<"62\n0\n";
	}
	for(int Point=0;Point<=LastPoint;Point++)
	{
		double theta=2.0*PI*double(Point)/double(LastPoint);
		double sDash=sqrt(xDash(theta,n)*xDash(theta,n)+yDash(theta,n)*yDash(theta,n));
		Curvature[Point]=(yDashDash(theta,n)*xDash(theta,n)-xDashDash(theta,n)*yDash(theta,n))/(sDash*sDash*sDash);
		
		
		xCentreOfCurvature[Point]=xCoordinate[Point]-yDash(theta,n)/(sDash*Curvature[Point]);
		yCentreOfCurvature[Point]=yCoordinate[Point]+xDash(theta,n)/(sDash*Curvature[Point]);
	}
				for(int Point=0;Point<=LastPoint;Point+=PointsPerCircle)
				{
					if(fabs(Curvature[Point])>=MinCurvature)
					{
						Chris<<"0\nLINE\n8\nRadiusOfCurvature\n";
						Chris<<"10\n"<<xCoordinate[Point]<<"\n";
						Chris<<"20\n"<<yCoordinate[Point]<<"\n";
						Chris<<"30\n"<<0.0<<"\n";
						Chris<<"11\n"<<xCentreOfCurvature[Point]<<"\n";
						Chris<<"21\n"<<yCentreOfCurvature[Point]<<"\n";
						Chris<<"31\n"<<0.0<<"\n";
						Chris<<"62\n0\n";
					}
				}
	
	for(int Point=0;Point<=LastPoint-1;Point++)
	{
		if(fabs(Curvature[Point])>=MinCurvature&&fabs(Curvature[Point+1])>=MinCurvature)
		{
			Chris<<"0\nLINE\n8\nEvolute\n";
			Chris<<"10\n"<<xCentreOfCurvature[Point]<<"\n";
			Chris<<"20\n"<<yCentreOfCurvature[Point]<<"\n";
			Chris<<"30\n"<<0.0<<"\n";
			Chris<<"11\n"<<xCentreOfCurvature[Point+1]<<"\n";
			Chris<<"21\n"<<yCentreOfCurvature[Point+1]<<"\n";
			Chris<<"31\n"<<0.0<<"\n";
			Chris<<"62\n3\n";
		}
    }
    
	for(int Point=0;Point<=LastPoint;Point+=PointsPerCircle)
	{
		if(fabs(Curvature[Point])>=MinCurvature)
		{
			Chris<<"0\nCIRCLE\n8\nCentreOfCurvature\n";
			Chris<<"10\n"<<xCentreOfCurvature[Point]<<"\n";
			Chris<<"20\n"<<yCentreOfCurvature[Point]<<"\n";
			Chris<<"30\n"<<0.0<<"\n";
			Chris<<"40\n"<<0.01<<"\n";
			Chris<<"62\n1\n";
		}
	}
	
	Chris<<"0\nENDSEC\n0\nEOF\n";
	Chris.close();
	cout<<"DXF file written, end of program\n";
	
	return 0;
}

double x(double beta,int Ellie)
{
	return (5.0+cos(double(Ellie)*beta))*cos(beta);
}

double y(double beta,int Ellie)
{
	return (5.0+cos(double(Ellie)*beta))*sin(beta);
}

double xDash(double beta,int Ellie)
{
	return -double(Ellie)*sin(double(Ellie)*beta)*cos(beta)-(5.0+cos(double(Ellie)*beta))*sin(beta);
}

double yDash(double beta,int Ellie)
{
	return -double(Ellie)*sin(double(Ellie)*beta)*sin(beta)+(5.0+cos(double(Ellie)*beta))*cos(beta);
}

double xDashDash(double beta,int Ellie)
{
	return -double(Ellie)*double(Ellie)*cos(double(Ellie)*beta)*cos(beta)+2.0*double(Ellie)*sin(double(Ellie)*beta)*sin(beta)
	-(5.0+cos(double(Ellie)*beta))*cos(beta);
}

double yDashDash(double beta,int Ellie)
{
	return -double(Ellie)*double(Ellie)*cos(double(Ellie)*beta)*sin(beta)-2.0*double(Ellie)*sin(double(Ellie)*beta)*cos(beta)
	-(5.0+cos(double(Ellie)*beta))*sin(beta);
}
