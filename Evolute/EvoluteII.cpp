#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   PointsPerCircle 50
#define   LastCircle 60

double xCoordinate[LastCircle*PointsPerCircle+1],yCoordinate[LastCircle*PointsPerCircle+1],
xCentreOfCurvature[LastCircle*PointsPerCircle+1],yCentreOfCurvature[LastCircle*PointsPerCircle+1],
Curvature[LastCircle*PointsPerCircle+1],PI;

double x(double,double);
double y(double,double);
double xDash(double,double);
double yDash(double,double);

ofstream Chris("Evolute.dxf");

int main(void)
{
	PI=4.0*atan(1.0);
	
	int LastPoint=LastCircle*PointsPerCircle;
	double a=1000.0;
	
	Chris<<"0\nSECTION\n2\nENTITIES\n";
	
	for(int Point=0;Point<=LastPoint;Point++)
	{
		double theta=2.0*PI*double(Point)/double(LastPoint);
		xCoordinate[Point]=x(a,theta);
		yCoordinate[Point]=y(a,theta);
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
		Chris<<"40\n"<<a/100.0<<"\n";
		Chris<<"62\n0\n";
	}
	for(int Point=0;Point<=LastPoint;Point++)
	{
		double theta=2.0*PI*double(Point)/double(LastPoint);
		
		if(theta!=0.0)
		{
			double sDash=sqrt(xDash(a,theta)*xDash(a,theta)+yDash(a,theta)*yDash(a,theta));
			xCentreOfCurvature[Point]=xCoordinate[Point]-yDash(a,theta)*a*theta/sDash;
			yCentreOfCurvature[Point]=yCoordinate[Point]+xDash(a,theta)*a*theta/sDash;
		}
		else
		{
			xCentreOfCurvature[Point]=xCoordinate[Point];
			yCentreOfCurvature[Point]=yCoordinate[Point];
		}
	}
				for(int Point=0;Point<=LastPoint;Point+=PointsPerCircle)
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
	
	for(int Point=0;Point<=LastPoint-1;Point++)
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
    
	for(int Point=0;Point<=LastPoint;Point+=PointsPerCircle)
	{
		Chris<<"0\nCIRCLE\n8\nCentreOfCurvature\n";
		Chris<<"10\n"<<xCentreOfCurvature[Point]<<"\n";
		Chris<<"20\n"<<yCentreOfCurvature[Point]<<"\n";
		Chris<<"30\n"<<0.0<<"\n";
		Chris<<"40\n"<<a/100.0<<"\n";
		Chris<<"62\n1\n";
	}
	
	Chris<<"0\nENDSEC\n0\nEOF\n";
	Chris.close();
	cout<<"DXF file written, end of program\n";
	
	return 0;
}

double x(double C,double beta)
{
	return C*(cos(beta)+beta*sin(beta));
}

double y(double C,double beta)
{
	return C*(sin(beta)-beta*cos(beta));
}

double xDash(double C,double beta)
{
	return C*beta*cos(beta);
}

double yDash(double C,double beta)
{
	return C*beta*sin(beta);
}
