#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   LastPoint 500

double x[LastPoint+1],y[LastPoint+1],z[LastPoint+1];

ofstream Chris("Ellipse.dxf");

int main(void)
{
	double PI=4.0*atan(1.0);
	
	double a=10000.0;
	double b=2000.0;
	
	for(int Point=0;Point<=LastPoint;Point++)
	{
		double theta=2.0*PI*double(Point)/double(LastPoint+1);
		x[Point]=a*cos(theta);
		y[Point]=b*sin(theta);
		z[Point]=0.0;
	}
	
	Chris<<"0\nSECTION\n2\nENTITIES\n";
	
	int Layer=0;int Colour=0;
	for(int Point=0;Point<=LastPoint;Point++)
	{
		Chris<<"0\nLINE\n8\n"<<Layer<<"\n";
		Chris<<"10\n"<<x[Point]<<"\n";
		Chris<<"20\n"<<y[Point]<<"\n";
		Chris<<"30\n"<<z[Point]<<"\n";
		int nextPoint=Point+1;
		if(nextPoint>LastPoint)nextPoint-=LastPoint+1;
		Chris<<"11\n"<<x[nextPoint]<<"\n";
		Chris<<"21\n"<<y[nextPoint]<<"\n";
		Chris<<"31\n"<<z[nextPoint]<<"\n";
		Chris<<"62\n"<<Colour<<"\n";
	}
	
	Chris<<"0\nENDSEC\n0\nEOF\n";
	Chris.close();
	cout<<"DXF file written, end of program\n";
	
	return 0;
}