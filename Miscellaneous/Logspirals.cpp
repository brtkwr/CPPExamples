#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define LastPoint 200
#define LastCurve 100
double x[LastCurve+1][LastPoint+1],y[LastCurve+1][LastPoint+1],z[LastCurve+1][LastPoint+1];
ofstream Chris("LogSpiralCurves.dxf");
ofstream Maud("LogSpiralSurface.dxf");
int main(void)
{
	double PI = 4.0*atan(1.0);
	double a = 1.0;
	Chris<<"0\nSECTION\n2\nENTITIES\n";
	Maud<<"0\nSECTION\n2\nENTITIES\n";
	for(int clock = -1;clock <= 1;clock += 2)
	{
		for(int Curve = 0;Curve <= LastCurve;Curve++)
		{
			double thetaStart = PI*double(Curve)/double(LastCurve);
			for(int Point = 0;Point <= LastPoint;Point++)
			{
				double theta = PI*double(Point)/double(LastPoint);
				double r = a*exp(theta);
				theta += thetaStart;
				theta *= double(clock);
				x[Curve][Point] = r*cos(theta);
				y[Curve][Point] = r*sin(theta);
				z[Curve][Point] = a*float(Curve);
			}
		}
		for(int Curve = 0;Curve <= LastCurve;Curve++)
		{
			int Layer = 0;int Colour = 0;
			if(Curve == 27)Colour = 1;
			for(int Point = 0;Point <= LastPoint - 1;Point++)
			{
				Chris<<"0\nLINE\n8\n"<<Layer<<"\n";
				Chris<<"10\n"<<x[Curve][Point  ]<<"\n20\n"<<y[Curve][Point  ]<<"\n30\n"<<z[Curve][Point  ]<<"\n";
				Chris<<"11\n"<<x[Curve][Point+1]<<"\n21\n"<<y[Curve][Point+1]<<"\n31\n"<<z[Curve][Point+1]<<"\n";
				Chris<<"62\n"<<Colour<<"\n";
			}
		}
		for(int Curve = 0;Curve <= LastCurve - 1;Curve++)
		{
			int Layer = 0;
			for(int Point = 0;Point <= LastPoint - 1;Point++)
			{
				Maud<<"0\n3DFACE\n8\n"<<Layer<<"\n";
				Maud<<"10\n"<<x[Curve  ][Point  ]<<"\n20\n"<<y[Curve  ][Point  ]<<"\n30\n"<<z[Curve  ][Point  ]<<"\n";
				Maud<<"11\n"<<x[Curve  ][Point+1]<<"\n21\n"<<y[Curve  ][Point+1]<<"\n31\n"<<z[Curve  ][Point+1]<<"\n";
				Maud<<"12\n"<<x[Curve+1][Point+1]<<"\n22\n"<<y[Curve+1][Point+1]<<"\n32\n"<<z[Curve+1][Point+1]<<"\n";
				Maud<<"13\n"<<x[Curve+1][Point  ]<<"\n23\n"<<y[Curve+1][Point  ]<<"\n33\n"<<z[Curve+1][Point  ]<<"\n";
			}
		}
	}
	Chris<<"0\nENDSEC\n0\nEOF\n";Chris.close();
	Maud<<"0\nENDSEC\n0\nEOF\n";Maud.close();
	cout<<"DXF files written, end of program\n";
	return 0;
}