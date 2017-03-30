#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   LastPoint 500

int   i,Layer,Colour;

float x[LastPoint+1],y[LastPoint+1],PI;

ofstream Chris("Hyperbola.dxf");

int main(void)
{
PI=4.0*atan(1.0);

double a=1000.0;
double thetaMax=1.5;
int LastLine=10;

Chris<<"0\nSECTION\n2\nENTITIES\n";

Chris<<"0\nLINE\n8\nLines\n";
Chris<<"10\n"<<0.0<<"\n";
Chris<<"20\n"<<0.0<<"\n";
Chris<<"30\n"<<0.0<<"\n";
Chris<<"11\n"<<a*exp(thetaMax)<<"\n";
Chris<<"21\n"<<0.0<<"\n";
Chris<<"31\n"<<0.0<<"\n";
Chris<<"62\n0\n";

Chris<<"0\nLINE\n8\nLines\n";
Chris<<"10\n"<<0.0<<"\n";
Chris<<"20\n"<<0.0<<"\n";
Chris<<"30\n"<<0.0<<"\n";
Chris<<"11\n"<<0.0<<"\n";
Chris<<"21\n"<<a*exp(thetaMax)<<"\n";
Chris<<"31\n"<<0.0<<"\n";
Chris<<"62\n0\n";

for(int Line=0;Line<=LastLine;Line++)
{
double theta=double(2*Line-LastLine)*thetaMax/double(LastLine);
double thisx=a*exp(+theta);
double thisy=a*exp(-theta);
double CircleRadius=a*exp(thetaMax)/500.0;

Chris<<"0\nLINE\n8\nLines\n";
Chris<<"10\n"<<0.0<<"\n";
Chris<<"20\n"<<thisy<<"\n";
Chris<<"30\n"<<0.0<<"\n";
Chris<<"11\n"<<thisx<<"\n";
Chris<<"21\n"<<0.0<<"\n";
Chris<<"31\n"<<0.0<<"\n";
Chris<<"62\n0\n";

Chris<<"0\nCIRCLE\n8\nLines\n";
Chris<<"10\n"<<0.0<<"\n";
Chris<<"20\n"<<thisy<<"\n";
Chris<<"30\n"<<0.0<<"\n";
Chris<<"40\n"<<CircleRadius<<"\n";
Chris<<"62\n0\n";

Chris<<"0\nCIRCLE\n8\nLines\n";
Chris<<"10\n"<<thisx<<"\n";
Chris<<"20\n"<<0.0<<"\n";
Chris<<"30\n"<<0.0<<"\n";
Chris<<"40\n"<<CircleRadius<<"\n";
Chris<<"62\n0\n";

Chris<<"0\nCIRCLE\n8\nLines\n";
Chris<<"10\n"<<thisx/2.0<<"\n";
Chris<<"20\n"<<thisy/2.0<<"\n";
Chris<<"30\n"<<0.0<<"\n";
Chris<<"40\n"<<CircleRadius<<"\n";
Chris<<"62\n1\n";
}

for(int Point=0;Point<=LastPoint;Point++)
{
double theta=double(2*Point-LastPoint)*thetaMax/double(LastPoint);
x[Point]=a*exp(+theta)/2.0;
y[Point]=a*exp(-theta)/2.0;
}

for(int Point=0;Point<=LastPoint-1;Point++)
{
Chris<<"0\nLINE\n8\nCurve\n";
Chris<<"10\n"<<x[Point]<<"\n";
Chris<<"20\n"<<y[Point]<<"\n";
Chris<<"30\n"<<0.0<<"\n";
Chris<<"11\n"<<x[Point+1]<<"\n";
Chris<<"21\n"<<y[Point+1]<<"\n";
Chris<<"31\n"<<0.0<<"\n";
Chris<<"62\n2\n";
}

Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();
cout<<"DXF file written, end of program\n";

return 0;
}