#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   NoSegments 500

double x[NoSegments+1],y[NoSegments+1],z[NoSegments+1],PI;

ofstream Chris("SuperEllipse.dxf");

int main(void)
{
PI=4.0*atan(1.0);

double a=10000.0;
double b=a/sqrt(2.0);

Chris<<"0\nSECTION\n2\nENTITIES\n";


for(int Curve=1;Curve<=50;Curve++)
{
double n=double(Curve)/2.0;
for(int i=0;i<=NoSegments;i++)
{
double theta=(2.0*PI*double(i))/double(NoSegments);
x[i]=a*pow(fabs(cos(theta)),2.0/n);if(cos(theta)<0.0)x[i]=-x[i];
y[i]=b*pow(fabs(sin(theta)),2.0/n);if(sin(theta)<0.0)y[i]=-y[i];
z[i]=0.0;
}

for(int i=0;i<=NoSegments-1;i++)
{
Chris<<"0\nLINE\n8\n0Curves\n";
Chris<<"10\n"<<x[i]<<"\n";
Chris<<"20\n"<<y[i]<<"\n";
Chris<<"30\n"<<z[i]<<"\n";
Chris<<"11\n"<<x[i+1]<<"\n";
Chris<<"21\n"<<y[i+1]<<"\n";
Chris<<"31\n"<<z[i+1]<<"\n";
Chris<<"62\n0\n";
}
}

Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();
cout<<"DXF file written, end of program\n";

return 0;
}