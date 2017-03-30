#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
ofstream Chris("Circles.dxf");
int main(void)
{
double PI=4.0*atan(1.0);
double R=1.0e3;
double RSq=R*R;
double x,y;
Chris<<"0\nSECTION\n2\nENTITIES\n";
for(int i=0;i<=100000;i+=1)
{
for(;;)
{
x=R*double(2*rand()-RAND_MAX)/double(RAND_MAX);
y=R*double(2*rand()-RAND_MAX)/double(RAND_MAX);
double rSq=x*x+y*y;
double test=double(rand())/double(RAND_MAX);
if(rSq<RSq&&exp(-4.0*rSq/RSq)*(1.0+cos(2.0*sqrt(rSq/R)*PI))>test)break;
}
Chris<<"0\nCIRCLE\n8\n0\n";
Chris<<"10\n"<<x<<"\n";
Chris<<"20\n"<<y<<"\n";
Chris<<"40\n1.0\n";
}
Chris<<"0\nENDSEC\n0\nEOF\n";
cout<<"DXF file written, end of program\n";
return 0;
}
