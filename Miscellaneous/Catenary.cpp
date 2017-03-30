#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   MaxNo 501
int   i,m;

float x[MaxNo],y[MaxNo],PI,
      a,b,c,H,r;

ofstream Chris("Catenary.dxf");

int main(void)
{
PI=4.0*atan(1.0);

m=MaxNo-1;
a=1.0;b=-0.3;c=1.5;r=0.01;H=0.3;

for(i=0;i<=m;i+=1)
{
x[i]=b+((a-b)*i)/(1.0*m);
y[i]=c*cosh(x[i]/c);
}

Chris<<"0\nSECTION\n2\nENTITIES\n";

for(i=1;i<=m;i+=1)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[i-1]<<"\n";
Chris<<"20\n"<<y[i-1]<<"\n";
Chris<<"11\n"<<x[i]<<"\n";
Chris<<"21\n"<<y[i]<<"\n";
}

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[0]<<"\n";
Chris<<"20\n"<<y[0]<<"\n";
Chris<<"11\n"<<x[m]<<"\n";
Chris<<"21\n"<<y[m]<<"\n";

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[0]<<"\n";
Chris<<"20\n"<<y[0]-H*sinh(x[0]/c)<<"\n";
Chris<<"11\n"<<x[0]-H<<"\n";
Chris<<"21\n"<<y[0]-H*sinh(x[0]/c)<<"\n";
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[0]<<"\n";
Chris<<"20\n"<<y[0]<<"\n";
Chris<<"11\n"<<x[0]<<"\n";
Chris<<"21\n"<<y[0]-H*sinh(x[0]/c)<<"\n";
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[0]<<"\n";
Chris<<"20\n"<<y[0]<<"\n";
Chris<<"11\n"<<x[0]-H<<"\n";
Chris<<"21\n"<<y[0]-H*sinh(x[0]/c)<<"\n";
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[0]<<"\n";
Chris<<"20\n"<<y[0]<<"\n";
Chris<<"11\n"<<x[0]-H<<"\n";
Chris<<"21\n"<<y[0]-H*(y[m]-y[0])/(x[m]-x[0])<<"\n";
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[0]-H<<"\n";
Chris<<"20\n"<<y[0]-H*(y[m]-y[0])/(x[m]-x[0])<<"\n";
Chris<<"11\n"<<x[0]-H<<"\n";
Chris<<"21\n"<<y[0]-H*sinh(x[0]/c)<<"\n";

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[m]<<"\n";
Chris<<"20\n"<<y[m]+H*sinh(x[m]/c)<<"\n";
Chris<<"11\n"<<x[m]+H<<"\n";
Chris<<"21\n"<<y[m]+H*sinh(x[m]/c)<<"\n";
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[m]<<"\n";
Chris<<"20\n"<<y[m]<<"\n";
Chris<<"11\n"<<x[m]<<"\n";
Chris<<"21\n"<<y[m]+H*sinh(x[m]/c)<<"\n";
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[m]<<"\n";
Chris<<"20\n"<<y[m]<<"\n";
Chris<<"11\n"<<x[m]+H<<"\n";
Chris<<"21\n"<<y[m]+H*sinh(x[m]/c)<<"\n";
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[m]<<"\n";
Chris<<"20\n"<<y[m]<<"\n";
Chris<<"11\n"<<x[m]+H<<"\n";
Chris<<"21\n"<<y[m]+H*(y[m]-y[0])/(x[m]-x[0])<<"\n";
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[m]+H<<"\n";
Chris<<"20\n"<<y[m]+H*(y[m]-y[0])/(x[m]-x[0])<<"\n";
Chris<<"11\n"<<x[m]+H<<"\n";
Chris<<"21\n"<<y[m]+H*sinh(x[m]/c)<<"\n";

Chris<<"0\nCIRCLE\n8\n0\n";
Chris<<"10\n"<<x[0]<<"\n";
Chris<<"20\n"<<y[0]<<"\n";
Chris<<"40\n"<<r<<"\n";
Chris<<"0\nCIRCLE\n8\n0\n";
Chris<<"10\n"<<x[m]<<"\n";
Chris<<"20\n"<<y[m]<<"\n";
Chris<<"40\n"<<r<<"\n";

Chris<<"0\nENDSEC\n0\nEOF\n";

cout<<"DXF file written, end of program\n";

return 0;
}