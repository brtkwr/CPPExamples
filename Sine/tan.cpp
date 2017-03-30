#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   MaxNo 501

int   i,m;

double x[MaxNo],y[MaxNo],PI,maxvalue;

ofstream Chris("tan.dxf");

int main(void)
{
PI=4.0*atan(1.0);
m=MaxNo-1;

for(i=0;i<=m;i+=1)
{
x[i]=(2.0*PI*i)/(1.0*m);
if(fabs(cos(x[i]))<1.0e12)y[i]=tan(x[i]);
}

Chris<<"0\nSECTION\n2\nENTITIES\n";

maxvalue=10.0;
for(i=0;i<=m-1;i+=1)
{
if(y[i]<maxvalue&&y[i+1]<maxvalue&&y[i]>-maxvalue&&y[i+1]>-maxvalue)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[i]  <<"\n20\n"<<y[i]  <<"\n30\n0.0\n";
Chris<<"11\n"<<x[i+1]<<"\n21\n"<<y[i+1]<<"\n31\n0.0\n";
}
}

Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();

cout<<"DXF file written, end of program\n";

return 0;
}