#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   MaxNo 501

int   i,m,harmonic;

double x[MaxNo],y[MaxNo],PI;

ofstream Chris("Sine.dxf");

int main(void)
{
PI=4.0*atan(1.0);
m=MaxNo-1;

for(i=0;i<=m;i+=1)y[i]=0.0;

Chris<<"0\nSECTION\n2\nENTITIES\n";

for(harmonic=1;harmonic<=21;harmonic+=2)
{
for(i=0;i<=m;i+=1)
{
x[i]=(1.0*i)/(1.0*m);
y[i]+=(2.0/PI)*sin(harmonic*2.0*PI*x[i])/(1.0*harmonic);
}

for(i=0;i<=m-1;i+=1)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[i]  <<"\n20\n"<<y[i]  <<"\n30\n0.0\n";
Chris<<"11\n"<<x[i+1]<<"\n21\n"<<y[i+1]<<"\n31\n0.0\n";
Chris<<"62\n"<<harmonic<<"\n";
}
}

Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();

cout<<"DXF file written, end of program\n";

return 0;
}