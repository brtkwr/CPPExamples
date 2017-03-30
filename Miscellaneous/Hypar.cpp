#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   m 50
#define   n 50

void DXFSetUp(void);
void DXFSurface(int,int);
void DXFFinishOff(void);

double x[m+1][n+1],y[m+1][n+1],z[m+1][n+1],
       PI,a,b,u,v;

ofstream Chris("Surface.dxf");

int main(void)
{
PI=4.0*atan(1.0);

a=1000.0;
b=500.0;

for(int i=0;i<=m;i++)
{
u=(1.0*i)/(1.0*m);
for(int j=0;j<=n;j++)
{
v=(1.0*j)/(1.0*n);
x[i][j]=u*a*cos(2.0*PI*v);
y[i][j]=u*a*sin(2.0*PI*v);
z[i][j]=b*x[i][j]*y[i][j]/(a*a);
}
}

DXFSetUp();

for(int i=0;i<=m-1;i++)
{
for(int j=0;j<=n-1;j++)
{
DXFSurface(i,j);
}
}

DXFFinishOff();
cout<<"DXF file written, end of program\n";
return 0;
}

void DXFSetUp(void)
{
Chris<<"0\nSECTION\n2\nENTITIES\n";
}

void DXFSurface(int i,int j)
{
Chris<<"0\n3DFACE\n8\nSurface\n";
Chris<<"10\n"<<x[i][j]<<"\n";
Chris<<"20\n"<<y[i][j]<<"\n";
Chris<<"30\n"<<z[i][j]<<"\n";
Chris<<"11\n"<<x[i][j+1]<<"\n";
Chris<<"21\n"<<y[i][j+1]<<"\n";
Chris<<"31\n"<<z[i][j+1]<<"\n";
Chris<<"12\n"<<x[i+1][j+1]<<"\n";
Chris<<"22\n"<<y[i+1][j+1]<<"\n";
Chris<<"32\n"<<z[i+1][j+1]<<"\n";
Chris<<"13\n"<<x[i+1][j]<<"\n";
Chris<<"23\n"<<y[i+1][j]<<"\n";
Chris<<"33\n"<<z[i+1][j]<<"\n";
}

void DXFFinishOff(void)
{
Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();
}