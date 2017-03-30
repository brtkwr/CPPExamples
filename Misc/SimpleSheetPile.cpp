#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   MaxmPlus1  41
#define   MaxnPlus1 101

void Calcdz(void);

int   i,j,m,n,step;

double x[MaxmPlus1][MaxnPlus1],y[MaxmPlus1][MaxnPlus1],
       PI,c,theta,psi;

ofstream Chris("SheetPile.dxf");

int main(void)
{
PI=4.0*atan(1.0);
c=1000.0;
step=5;

m=MaxmPlus1-1;n=MaxnPlus1-1;
m=m/step;m=m*step;
n=n/step;n=n*step;

for(i=0;i<=m;i++)
{
for(j=0;j<=n;j++)
{
theta=(1.0*PI*i)/(1.0*n);
psi=(1.0*PI*j)/(1.0*n)-PI/2.0;
x[i][j]=+c*sinh(theta)*sin(psi);
y[i][j]=-c*cosh(theta)*cos(psi);
}
}
Chris<<"0\nSECTION\n2\nENTITIES\n";

for(j=0;j<=n;j+=step)
{
for(i=0;i<=m-1;i++)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[i][j]<<"\n";
Chris<<"20\n"<<y[i][j]<<"\n";
Chris<<"11\n"<<x[i+1][j]<<"\n";
Chris<<"21\n"<<y[i+1][j]<<"\n";
}
}

for(i=0;i<=m;i+=step)
{
for(j=0;j<=n-1;j++)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[i][j]<<"\n";
Chris<<"20\n"<<y[i][j]<<"\n";
Chris<<"11\n"<<x[i][j+1]<<"\n";
Chris<<"21\n"<<y[i][j+1]<<"\n";
}
}


Chris<<"0\nENDSEC\n0\nEOF\n";

cout<<"DXF file written, end of program\n";

return 0;
}
