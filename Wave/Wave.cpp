#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define   m 100
#define   n  50
int   i,j;
double x[m+1][n+1],y[m+1][n+1],PI,a,d,amp,lambda;
ofstream Chris("Wave.dxf");
int main(void)
{
PI=4.0*atan(1.0);
a=10000.0;d=a/(1.0*m);
lambda=a/3.0;amp=(0.04*d*n)/sinh((2.0*PI*d*n)/lambda);
for(i=0;i<=m;i++)
{
for(j=0;j<=n;j++)
{
x[i][j]=d*i+amp*cosh((2.0*PI*d*j)/lambda)*cos((2.0*PI*d*i)/lambda);
y[i][j]=d*j-amp*sinh((2.0*PI*d*j)/lambda)*sin((2.0*PI*d*i)/lambda);
}
}
Chris<<"0\nSECTION\n2\nENTITIES\n";
for(i=0;i<=m;i++)
{
for(j=0;j<=n-1;j++)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[i][j+0]<<"\n20\n"<<y[i][j+0]<<"\n";
Chris<<"11\n"<<x[i][j+1]<<"\n21\n"<<y[i][j+1]<<"\n";
}
}
for(i=0;i<=m-1;i++)
{
for(j=0;j<=n;j++)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[i+0][j]<<"\n20\n"<<y[i+0][j]<<"\n";
Chris<<"11\n"<<x[i+1][j]<<"\n21\n"<<y[i+1][j]<<"\n";
}
}
Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();
cout<<"DXF file written, end of program\n";
return 0;
}