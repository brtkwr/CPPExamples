#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
ofstream Julia("Map.dxf");

#define mplus1   51
#define nplus1   51

double x[mplus1][nplus1],y[mplus1][nplus1],z[mplus1][nplus1],
       PI,phi,psi,real,imag,a;

int i,j,m,n,TopOrBot;

void SquareRoot(void);
void DrawMap(void);

int main()
{
PI=4.0*atan(1.0);

a=1000.0;

m=mplus1-1;
n=nplus1-1;

Julia<<"0\nSECTION\n2\nENTITIES\n";

for(TopOrBot=-1;TopOrBot<=1;TopOrBot+=2)
{
for(i=0;i<=m;i++)
{
phi=2.0*PI*(i-25)/(1.0*n);
for(j=0;j<=n;j++)
{
z[i][j]=0.0;
psi=2.0*PI*(0.5*(1+TopOrBot)+(1.0*j)/(1.0*n));
real=exp(phi)*cos(psi)-1.0;
imag=exp(phi)*sin(psi);
SquareRoot();
x[i][j]=a*x[i][j];
y[i][j]=a*y[i][j];
if(x[i][j]*cos(psi/2.0)<0.0||y[i][j]*TopOrBot>0.0)
{
x[i][j]=-x[i][j];
y[i][j]=-y[i][j];
}
}
}

DrawMap();
}

Julia<<"0\nENDSEC\n0\nEOF\n";
Julia.close();
cout<<"Finished\n";
return 0;
}

void SquareRoot(void)
{
double tempdouble;

tempdouble=sqrt(real*real+imag*imag);
if(tempdouble>1.0e-12)
{
if(real>=0.0)
{
x[i][j]=sqrt((tempdouble+real)/2.0);
y[i][j]=imag/(2.0*x[i][j]);
}
else
{
y[i][j]=sqrt((tempdouble-real)/2.0);
x[i][j]=imag/(2.0*y[i][j]);
}
}
else
{
y[i][j]=0.0;
x[i][j]=0.0;
}
}

void DrawMap(void)
{
for(i=0;i<=m-1;i++)
{
for(j=0;j<=n-1;j++)
{
Julia<<"0\n3DFACE\n8\n0\n";
Julia<<"10\n"<<x[i][j]<<"\n";
Julia<<"20\n"<<y[i][j]<<"\n";
Julia<<"30\n"<<z[i][j]<<"\n";
Julia<<"11\n"<<x[i+1][j]<<"\n";
Julia<<"21\n"<<y[i+1][j]<<"\n";
Julia<<"31\n"<<z[i+1][j]<<"\n";
Julia<<"12\n"<<x[i+1][j+1]<<"\n";
Julia<<"22\n"<<y[i+1][j+1]<<"\n";
Julia<<"32\n"<<z[i+1][j+1]<<"\n";
Julia<<"13\n"<<x[i][j+1]<<"\n";
Julia<<"23\n"<<y[i][j+1]<<"\n";
Julia<<"33\n"<<z[i][j+1]<<"\n";
}
}
}