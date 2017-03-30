#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
ofstream Julia("Bridge.dxf");

#define mplus1  101
#define nplus1   26

double x[mplus1][nplus1],y[mplus1][nplus1],z[mplus1][nplus1],
       PI,phi,psi,real,imag,a;

int i,j,m,n,TopOrBot,zeropoint;

void SquareRoot(void);
double MyFunction(int thisvalue);
void DrawMap(void);

int main()
{
PI=4.0*atan(1.0);

m=mplus1-1;
n=nplus1-1;
zeropoint=int(m/2.5);

Julia<<"0\nSECTION\n2\nENTITIES\n";

for(TopOrBot=-1;TopOrBot<=1;TopOrBot+=2)
{
for(i=0;i<=m;i++)
{
a=1000.0*(1.0+MyFunction(2*i-m)*MyFunction(2*i-m)/10.0);
phi=2.0*PI*(MyFunction(2*zeropoint-m)*MyFunction(2*zeropoint-m)
        -MyFunction(2*i-m)*MyFunction(2*i-m))/10.0;

for(j=0;j<=n;j++)
{
z[i][j]=20.0*a*MyFunction(2*i-m);
psi=2.0*PI*(0.5*(1+TopOrBot)+(1.0*j)/(1.0*n));
psi=psi-sin(psi);
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
y[i][j]-=1.2*a;
y[i][j]+=y[i][j]*z[i][j]*z[i][j]/(1000.0*a*a);
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

double MyFunction(int thisvalue)
{
return (1.0*thisvalue)/(1.0*m)+((2.0*zeropoint-1.0*m)/(PI*m))*sin(PI*thisvalue/(2.0*zeropoint-1.0*m));
}

void DrawMap(void)
{
for(i=0;i<=m-1;i++)
{
for(j=0;j<=n-1;j++)
{
Julia<<"0\n3DFACE\n8\n0\n";
Julia<<"10\n"<<z[i][j]<<"\n";
Julia<<"20\n"<<x[i][j]<<"\n";
Julia<<"30\n"<<y[i][j]<<"\n";
Julia<<"11\n"<<z[i+1][j]<<"\n";
Julia<<"21\n"<<x[i+1][j]<<"\n";
Julia<<"31\n"<<y[i+1][j]<<"\n";
Julia<<"12\n"<<z[i+1][j+1]<<"\n";
Julia<<"22\n"<<x[i+1][j+1]<<"\n";
Julia<<"32\n"<<y[i+1][j+1]<<"\n";
Julia<<"13\n"<<z[i][j+1]<<"\n";
Julia<<"23\n"<<x[i][j+1]<<"\n";
Julia<<"33\n"<<y[i][j+1]<<"\n";
}
}
}
