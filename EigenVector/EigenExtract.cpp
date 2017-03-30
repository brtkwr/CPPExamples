#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include "size.h"

void Eigen(int n,double A[][nplus1],double x[][nplus1],double Lambda[])
{
//Note Stodola's method. Matrix must be symmetric
int i,j,eigen,cycle,maxcycles,previous;
double y[nplus1],
       sumsq,rootsumsq,top,bottom,oldvalue,mult;
       
maxcycles=n*n*n*n;
for(eigen=0;eigen<=n;eigen++)
{
sumsq=0.0;
for(i=0;i<=n;i++)
{
x[i][eigen]=rand()/32768.0-0.5;
sumsq+=x[i][eigen]*x[i][eigen];
}
rootsumsq=sqrt(sumsq);for(i=0;i<=n;i++)x[i][eigen]=x[i][eigen]/rootsumsq;
oldvalue=0.0;
for(cycle=1;cycle<=maxcycles;cycle++)
{
if(eigen!=0)
{
for(previous=0;previous<=eigen-1;previous++)
{
mult=0.0;for(i=0;i<=n;i++)mult+=x[i][previous]*x[i][eigen];
for(i=0;i<=n;i++)x[i][eigen]-=mult*x[i][previous];
}
}
for(i=0;i<=n;i++)
{
y[i]=0.0;
for(j=0;j<=n;j++)y[i]+=A[i][j]*x[j][eigen];
}
top=0.0;bottom=0.0;
for(i=0;i<=n;i++)
{
top+=x[i][eigen]*y[i];
bottom+=x[i][eigen]*x[i][eigen];
}
Lambda[eigen]=top/bottom;
if(fabs(oldvalue/Lambda[eigen]-1.0)<1.0e-6)break;
oldvalue=Lambda[eigen];
sumsq=0.0;
for(i=0;i<=n;i++)
{
x[i][eigen]=y[i];
sumsq+=x[i][eigen]*x[i][eigen];
}
rootsumsq=sqrt(sumsq);for(i=0;i<=n;i++)x[i][eigen]=x[i][eigen]/rootsumsq;
if(cycle==maxcycles)cout<<"Has not converged\n";
}
}
}