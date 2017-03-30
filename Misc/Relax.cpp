#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

using namespace std;
ofstream Julia("Random.dxf");

#define mplus1   101
#define nplus1   101

double x[mplus1][nplus1],y[mplus1][nplus1],z[mplus1][nplus1],
       PI,a,b;

int i,j,m,n,cycle;

int main()
{
PI=4.0*atan(1.0);
m=mplus1-1;
n=nplus1-1;

a=10000.0/(1.0*m);
b=1.0e-5;

for(i=0;i<=m;i++)
{
x[i][0]=a*(2*i-m)+b*rand();
y[i][0]=-a*n+b*rand();
z[i][0]=+b*rand();

x[i][n]=a*(2*i-m)+b*rand();
y[i][n]=+a*n+b*rand();;
z[i][n]=+b*rand();
}

for(j=1;j<=n-1;j++)
{
x[0][j]=-a*m+b*rand();
y[0][j]=a*(2*j-n)+b*rand();
z[0][j]=b*rand();

x[m][j]=+a*m+b*rand();
y[m][j]=a*(2*j-n)+b*rand();
z[m][j]=b*rand();
}

for(i=1;i<=m-1;i++)
{
for(j=1;j<=n-1;j++)
{
x[i][j]=a*(2*i-m);
y[i][j]=a*(2*j-n);
z[i][j]=0.0;
}
}

for(cycle=0;cycle<=1000;cycle++)
{
if(100*(cycle/100)==cycle)cout<<cycle<<"\n";

for(i=1;i<=m-1;i++)
{
for(j=1;j<=n-1;j++)
{
x[i][j]=(x[i-1][j]+x[i+1][j]+x[i][j-1]+x[i][j+1])/4.0;
y[i][j]=(y[i-1][j]+y[i+1][j]+y[i][j-1]+y[i][j+1])/4.0;
z[i][j]=(z[i-1][j]+z[i+1][j]+z[i][j-1]+z[i][j+1])/4.0;
}
}
}

Julia<<"0\nSECTION\n2\nENTITIES\n";

for(j=0;j<=n-1;j++)
{
for(i=0;i<=m-1;i++)
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

Julia<<"0\nENDSEC\n0\nEOF\n";
Julia.close();
cout<<"Finished\n";
return 0;
}
