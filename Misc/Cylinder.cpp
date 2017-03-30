#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   m 100
#define   n 200

void DXFSetUp(void);
void DXFSurface(void);
void DXFFinishOff(void);

int   i,j,Layer,cycle,jplus1,jminus1;

float x[m+1][n+1],y[m+1][n+1],psi[m+1][n+1],error[m+1][n+1],move[m+1][n+1],
	  PI,a,r,theta,errorSq;

ofstream Chris("Cylinder.dxf");

int main(void)
{
PI=4.0*atan(1.0);

a=1000.0;

for(j=0;j<=n;j+=1)
{
theta=(2.0*j*PI)/(1.0*(n+1));
for(i=0;i<=m;i+=1)
{
r=a*exp((2.0*i*PI)/(1.0*(n+1)));
x[i][j]=r*cos(theta);
y[i][j]=r*sin(theta);
if(i!=m)psi[i][j]=0.0;else psi[i][j]=y[i][j];
move[i][j]=0.0;
}
}

for(cycle=0;cycle<=2000;cycle++)
{
errorSq=0.0;
for(i=1;i<=m-1;i+=1)
{
for(j=0;j<=n;j+=1)
{
jplus1=j+1;if(j==n)jplus1=0;
jminus1=j-1;if(j==0)jminus1=n;
error[i][j]=(psi[i-1][j]+psi[i+1][j]+psi[i][jminus1]+psi[i][jplus1])/4.0-psi[i][j];
errorSq+=error[i][j]*error[i][j];
}
}
for(i=1;i<=m-1;i+=1)
{
for(j=0;j<=n;j+=1)
{
move[i][j]=0.99*move[i][j]+error[i][j];
psi[i][j]+=move[i][j];
}
}
if(100*int(cycle/100)==cycle)cout<<cycle<<"  "<<errorSq<<"\n";
}

DXFSetUp();

Layer=0;
for(i=0;i<=m-1;i+=1)
{
for(j=0;j<=n;j+=1)
{
jplus1=j+1;if(j==n)jplus1=0;
DXFSurface();
}
}

DXFFinishOff();
cout<<"DXF file written, end of program\n";
return 0;
}

void DXFSetUp(void)
{
Chris<<"0\n";
Chris<<"SECTION\n";
Chris<<"2\n";
Chris<<"ENTITIES\n";
}

void DXFSurface(void)
{
Chris<<"0\n3DFACE\n8\npsi\n";
Chris<<"10\n"<<x[i][j]<<"\n";
Chris<<"20\n"<<y[i][j]<<"\n";
Chris<<"30\n"<<psi[i][j]<<"\n";
Chris<<"11\n"<<x[i][jplus1]<<"\n";
Chris<<"21\n"<<y[i][jplus1]<<"\n";
Chris<<"31\n"<<psi[i][jplus1]<<"\n";
Chris<<"12\n"<<x[i+1][jplus1]<<"\n";
Chris<<"22\n"<<y[i+1][jplus1]<<"\n";
Chris<<"32\n"<<psi[i+1][jplus1]<<"\n";
Chris<<"13\n"<<x[i+1][j]<<"\n";
Chris<<"23\n"<<y[i+1][j]<<"\n";
Chris<<"33\n"<<psi[i+1][j]<<"\n";
}

void DXFFinishOff(void)
{
Chris<<"0\n";
Chris<<"ENDSEC\n";
Chris<<"0\n";
Chris<<"EOF\n";
Chris.close();
}