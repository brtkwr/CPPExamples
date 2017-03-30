#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
ofstream Julia("GeodesicCoordinates.dxf");

#define mplus1   701
#define nplus1   801

double theta1[mplus1][nplus1],theta2[mplus1][nplus1],theta11[mplus1][nplus1],theta21[mplus1][nplus1],
       PI,x,x1,x2,y,myy1,y2,z,z1,z2,
       g11,g12,g22,g,theta1temp,theta2temp,theta12,theta22,L,
       temp;

int i,j,m,n,istep,jstep,igrid,jgrid,itemp;

void curve(void);
void surface(void);

int main()
{
PI=4.0*atan(1.0);

istep=5;igrid=2*istep;m=mplus1-1;m=m/igrid;m=m*igrid;
jstep=5;jgrid=2*jstep;n=nplus1-1;n=n/jgrid;n=n*jgrid;

for(j=0;j<=n;j++)
{
curve();
theta1[0][j]=theta1temp;
theta2[0][j]=theta2temp;
}

for(i=0;i<=m-1;i++)
{
itemp=(i+1)/igrid;if(itemp*igrid==i+1)cout<<i+1<<"\n";
for(j=0;j<=n;j++)
{
if(j==0)
{
theta12=theta1[i][1]-theta1[i][0]-(1.0/2.0)*(theta1[i][2]-2.0*theta1[i][1]+theta1[i][0])/2.0;
theta22=theta2[i][1]-theta2[i][0]-(1.0/2.0)*(theta2[i][2]-2.0*theta2[i][1]+theta2[i][0])/2.0;
}
if(j!=0&&j!=n)
{
theta12=(theta1[i][j+1]-theta1[i][j-1])/2.0;
theta22=(theta2[i][j+1]-theta2[i][j-1])/2.0;
}
if(j==n)
{
theta12=theta1[i][n]-theta1[i][n-1]+(1.0/2.0)*(theta1[i][n]-2.0*theta1[i][n-1]+theta1[i][n-2])/2.0;
theta22=theta2[i][n]-theta2[i][n-1]+(1.0/2.0)*(theta2[i][n]-2.0*theta2[i][n-1]+theta2[i][n-2])/2.0;
}

theta1temp=theta1[i][j];
theta2temp=theta2[i][j];
surface();

g11=x1*x1+myy1*myy1+z1*z1;
g12=x1*x2+myy1*y2+z1*z2;
g22=x2*x2+y2*y2+z2*z2;

g=g11*g22-g12*g12;

temp=sqrt(g*(g11*theta12*theta12+2.0*g12*theta12*theta22+g22*theta22*theta22));

theta11[i][j]=+L*(g12*theta12+g22*theta22)/temp;
theta21[i][j]=-L*(g11*theta12+g12*theta22)/temp;

theta1[i+1][j]=theta1[i][j]+theta11[i][j];
theta2[i+1][j]=theta2[i][j]+theta21[i][j];
if(i!=0)
{
theta1[i+1][j]+=(1.0/2.0)*(theta11[i][j]-theta11[i-1][j]);
theta2[i+1][j]+=(1.0/2.0)*(theta21[i][j]-theta21[i-1][j]);
}
}
}

Julia<<"0\nSECTION\n2\nENTITIES\n";

for(j=0;j<=n;j+=jgrid)
{
for(i=0;i<=m-1;i+=istep)
{
theta1temp=theta1[i][j];
theta2temp=theta2[i][j];
surface();
Julia<<"0\nLINE\n8\n0\n";
Julia<<"10\n"<<x<<"\n";
Julia<<"20\n"<<y<<"\n";
Julia<<"30\n"<<z<<"\n";
theta1temp=theta1[i+istep][j];
theta2temp=theta2[i+istep][j];
surface();
Julia<<"11\n"<<x<<"\n";
Julia<<"21\n"<<y<<"\n";
Julia<<"31\n"<<z<<"\n";
Julia<<"62\n0\n";
}
}

for(i=0;i<=m;i+=igrid)
{
for(j=0;j<=n-1;j+=jstep)
{
theta1temp=theta1[i][j];
theta2temp=theta2[i][j];
surface();
Julia<<"0\nLINE\n8\n0\n";
Julia<<"10\n"<<x<<"\n";
Julia<<"20\n"<<y<<"\n";
Julia<<"30\n"<<z<<"\n";
theta1temp=theta1[i][j+jstep];
theta2temp=theta2[i][j+jstep];
surface();
Julia<<"11\n"<<x<<"\n";
Julia<<"21\n"<<y<<"\n";
Julia<<"31\n"<<z<<"\n";
if(i!=0)Julia<<"62\n0\n";else Julia<<"62\n1\n";
}
}

Julia<<"0\nENDSEC\n0\nEOF\n";
Julia.close();
cout<<"Finished\n";
return 0;
}

void surface(void)
{
double a1,a2,a3,a4;

a1=1000.0;a2=1.0;a3=0.01;a4=1.0;

x=a1*theta1temp*cos(theta2temp);

x1=a1*cos(theta2temp);
x2=-a1*theta1temp*sin(theta2temp);

y=a1*theta1temp*sin(theta2temp);
myy1=a1*sin(theta2temp);
y2=a1*theta1temp*cos(theta2temp);

z=a4*a1*cos(a2*theta1temp)/cosh(a3*theta1temp);
z1=-a4*a1*a2*sin(a2*theta1temp)/cosh(a3*theta1temp)
   -a4*a1*cos(a2*theta1temp)*a3*sinh(a3*theta1temp)/(cosh(a3*theta1temp)*cosh(a3*theta1temp));
z2=0.0;
}

void curve(void)
{
double b1,b2;

b1=0.01;
b2=1.5;
L=10.0;

theta1temp=b1*j+b2;
theta2temp=PI/2.0;
}