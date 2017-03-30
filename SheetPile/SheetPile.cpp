#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   MaxmPlus1 500
#define   MaxnPlus1 500

void Calcdz(void);

int   i,j,m,n,imax,jmax,step,mirror,reflect;

double x[MaxmPlus1][MaxnPlus1],y[MaxmPlus1][MaxnPlus1],
       PI,a,myh,p,q,dphi,dpsi,dx,dy,cospiaoverh,rootthing,
       delta,xtemp,ytemp,myu,myv,realdzbydeta,imagdzbydeta;

ofstream Chris("SheetPile.dxf");

int main(void)
{
PI=4.0*atan(1.0);
step=10;

myh=1.0;a=0.8;

m=MaxmPlus1-1;n=MaxnPlus1-1;

cospiaoverh=cos(PI*a/myh);
rootthing=sqrt(1.0-cospiaoverh);

delta=PI/(1.5*m);

i=0;
x[i][0]=0.0;y[i][0]=-myh;
dphi=0.0;dpsi=delta;
for(j=1;j<=n;j++)
{
xtemp=x[i][j-1];ytemp=y[i][j-1];Calcdz();
xtemp+=dx/2.0;ytemp+=dy/2.0;Calcdz();
x[i][j]=x[i][j-1]+dx;
y[i][j]=y[i][j-1]+dy;
jmax=j;if(y[i][j]>-a*(1.00001))break;
}
dphi=delta;dpsi=0.0;
for(j=jmax;j>=0;j-=1)
{
for(i=1;i<=m;i++)
{
xtemp=x[i-1][j];ytemp=y[i-1][j];Calcdz();
xtemp+=dx/2.0;ytemp+=dy/2.0;Calcdz();
x[i][j]=x[i-1][j]+dx;
y[i][j]=y[i-1][j]+dy;
if(j==jmax)imax=i-1;
if(j==jmax&&y[i][j]>0.0)break;
if(j!=jmax&&i==imax)break;
}
}
Chris<<"0\nSECTION\n2\nENTITIES\n";

for(reflect=-1;reflect<=1;reflect+=2)
{
for(j=0;j<=jmax;j+=step)
{
for(i=1;i<=imax;i++)
{
Chris<<"0\nLINE\n8\n0\n";
if(i-1!=0)mirror=reflect;else mirror=1;
Chris<<"10\n"<<mirror*x[i-1][j]<<"\n";
Chris<<"20\n"<<y[i-1][j]<<"\n";
mirror=reflect;
Chris<<"11\n"<<mirror*x[i][j]<<"\n";
Chris<<"21\n"<<y[i][j]<<"\n";
}
}

mirror=reflect;
for(i=0;i<=imax;i+=step)
{
if(i!=0||reflect==1)
{
for(j=1;j<=jmax;j++)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<mirror*x[i][j-1]<<"\n";
Chris<<"20\n"<<y[i][j-1]<<"\n";
Chris<<"11\n"<<mirror*x[i][j]<<"\n";
Chris<<"21\n"<<y[i][j]<<"\n";
}
}
}
}

Chris<<"0\nENDSEC\n0\nEOF\n";

cout<<"DXF file written, end of program\n";

return 0;
}

void Calcdz(void)
{
p=cosh(PI*xtemp/myh)*cos(PI*ytemp/myh)-cospiaoverh;
q=sinh(PI*xtemp/myh)*sin(PI*ytemp/myh);

myv=sqrt((-p+sqrt(p*p+q*q))/2.0);
myu=q/(2.0*myv);

realdzbydeta=+a*myv/rootthing;
imagdzbydeta=-a*myu/rootthing;

dx=realdzbydeta*dphi-imagdzbydeta*dpsi;
dy=imagdzbydeta*dphi+realdzbydeta*dpsi;
}
