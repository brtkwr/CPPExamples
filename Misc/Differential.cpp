
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

int   i,j,m,n,flip,Colour;

double x[2],y[2],PI,c,L,curve,dx,dy,ds;

ofstream Chris("Curve.dxf");

int main(void)
{
PI=4.0*atan(1.0);
c=1.0;L=1.0;

Chris<<"0\nSECTION\n2\nENTITIES\n";

m=50;
n=40;
for(i=0;i<=m;i+=1)
{
x[0]=(4.0*L*i)/(1.0*m)-2.0*L;
for(j=0;j<=n;j++)
{
y[0]=(3.0*L*j)/(1.0*n)-1.5*L;
dx=y[0];
dy=c*sin(PI*x[0]/L);
ds=sqrt(dx*dx+dy*dy);
if(ds>1.0e-12)
{
dx=0.05*L*dx/ds;
dy=0.05*L*dy/ds;
Chris<<"0\nLINE\n8\nDirections\n";
Chris<<"10\n"<<x[0]+dx/2.0<<"\n20\n"<<y[0]+dy/2.0<<"\n30\n0.0\n";
Chris<<"11\n"<<x[0]-dx/2.0<<"\n21\n"<<y[0]-dy/2.0<<"\n31\n0.0\n";
Chris<<"62\n0\n";
}
}
}

m=250;
Colour=0;
for(curve=-0.8;curve<=2.0;curve+=0.2)
{
Colour+=1;
for(i=0;i<=m-1;i+=1)
{
x[0]=(4.0*L*i)/(1.0*m)-2.0*L;
x[1]=(4.0*L*(i+1))/(1.0*m)-2.0*L;
for(j=0;j<=1;j++)y[j]=2.0*(c*L/PI)*(curve-cos(PI*x[j]/L));

if(y[0]>0.0&&y[1]<0.0){x[1]=x[0]-(x[1]-x[0])*y[0]/(y[1]-y[0]);y[1]=0.0;}
if(y[1]>0.0&&y[0]<0.0){x[0]=x[1]-(x[1]-x[0])*y[1]/(y[1]-y[0]);y[0]=0.0;}

if(y[0]>=0.0&&y[1]>=0.0)
{
for(flip=-1;flip<=1;flip+=2)
{
Chris<<"0\nLINE\n8\nCurves\n";
Chris<<"10\n"<<x[0]<<"\n20\n"<<flip*sqrt(y[0])<<"\n30\n0.0\n";
Chris<<"11\n"<<x[1]<<"\n21\n"<<flip*sqrt(y[1])<<"\n31\n0.0\n";
Chris<<"62\n"<<Colour<<"\n";
}
}
}
}

Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();

cout<<"DXF file written, end of program\n";

return 0;
}
