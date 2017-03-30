#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   Maxm 101
#define   Maxn 51

int    i,m,j,n,step,integern,flip,slice;
double x[Maxm][Maxn],y[Maxm][Maxn],PI,a,phi,psi,doublen,
       real,imag,argument,modulus;

ofstream Chris("Sources.dxf");

int main(void)
{
PI=4.0*atan(1.0);
step=2*4;//Must be even
m=Maxm-1;m=m/step;m=m/2;m=2*m-1;m=m*step;
n=Maxn-1;n=n/step;n=n*step;
a=10000.0;
integern=5;
doublen=1.0*integern;

Chris<<"0\nSECTION\n2\nENTITIES\n";

for(slice=0;slice<=integern-1;slice++)
{
for(flip=-1;flip<=1;flip+=2)
{
for(i=0;i<=m;i+=1)
{
for(j=0;j<=n;j+=1)
{
phi=(PI*(2.0*i-1.0*m))/(2.0*n);
psi=(PI*j)/(1.0*n);
real=1.0-exp(phi)*cos(psi);
imag=exp(phi)*sin(psi);
modulus=sqrt(real*real+imag*imag);
if(fabs(real)>fabs(imag))
{
argument=asin(imag/modulus);
if(real<0.0)argument=PI-argument;
}
else
{
argument=acos(real/modulus);
if(imag<0.0)argument=-argument;
}
if(argument<0.0)argument+=2.0*PI;

x[i][j]=+a*cos(flip*argument/doublen+slice*2.0*PI/doublen)/pow(modulus,1.0/doublen);
y[i][j]=-a*sin(flip*argument/doublen+slice*2.0*PI/doublen)/pow(modulus,1.0/doublen);
}
}

for(i=0;i<=m;i+=step)
{
for(j=0;j<=n-1;j+=1)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[i][j]  <<"\n20\n"<<y[i][j]  <<"\n30\n0.0\n";
Chris<<"11\n"<<x[i][j+1]<<"\n21\n"<<y[i][j+1]<<"\n31\n0.0\n";
}
}
for(i=0;i<=m-1;i+=1)
{
for(j=step/2;j<=n-step/2;j+=step)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[i][j]  <<"\n20\n"<<y[i][j]  <<"\n30\n0.0\n";
Chris<<"11\n"<<x[i+1][j]<<"\n21\n"<<y[i+1][j]<<"\n31\n0.0\n";
}
}
}
}

Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();

cout<<"DXF file written, end of program\n";

return 0;
}
