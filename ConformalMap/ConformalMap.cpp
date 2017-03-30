#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
ofstream Julia("Map.dxf");

#define mplus1   501
#define nplus1   151

double x[mplus1][nplus1],y[mplus1][nplus1],
       c,lambda,myu,myv,uarg,varg,
       tempdouble,r,theta,PI,
       realresult,imagresult,realtop,imagtop,realbot,imagbot;

int i,j,m,n,ymirror,plusorminus,step;

void myfraction(void);
void Lens(void);
void DrawMap(void);

int main()
{
PI=4.0*atan(1.0);

c=100.0;
lambda=1.5;

m=mplus1-1;
n=nplus1-1;
step=2;
m=m/(2*step);m=m*2*step;
n=n/step;n=n*step;

for(i=0;i<=m;i++)
{
for(j=0;j<=n;j++)
{
myu=lambda*c*((2*i-m))/(2.0*int(0.4*m));
myv=lambda*c*((2*j-0))/(2.0*int(0.4*m));
Lens();
}
}

DrawMap();

cout<<"Finished\n";
return 0;
}

void myfraction(void)
{
tempdouble=realbot*realbot+imagbot*imagbot;

realresult=(realtop*realbot+imagtop*imagbot)/tempdouble;
imagresult=(imagtop*realbot-realtop*imagbot)/tempdouble;
}

void Lens(void)
{
realtop=0.0;
imagtop=0.0;

realbot=0.0;
imagbot=0.0;

for(plusorminus=-1;plusorminus<=1;plusorminus+=2)
{
uarg=myu+plusorminus*lambda*c;
varg=myv;

if(uarg!=0.0||varg!=0.0)
{
r=sqrt(uarg*uarg+varg*varg);
if(fabs(uarg)<varg)theta=acos(uarg/r);
else {theta=asin(varg/r);if(uarg<0.0)theta=PI-theta;}
r=pow(r,1.0/lambda);
theta=theta/lambda;

realtop+=r*cos(theta);
imagtop+=r*sin(theta);

realbot+=plusorminus*r*cos(theta);
imagbot+=plusorminus*r*sin(theta);
}
}

myfraction();
x[i][j]=c*realresult;
y[i][j]=c*imagresult;
}

void DrawMap(void)
{
Julia<<"0\nSECTION\n2\nENTITIES\n";
for(ymirror=-1;ymirror<=1;ymirror+=2)
{
for(i=0;i<=m;i+=step)
{
for(j=0;j<=n-1;j++)
{
Julia<<"0\nLINE\n8\n0\n";
Julia<<"10\n"<<x[i][j]<<"\n";
Julia<<"20\n"<<ymirror*y[i][j]<<"\n";
Julia<<"11\n"<<x[i][j+1]<<"\n";
Julia<<"21\n"<<ymirror*y[i][j+1]<<"\n";
}
}
for(i=0;i<=m-1;i++)
{
for(j=0;j<=n;j+=step)
{
Julia<<"0\nLINE\n8\n0\n";
Julia<<"10\n"<<x[i][j]<<"\n";
Julia<<"20\n"<<ymirror*y[i][j]<<"\n";
Julia<<"11\n"<<x[i+1][j]<<"\n";
Julia<<"21\n"<<ymirror*y[i+1][j]<<"\n";
}
}
}
Julia<<"0\nENDSEC\n0\nEOF\n";
Julia.close();
}