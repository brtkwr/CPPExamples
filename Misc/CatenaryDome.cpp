#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   m 50
#define   n  2

void DrawLine(void);

int   i,j,otheri,otherj;

float x[m*n+1][m*n+1],y[m*n+1][m*n+1],z[m*n+1][m*n+1],
      a,b,c,sinalpha,cosalpha,u,v;

ofstream Maud("Dome.dxf");

int main(void)
{
a=500.0;
b=5.0e3;
c=5.0e3;
cosalpha=10.0/sqrt(10.0*10.0+7.5*7.5);
sinalpha=sqrt(1.0-cosalpha*cosalpha);

for(i=0;i<=m*n;i++)
{
u=a*((1.0*i)/(1.0*n)-0.5*m);
for(j=0;j<=m*n;j++)
{
v=a*((1.0*j)/(1.0*n)-0.5*m);
x[i][j]=c*cosalpha*(asinh(u/c)-asinh(v/c));
y[i][j]=c*sinalpha*(asinh(u/c)+asinh(v/c));
z[i][j]=b-c*(sqrt(1.0+(u/c)*(u/c))+sqrt(1.0+(v/c)*(v/c))-2.0);
}
}

Maud<<"0\nSECTION\n2\nENTITIES\n";

for(i=0;i<=m*n;i+=n)
{
for(j=0;j<=m*n-1;j++)
{
otheri=i;
otherj=j+1;
DrawLine();
}
}

for(j=0;j<=m*n;j+=n)
{
for(i=0;i<=m*n-1;i++)
{
otheri=i+1;
otherj=j;
DrawLine();
}
}

Maud<<"0\nENDSEC\n0\nEOF\n";
Maud.close();
cout<<"DXF file written, end of program\n";

return 0;
}

void DrawLine(void)
{
if(z[i][j]>0.0||z[otheri][otherj]>0.0)Maud<<"0\nLINE\n8\nShell\n";
else Maud<<"0\nLINE\n8\nExcess";
if(z[i][j]<0.0&&z[otheri][otherj]>0.0)
{
Maud<<"\n10\n"<<x[otheri][otherj]-(x[i][j]-x[otheri][otherj])*z[otheri][otherj]/(z[i][j]-z[otheri][otherj]);
Maud<<"\n20\n"<<y[otheri][otherj]-(y[i][j]-y[otheri][otherj])*z[otheri][otherj]/(z[i][j]-z[otheri][otherj]);
Maud<<"\n30\n"<<0.0<<"\n";
}
else
{Maud<<"\n10\n"<<x[i][j]<<"\n20\n"<<y[i][j]<<"\n30\n"<<z[i][j]<<"\n";}
if(z[i][j]>0.0&&z[otheri][otherj]<0.0)
{
Maud<<"\n11\n"<<x[i][j]-(x[otheri][otherj]-x[i][j])*z[i][j]/(z[otheri][otherj]-z[i][j]);
Maud<<"\n21\n"<<y[i][j]-(y[otheri][otherj]-y[i][j])*z[i][j]/(z[otheri][otherj]-z[i][j]);
Maud<<"\n31\n"<<0.0<<"\n";
}
else
{Maud<<"11\n"<<x[otheri][otherj]<<"\n21\n"<<y[otheri][otherj]<<"\n31\n"<<z[otheri][otherj]<<"\n";}
if(z[i][j]>0.0||z[otheri][otherj]>0.0)Maud<<"62\n0\n";
else Maud<<"62\n1\n";
}