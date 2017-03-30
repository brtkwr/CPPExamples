#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define   m 200
#define   n 25
int   LastSurface;
double PI,a,u,v,lambda,
      x[m+1][n+1],y[m+1][n+1],z[m+1][n+1];
ofstream Maud("Surfaces.dxf");
int main(void)
{
PI=4.0*atan(1.0);a=10000.0;LastSurface=5;
Maud<<"0\nSECTION\n2\nENTITIES\n";
for(int Surface=0;Surface<=LastSurface;Surface++)
{
lambda=PI*double(Surface)/double(2*LastSurface);
for(int i=0;i<=m;i+=1)
{
u=double(2*i-m)*PI/double(m);
for(int j=0;j<=n;j+=1)
{
v=double(2*j-n)*PI/double(m);
x[i][j]=a*(cos(lambda)*cosh(v)*cos(u)+sin(lambda)*sinh(v)*sin(u));
y[i][j]=a*(cos(lambda)*cosh(v)*sin(u)-sin(lambda)*sinh(v)*cos(u));
z[i][j]=a*(cos(lambda)*v-sin(lambda)*u);
if(Surface>=3)
{
y[i][j]-=5.0*a*double(Surface-4);
z[i][j]-=3.3*a;
}
else
{
y[i][j]-=5.0*a*double(Surface-1);
z[i][j]+=3.3*a;
}
}
}
for(int i=0;i<=m-1;i+=1)
{
for(int j=0;j<=n-1;j+=1)
{
Maud<<"0\n3DFACE\n8\n0\n";
Maud<<"10\n"<<x[i+0][j+0]<<"\n20\n"<<y[i+0][j+0]<<"\n30\n"<<z[i+0][j+0]<<"\n";
Maud<<"11\n"<<x[i+0][j+1]<<"\n21\n"<<y[i+0][j+1]<<"\n31\n"<<z[i+0][j+1]<<"\n";
Maud<<"12\n"<<x[i+1][j+1]<<"\n22\n"<<y[i+1][j+1]<<"\n32\n"<<z[i+1][j+1]<<"\n";
Maud<<"13\n"<<x[i+1][j+0]<<"\n23\n"<<y[i+1][j+0]<<"\n33\n"<<z[i+1][j+0]<<"\n";
}
}
}
Maud<<"0\nENDSEC\n0\nEOF\n";Maud.close();
cout<<"DXF file written, end of program\n";
return 0;
}
