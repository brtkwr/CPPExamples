#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define  m 50
#define  n 50
int   i,j;
double x[m+1][n+1],y[m+1][n+1],z[m+1][n+1],
       PI,R,r,theta,phi,phimax,ratio;
ofstream Maud("Torus lines.dxf");
ofstream Madeleine("Torus faces.dxf");
int main(void)
{
PI=4.0*atan(1.0);
r=1000.0;
cout<<"Enter ratio R over r and then press return.\n";
cout<<"Any value greater than -1 is possible.\n";
cout<<"Negative values give cigar.\n";
cout<<"0 gives sphere.\n";
cout<<"Positive values greater than 1 give complete torus.\n";
for(;;)
{
cin>>ratio;
if(ratio>-1.0)break;else cout<<"Must be greater than -1\n";
}
R=ratio*r;
if(ratio<1.0)phimax=acos(-ratio);
else phimax=PI;
for(i=0;i<=m;i++)
{
theta=(2.0*PI*i)/(1.0*m);
for(j=0;j<=n;j++)
{
phi=phimax*(2.0*j-1.0*n)/(1.0*n);
x[i][j]=(R+r*cos(phi))*cos(theta);
y[i][j]=(R+r*cos(phi))*sin(theta);
z[i][j]=r*sin(phi);
}
}
Maud<<"0\nSECTION\n2\nENTITIES\n";
for(i=0;i<=m-1;i++)
{
for(j=0;j<=n-1;j++)
{
Maud<<"0\nLINE\n8\nLines 0\n";
Maud<<"10\n"<<x[i][j]<<"\n";
Maud<<"20\n"<<y[i][j]<<"\n";
Maud<<"30\n"<<z[i][j]<<"\n";
Maud<<"11\n"<<x[i][j+1]<<"\n";
Maud<<"21\n"<<y[i][j+1]<<"\n";
Maud<<"31\n"<<z[i][j+1]<<"\n";
Maud<<"62\n0\n";

Maud<<"0\nLINE\n8\nLines 1\n";
Maud<<"10\n"<<x[i][j]<<"\n";
Maud<<"20\n"<<y[i][j]<<"\n";
Maud<<"30\n"<<z[i][j]<<"\n";
Maud<<"11\n"<<x[i+1][j]<<"\n";
Maud<<"21\n"<<y[i+1][j]<<"\n";
Maud<<"31\n"<<z[i+1][j]<<"\n";
Maud<<"62\n1\n";
}
}
Maud<<"0\nENDSEC\n0\nEOF\n";Maud.close();
Madeleine<<"0\nSECTION\n2\nENTITIES\n";
for(i=0;i<=m-1;i++)
{
for(j=0;j<=n-1;j++)
{
Madeleine<<"0\n3DFACE\n8\nFaces\n";
Madeleine<<"10\n"<<x[i][j]<<"\n";
Madeleine<<"20\n"<<y[i][j]<<"\n";
Madeleine<<"30\n"<<z[i][j]<<"\n";
Madeleine<<"11\n"<<x[i+1][j]<<"\n";
Madeleine<<"21\n"<<y[i+1][j]<<"\n";
Madeleine<<"31\n"<<z[i+1][j]<<"\n";
Madeleine<<"12\n"<<x[i+1][j+1]<<"\n";
Madeleine<<"22\n"<<y[i+1][j+1]<<"\n";
Madeleine<<"32\n"<<z[i+1][j+1]<<"\n";
Madeleine<<"13\n"<<x[i][j+1]<<"\n";
Madeleine<<"23\n"<<y[i][j+1]<<"\n";
Madeleine<<"33\n"<<z[i][j+1]<<"\n";
Madeleine<<"62\n0\n";
}
}
Madeleine<<"0\nENDSEC\n0\nEOF\n";Madeleine.close();
cout<<"DXF files written, end of program\n";
return 0;
}