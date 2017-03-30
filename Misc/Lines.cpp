#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

using namespace std;

ofstream Julia("Lines.dxf");

void StartDrawing(void);
void FindParameters(double&,double&,double&,double&);
void DrawLine(double[],double[],double[]);
void FinishDrawing(void);

int main()
{
double PI,a,b,c,phi,x[2],y[2],z[2],mycos,mytan;
int i,m,mirror;

PI=4.0*atan(1.0);
m=200;

StartDrawing();
FindParameters(a,b,c,phi);

mycos=cos(phi/2.0);
mytan=tan(phi/2.0);

for(mirror=-1;mirror<=1;mirror+=2)
{
for(i=0;i<=m;i++)
{
x[0]=+(a/mycos)*cos((2.0*PI*i)/(1.0*m));
y[0]=+(b/mycos)*sin((2.0*PI*i)/(1.0*m));
z[0]=-(c*mytan);
x[1]=+(a/mycos)*cos((2.0*PI*i)/(1.0*m)+mirror*phi);
y[1]=+(b/mycos)*sin((2.0*PI*i)/(1.0*m)+mirror*phi);
z[1]=+(c*mytan);

DrawLine(x,y,z);
}
}

FinishDrawing();
return 0;
}

void StartDrawing(void)
{
Julia<<"0\nSECTION\n2\nENTITIES\n";
}

void FindParameters(double& aValue,double& bValue,double& cValue,double& phiValue)
{
double halfheight;

cout<<"phi must be greater or equal to 0 and less than 180 degrees\nIn degrees phi = \n";cin>>phiValue;
phiValue=phiValue*atan(1.0)/45.0;
cout<<"a = \n";cin>>aValue;
cout<<"b = \n";cin>>bValue;
cout<<"Half height = c times tan(phi/2)\nHalf height = \n";cin>>halfheight;

cValue=halfheight/tan(phiValue/2.0);
}

void DrawLine(double xCoords[],double yCoords[],double zCoords[])
{
Julia<<"0\nLINE\n8\n0\n";
Julia<<"10\n"<<xCoords[0]<<"\n20\n"<<yCoords[0]<<"\n30\n"<<zCoords[0]<<"\n";
Julia<<"11\n"<<xCoords[1]<<"\n21\n"<<yCoords[1]<<"\n31\n"<<zCoords[1]<<"\n";
}

void FinishDrawing(void)
{
Julia<<"0\nENDSEC\n0\nEOF\n";
Julia.close();
}