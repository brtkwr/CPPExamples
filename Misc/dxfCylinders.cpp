#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
using namespace std;
#define   LastNode 50

void MakeLine(double radius,double x1,double y1,double z1,double x2,double y2,double z2);
double x[LastNode+1],y[LastNode+1],z[LastNode+1];
int    StrutAtNode[LastNode+1],TieAtNode[LastNode+1],Node,otherNode,Done,TieGo;

ofstream Julia("Cylinders.dxf");
int main()
{
Julia<<"0\nSECTION\n2\nENTITIES\n";
for(Node=0;Node<=LastNode;Node++)
{
x[Node]=10000.0*rand()/(1.0*RAND_MAX);
y[Node]=10000.0*rand()/(1.0*RAND_MAX);
z[Node]=10000.0*rand()/(1.0*RAND_MAX);
StrutAtNode[Node]=0;
TieAtNode[Node]=0;
}
for(Node=0;Node<=LastNode-1;Node++)
{
if(StrutAtNode[Node]==0)
{
Done=0;
for(otherNode=Node+1;otherNode<=LastNode;otherNode++)
{
if(StrutAtNode[otherNode]==0)
{
MakeLine(50.0,x[Node],y[Node],z[Node],x[otherNode],y[otherNode],z[otherNode]);
StrutAtNode[Node]++;
StrutAtNode[otherNode]++;
Done++;
}
if(Done==1)break;
}
}
}

for(TieGo=0;TieGo<=2;TieGo++)
{
for(Node=0;Node<=LastNode-1;Node++)
{
if(TieAtNode[Node]<2)
{
Done=0;
for(otherNode=Node+1;otherNode<=LastNode;otherNode++)
{
if(TieAtNode[otherNode]<2)
{
MakeLine(5.0,x[Node],y[Node],z[Node],x[otherNode],y[otherNode],z[otherNode]);
TieAtNode[Node]++;
TieAtNode[otherNode]++;
Done++;
}
if(Done==3)break;
}
}
}
}

Julia<<"0\nENDSEC\n0\nEOF\n";
cout<<"Finished\n";
return 0;
}

void MakeLine(double radius,double x1,double y1,double z1,double x2,double y2,double z2)
{
double memberlength,Xx,Xy,Xz,Yx,Yy,Yz,Zx,Zy,Zz,thislength;

memberlength=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
Zx=(x1-x2)/memberlength;
Zy=(y1-y2)/memberlength;
Zz=(z1-z2)/memberlength;
if(fabs(Zx)<1.0/64.0&&fabs(Zy)<1.0/64.0)
{
Xx=Zz;Xy=0.0;Xz=-Zx;
}
else
{
Xx=-Zy;Xy=+Zx;Xz=0.0;
}
thislength=sqrt(Xx*Xx+Xy*Xy+Xz*Xz);
Xx=Xx/thislength;Xy=Xy/thislength;Xz=Xz/thislength;
Yx=Zy*Xz-Zz*Xy;
Yy=Zz*Xx-Zx*Xz;
Yz=Zx*Xy-Zy*Xx;
if(radius==5.0)Julia<<"0\nCIRCLE\n8\nTies\n";
else Julia<<"0\nCIRCLE\n8\nStruts\n";
Julia<<"39\n"<<memberlength<<"\n";
Julia<<"10\n"<<x2*Xx+y2*Xy+z2*Xz<<"\n";
Julia<<"20\n"<<x2*Yx+y2*Yy+z2*Yz<<"\n";
Julia<<"30\n"<<x2*Zx+y2*Zy+z2*Zz<<"\n";
Julia<<"40\n"<<radius<<"\n";
Julia<<"210\n"<<Zx<<"\n";
Julia<<"220\n"<<Zy<<"\n";
Julia<<"230\n"<<Zz<<"\n";
}