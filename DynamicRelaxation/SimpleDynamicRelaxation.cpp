#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#define   MaxLastNode    17001
#define   MaxLastElement 100001
int    NodeType[MaxLastNode+1],End0[MaxLastElement+1],End1[MaxLastElement+1],
       i,j,m,n,Node,LastNode,Element,LastElement,
	   Cycle,LastCycle,PrintCycle;
double x[MaxLastNode+1],y[MaxLastNode+1],z[MaxLastNode+1],
       xLoad[MaxLastNode+1],yLoad[MaxLastNode+1],zLoad[MaxLastNode+1],
	   xVely[MaxLastNode+1],yVely[MaxLastNode+1],zVely[MaxLastNode+1],
	   xForce[MaxLastNode+1],yForce[MaxLastNode+1],zForce[MaxLastNode+1],
	   Stiffness[MaxLastNode+1],
	   L[MaxLastElement+1],EA[MaxLastElement+1],
	   a,LConst,EAConst,EdgeStiff,deltax,deltay,deltaz,ActualLength,TensionCoefficient,
	   CarryOver,Relax,SumSq,DrawRadius;

ofstream Julia("Net.dxf");
int main(void)
{
LastCycle=10000;PrintCycle=100;
m=50;n=50;
a=1000.0;DrawRadius=a/100.0;
LConst=1.0*a;
EAConst=1.0;EdgeStiff=100.0*EAConst;
Node=-1;
CarryOver=0.9;
Relax=0.5;
for(j=0;j<=n;j++)
{
for(i=0;i<=m;i++)
{
Node++;if(Node>MaxLastNode){cout<<"Too many Nodes\n";return 0;}
x[Node]=i*a;y[Node]=a*j;z[Node]=0.0;
if((i==0&&j==0)
 ||(i==0&&j==n)
 ||(i==m&&j==0)
 ||(i==m&&j==n))
NodeType[Node]=1;else NodeType[Node]=0;
xLoad[Node]=0.0;yLoad[Node]=0.0;zLoad[Node]=0.01;
xVely[Node]=0.0;yVely[Node]=0.0;zVely[Node]=0.0;
Stiffness[Node];
}
}
LastNode=Node;
Element=-1;
for(j=0;j<=n;j++)
{
for(i=0;i<=m-1;i++)
{
Element++;if(Element>MaxLastElement){cout<<"Too many Elements\n";return 0;}
End0[Element]=(m+1)*j+i;End1[Element]=(m+1)*j+i+1;
if(j!=0&&j!=n)EA[Element]=EAConst;else EA[Element]=EdgeStiff;
L[Element]=LConst;
}
}
for(i=0;i<=m;i++)
{
for(j=0;j<=n-1;j++)
{
Element++;if(Element>MaxLastElement){cout<<"Too many Elements\n";return 0;}
End0[Element]=(m+1)*j+i;End1[Element]=(m+1)*(j+1)+i;
if(i!=0&&i!=m)EA[Element]=EAConst;else EA[Element]=EdgeStiff;
L[Element]=LConst;
}
}
LastElement=Element;
for(Element=0;Element<=LastElement;Element++)
{
Stiffness[End0[Element]]+=EA[Element]/L[Element];
Stiffness[End1[Element]]+=EA[Element]/L[Element];
}
Julia<<"0\nSECTION\n2\nENTITIES\n";
for(Node=0;Node<=LastNode;Node++)
{
Julia<<"0\nCIRCLE\n8\nOriginalNodes\n";
Julia<<"10\n"<<x[Node]<<"\n";
Julia<<"20\n"<<y[Node]<<"\n";
Julia<<"30\n"<<z[Node]<<"\n";
Julia<<"40\n"<<DrawRadius<<"\n";
Julia<<"62\n"<<NodeType[Node]<<"\n";
}

for(Element=0;Element<=LastElement;Element++)
{
Julia<<"0\nLINE\n8\nUnrelaxed\n";
Julia<<"10\n"<<x[End0[Element]]<<"\n";
Julia<<"20\n"<<y[End0[Element]]<<"\n";
Julia<<"30\n"<<z[End0[Element]]<<"\n";
Julia<<"11\n"<<x[End1[Element]]<<"\n";
Julia<<"21\n"<<y[End1[Element]]<<"\n";
Julia<<"31\n"<<z[End1[Element]]<<"\n";
Julia<<"62\n1\n";
}

for(Cycle=0;Cycle<=LastCycle;Cycle++)
{
for(Node=0;Node<=LastNode;Node++)
{
xForce[Node]=xLoad[Node];
yForce[Node]=yLoad[Node];
zForce[Node]=zLoad[Node];
}
for(Element=0;Element<=LastElement;Element++)
{
deltax=x[End1[Element]]-x[End0[Element]];
deltay=y[End1[Element]]-y[End0[Element]];
deltaz=z[End1[Element]]-z[End0[Element]];
ActualLength=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
if(ActualLength>1.0e-12*a)
{
TensionCoefficient=EA[Element]*(ActualLength-L[Element])/(ActualLength*L[Element]);
xForce[End0[Element]]+=TensionCoefficient*deltax;
yForce[End0[Element]]+=TensionCoefficient*deltay;
zForce[End0[Element]]+=TensionCoefficient*deltaz;
xForce[End1[Element]]-=TensionCoefficient*deltax;
yForce[End1[Element]]-=TensionCoefficient*deltay;
zForce[End1[Element]]-=TensionCoefficient*deltaz;
}
}
SumSq=0.0;
for(Node=0;Node<=LastNode;Node++)
{
if(NodeType[Node]==0)
{
SumSq+=xForce[Node]*xForce[Node]+yForce[Node]*yForce[Node]+zForce[Node]*zForce[Node];
xVely[Node]=CarryOver*xVely[Node]+Relax*xForce[Node]/Stiffness[Node];
yVely[Node]=CarryOver*yVely[Node]+Relax*yForce[Node]/Stiffness[Node];
zVely[Node]=CarryOver*zVely[Node]+Relax*zForce[Node]/Stiffness[Node];
x[Node]+=xVely[Node];
y[Node]+=yVely[Node];
z[Node]+=zVely[Node];
}
}
if(PrintCycle*int(Cycle/PrintCycle)==Cycle)cout<<Cycle<<"  "<<SumSq<<"\n";
}
for(Node=0;Node<=LastNode;Node++)
{
Julia<<"0\nCIRCLE\n8\nMovedNodes\n";
Julia<<"10\n"<<x[Node]<<"\n";
Julia<<"20\n"<<y[Node]<<"\n";
Julia<<"30\n"<<z[Node]<<"\n";
Julia<<"40\n"<<DrawRadius<<"\n";
Julia<<"62\n"<<NodeType[Node]<<"\n";
}
for(Element=0;Element<=LastElement;Element++)
{
Julia<<"0\nLINE\n8\nRelaxed\n";
Julia<<"10\n"<<x[End0[Element]]<<"\n";
Julia<<"20\n"<<y[End0[Element]]<<"\n";
Julia<<"30\n"<<z[End0[Element]]<<"\n";
Julia<<"11\n"<<x[End1[Element]]<<"\n";
Julia<<"21\n"<<y[End1[Element]]<<"\n";
Julia<<"31\n"<<z[End1[Element]]<<"\n";
Julia<<"62\n0\n";
}
Julia<<"0\nENDSEC\n0\nEOF\n";
Julia.close();
cout<<"Finished\n";return 0;
}