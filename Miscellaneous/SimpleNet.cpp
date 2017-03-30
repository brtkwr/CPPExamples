#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

void DefineStructure(void);
void DrawStructure(void);
void dxfSetUp(void);
void dxfMember(void);
void dxfMemberNo(void);
void dxfNodeNo(void);
void dxfFinishOff(void);

#define MaxLastNode 100000
#define MaxLastMember 200000
#define MaxDim  20

int	NodeNumber[2*MaxDim+1][2*MaxDim+1],End[2][MaxLastMember+1],
    NodeType[MaxLastNode+1],
	i,j,k,
    Node,LastNode,Member,LastMember,cycle;

double x[3][MaxLastNode+1],Force[3][MaxLastNode+1],Movement[3][MaxLastNode+1],Stiffness[MaxLastNode+1],
       EA[MaxLastMember+1],L[MaxLastMember+1],SigmaForceSq,Tension,Length,
       textcircsize;

ofstream Picy("Net.dxf");

int main(void)
{
DefineStructure();

for(Node=0;Node<LastNode;Node++)
{
Movement[0][Node]=0.0;
Movement[1][Node]=0.0;
Movement[2][Node]=0.0;
}

for(cycle=0;cycle<=2000;cycle++)
{
for(Node=0;Node<LastNode;Node++)
{
Force[0][Node]=0.0;
Force[1][Node]=0.0;
Force[2][Node]=-0.1;

Stiffness[Node]=0.0;
}

for(Member=0;Member<=LastMember;Member++)
{
Length=sqrt((x[0][End[1][Member]]-x[0][End[0][Member]])
           *(x[0][End[1][Member]]-x[0][End[0][Member]])
		   +(x[1][End[1][Member]]-x[1][End[0][Member]])
           *(x[1][End[1][Member]]-x[1][End[0][Member]])
		   +(x[2][End[1][Member]]-x[2][End[0][Member]])
           *(x[2][End[1][Member]]-x[2][End[0][Member]]));
Tension=EA[Member]*(Length-L[Member])/L[Member];//cout<<Member<<"  "<<Length<<"\n";
if(Tension>0.0)
{
for(k=0;k<=2;k++)
{
Force[k][End[0][Member]]+=Tension*(x[k][End[1][Member]]-x[k][End[0][Member]])/Length;
Force[k][End[1][Member]]-=Tension*(x[k][End[1][Member]]-x[k][End[0][Member]])/Length;
}
}
Stiffness[End[0][Member]]+=EA[Member]/L[Member];
Stiffness[End[1][Member]]+=EA[Member]/L[Member];
}

SigmaForceSq=0.0;
for(Node=0;Node<LastNode;Node++)
{
if(NodeType[Node]==0)
{
for(k=0;k<=2;k++)
{
//cout<<Node<<"  "<<Stiffness[Node]<<"\n";
Movement[k][Node]=0.2*Force[k][Node]/Stiffness[Node]+0.99*Movement[k][Node];
x[k][Node]+=Movement[k][Node];
SigmaForceSq=Force[k][Node]*Force[k][Node];
}
}
}
if(50*int(cycle/50)==cycle)cout<<cycle<<"  "<<SigmaForceSq<<"\n";
}

DrawStructure();
}

void DefineStructure(void)
{
Node=-1;
for(i=0;i<=2*MaxDim;i++)
{
for(j=0;j<=2*MaxDim;j++)
{
if(i+j>=MaxDim&&i+j<=3*MaxDim&&i-j<=MaxDim&&i-j>=-MaxDim)
{
Node++;if(Node>MaxLastNode)cout<<"Too many nodes\n";
NodeNumber[i][j]=Node;
x[0][Node]=1.5*i;x[1][Node]=1.5*j;x[2][Node]=0.0;
if(i==0||i==2*MaxDim||j==0||j==2*MaxDim)NodeType[Node]=1;else NodeType[Node]=0;
}
else NodeNumber[i][j]=-1;
}
}

Member=-1;
for(i=0;i<=2*MaxDim-1;i++)
{
for(j=0;j<=2*MaxDim-1;j++)
{
if(NodeNumber[i][j]>=0&&NodeNumber[i+1][j]>=0)
{
Member++;if(Member>MaxLastMember)cout<<"Too many members\n";
End[0][Member]=NodeNumber[i][j];
End[1][Member]=NodeNumber[i+1][j];
EA[Member]=1.0;L[Member]=1.0;
}
if(NodeNumber[i][j]>=0&&NodeNumber[i][j+1]>=0)
{
Member++;if(Member>MaxLastMember)cout<<"Too many members\n";
End[0][Member]=NodeNumber[i][j];
End[1][Member]=NodeNumber[i][j+1];
EA[Member]=1.0;L[Member]=1.0;
}
}
}

for(i=0;i<=MaxDim-1;i++)
{
Member++;if(Member>MaxLastMember)cout<<"Too many members\n";
End[0][Member]=NodeNumber[i][MaxDim+i];
End[1][Member]=NodeNumber[i+1][MaxDim+i+1];
EA[Member]=10.0;L[Member]=sqrt(2.0);
Member++;if(Member>MaxLastMember)cout<<"Too many members\n";
End[0][Member]=NodeNumber[i][MaxDim-i];
End[1][Member]=NodeNumber[i+1][MaxDim-i-1];
EA[Member]=10.0;L[Member]=sqrt(2.0);
}

for(i=MaxDim;i<=2*MaxDim-1;i++)
{
Member++;if(Member>MaxLastMember)cout<<"Too many members\n";
End[0][Member]=NodeNumber[i][i-MaxDim];
End[1][Member]=NodeNumber[i+1][i-MaxDim+1];
EA[Member]=10.0;L[Member]=sqrt(2.0);
Member++;if(Member>MaxLastMember)cout<<"Too many members\n";
End[0][Member]=NodeNumber[i][3*MaxDim-i];
End[1][Member]=NodeNumber[i+1][3*MaxDim-i-1];
EA[Member]=10.0;L[Member]=sqrt(2.0);
}

LastNode=Node;
LastMember=Member;
}

void DrawStructure(void)
{
dxfSetUp();

textcircsize=0.05;
for(Member=0;Member<=LastMember;Member+=1)dxfMember();

for(Member=0;Member<=LastMember;Member+=1)dxfMemberNo();

for(Node=0;Node<=LastNode;Node+=1)dxfNodeNo();

dxfFinishOff();
}

void dxfSetUp(void)
{
Picy<<"0\n"<<"SECTION\n"<<"2\n"
<<"ENTITIES\n";
}

void dxfMember(void)
{
Picy<<"0\nLINE\n8\nStructure\n"
<<"10\n"<<x[0][End[0][Member]]<<"\n"
<<"20\n"<<x[1][End[0][Member]]<<"\n"
<<"30\n"<<x[2][End[0][Member]]<<"\n"
<<"11\n"<<x[0][End[1][Member]]<<"\n"
<<"21\n"<<x[1][End[1][Member]]<<"\n"
<<"31\n"<<x[2][End[1][Member]]<<"\n"
<<"62\n0\n";
}
void dxfMemberNo(void)
{
Picy<<"0\nTEXT\n8\nInformation\n"
<<"10\n"<<(x[0][End[0][Member]]+x[0][End[1][Member]])/2.0<<"\n"
<<"20\n"<<(x[1][End[0][Member]]+x[1][End[1][Member]])/2.0<<"\n"
<<"40\n"<<textcircsize<<"\n"<<"62\n0\n"
<<"1\n"<<Member<<"\n";
}

void dxfNodeNo(void)
{
Picy<<"0\nTEXT\n8\nInformation\n"
<<"10\n"<<x[0][Node]+textcircsize<<"\n"
<<"20\n"<<x[1][Node]+textcircsize<<"\n"
<<"30\n"<<x[2][Node]+textcircsize<<"\n"
<<"40\n"<<textcircsize<<"\n"<<"62\n1\n"
<<"1\n"<<Node<<"\n";
}

void dxfFinishOff(void)
{
Picy<<"0\n"<<"ENDSEC\n"<<"0\n"<<"EOF\n";
Picy.close();
}

