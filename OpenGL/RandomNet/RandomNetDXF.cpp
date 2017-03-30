#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define maxpts 5000

double x[3][maxpts],force[3][maxpts],movement[3][maxpts],a,sumsq;

int    end[2][2*maxpts],weighting[maxpts],
       node,member,lastnode,lastmember,lastcontrolnode,coord,cycle,printcycle;

ofstream Julia("Random.dxf");

int main(void)
{
lastnode=maxpts-1;
lastmember=2*maxpts-1;
lastcontrolnode=250;
a=10000.0;
printcycle=100;

cout<<"Maximum random number with this compiler is "<<RAND_MAX<<"\n";

for(node=0;node<=lastnode;node++)
{
if(node<=lastcontrolnode)
{
for(;;)
{
for(coord=0;coord<=2;coord++)x[coord][node]=a*((2.0*rand())/(1.0*RAND_MAX)-1.0);
if(x[0][node]*x[0][node]+x[1][node]*x[1][node]+x[2][node]*x[2][node]<=a*a)break;
}
}
else
{
for(coord=0;coord<=2;coord++)x[coord][node]=0.0;
}
for(coord=0;coord<=2;coord++)movement[coord][node]=0.0;

end[0][2*node+0]=node;end[1][2*node+0]=rand()%lastnode;
end[0][2*node+1]=node;end[1][2*node+1]=rand()%lastnode;
weighting[node]=0;
}

for(member=0;member<=lastmember;member++)
{
weighting[end[0][member]]+=1;
weighting[end[1][member]]+=1;
}

for(cycle=0;cycle<=2000;cycle++)
{
sumsq=0.0;
for(node=0;node<=lastnode;node++)
{
for(coord=0;coord<=2;coord++)force[coord][node]=0.0;
}

for(member=0;member<=lastmember;member++)
{
for(coord=0;coord<=2;coord++)
{
force[coord][end[0][member]]+=x[coord][end[1][member]]-x[coord][end[0][member]];
force[coord][end[1][member]]+=x[coord][end[0][member]]-x[coord][end[1][member]];
}
}

for(node=lastcontrolnode+1;node<=lastnode;node++)
{
for(coord=0;coord<=2;coord++)
{
movement[coord][node]=force[coord][node]/(1.0*weighting[node])+0.99*movement[coord][node];
sumsq+=movement[coord][node]*movement[coord][node];
x[coord][node]+=movement[coord][node];
}
}
if(printcycle*int(cycle/printcycle)==cycle)cout<<cycle<<"   "<<sumsq<<"\n";
}

Julia<<"0\nSECTION\n2\nENTITIES\n";

for(member=0;member<=lastmember;member++)
{
Julia<<"0\nLINE\n8\n0\n";
Julia<<"10\n"<<x[0][end[0][member]]<<"\n";
Julia<<"20\n"<<x[1][end[0][member]]<<"\n";
Julia<<"30\n"<<x[2][end[0][member]]<<"\n";
Julia<<"11\n"<<x[0][end[1][member]]<<"\n";
Julia<<"21\n"<<x[1][end[1][member]]<<"\n";
Julia<<"31\n"<<x[2][end[1][member]]<<"\n";
}

for(node=0;node<=lastnode;node++)
{
if(node<=lastcontrolnode)Julia<<"0\nCIRCLE\n8\nFixed\n";
else Julia<<"0\nCIRCLE\n8\nMovable\n";
Julia<<"10\n"<<x[0][node]<<"\n";
Julia<<"20\n"<<x[1][node]<<"\n";
Julia<<"30\n"<<x[2][node]<<"\n";
Julia<<"40\n"<<a/500.0<<"\n";
if(node<=lastcontrolnode)Julia<<"62\n5\n";else Julia<<"62\n1\n";
}

Julia<<"0\nENDSEC\n0\nEOF\n";Julia.close();
cout<<"Finished\n";
return 0;
}