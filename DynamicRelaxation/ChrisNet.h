#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#define   MaxNodesplus1    17001
#define   MaxElementsplus1  8001
int WriteNet(void);
int ReadNet(void);
int    nodeType[MaxNodesplus1],end[MaxElementsplus1][2],node,lastnode,element,lastelement,direction;
double x[MaxNodesplus1][3],load[MaxNodesplus1][3],L[MaxElementsplus1],EA[MaxElementsplus1];

int WriteNet(void)
{
ofstream Maud("NetData.txt");
for(node=0;node<=lastnode;node++)
Maud<<node<<"  "<<x[node][0]<<"  "<<x[node][1]<<"  "<<x[node][2]<<"  "<<nodeType[node]<<"\n";
Maud<<"-1\n";
for(element=0;element<=lastelement;element++)
Maud<<element<<"  "<<end[element][0]<<"  "<<end[element][1]<<"  "<<L[element]<<"  "<<EA[element]<<"\n";
Maud<<"-1\n";
for(node=0;node<=lastnode;node++)
{
for(direction=0;direction<=2;direction++)
{
if(load[node][direction]!=0.0)Maud<<node<<"  "<<direction<<"  "<<load[node][direction]<<"\n";
}
}
Maud<<"-1\n";
Maud.close();return 0;
}

int ReadNet(void)
{
ifstream Julia("NetData.txt");lastnode=-1;
for(;;)
{
Julia>>node;if(node<0)break;if(node>MaxNodesplus1-1){cout<<"Too many nodes\n";return 0;}
Julia>>x[node][0]>>x[node][1]>>x[node][2]>>nodeType[node];
lastnode+=1;
}
lastelement=-1;
for(;;)
{
Julia>>element;if(element<0)break;if(element>MaxElementsplus1-1){cout<<"Too many elements\n";return 0;}
Julia>>end[element][0]>>end[element][1]>>L[element]>>EA[element];
lastelement+=1;
}
for(node=0;node<=lastnode;node++)
{
load[node][0]=0.0;load[node][1]=0.0;load[node][2]=0.0;
}
for(;;)
{
Julia>>node;if(node<0)break;Julia>>direction;Julia>>load[node][direction];
}
Julia.close();return 0;
}
