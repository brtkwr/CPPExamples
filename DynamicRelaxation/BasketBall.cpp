#include "ChrisNet.h"
int    i,j,m,n,even;
double PI,LConst,EAConst,R,angle,mysin,mycos;
int main(void)
{
PI=4.0*atan(1.0);
mysin=0.9;mycos=sqrt(1.0-mysin*mysin);
m=30;n=20;
EAConst=1.0;LConst=1000.0;R=2.0*mysin*LConst*(1.0*m+1.0)/(2.0*PI);
node=-1;even=-1;
for(j=0;j<=n;j++)
{
even=-even;
for(i=0;i<=m;i++)
{
node+=1;
angle=2.0*PI*(1.0*i+(1.0-1.0*even)/4.0)/(1.0*m+1.0);
x[node][0]=R*cos(angle);x[node][1]=R*sin(angle);x[node][2]=-j*LConst*mycos;
if(j==0)nodeType[node]=1;else nodeType[node]=0;
load[node][0]=0.0;load[node][1]=0.0;load[node][2]=-0.001/(1.0*n);
}
}
lastnode=node;
element=-1;
for(j=0;j<=n;j+=2)
{
for(i=0;i<=m;i++)
{
if(j>0)
{
element+=1;end[element][0]=(m+1)*j+i;if(i>0)end[element][1]=end[element][0]-(m+2);else end[element][1]=end[element][0]-1;
element+=1;end[element][0]=(m+1)*j+i;end[element][1]=end[element][0]-(m+1);
}
if(j<n)
{
element+=1;end[element][0]=(m+1)*j+i;if(i>0)end[element][1]=(m+1)*(j+1)+(i-1);else end[element][1]=(m+1)*(j+2)-1;
element+=1;end[element][0]=(m+1)*j+i;end[element][1]=(m+1)*(j+1)+i;
}
}
}
lastelement=element;
for(element=0;element<=lastelement;element++)
{
EA[element]=EAConst;L[element]=LConst;
}
WriteNet();cout<<"Data file written\n";return 0;
}