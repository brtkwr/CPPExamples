#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   m 100
#define   n 200

void Contour(float ContourInterval,float xTriangle[3],float yTriangle[3],float origpsiTraingle[3]);

int   i,j,Layer,cycle,jplus1,jminus1;

float x[m+1][n+1],y[m+1][n+1],psi[m+1][n+1],error[m+1][n+1],move[m+1][n+1],
	  xplot[3],yplot[3],psiplot[3],PI,a,r,theta,errorSq,ContInt;

ofstream Chris("Cylinder.dxf");
ofstream Maud("Contours.dxf");

int main(void)
{
PI=4.0*atan(1.0);

a=1000.0;
ContInt=200.0;

for(j=0;j<=n;j+=1)
{
theta=(2.0*j*PI)/(1.0*(n+1));
for(i=0;i<=m;i+=1)
{
r=a*exp((2.0*i*PI)/(1.0*(n+1)));
x[i][j]=r*cos(theta);
y[i][j]=r*sin(theta);
if(i!=m)psi[i][j]=0.0;else psi[i][j]=y[i][j];
move[i][j]=0.0;
}
}

for(cycle=0;cycle<=2000;cycle++)
{
errorSq=0.0;
for(i=1;i<=m-1;i+=1)
{
for(j=0;j<=n;j+=1)
{
jplus1=j+1;if(j==n)jplus1=0;
jminus1=j-1;if(j==0)jminus1=n;
error[i][j]=(psi[i-1][j]+psi[i+1][j]+psi[i][jminus1]+psi[i][jplus1])/4.0-psi[i][j];
errorSq+=error[i][j]*error[i][j];
}
}
for(i=1;i<=m-1;i+=1)
{
for(j=0;j<=n;j+=1)
{
move[i][j]=0.99*move[i][j]+error[i][j];
psi[i][j]+=move[i][j];
}
}
if(100*int(cycle/100)==cycle)cout<<cycle<<"  "<<errorSq<<"\n";
}

Chris<<"0\nSECTION\n2\nENTITIES\n";
Maud<<"0\nSECTION\n2\nENTITIES\n";

Layer=0;
for(i=0;i<=m-1;i+=1)
{
for(j=0;j<=n;j+=1)
{
jplus1=j+1;if(j==n)jplus1=0;
Chris<<"0\n3DFACE\n8\npsi\n";
Chris<<"10\n"<<x[i][j]<<"\n";
Chris<<"20\n"<<y[i][j]<<"\n";
Chris<<"30\n"<<psi[i][j]<<"\n";
Chris<<"11\n"<<x[i][jplus1]<<"\n";
Chris<<"21\n"<<y[i][jplus1]<<"\n";
Chris<<"31\n"<<psi[i][jplus1]<<"\n";
Chris<<"12\n"<<x[i+1][jplus1]<<"\n";
Chris<<"22\n"<<y[i+1][jplus1]<<"\n";
Chris<<"32\n"<<psi[i+1][jplus1]<<"\n";
Chris<<"13\n"<<x[i+1][j]<<"\n";
Chris<<"23\n"<<y[i+1][j]<<"\n";
Chris<<"33\n"<<psi[i+1][j]<<"\n";
}
}

for(j=0;j<=n;j+=1)
{
jplus1=j+1;if(j==n)jplus1=0;
Maud<<"0\nLINE\n8\nCylinder\n";
Maud<<"10\n"<<x[0][j]<<"\n";
Maud<<"20\n"<<y[0][j]<<"\n";
Maud<<"11\n"<<x[0][jplus1]<<"\n";
Maud<<"21\n"<<y[0][jplus1]<<"\n";
}

for(i=0;i<=m-1;i+=1)
{
for(j=0;j<=n;j+=1)
{
jplus1=j+1;if(j==n)jplus1=0;
xplot[0]=x[i][j];yplot[0]=y[i][j];psiplot[0]=psi[i][j];
xplot[1]=x[i+1][j];yplot[1]=y[i+1][j];psiplot[1]=psi[i+1][j];
xplot[2]=x[i+1][jplus1];yplot[2]=y[i+1][jplus1];psiplot[2]=psi[i+1][jplus1];
Contour(ContInt,xplot,yplot,psiplot);
xplot[0]=x[i][j];yplot[0]=y[i][j];psiplot[0]=psi[i][j];
xplot[1]=x[i][jplus1];yplot[1]=y[i][jplus1];psiplot[1]=psi[i][jplus1];
xplot[2]=x[i+1][jplus1];yplot[2]=y[i+1][jplus1];psiplot[2]=psi[i+1][jplus1];
Contour(ContInt,xplot,yplot,psiplot);
}
}
Chris<<"0\nENDSEC\n0\nEOF\n";Chris.close();
Maud<<"0\nENDSEC\n0\nEOF\n";Maud.close();
cout<<"DXF files written, end of program\n";
return 0;
}

void Contour(float ContourInterval,float xTriangle[3],float yTriangle[3],float origpsiTraingle[3])
{
int Corner,CornerPlus1,CornerPlus2,intpsi,intpsiStart,intpsiStop;
float psiTriangle[3],Pointxy[2],psiStart,psiStop,thispsi;
for(Corner=0;Corner<=2;Corner++)psiTriangle[Corner]=origpsiTraingle[Corner]+0.5*ContourInterval;
for(Corner=0;Corner<=2;Corner++)
{
CornerPlus1=Corner+1;if(Corner+1>2)CornerPlus1-=3;
CornerPlus2=Corner+2;if(Corner+2>2)CornerPlus2-=3;
if((psiTriangle[CornerPlus1]-psiTriangle[Corner])*(psiTriangle[CornerPlus2]-psiTriangle[Corner])>
1.0e-6*ContourInterval*ContourInterval);
{
if(psiTriangle[CornerPlus1]-psiTriangle[Corner]>0.0)
{
psiStart=psiTriangle[Corner];
if(psiTriangle[CornerPlus1]>psiTriangle[CornerPlus2])psiStop=psiTriangle[CornerPlus2];
else psiStop=psiTriangle[CornerPlus1];
}
else
{
psiStop=psiTriangle[Corner];
if(psiTriangle[CornerPlus1]<psiTriangle[CornerPlus2])psiStart=psiTriangle[CornerPlus2];
else psiStart=psiTriangle[CornerPlus1];
}
}
intpsiStart=int(psiStart/ContourInterval);if(psiStart<0.0)intpsiStart-=1;
intpsiStop=int(psiStop/ContourInterval);  if(psiStop<0.0)intpsiStop-=1;
if(intpsiStart<intpsiStop)
{
for(intpsi=intpsiStart+1;intpsi<=intpsiStop;intpsi++)
{
thispsi=intpsi*ContourInterval;
Maud<<"0\nLINE\n8\npsi\n";
Pointxy[0]=xTriangle[Corner]+(xTriangle[CornerPlus1]-xTriangle[Corner])*(thispsi-psiTriangle[Corner])/(psiTriangle[CornerPlus1]-psiTriangle[Corner]);
Pointxy[1]=yTriangle[Corner]+(yTriangle[CornerPlus1]-yTriangle[Corner])*(thispsi-psiTriangle[Corner])/(psiTriangle[CornerPlus1]-psiTriangle[Corner]);
Maud<<"10\n"<<Pointxy[0]<<"\n";
Maud<<"20\n"<<Pointxy[1]<<"\n";
Pointxy[0]=xTriangle[Corner]+(xTriangle[CornerPlus2]-xTriangle[Corner])*(thispsi-psiTriangle[Corner])/(psiTriangle[CornerPlus2]-psiTriangle[Corner]);
Pointxy[1]=yTriangle[Corner]+(yTriangle[CornerPlus2]-yTriangle[Corner])*(thispsi-psiTriangle[Corner])/(psiTriangle[CornerPlus2]-psiTriangle[Corner]);
Maud<<"11\n"<<Pointxy[0]<<"\n";
Maud<<"21\n"<<Pointxy[1]<<"\n";
}
}
}
}