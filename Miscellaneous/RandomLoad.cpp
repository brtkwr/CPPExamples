#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   maxmplus1 2001
#define   maxnplus1 2001

int   i,m,j,n;

float t[maxmplus1],p[maxmplus1],x[maxmplus1],
      omega[maxnplus1],amp[maxnplus1],phi[maxnplus1],ReX[maxnplus1],ImX[maxnplus1],
      PI,omegan,c,omegamin,omegamax,omegaav,omegahalfdif,pmean,pstd,xmean,xpow,xstd,sumsq,sum,
      totaltime,realbot,imagbot,magbotsq,top,bot;

ofstream Chris("RandomLoad.dxf");
void DrawLines(void);

int main(void)
{
PI=4.0*atan(1.0);
m=maxmplus1-1;
n=maxnplus1-1;
omegamin=0.0;omegamax=10.0;
omegaav=(omegamax+omegamin)/2.0;
omegahalfdif=(omegamax-omegamin)/2.0;
pmean=5.0;
pstd=2.0;
totaltime=50.0;
omegan=3.0;
c=0.05;

sumsq=0.0;
for(j=0;j<=n;j++)
{
omega[j]=omegamin+2.0*omegahalfdif*(1.0*j)/(1.0*n);
amp[j]=exp(-2.0*(omega[j]-omegaav)*(omega[j]-omegaav)/(omegahalfdif*omegahalfdif));
sumsq+=(1.0/2.0)*amp[j]*amp[j];
top=1.0*rand();bot=1.0*rand();top-=1.0*rand();bot-=1.0*rand();
if(bot!=0.0)phi[j]=2.0*atan(top/bot);else phi[j]=PI;
realbot=1.0-omega[j]*omega[j]/(omegan*omegan);
imagbot=2.0*c*omega[j]/omegan;
magbotsq=realbot*realbot+imagbot*imagbot;
ReX[j]=+realbot/magbotsq;
ImX[j]=-imagbot/magbotsq;
}
sum=sqrt(sumsq);
xpow=0.0;
for(j=0;j<=n;j++)
{
amp[j]=pstd*amp[j]/sum;
ReX[j]=pstd*ReX[j]/sum;
ImX[j]=pstd*ImX[j]/sum;
xpow+=(1.0/2.0)*(ReX[j]*ReX[j]+ImX[j]*ImX[j]);
}
xstd=sqrt(xpow);
xmean=pmean;
for(i=0;i<=m;i+=1)
{
t[i]=(totaltime*i)/(1.0*m);
p[i]=pmean;
x[i]=pmean;
for(j=0;j<=n;j++)
{
p[i]+=amp[j]*cos(omega[j]*t[i]+phi[j]);
x[i]+=ReX[j]*cos(omega[j]*t[i]+phi[j])-ImX[j]*sin(omega[j]*t[i]+phi[j]);
}
}

Chris<<"0\nSECTION\n2\nENTITIES\n";
DrawLines();
for(i=0;i<=m-1;i+=1)
{
Chris<<"0\nLINE\n8\nLoad\n";
Chris<<"10\n"<<t[i]  <<"\n20\n"<<p[i]  <<"\n30\n0.0\n";
Chris<<"11\n"<<t[i+1]<<"\n21\n"<<p[i+1]<<"\n31\n0.0\n";
Chris<<"62\n3\n";
Chris<<"0\nLINE\n8\nDisplacement\n";
Chris<<"10\n"<<t[i]  <<"\n20\n"<<x[i]  <<"\n30\n0.0\n";
Chris<<"11\n"<<t[i+1]<<"\n21\n"<<x[i+1]<<"\n31\n0.0\n";
Chris<<"62\n0\n";
}

Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();

cout<<"DXF file written, end of program\n";

return 0;
}

void DrawLines(void)
{
Chris<<"0\nLINE\n8\nLoad\n";
Chris<<"10\n0.0\n20\n0.0\n30\n0.0\n";
Chris<<"11\n"<<totaltime<<"\n21\n0.0\n31\n0.0\n";
Chris<<"0\nLINE\n8\nLoad\n";
Chris<<"10\n0.0\n20\n"<<pmean<<"\n30\n0.0\n";
Chris<<"11\n"<<totaltime<<"\n21\n"<<pmean<<"\n31\n0.0\n";
Chris<<"0\nLINE\n8\nLoad\n";
Chris<<"10\n0.0\n20\n"<<pmean+pstd<<"\n30\n0.0\n";
Chris<<"11\n"<<totaltime<<"\n21\n"<<pmean+pstd<<"\n31\n0.0\n";
Chris<<"0\nLINE\n8\nLoad\n";
Chris<<"10\n0.0\n20\n"<<pmean-pstd<<"\n30\n0.0\n";
Chris<<"11\n"<<totaltime<<"\n21\n"<<pmean-pstd<<"\n31\n0.0\n";
Chris<<"0\nLINE\n8\nLoad\n";
Chris<<"10\n0.0\n20\n0.0\n30\n0.0\n";
Chris<<"11\n0.0\n21\n"<<pmean+4.0*pstd<<"\n31\n0.0\n";

Chris<<"0\nLINE\n8\nDisplacement\n";
Chris<<"10\n0.0\n20\n0.0\n30\n0.0\n";
Chris<<"11\n"<<totaltime<<"\n21\n0.0\n31\n0.0\n";
Chris<<"0\nLINE\n8\nDisplacement\n";
Chris<<"10\n0.0\n20\n"<<xmean<<"\n30\n0.0\n";
Chris<<"11\n"<<totaltime<<"\n21\n"<<xmean<<"\n31\n0.0\n";
Chris<<"0\nLINE\n8\nDisplacement\n";
Chris<<"10\n0.0\n20\n"<<xmean+xstd<<"\n30\n0.0\n";
Chris<<"11\n"<<totaltime<<"\n21\n"<<xmean+xstd<<"\n31\n0.0\n";
Chris<<"0\nLINE\n8\nDisplacement\n";
Chris<<"10\n0.0\n20\n"<<xmean-xstd<<"\n30\n0.0\n";
Chris<<"11\n"<<totaltime<<"\n21\n"<<xmean-xstd<<"\n31\n0.0\n";
Chris<<"0\nLINE\n8\nDisplacement\n";
Chris<<"10\n0.0\n20\n0.0\n30\n0.0\n";
Chris<<"11\n0.0\n21\n"<<xmean+4.0*xstd<<"\n31\n0.0\n";
}