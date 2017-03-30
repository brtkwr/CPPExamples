#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define   m 500
int   i,NoTerms,n;
double t[m+1],f[m+1],
       PI,T,Amp,plusorminus;
ofstream Chris("Fourier.dxf");
int main(void)
{
PI=4.0*atan(1.0);
T=2.0;
Amp=1.0;
for(i=0;i<=m;i+=1)
{
t[i]=(2.0*T*i)/(1.0*m)-T;
f[i]=0.0;
}
NoTerms=10;
plusorminus=-1.0;
Chris<<"0\nSECTION\n2\nENTITIES\n";
for(n=1;n<=NoTerms;n+=1)
{
plusorminus=-plusorminus;
for(i=0;i<=m;i+=1)
f[i]+=plusorminus*8.0*Amp*sin(2.0*PI*(2.0*n-1.0)*t[i]/T)/(PI*PI*(2.0*n-1.0)*(2.0*n-1.0));
for(i=1;i<=m-1;i+=1)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<t[i]<<"\n";
Chris<<"20\n"<<f[i]<<"\n";
Chris<<"11\n"<<t[i+1]<<"\n";
Chris<<"21\n"<<f[i+1]<<"\n";
Chris<<"62\n"<<n<<"\n";
}
}
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<-T<<"\n";
Chris<<"20\n0.0\n";
Chris<<"11\n"<<+T<<"\n";
Chris<<"21\n0.0\n";
Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();
cout<<"DXF file written, end of program\n";
return 0;
}
