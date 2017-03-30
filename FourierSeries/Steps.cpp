#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define   m 5000
double t[m+1],f[m+1],
       PI,tlimit;
ofstream Chris("Fourier.dxf");
int main(void)
{
PI=4.0*atan(1.0);
tlimit=2.2;
for(int i=0;i<=m;i++)
{
t[i]=(2.0*tlimit*i)/(1.0*m)-tlimit;
f[i]=0.0;
}
int NoTerms=50;
for(int i=0;i<=m;i++)
{
for(int n=-NoTerms;n<=NoTerms;n++)
{
if(n!=0)
f[i]+=(1.0/(2.0*PI*double(n)))*sin(2.0*PI*double(n)*t[i]);
else
f[i]+=t[i];
}
}
Chris<<"0\nSECTION\n2\nENTITIES\n";
for(int i=0;i<=m-1;i++)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<t[i]<<"\n";
Chris<<"20\n"<<f[i]<<"\n";
Chris<<"11\n"<<t[i+1]<<"\n";
Chris<<"21\n"<<f[i+1]<<"\n";
}
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<-tlimit<<"\n";
Chris<<"20\n0.0\n";
Chris<<"11\n"<<+tlimit<<"\n";
Chris<<"21\n0.0\n";

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n0.0\n";
Chris<<"20\n"<<-tlimit-0.5<<"\n";
Chris<<"11\n0.0\n";
Chris<<"21\n"<<tlimit+0.5<<"\n";

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n1.0\n";
Chris<<"20\n-0.2\n";
Chris<<"11\n1.0\n";
Chris<<"21\n0.2\n";

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n-0.2\n";
Chris<<"20\n0.5\n";
Chris<<"11\n0.2\n";
Chris<<"21\n0.5\n";

Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();
cout<<"DXF file written, end of program\n";
return 0;
}
