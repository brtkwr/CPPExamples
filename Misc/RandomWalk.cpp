#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
ofstream Chris("RandomWalk.dxf");
double a[3],V[3],x[3],oldx[3];
int main(void)
{
double weighting;
cout<<"What is the weighting for previous velocity?\n";
cin>>weighting;
double asq;
for(int xyz=0;xyz<=2;xyz++)
{
V[xyz]=0.0;
x[xyz]=0.0;
oldx[xyz]=0.0;
}
Chris<<"0\nSECTION\n2\nENTITIES\n";
for(int i=0;i<=100000;i++)
{
div_t divresult;
divresult = div (i,10000);
if(divresult.rem==0)cout<<divresult.quot<<"\n";
for(;;)
{
for(int xyz=0;xyz<=2;xyz++)a[xyz]=2.0*double(rand())/double(RAND_MAX)-1.0;
asq=0.0;
for(int xyz=0;xyz<=2;xyz++)asq+=a[xyz]*a[xyz];
if(asq<=1.0)break;
}
double moda=sqrt(asq);
for(int xyz=0;xyz<=2;xyz++)V[xyz]=weighting*V[xyz]+a[xyz]/moda;
double Vsq=0.0;
for(int xyz=0;xyz<=2;xyz++)Vsq+=V[xyz]*V[xyz];
double modV=sqrt(Vsq);
for(int xyz=0;xyz<=2;xyz++)
{
V[xyz]=V[xyz]/modV;
x[xyz]+=V[xyz];
}
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<oldx[0]<<"\n";
Chris<<"20\n"<<oldx[1]<<"\n";
Chris<<"30\n"<<oldx[2]<<"\n";
Chris<<"11\n"<<x[0]<<"\n";
Chris<<"21\n"<<x[1]<<"\n";
Chris<<"31\n"<<x[2]<<"\n";
for(int xyz=0;xyz<=2;xyz++)oldx[xyz]=x[xyz];
}
Chris<<"0\nENDSEC\n0\nEOF\n";
cout<<"DXF file written, end of program\n";
return 0;
}
