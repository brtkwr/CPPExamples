#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define   NoSegments 200
int   i,j,NoTerms,m,n;
double x[NoSegments+1][NoSegments+1],y[NoSegments+1][NoSegments+1],z[NoSegments+1][NoSegments+1],PI,
       w,amn,a,L;
ofstream Chris("PlateLoad.dxf");
int main(void)
{
PI=4.0*atan(1.0);
a=0.08;
L=1.0;
w=-0.1*(-1.0/(1.0-L*L/(a*a)));
for(;;)
{
cout<<"How many terms do you want in the Fourier Series?\n";
cin>>NoTerms;
if(NoTerms>0)break;
cout<<"The number of terms should be positive\n";
}

for(i=0;i<=NoSegments;i+=1)
{
for(j=0;j<=NoSegments;j+=1)
{
x[i][j]=(3.0*i*L)/(1.0*NoSegments);
y[i][j]=(3.0*j*L)/(1.0*NoSegments);
z[i][j]=0.0;
for(m=1;m<=NoTerms;m++)
{
z[i][j]-=(2.0*w*L/(1.0*m*PI*a))*sin(1.0*m*PI*a/L)*(cos(2.0*m*PI*x[i][j]/L)+cos(2.0*m*PI*y[i][j]/L));
for(n=1;n<=NoTerms;n++)
z[i][j]-=(4.0*w*L*L/(1.0*m*n*PI*PI*a*a))*sin(1.0*m*PI*a/L)*sin(1.0*n*PI*a/L)
                     *cos(2.0*m*PI*x[i][j]/L)*cos(2.0*n*PI*y[i][j]/L);
}
}
}

Chris<<"0\n";
Chris<<"SECTION\n";
Chris<<"2\n";
Chris<<"ENTITIES\n";
for(i=1;i<=NoSegments;i+=1)
{
for(j=1;j<=NoSegments;j+=1)
{
Chris<<"0\n3DFACE\n8\nPlate Load\n";
Chris<<"10\n"<<x[i-1][j-1]<<"\n";
Chris<<"20\n"<<y[j-1][j-1]<<"\n";
Chris<<"30\n"<<z[i-1][j-1]<<"\n";
Chris<<"11\n"<<x[i][j-1]<<"\n";
Chris<<"21\n"<<y[i][j-1]<<"\n";
Chris<<"31\n"<<z[i][j-1]<<"\n";
Chris<<"12\n"<<x[i][j]<<"\n";
Chris<<"22\n"<<y[i][j]<<"\n";
Chris<<"32\n"<<z[i][j]<<"\n";
Chris<<"13\n"<<x[i-1][j]<<"\n";
Chris<<"23\n"<<y[i-1][j]<<"\n";
Chris<<"33\n"<<z[i-1][j]<<"\n";
}
}
Chris<<"0\n";
Chris<<"ENDSEC\n";
Chris<<"0\n";
Chris<<"EOF\n";
Chris.close();
cout<<"DXF file written, end of program\n";
return 0;
}