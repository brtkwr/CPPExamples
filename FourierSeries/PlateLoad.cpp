#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   NoSegments 100

void DXFSetUp(void);
void DXF3DFACE(void);
void DXFFinishOff(void);

int   i,j,NoTerms,n;

double x[NoSegments+1],y[NoSegments+1],z[NoSegments+1],PI,
       DXFx1,DXFy1,DXFz1,
       DXFx2,DXFy2,DXFz2,
       DXFx3,DXFy3,DXFz3,
       DXFx4,DXFy4,DXFz4,
       DXFtextx,DXFtexty,DXFtextz,
       textsize;

ofstream Chris("PlateLoad.dxf");

int main(void)
{
PI=4.0*atan(1.0);

for(;;)
{
cout<<"How many terms do you want in the Fourier Series?\n";
cin>>NoTerms;
if(NoTerms>0)break;
cout<<"The number of terms should be positive\n";
}

for(i=0;i<=NoSegments;i+=1)
{
x[i]=(1.0*i)/(1.0*NoSegments);
y[i]=x[i];
z[i]=0.0;
for(n=1;n<=NoTerms;n+=1)
z[i]+=(2.0/PI)*sin(2.0*PI*(2.0*n-1.0)*x[i])/(2.0*n-1.0);
}

DXFSetUp();
for(i=1;i<=NoSegments;i+=1)
{
for(j=1;j<=NoSegments;j+=1)
{
DXFx1=x[i-1];
DXFy1=y[j-1];
DXFz1=z[i-1]*z[j-1];
DXFx2=x[i];
DXFy2=y[j-1];
DXFz2=z[i]*z[j-1];
DXFx3=x[i];
DXFy3=y[j];
DXFz3=z[i]*z[j];
DXFx4=x[i-1];
DXFy4=y[j];
DXFz4=z[i-1]*z[j];
DXF3DFACE();
}
}
DXFFinishOff();
cout<<"DXF file written, end of program\n";

return 0;
}

void DXFSetUp(void)
{
Chris<<"0\n";
Chris<<"SECTION\n";
Chris<<"2\n";
Chris<<"ENTITIES\n";
}

void DXF3DFACE(void)
{
Chris<<"0\n3DFACE\n8\nPlate Load\n";
Chris<<"10\n"<<DXFx1<<"\n";
Chris<<"20\n"<<DXFy1<<"\n";
Chris<<"30\n"<<DXFz1<<"\n";
Chris<<"11\n"<<DXFx2<<"\n";
Chris<<"21\n"<<DXFy2<<"\n";
Chris<<"31\n"<<DXFz2<<"\n";
Chris<<"12\n"<<DXFx3<<"\n";
Chris<<"22\n"<<DXFy3<<"\n";
Chris<<"32\n"<<DXFz3<<"\n";
Chris<<"13\n"<<DXFx4<<"\n";
Chris<<"23\n"<<DXFy4<<"\n";
Chris<<"33\n"<<DXFz4<<"\n";
}

void DXFFinishOff(void)
{
Chris<<"0\n";
Chris<<"ENDSEC\n";
Chris<<"0\n";
Chris<<"EOF\n";
Chris.close();
}
