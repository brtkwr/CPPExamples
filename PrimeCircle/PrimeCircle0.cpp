#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   MaxNo 100
#define NumPrimesPlus1 13

void PrimeNumbers(void);
void CalculateDraw(void);

int   Prime[NumPrimesPlus1],
      whichPrime,NumPrimes,
	  i,j;

double x[MaxNo],y[MaxNo],PI,
       r,theta,thetaStart;

ofstream Chris("Circle.dxf");

int main(void)
{
NumPrimes=NumPrimesPlus1-1;
PrimeNumbers();

PI=4.0*atan(1.0);

Chris<<"0\nSECTION\n2\nENTITIES\n";
r=100000.0;
thetaStart=-PI/2.0;
for(whichPrime=2;whichPrime<=NumPrimes;whichPrime++)
{
if(Prime[whichPrime]>MaxNo){cout<<"Prime too big\n";return 0;}
CalculateDraw();
}

r=1.2*r;

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<+r<<"\n";Chris<<"20\n"<<+r<<"\n";
Chris<<"11\n"<<-r<<"\n";Chris<<"21\n"<<+r<<"\n";
Chris<<"62\n"<<0<<"\n";

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<-r<<"\n";Chris<<"20\n"<<+r<<"\n";
Chris<<"11\n"<<-r<<"\n";Chris<<"21\n"<<-r<<"\n";
Chris<<"62\n"<<0<<"\n";

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<-r<<"\n";Chris<<"20\n"<<-r<<"\n";
Chris<<"11\n"<<+r<<"\n";Chris<<"21\n"<<-r<<"\n";
Chris<<"62\n"<<0<<"\n";

Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<+r<<"\n";Chris<<"20\n"<<-r<<"\n";
Chris<<"11\n"<<+r<<"\n";Chris<<"21\n"<<+r<<"\n";
Chris<<"62\n"<<0<<"\n";

Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();
cout<<"DXF file written, end of program\n";

return 0;
}

void PrimeNumbers(void)
{
int Number,integerFraction,PrimesSoFar,TempPrimesSoFar;
double  Fraction;

Prime[0]=1;
PrimesSoFar=1;Prime[1]=2;
cout<<"1  2\n";
for(Number=3;;Number++)
{
for(whichPrime=1;whichPrime<=PrimesSoFar;whichPrime+=1)
{
Fraction=(1.0*Number)/(1.0*Prime[whichPrime]);//cout<<Fraction<<"\n";
integerFraction=int(Fraction);
if(Fraction==1.0*integerFraction)break;
if(whichPrime==PrimesSoFar||Prime[whichPrime]*Prime[whichPrime]>Number)
{
TempPrimesSoFar=PrimesSoFar+1;
Prime[TempPrimesSoFar]=Number;
cout<<TempPrimesSoFar<<"  "<<Number<<"\n";
}
if(Prime[whichPrime]*Prime[whichPrime]>Number)break;
}
PrimesSoFar=TempPrimesSoFar;
if(PrimesSoFar>=NumPrimes)break;
}
}

void CalculateDraw(void)
{
for(i=0;i<=Prime[whichPrime]-1;i++)
{
theta=(2.0*PI*i)/(1.0*Prime[whichPrime])+thetaStart;
x[i]=r*cos(theta);
y[i]=r*sin(theta);
}

for(i=0;i<=Prime[whichPrime]-2;i++)
{
for(j=i+1;j<=Prime[whichPrime]-1;j++)
{
Chris<<"0\nLINE\n8\n"<<Prime[whichPrime]<<"\n";
Chris<<"10\n"<<x[i]<<"\n";
Chris<<"20\n"<<y[i]<<"\n";
Chris<<"11\n"<<x[j]<<"\n";
Chris<<"21\n"<<y[j]<<"\n";
Chris<<"62\n"<<0<<"\n";
}
}
}