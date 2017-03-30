#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define NumPrimes 1001

int Prime[NumPrimes],Number,integerFraction,whichPrime,PrimesSoFar,TempPrimesSoFar;
double  Fraction;

int main(void)
{
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
if(PrimesSoFar>=NumPrimes-1)break;
}

return 0;
}