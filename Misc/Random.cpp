#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;

#define   MaxNo 1000

int   i;
double x,xMin,xMax;

ofstream Madeleine("Numbers.txt");

int main(void)
{
srand(int(time(NULL)));
for(i=0;i<=MaxNo;i++)
{
x=(1.0*rand())/(1.0*RAND_MAX);
if(xMax<x||i==0)xMax=x;
if(xMin>x||i==0)xMin=x;
Madeleine<<x<<"\n";
}

Madeleine.close();

cout<<"Largest value = "<<xMax<<"\n";
cout<<"Smallest value = "<<xMin<<"\n";
return 0;
}