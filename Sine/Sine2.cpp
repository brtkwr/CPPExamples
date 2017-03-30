#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

int angle;
double PI,ChrisAngle,value;

ofstream Chris("Sine.txt");

int main()
{
PI=4.0*atan(1.0);
for(angle=0;angle<=90;angle++)
{
ChrisAngle=PI*angle/180.0;
value=sin(ChrisAngle);
Chris<<angle<<"  "<<value<<"\n";
}
cout<<"Values writen to file Sine.txt\n";
return 0;
}
