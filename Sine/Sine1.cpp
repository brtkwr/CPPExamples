#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

double angle,PI;

int main()
{
PI=4.0*atan(1.0);
cout<<"Type an angle in degrees \n";
cin>>angle;
cout<<"Sine of angle is "<<sin(angle*PI/180.0)<<"\n";
return 0;
}
