#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
ofstream Chris("MyLines.dxf");
double ChrisDXFline(double xThisEnd,double yThisEnd,double xThatEnd,double yThatEnd);
int main(void)
{
int   Line,LastLine;
double theta,a,thetamax,x,y,length,TotalLength;
a=1000.0;
thetamax=1.0;
cout<<"How many lines do you want?\n";
cin>>LastLine;LastLine-=1;
TotalLength=0.0;
Chris<<"0\nSECTION\n2\nENTITIES\n";
for(Line=0;Line<=LastLine;Line++)
{
theta=thetamax*(2.0*Line-1.0*LastLine)/(1.0*LastLine);
x=a*exp(theta);
y=a*exp(-theta);
length=ChrisDXFline(x,0.0,0.0,y);
TotalLength+=length;
}
Chris<<"0\nENDSEC\n0\nEOF\n";
Chris.close();
cout<<"Total length of lines = "<<TotalLength<<"\n";
cout<<"DXF file written, end of program\n";
system("PAUSE");//Only for a PC!
return 0;
}
double ChrisDXFline(double xThisEnd,double yThisEnd,double xThatEnd,double yThatEnd)
{
Chris<<"0\nLINE\n8\nChris's lines\n";
Chris<<"10\n"<<xThisEnd<<"\n";
Chris<<"20\n"<<yThatEnd<<"\n";
Chris<<"11\n"<<xThatEnd<<"\n";
Chris<<"21\n"<<yThatEnd<<"\n";
return sqrt((xThatEnd-xThisEnd)*(xThatEnd-xThisEnd)+(yThatEnd-yThisEnd)*(yThatEnd-yThisEnd));
}
