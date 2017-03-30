#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include <math.h>
#define   m 500
ofstream Julia("Parabola.dxf");
int main(void)
{
double xmax=130000.0;
double ymax=75000.0;
double tmax=2.0*xmax/ymax;
double a=xmax/(tmax*tmax);
Julia<<"0\nSECTION\n2\nENTITIES";
Julia<<"\n0\nPOLYLINE";
Julia<<"\n100\nAcDbEntity\n8\nParabola";
Julia<<"\n100\nAcDb3dPolyline\n66\n1";
Julia<<"\n10\n0\n20\n0\n30\n0.0\n70\n8";//For open polyline
//Julia<<"\n10\n0\n20\n0\n30\n0.0\n70\n9";//For closed polyline
for(int i=0;i<=m;i+=1){
double t=tmax*(2.0*i-1.0*m)/(1.0*m);
Julia<<"\n0\nVERTEX";
Julia<<"\n100\nAcDbEntity\n8\nParabola";
Julia<<"\n100\nAcDb3dPolylineVertex";
Julia<<"\n10\n"<<a*t*t;
Julia<<"\n20\n"<<2.0*a*t;
Julia<<"\n30\n"<<0.0;}
Julia<<"\n0\nSEQEND";
Julia<<"\n0\nENDSEC\n0\nEOF\n";
cout<<"DXF file written, end of program\n";
return 0;
}