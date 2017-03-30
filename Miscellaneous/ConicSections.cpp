#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include <math.h>
ofstream Julia("ConicSections.dxf");
int main(void)
{
	double PI=4.0*atan(1.0);
	int m=500;int n=50;
	double f=1000.0;
	double small=1.0e-3;
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	for(int j=0;j<=n;j++)
	{
		double e=double(2*j)/double(n);
		int type=2;
		if(e<1.0-small)type=0;//Ellipse
			if(e>1.0+small)type=1;//Hyperbola
				Julia<<"0\nPOLYLINE\n";
				Julia<<"100\nAcDbEntity\n8\nEllipses\n";
				Julia<<"100\nAcDb3dPolyline\n66\n1\n";
				Julia<<"10\n0\n20\n0\n30\n0.0\n70\n8\n";
				Julia<<"62\n"<<type<<"\n";
				for(int i=0;i<=m;i+=1)
				{
					Julia<<"0\nVERTEX\n";
					Julia<<"100\nAcDbEntity\n8\nHyperbola\n";
					Julia<<"100\nAcDb3dPolylineVertex\n";
					if(type==0)
					{
						double theta=2.0*PI*double(i)/double(m);
						Julia<<"10\n"<<(f/(1.0-e))*(cos(theta)+1.0)<<"\n";
						Julia<<"20\n"<<(f*sqrt(1.0-e*e)/(1.0-e))*sin(theta)<<"\n";
					}
					else
					{
						double zeta=12.0*(double(i)/double(m)-0.5)/(e*e);
						if(type==1)
						{
							Julia<<"10\n"<<(f/(e-1.0))*(cosh(zeta*sqrt(e*e-1.0))-1.0)<<"\n";
							Julia<<"20\n"<<(f*sqrt(e*e-1.0)/(e-1.0))*sinh(zeta*sqrt(e*e-1.0))<<"\n";
						}
						else
						{
							Julia<<"10\n"<<f*zeta*zeta<<"\n";
							Julia<<"20\n"<<2*f*zeta<<"\n";
						}
					}
					Julia<<"30\n0.0\n";
					Julia<<"62\n"<<type<<"\n";
				}
				Julia<<"0\nSEQEND\n";
	}
	Julia<<"0\nENDSEC\n0\nEOF\n\n";
	cout<<"DXF file written, end of program\n\n";
	return 0;
}

