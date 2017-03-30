#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   LastPoint 5000

double t[LastPoint+1],i[LastPoint+1];

ofstream Chris("Fourier.dxf");

int main(void)
{
	double PI=4.0*atan(1.0);
	double sigma = 0.4;
	double T = 1.0;
	
	int n_min = 4;
	int n_max = 40;
	
	int LastTerm = n_max - n_min;
	double c = sigma * sqrt(2.0 /(double(1 + LastTerm)));
	
	Chris<<"0\nSECTION\n2\nENTITIES\n";
	
	for(int Point=0;Point<=LastPoint;Point++)
	{
		t[Point]=double(2*Point)/double(LastPoint);
		i[Point]=0.0;
	}
	for(int Term=0;Term<=LastTerm;Term++)
	{
	double phi = 2.0 * PI * double(rand())/double(RAND_MAX);
	int n = n_min + Term;
	for(int Point=0;Point<=LastPoint;Point++)i[Point]+=c*cos(2.0*PI*double(n)*t[Point]/T + phi);
		}
		
		for(int Point=0;Point<=LastPoint-1;Point++)
		{
			Chris<<"0\nLINE\n8\nSignal\n";
			Chris<<"10\n"<<t[Point]<<"\n";
			Chris<<"20\n"<<i[Point]<<"\n";
			Chris<<"11\n"<<t[Point+1]<<"\n";
			Chris<<"21\n"<<i[Point+1]<<"\n";
		}
		
		Chris<<"0\nLINE\n8\nrms\n";
			Chris<<"10\n"<<0.0<<"\n";
			Chris<<"20\n"<<sigma<<"\n";
			Chris<<"11\n"<<2.0<<"\n";
			Chris<<"21\n"<<sigma<<"\n";
			Chris<<"62\n1\n";
			
			Chris<<"0\nLINE\n8\nrms\n";
			Chris<<"10\n"<<0.0<<"\n";
			Chris<<"20\n"<<-sigma<<"\n";
			Chris<<"11\n"<<2.0<<"\n";
			Chris<<"21\n"<<-sigma<<"\n";
			Chris<<"62\n1\n";
	
	Chris<<"0\nLINE\n8\nAxes\n";
	Chris<<"10\n"<<0.0<<"\n";
	Chris<<"20\n"<<0.0<<"\n";
	Chris<<"11\n"<<2.0<<"\n";
	Chris<<"21\n"<<0.0<<"\n";
	Chris<<"62\n"<<0<<"\n";
	
	Chris<<"0\nLINE\n8\nAxes\n";
	Chris<<"10\n"<<0.0<<"\n";
	Chris<<"20\n"<<-3.0*sigma<<"\n";
	Chris<<"11\n"<<0.0<<"\n";
	Chris<<"21\n"<<3.0*sigma<<"\n";
	Chris<<"62\n"<<0<<"\n";
	
	Chris<<"0\nENDSEC\n0\nEOF\n";
	Chris.close();
	cout<<"DXF file written, end of program\n";
	
	return 0;
}