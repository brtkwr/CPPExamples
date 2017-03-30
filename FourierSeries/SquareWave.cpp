#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   LastPoint 5000

double x[LastPoint+1],y[LastPoint+1];

ofstream Chris("Fourier.dxf");

int main(void)
{
	double PI=4.0*atan(1.0);
	double c = 0.5;
	int LastTerm = 0;
	
	while(LastTerm < 1)
	{
		cout<<"How many terms do you want in the Fourier Series?\n";
		cin>>LastTerm;
		if(LastTerm<1)cout<<"The number of terms should be greater or equal to 1\n";
	}
	LastTerm-=1;
	
	Chris<<"0\nSECTION\n2\nENTITIES\n";
	
	for(int Point=0;Point<=LastPoint;Point++)
	{
		x[Point]=double(2*Point)/double(LastPoint);
		y[Point]=0.0;
	}
	double PlusOrMinus=1.0;
	for(int Term=0;Term<=LastTerm;Term++)
	{
		int n=2*Term+1;
		for(int Point=0;Point<=LastPoint;Point++)y[Point]+=PlusOrMinus*(4.0*c/PI)*cos(2.0*PI*double(n)*x[Point])/double(n);
		
		for(int Point=0;Point<=LastPoint-1;Point++)
		{
			Chris<<"0\nLINE\n8\n"<<Term<<"\n";
			Chris<<"10\n"<<x[Point]<<"\n";
			Chris<<"20\n"<<y[Point]<<"\n";
			Chris<<"11\n"<<x[Point+1]<<"\n";
			Chris<<"21\n"<<y[Point+1]<<"\n";
			Chris<<"62\n"<<Term<<"\n";
		}
		PlusOrMinus=-PlusOrMinus;
	}
	
	
	Chris<<"0\nLINE\n8\nAxes\n";
	Chris<<"10\n"<<0.0<<"\n";
	Chris<<"20\n"<<0.0<<"\n";
	Chris<<"11\n"<<2.0<<"\n";
	Chris<<"21\n"<<0.0<<"\n";
	Chris<<"62\n"<<0<<"\n";
	
	Chris<<"0\nLINE\n8\nAxes\n";
	Chris<<"10\n"<<0.0<<"\n";
	Chris<<"20\n"<<-0.6<<"\n";
	Chris<<"11\n"<<0.0<<"\n";
	Chris<<"21\n"<<0.6<<"\n";
	Chris<<"62\n"<<0<<"\n";
	
	for(int i=1;i<=20;i++)
	{
		Chris<<"0\nLINE\n8\nAxes\n";
		Chris<<"10\n"<<0.1*i<<"\n";
		Chris<<"20\n"<<0.0<<"\n";
		Chris<<"11\n"<<0.1*i<<"\n";
		Chris<<"21\n"<<-0.02<<"\n";
		Chris<<"62\n"<<0<<"\n";
	}
	
	for(int i=-6;i<=6;i++)
	{
		Chris<<"0\nLINE\n8\nAxes\n";
		Chris<<"10\n"<<0.0<<"\n";
		Chris<<"20\n"<<0.1*i<<"\n";
		Chris<<"11\n"<<-0.02<<"\n";
		Chris<<"21\n"<<0.1*i<<"\n";
		Chris<<"62\n"<<0<<"\n";
	}
	
	Chris<<"0\nENDSEC\n0\nEOF\n";
	Chris.close();
	cout<<"DXF file written, end of program\n";
	
	return 0;
}