#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   NoSegments 50000

void DXFSetUp(void);
void DXFLines(void);
void DXFText(void);
void DXFFinishOff(void);

double t[NoSegments+1],f[NoSegments+1],PI,
DXFx1,DXFy1,DXFz1,
DXFx2,DXFy2,DXFz2,
DXFtextx,DXFtexty,DXFtextz,
textsize;

int Layer=0,Colour=0;

ofstream Chris("Fourier.dxf");

int main(void)
{
	PI=4.0*atan(1.0);
	double Beta=0.1;
	double T=1.0;
	int NumberOfTerms;
	for(;;)
	{
		cout<<"How many terms do you want in the Fourier Series?\n";
		cin>>NumberOfTerms;
		if(NumberOfTerms>0)break;
		cout<<"The number of terms should be positive\n";
	}
	
	for(int i=0;i<=NoSegments;i++)
	{
		t[i]=double(2*i)*T/double(NoSegments);
		f[i]=4.0*Beta/3.0;
		for(int n=1;n<=NumberOfTerms;n++)
		{
			double ro=2.0*PI*double(n)*Beta;
            f[i]+=2.0*(4.0*Beta/(ro*ro))*(sin(ro)/ro-cos(ro))*cos(2.0*PI*double(n)*t[i]/T);
		}
	}
	
	DXFSetUp();
	
	for(int i=1;i<=NoSegments;i++)
	{
		DXFx1=t[i-1];
		DXFy1=f[i-1];
		DXFz1=0.0;
		DXFx2=t[i];
		DXFy2=f[i];
		DXFz2=0.0;
		DXFLines();
	}
	
	DXFx1=0.0;DXFy1=0.0;DXFz1=0.0;
	DXFx2=2.1;DXFy2=0.0;DXFz2=0.0;
	DXFLines();
	
	DXFy1=0.0;
	DXFy2=-0.02;
	for(int i=1;i<=21;i++){DXFx1=0.1*i;DXFx2=DXFx1;DXFLines();}
	
	DXFx1=0.0;DXFy1=0.0;
	DXFx2=0.0;DXFy2=1.2;
	DXFLines();
	
	DXFx1=0.0;
	DXFx2=-0.02;
	for(int i=0;i<=12;i++){DXFy1=0.1*i;DXFy2=DXFy1;DXFLines();}
	
	Layer=1;Colour=1;textsize=0.1;
	DXFtextx=0.0;
	DXFtexty=-0.8;
	DXFtextz=0.0;
	DXFText();
	Chris<<"Fourier Series\n";
	
	DXFFinishOff();
	cout<<"DXF file written, end of program\n";
	
	return 0;
}

void DXFSetUp(void)
{
	Chris<<"0\n";
	Chris<<"SECTION\n";
	Chris<<"2\n";
	Chris<<"ENTITIES\n";
}

void DXFLines(void)
{
	Chris<<"0\nLINE\n8\n"<<Layer<<"\n";
	Chris<<"10\n"<<DXFx1<<"\n";
	Chris<<"20\n"<<DXFy1<<"\n";
	Chris<<"30\n"<<DXFz1<<"\n";
	Chris<<"11\n"<<DXFx2<<"\n";
	Chris<<"21\n"<<DXFy2<<"\n";
	Chris<<"31\n"<<DXFz2<<"\n";
	Chris<<"62\n"<<Colour<<"\n";
}

void DXFText(void)
{
	Chris<<"0\nTEXT\n8\n"<<Layer<<"\n";
	Chris<<"10\n"<<DXFtextx<<"\n";
	Chris<<"20\n"<<DXFtexty<<"\n";
	Chris<<"30\n"<<DXFtextz<<"\n";
	Chris<<"40\n"<<textsize<<"\n";
	Chris<<"62\n"<<Colour<<"\n";
	Chris<<"1\n";
}

void DXFFinishOff(void)
{
	Chris<<"0\n";
	Chris<<"ENDSEC\n";
	Chris<<"0\n";
	Chris<<"EOF\n";
	Chris.close();
}
