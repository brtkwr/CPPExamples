//Note on Macintosh run program once from 2DPinned.xcodeproj.
//This will create a folder 'build' and a folder 'Release' withing 'build'.
//Put the data file(s) in the folder Release and run again.
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include "size.h"

double A[nplus1][nplus1],x[nplus1][nplus1],Lambda[nplus1];
int i,j,n,eigen;

void Eigen(int n,double A[][nplus1],double x[][nplus1],double Lambda[]);
ifstream Chris("MyMatrix.txt");

int main (void)
{

Chris>>n;n=n-1;
if(n>nplus1-1){cout<<"Problem too big\n";return 0;}


for(i=0;i<=n;i++){for(j=i;j<=n;j++)Chris>>A[i][j];}
for(i=1;i<=n;i++){for(j=0;j<=i-1;j++)A[i][j]=A[j][i];}
for(i=0;i<=n;i++){for(j=0;j<=n;j++)cout<<A[i][j]<<"  ";cout<<"\n";}

Eigen(n,A,x,Lambda);

for(eigen=0;eigen<=n;eigen++)
{
cout<<Lambda[eigen]<<"\n";for(i=0;i<=n;i++)cout<<x[i][eigen]<<"  ";cout<<"\n";
}
cout<<"type any number and then press enter to finish\n";//Only required for PC
    float something;cin>>something;//Only required for PC
   //system("PAUSE");//Could use this for PC, but not part of C++ language
return 0;
}

