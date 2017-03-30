//Note on Macintosh run program once from 2DPinned.xcodeproj.
//This will create a folder 'build' and a folder 'Release' withing 'build'.
//Put the data file(s) in the folder Release and run again.
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include "size.h"

void Invert(int n,double A[][nplus1],double B[][nplus1]);
void Mult(int n,double A[][nplus1],double B[][nplus1],double C[][nplus1]);
ifstream Chris("ArrayP.txt");

int main (void)
{
double P[nplus1][nplus1],Q[nplus1][nplus1],R[nplus1][nplus1];
int i,j,n;

Chris>>n;n=n-1;
if(n>nplus1-1){cout<<"Too big\n";return 0;}

for(i=0;i<=n;i++){for(j=0;j<=n;j++)Chris>>P[i][j];}
Invert(n,P,Q);
Mult(n,P,Q,R);
for(i=0;i<=n;i++){for(j=0;j<=n;j++)cout<<R[i][j]<<"  ";cout<<"\n";}cout<<"\n";
    cout<<"type any number and then press enter to finish\n";//Only required for PC
    float something;cin>>something;//Only required for PC
   //system("PAUSE");//Could use this for PC, but not part of C++ language
   return 0;
}
