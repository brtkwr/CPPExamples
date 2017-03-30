#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   MaxNo 10000


int   i,n,Number[MaxNo],SortedNumber[MaxNo],
      Crossout[MaxNo],SoFar,CurrentBiggest;

ifstream Julia("Numbers.txt");
ofstream Maud("SortedNumbers.txt");

int main(void)
{
Julia>>n;
if(n>MaxNo){cout<<"n too big\n";return 0;}

for(i=0;i<=n-1;i++){Julia>>Number[i];Crossout[i]=0;}
Julia.close();

SoFar=0;
for(;;)
{
for(i=0;i<=n-1;i++)
{
CurrentBiggest=i;
if(Crossout[i]==0)break;
}
for(i=0;i<=n-1;i++)
{
if(Number[i]>Number[CurrentBiggest]&&Crossout[i]==0)
CurrentBiggest=i;
}
Crossout[CurrentBiggest]=1;
SortedNumber[SoFar]=Number[CurrentBiggest];
SoFar+=1;
if(SoFar>n-1)break;
}

Maud<<n<<"\n";
for(i=0;i<=n-1;i++)Maud<<SortedNumber[i]<<"\n";

Maud.close();

cout<<"Numbers have been sorted\n";

return 0;
}