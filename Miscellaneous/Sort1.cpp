#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define   MaxNo 10000


int   i,n,Number[MaxNo],temp,go;

ifstream Julia("Numbers.txt");
ofstream Maud("SortedNumbers.txt");

int main(void)
{
Julia>>n;
if(n>MaxNo){cout<<"n too big\n";return 0;}

for(i=0;i<=n-1;i++)Julia>>Number[i];
Julia.close();

for(go=0;go<=n-2;go++)
{
for(i=0;i<=n-2;i++)
{
if(Number[i]<Number[i+1])
{
temp=Number[i];Number[i]=Number[i+1];Number[i+1]=temp;
}
}
}

Maud<<n<<"\n";
for(i=0;i<=n-1;i++)Maud<<Number[i]<<"\n";

Maud.close();

cout<<"Numbers have been sorted\n";

return 0;
}