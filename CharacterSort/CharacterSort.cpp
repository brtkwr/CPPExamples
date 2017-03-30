#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define MaxPeoplePlus1 100
#define MaxCharactersPlus1 50

char   Name[MaxPeoplePlus1][MaxCharactersPlus1],temp[MaxCharactersPlus1];
int i,n,go;

ifstream Julia("Names.txt");

int main(void)
{
Julia>>n;if(n>MaxPeoplePlus1-1){cout<<"n too big\n";return 0;}

//for(i=0;i<=n;i++)Julia.getline(Name[i],MaxCharactersPlus1,'\n');//USE FOR MOST COMPUTERS
for(i=0;i<=n;i++)Julia.getline(Name[i],MaxCharactersPlus1,'\r');//USE THIS FOR OS X ON MAC
Julia.close();

cout<<"\n";
for(i=1;i<=n;i++)cout<<i<<"  "<<Name[i]<<"\n";
cout<<"\n";

for(go=0;go<=n-1;go++)
{
for(i=1;i<=n-1;i++)
{
if(strcmp(Name[i],Name[i+1])>0)
{
strcpy(temp,Name[i]);
strcpy(Name[i],Name[i+1]);
strcpy(Name[i+1],temp);
}
}
}
for(i=1;i<=n;i++)cout<<i<<"  "<<Name[i]<<"\n";
cout<<"\n";
cout<<"type any number and then press enter to finish\n";//Only required for PC
    float something;cin>>something;//Only required for PC
   //system("PAUSE");//Could use this for PC, but not part of C++ language
return 0;
}
