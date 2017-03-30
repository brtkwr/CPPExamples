#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define LastPerson 7
#define MaxCharactersPlus1 50
int  Age[LastPerson+1];
char Name[LastPerson+1][MaxCharactersPlus1];
int Swap(int);
void PrintThePeople(void);
int main(void)
{
	strcpy(Name[0],"Ludwig   ");Age[0]=62;strcpy(Name[1],"Bertrand ");Age[1]=97;
	strcpy(Name[2],"S¿ren    ");Age[2]=42;strcpy(Name[3],"Martin   ");Age[3]=86;
	strcpy(Name[4],"Jean-Paul");Age[4]=74;strcpy(Name[5],"Immanuel ");Age[5]=79;
	strcpy(Name[6],"Friedrich");Age[6]=55;strcpy(Name[7],"Blaise   ");Age[7]=39;
	int TotalAge=0;
	for(int Person=0;Person<=LastPerson;Person++)TotalAge+=Age[Person];
	cout<<"\nSum of the ages of all the people = "<<TotalAge<<"\n";
	PrintThePeople();
	for(int go=0;go<=LastPerson-1;go++)
	{
		for(int Person=0;Person<=LastPerson-1;Person++)
		{
			if(strcmp(Name[Person],Name[Person+1])>0)Swap(Person);
		}
	}
	PrintThePeople();
	for(int go=0;go<=LastPerson-1;go++)
	{
		for(int Person=0;Person<=LastPerson-1;Person++)
		{
			if(Age[Person]>Age[Person+1])Swap(Person);
		}
	}
	PrintThePeople();
	return 0;
}
int Swap(int ThisPerson)
{
	char TemporaryName[MaxCharactersPlus1];
	strcpy(TemporaryName,Name[ThisPerson]);
	strcpy(Name[ThisPerson],Name[ThisPerson+1]);
	strcpy(Name[ThisPerson+1],TemporaryName);
	int TemporaryAge=Age[ThisPerson];
	Age[ThisPerson]=Age[ThisPerson+1];
	Age[ThisPerson+1]=TemporaryAge;
}
void PrintThePeople(void)
{
	cout<<"\n";
	for(int Person=0;Person<=LastPerson;Person++)cout<<Name[Person]<<"\t"<<Age[Person]<<"\n";
}
