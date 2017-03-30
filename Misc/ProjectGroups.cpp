#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define MaxLastPerson 100
#define MaxLastCharacter 256
#define MaxLastGroup 4

char   Text[MaxLastPerson+1][MaxLastCharacter],LineStart[4];
int PersonNumber[MaxLastPerson+1],LastPersonInGroup[MaxLastGroup+1],WhichGroup[MaxLastPerson+1];

ifstream Julia("Persons.txt");ofstream Madeleine("Groups.txt");
// Note this programm expects a file, Persons.txt, with a group of students, one per line. Students whose names begine with "Full" or "Sand"
// (short for Fulltime and Sandwich) are ignored. The list must end with a line "End".
int main(void)
{
	int LastGroup=4;if(LastGroup>MaxLastGroup){cout<<"Too many groups\n";return 0;}
	cout<<"Total number of groups is "<<LastGroup+1<<"\n";
	int LineNumber=-1;
	int PersonCounter=-1;
	do{
		LineNumber++;
		Julia.getline(Text[LineNumber],MaxLastCharacter,'\r');//USE '\n' for Windows Machine and '\r' for Mac OSX
			strncpy(LineStart,Text[LineNumber],4);
			if((strcmp(LineStart,"Full")!=0)&&(strcmp(LineStart,"Sand")!=0))
			{
				PersonCounter++;
				PersonNumber[LineNumber]=PersonCounter;
			}
			else PersonNumber[LineNumber]=-1;
	}
	while (strcmp(Text[LineNumber],"End")!=0);
	Julia.close();
	int LastLine=LineNumber-1;
	int LastPerson=PersonCounter;
	cout<<"The total number of people is "<<LastPerson+1<<"\n";
	for(int Group=0;Group<=LastGroup;Group++)LastPersonInGroup[Group]=-1;
	int CurrentLastNumberOfPersons=0;
	for(int Person=0;Person<=LastPerson;Person++)
	{
		int Tested=-1;
		for(int Group=0;Group<=LastGroup;Group++)
		{
			if(LastPersonInGroup[Group]==CurrentLastNumberOfPersons)Tested++;
		}
		if(Tested>=LastGroup)
		{
			CurrentLastNumberOfPersons++;
			cout<<"Maximum number of people in a group is now "<<CurrentLastNumberOfPersons+1<<"\n";
		}
		int Success=0;
		do
		{
			int TrialGroup=rand()%(LastGroup+1);
			if(LastPersonInGroup[TrialGroup]<CurrentLastNumberOfPersons)
			{
				WhichGroup[Person]=TrialGroup;
				LastPersonInGroup[TrialGroup]=CurrentLastNumberOfPersons;
				Success=1;
			}
		}while(Success==0);
	}
	
	for(int LineNumber=0;LineNumber<=LastLine;LineNumber++)
	{
		if(PersonNumber[LineNumber]!=-1)
			Madeleine<<Text[LineNumber]<<"\t"<<WhichGroup[PersonNumber[LineNumber]]<<"\n";
		else
			Madeleine<<Text[LineNumber]<<"\n";
	}
	Madeleine.close();
	return 0;
}
