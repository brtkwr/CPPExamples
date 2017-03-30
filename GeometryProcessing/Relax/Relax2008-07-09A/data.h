//Written by Chris J K Williams, University of Bath, UK
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include <stdio.h>
#include <string.h>

#define   MaxLastMember    50808
#define   MaxLastNode       8573
#define   MaxLastMemberType  500
#define   MaxCharacters 50
#define   FunctionLastInteger  1000

int ReadDescription(void);
int ReadPostiveOrZeroInteger(void);
void FindNextLine(void);
void WriteResults(void);

double    Coord[MaxLastNode+1][3],InitialCoord[MaxLastNode+1][3],a[MaxLastNode+1][3],
AppliedNodeLoad[MaxLastNode+1][3],NodeForce[MaxLastNode+1][3],
AppliedNodeMom[MaxLastNode+1][3],NodeMom[MaxLastNode+1][3],
Stiff[MaxLastNode+1][3][3],RotStiff[MaxLastNode+1],
A[MaxLastNode+1][3][3],
Displacement[MaxLastNode+1][3]={0.0},Rotation[MaxLastNode+1][3]={0.0},

OrigLength[MaxLastMember+1],
x[2][MaxLastMember+1][3],
y[2][MaxLastMember+1][3],
Rotatedx[2][MaxLastMember+1][3],
Rotatedy[2][MaxLastMember+1][3],

EA[MaxLastMemberType+1],rogA[MaxLastMemberType+1],
EIxx[MaxLastMemberType+1],EIyy[MaxLastMemberType+1],
GJ[MaxLastMemberType+1],
AxialPeak[MaxLastMemberType+1],
MPeakx[MaxLastMemberType+1],MPeaky[MaxLastMemberType+1],MPeakTorsion[MaxLastMemberType+1],

FunctionValue[MaxLastMemberType+1][FunctionLastInteger+1],
TotalLengthofMembers[MaxLastMemberType+1],

p[3],MaxCoord[3],MinCoord[3],AverageCoord[3],

DispCarryOver,RotCarryOver,
StartingAppliedSafetyFact,StartingOwnWtSafetyFact,
AppliedSafetyFact,OwnWtSafetyFact,SumForce,TotalFactoredLoad[3],TotalReaction[3],

DispMassMult,RotInertiaMult,
MaxHalfDimension,
LoadScale,TypicalLoad=0.0;

const double FunctionIncrement=5.0/double(FunctionLastInteger);

int       LastNode,NumberOfActiveNodes,LastMember,NonLinear,AxStiffType,NodalMatrixControl,
FirstCycle,PictureCycle,CalculationCycle,LastCalculationCycle;

short int MemberExists[MaxLastMember+1]={0};
short int MemberReachedPeak[MaxLastMember+1]={0};
short int TensionOnlyMember[MaxLastMemberType+1]={0};
short int NodeExists[MaxLastNode+1]={0};
short int MemberAtNode[MaxLastNode+1]={0};
short int MemberTypePropsDef[MaxLastMemberType]={0};
short int NodeDispType[MaxLastNode+1],NodeRotType[MaxLastNode+1];
short int End[2][MaxLastMember+1],MemberType[MaxLastMember+1];

char NumericalValue[101];
char *MyText;
char ExpectedDescription[MaxCharacters+1];

int ReadBasicData(void);

ifstream Control,Data;
ofstream WrittenData,Results;

char DataFile[MaxCharacters+1];

int ReadBasicData(void)
{	
	char WrittenDataFile[MaxCharacters+8]="Written";
	char  text[5]=".txt";
	Control.open("RelaxDataFile.txt");
	if(Control==NULL){cout<<"\nCannot find file \"RelaxDataFile.txt\".\n\n";
	cout<<"On a Macintosh using the Xcode compiler it should be in the folder \"Release\" inside the \"build\" folder which Xcode creates when running.\n\n";
	cout<<"On a Macintosh running the Unix excecutable it should be in your user's root directory, for example /Users/Chris/\n\n";
	cout<<"Using Dev C++ on a Windows machine, everything should be in the same folder.\n\n";return 0;}
	Control.getline(DataFile,MaxCharacters+1,'.');
	strcat(DataFile,text);
	cout<<"\nData file name is \""<<DataFile<<"\"\n";
	
	strcpy(ExpectedDescription,"Starting own weight safety factor ");if(ReadDescription()==1)return 0;
	Control>>StartingOwnWtSafetyFact;cout<<ExpectedDescription<<" = "<<StartingOwnWtSafetyFact<<"\n";
	OwnWtSafetyFact=StartingOwnWtSafetyFact;
	
	strcpy(ExpectedDescription,"Starting applied load safety factor ");if(ReadDescription()==1)return 0;
	Control>>StartingAppliedSafetyFact;cout<<ExpectedDescription<<" = "<<StartingAppliedSafetyFact<<"\n";
	AppliedSafetyFact=StartingAppliedSafetyFact;
	
	strcpy(ExpectedDescription,"Non-linear control ");if(ReadDescription()==1)return 0;
	Control>>NonLinear;cout<<ExpectedDescription<<" = "<<NonLinear<<"\n";
	
	strcpy(ExpectedDescription,"Axial stiffness control ");if(ReadDescription()==1)return 0;
	Control>>AxStiffType;cout<<ExpectedDescription<<" = "<<AxStiffType<<"\n";
	
	strcpy(ExpectedDescription,"Nodal matrix control ");if(ReadDescription()==1)return 0;
	Control>>NodalMatrixControl;cout<<ExpectedDescription<<" = "<<NodalMatrixControl<<"\n";
	
	strcpy(ExpectedDescription,"Displacement mass multiplier ");if(ReadDescription()==1)return 0;
	Control>>DispMassMult;cout<<ExpectedDescription<<" = "<<DispMassMult<<"\n";
	
	strcpy(ExpectedDescription,"Rotation inertia multiplier ");if(ReadDescription()==1)return 0;
	Control>>RotInertiaMult;cout<<ExpectedDescription<<" = "<<RotInertiaMult<<"\n";
	
	strcpy(ExpectedDescription,"Displacement carry-over factor ");if(ReadDescription()==1)return 0;
	Control>>DispCarryOver;cout<<ExpectedDescription<<" = "<<DispCarryOver<<"\n";
	
	strcpy(ExpectedDescription,"Rotation carry-over factor ");if(ReadDescription()==1)return 0;
	Control>>RotCarryOver;cout<<ExpectedDescription<<" = "<<RotCarryOver<<"\n";
	
	Control.close();
	
	Data.open(DataFile);
	if(Data==NULL){cout<<"\nCannot find data file refered to in \"RelaxDataFile.txt\". It should be in the same folder as\"RelaxDataFile.txt\".\n\n";return 0;}
	
	strcat (WrittenDataFile,DataFile);
	WrittenData.open(WrittenDataFile);
	WrittenData<<"Name of data file: "<<DataFile<<"\n\n\nThe following section contains:\n\n";
	WrittenData<<"Node number, x, y, z coordinates and initial rotation vector, often 0.0 0.0 0.0\n\n";
	
	{
		int Node;
		LastNode=0;
		int FirstNode=1;
		
		do{Node=ReadPostiveOrZeroInteger();}while(Node<0);
		if(Node>MaxLastNode){cout<<"Too many nodes\n";return 0;}
		while(Node>=0)
		{
			if(LastNode<Node)LastNode=Node;
			WrittenData<<Node;
			NodeExists[Node]=1;
			for(int i=0;i<=2;i++){Data>>Coord[Node][i];InitialCoord[Node][i]=Coord[Node][i];WrittenData<<"\t"<<Coord[Node][i];}
			WrittenData<<"\n";
			for(int i=0;i<=2;i++){Data>>a[Node][i];WrittenData<<"\t"<<a[Node][i];}
			WrittenData<<"\n";
			if(FirstNode==1)
			{
				for(int i=0;i<=2;i++)
				{
					MaxCoord[i]=Coord[Node][i];
					MinCoord[i]=Coord[Node][i];
				}
			}
			else
			{
				for(int i=0;i<=2;i++)
				{
					if(MaxCoord[i]<Coord[Node][i])MaxCoord[i]=Coord[Node][i];
					if(MinCoord[i]>Coord[Node][i])MinCoord[i]=Coord[Node][i];
				}
			}
			FirstNode=0;
			Node=ReadPostiveOrZeroInteger();
			if(Node>MaxLastNode){cout<<"Too many nodes\n";return 0;}
		}
		WrittenData<<"\n\nThe following section contains:\n\n";
		WrittenData<<"Member number,  nodes at ends, member type, slack length\n";
		WrittenData<<"and components of the local x and y axes at the two member ends\n\n";
		cout<<"\nNodes read, last node = "<<LastNode<<"\n";
	}
	
	double LongestMember=0.0;
	
	{
		int Member;
		LastMember=0;
		do{Member=ReadPostiveOrZeroInteger();}while(Member<0);
		if(Member>MaxLastMember){cout<<"Too many members\n";return 0;}
		while(Member>=0)
		{
			if(LastMember<Member)LastMember=Member;
			if(MemberExists[Member]==1)cout<<"Member number "<<Member<<" already exists\n";
			MemberExists[Member]=1;
			Data>>End[0][Member]>>End[1][Member]>>MemberType[Member];
			Data>>OrigLength[Member];if(LongestMember<OrigLength[Member])LongestMember=OrigLength[Member];
			WrittenData<<Member<<"\t"<<End[0][Member]<<"\t"<<End[1][Member]<<"\t"<<MemberType[Member]<<"\t"<<OrigLength[Member]<<"\n";
			if(End[0][Member]==End[1][Member])
			{
				MemberExists[Member]=0;
				cout<<"Member "<<Member<<" is connected the same node at both ends.\nThe remainder of the member data is read, but the member is removed.\n";
			}
			if(OrigLength[Member]==0.0)
			{
				MemberExists[Member]=0;
				cout<<"Member "<<Member<<" is has zero original length.\nThe remainder of the member data is read, but the member is removed.\n";
			}
			if(NodeExists[End[0][Member]]==0||NodeExists[End[1][Member]]==0)
			{
				MemberExists[Member]=0;
				cout<<"Member "<<Member<<" is connected to nodes, at least one of which which does not exist.\nThe remainder of the member data is read, but the member is removed.\n";
			}
			for(int MemberEnd=0;MemberEnd<=1;MemberEnd++)
			{
				for(int i=0;i<=2;i++)Data>>x[MemberEnd][Member][i];
				for(int i=0;i<=2;i++)Data>>y[MemberEnd][Member][i];
				double tempdouble=0.0;
				for(int i=0;i<=2;i++)tempdouble+=x[MemberEnd][Member][i]*x[MemberEnd][Member][i];
				tempdouble=sqrt(tempdouble);
				for(int i=0;i<=2;i++)x[MemberEnd][Member][i]=x[MemberEnd][Member][i]/tempdouble;
				tempdouble=0.0;
				for(int i=0;i<=2;i++)tempdouble+=x[MemberEnd][Member][i]*y[MemberEnd][Member][i];
				for(int i=0;i<=2;i++)y[MemberEnd][Member][i]-=x[MemberEnd][Member][i]*tempdouble;
				tempdouble=0.0;
				for(int i=0;i<=2;i++)tempdouble+=y[MemberEnd][Member][i]*y[MemberEnd][Member][i];
				tempdouble=sqrt(tempdouble);
				for(int i=0;i<=2;i++)y[MemberEnd][Member][i]=y[MemberEnd][Member][i]/tempdouble;
				
				double ScalarProduct=0.0;
				for(int i=0;i<=2;i++)ScalarProduct+=x[MemberEnd][Member][i]*(Coord[End[1][Member]][i]-Coord[End[0][Member]][i])/OrigLength[Member];
				if(fabs(ScalarProduct)>1.0e-3)cout<<"Member "<<Member<<" is initially curved, scalar product = "<<ScalarProduct<<"\n";
				ScalarProduct=0.0;
				for(int i=0;i<=2;i++)ScalarProduct+=y[MemberEnd][Member][i]*(Coord[End[1][Member]][i]-Coord[End[0][Member]][i])/OrigLength[Member];
				if(fabs(ScalarProduct)>1.0e-3)cout<<"Member "<<Member<<" is initially curved, scalar product = "<<ScalarProduct<<"\n";
				
				for(int i=0;i<=2;i++)WrittenData<<"\t"<<x[MemberEnd][Member][i];
				WrittenData<<"\n";
				for(int i=0;i<=2;i++)WrittenData<<"\t"<<y[MemberEnd][Member][i];
				WrittenData<<"\n";
			}
			
			double ActualLengthSq=0.0;
			for(int i=0;i<=2;i++)ActualLengthSq+=(Coord[End[1][Member]][i]-Coord[End[0][Member]][i])*(Coord[End[1][Member]][i]-Coord[End[0][Member]][i]);
			double strain=(sqrt(ActualLengthSq)-OrigLength[Member])/OrigLength[Member];
			if(fabs(strain)>1.0e-3)cout<<"Strain in member "<<Member<<" is "<<strain<<"\n";
			
			Member=ReadPostiveOrZeroInteger();if(Member>MaxLastMember){cout<<"Too many members\n";return 0;}
		}
		WrittenData<<"\n\nThe following section contains:\n\n";
		WrittenData<<"Member type, Young's modulus, Poisson's rato.\n";
		WrittenData<<"These next 2 factors control yielding, or to be more precise, non-linear elastic behaviour:\n";
		WrittenData<<"yield stress, factor controlling stress/strain curve (typically 5.0),\n";
		WrittenData<<"density times acceleration due to gravity,\n";
		WrittenData<<"cross-sectional area, Ixx, Iyy, J (as in GJ = torsional stiffness),\n";
		WrittenData<<"factor mutiplying area times yield stress to give yield tension.\n";
		WrittenData<<"The remaining factors are lengths with which to divide Ixx, Iyy and J to give the corresponding plastic section modulii.\n\n";
		cout<<"Members read, last member = "<<LastMember<<"\n";
	}
	
	{
		int ThisMemberType;
		do{ThisMemberType=ReadPostiveOrZeroInteger();}while(ThisMemberType<0);
		if(ThisMemberType>MaxLastMemberType){cout<<"Too many member properties\n";return 0;}
		while(ThisMemberType>=0)
		{
			double E,G,Poisson,rog,CrossSectArea,Ixx,Iyy,J,sigmay,lambda,areafactor,benlenx,benleny,twistlen;
			
			MemberTypePropsDef[ThisMemberType]=1;
			
			Data>>E>>Poisson>>sigmay>>lambda>>rog;
			Data>>CrossSectArea>>Ixx>>Iyy>>J;
			Data>>areafactor>>benlenx>>benleny>>twistlen;
			
			WrittenData<<ThisMemberType<<"\t"<<E<<"\t"<<Poisson<<"\t"<<sigmay<<"\t"<<lambda<<"\t"<<rog<<"\n";
			WrittenData<<"\t"<<CrossSectArea<<"\t"<<Ixx<<"\t"<<Iyy<<"\t"<<J<<"\n";
			WrittenData<<"\t"<<areafactor<<"\t"<<benlenx<<"\t"<<benleny<<"\t"<<twistlen<<"\n";
			
			if(Ixx==0.0||Iyy==0.0||J==0.0)TensionOnlyMember[ThisMemberType]=1;
			
			rogA[ThisMemberType]=rog*CrossSectArea;if(TypicalLoad<rogA[ThisMemberType]*LongestMember)TypicalLoad=rogA[ThisMemberType]*LongestMember;
			
			EA[ThisMemberType]=E*CrossSectArea;
			EIxx[ThisMemberType]=E*Ixx;
			EIyy[ThisMemberType]=E*Iyy;
			G=E/(2.0*(1.0+Poisson));
			GJ[ThisMemberType]=G*J;
			
			cout<<"Member Type "<<ThisMemberType<<" is defined.";
			if(TensionOnlyMember[ThisMemberType]==1)cout<<" They are tension only members.";
			cout<<"\n";
			
			if(EA[ThisMemberType]==0.0&&EIxx[ThisMemberType]==0.0&&EIyy[ThisMemberType]==0.0&&GJ[ThisMemberType])
			{
				cout<<"\nMember Type "<<ThisMemberType<<" has no stiffness contribution.\n\n";
				return 0;
			}
			
			AxialPeak[ThisMemberType]=areafactor*sigmay*CrossSectArea;
			MPeakx[ThisMemberType]=sigmay*Ixx/benlenx;
			MPeaky[ThisMemberType]=sigmay*Iyy/benleny;
			MPeakTorsion[ThisMemberType]=(sigmay/sqrt(3.0))*J/twistlen;
			
			for(int FunctionInteger=0;FunctionInteger<=FunctionLastInteger;FunctionInteger++)
			{
				double argument=double(FunctionInteger)*FunctionIncrement;
				double newargument=tanh(lambda)*tanh(lambda*argument/tanh(lambda));
				
				FunctionValue[ThisMemberType][FunctionInteger]=(1.0/lambda)*0.5*log((1.0+newargument)/(1.0-newargument));
				//cout<<argument<<"\t"<<FunctionValue[ThisMemberType][FunctionInteger]<<"\n";
			}
			
			ThisMemberType=ReadPostiveOrZeroInteger();if(ThisMemberType>MaxLastMemberType){cout<<"Too many member properties\n";return 0;}
		}
		WrittenData<<"\n\nThe following section contains:\n\n";
		WrittenData<<"Node number, restraint in direction and rotation according to adding:\n";
		WrittenData<<"The number 1 = x direction prevented\n";
		WrittenData<<"The number 2 = y direction prevented\n";
		WrittenData<<"The number 4 = z direction prevented\n";
		WrittenData<<"All three = 1 + 2 + 4 = 7\n";
		WrittenData<<"The number 1 = x rotation prevented\n";
		WrittenData<<"The number 2 = y rotation prevented\n";
		WrittenData<<"The number 4 = z rotation prevented\n";
		WrittenData<<"All three = 1 + 2 + 4 = 7\n\n";
		cout<<"Member properties read\n";
	}
	
	for(int Member=0;Member<=LastMember;Member++)
	{
		if(MemberExists[Member]==1&&MemberTypePropsDef[MemberType[Member]]!=1)
		{
			MemberExists[Member]==0;
			cout<<"Member "<<Member<<" is of type "<<MemberType[Member]<<" which has no properties defined. The member has now been removed.\n";
		}
		
		if(MemberExists[Member]==1)
		{
			MemberAtNode[End[0][Member]]=1;
			MemberAtNode[End[1][Member]]=1;
		}
	}
	
	for(int Node=0;Node<=LastNode;Node++)
	{
		if(NodeExists[Node]==1&&MemberAtNode[Node]==0)
		{
			cout<<"Node "<<Node<<" exists but has no member attached to it. The node has now been removed.\n";
			NodeExists[Node]=0;
		}
	}
	
	for(int Node=0;Node<=LastNode;Node++)
	{
		NodeDispType[Node]=0;
		NodeRotType[Node]=0;
	}
	
	{
		int Node;
		do{Node=ReadPostiveOrZeroInteger();}while(Node<0);
		while(Node>=0)
		{
			Data>>NodeDispType[Node];
			
			// 1 = x direction prevented
			// 2 = y direction prevented
			// 4 = z direction prevented
			
			Data>>NodeRotType[Node];
			
			// 1 = x rotation prevented
			// 2 = y rotation prevented
			// 4 = z rotation prevented
			
			WrittenData<<Node<<"\t"<<NodeDispType[Node]<<"\t"<<NodeRotType[Node]<<"\n";
			
			Node=ReadPostiveOrZeroInteger();
		}
		WrittenData<<"\n\nThe following section contains:\n\n";
		WrittenData<<"Node number and 3 components of applied load and applied moment\n\n";
		cout<<"Supports read\n";
	}
	
	for(int Node=0;Node<=LastNode;Node++)
	{
		for(int i=0;i<=2;i++)
		{
			AppliedNodeLoad[Node][i]=0.0;if(TypicalLoad<fabs(AppliedNodeLoad[Node][i]))TypicalLoad=fabs(AppliedNodeLoad[Node][i]);
			AppliedNodeMom[Node][i]=0.0;
		}
	}
	
	{
		int Node;
		do{Node=ReadPostiveOrZeroInteger();}while(Node<0);
		while(Node>=0)
		{
			for(int i=0;i<=2;i++)Data>>AppliedNodeLoad[Node][i];
			for(int i=0;i<=2;i++)Data>>AppliedNodeMom[Node][i];
			
			WrittenData<<Node;
			for(int i=0;i<=2;i++)WrittenData<<"\t"<<AppliedNodeLoad[Node][i];
			WrittenData<<"\n";
			for(int i=0;i<=2;i++)WrittenData<<"\t"<<AppliedNodeMom[Node][i];
			WrittenData<<"\n";
			
			Node=ReadPostiveOrZeroInteger();
		}
		WrittenData<<"\nEnd of data\n\n";
		cout<<"Loads read\n";
	}
	
	Data.close();
	WrittenData.close();
	
	NumberOfActiveNodes=0;
	for(int Node=0;Node<=LastNode;Node++)
	{
		if(NodeExists[Node]==1)NumberOfActiveNodes++;
	}
	
	for(int ThisMemberType=0;ThisMemberType<=MaxLastMemberType;ThisMemberType++)TotalLengthofMembers[ThisMemberType]=0.0;
	
	for(int Member=0;Member<=LastMember;Member++)
	{
		if(MemberExists[Member]==1)TotalLengthofMembers[MemberType[Member]]+=OrigLength[Member];
	}		
	
	return 1;
}

int ReadDescription(void)
{
	char Description[MaxCharacters+1];
	char EndOfLine=' ';
	while(EndOfLine!='\r'&&EndOfLine!='\n')EndOfLine=Control.get();
	EndOfLine=Control.peek();if(EndOfLine=='\r'||EndOfLine=='\n')EndOfLine=Control.get();
	Control.getline(Description,MaxCharacters+1,'=');
	if(strcmp(Description,ExpectedDescription)!=0)
	{
		cout<<"Was expecting the text \""<<ExpectedDescription<<"\" but found \""<<Description<<"=\"\n";
		return 1;
	}
	return 0;
}

int ReadPostiveOrZeroInteger(void)
{
	char c;int Number;
	c=Data.peek();
	while(c==' '||c=='\r'||c=='\n'||c=='\t')
	{c=Data.get();c=Data.peek();}
	if(c!='+'&&c!='0'&&c!='1'&&c!='2'&&c!='3'&&c!='4'&&c!='5'&&c!='6'&&c!='7'&&c!='8'&&c!='9')
	{Number=-1;FindNextLine();}else Data>>Number;
	return Number;
}

void FindNextLine(void)
{
	char c;
	do{c=Data.get();}
	while(c!='\r'&&c!='\n'&&c!=EOF);
}

void WriteResults(void)
{
	char ResultFile[MaxCharacters+8]="Results";
	strcat (ResultFile,DataFile);
	Results.open(ResultFile);
	
	double SumTotalWeights=0.0;
	Results<<"Member type\tTotal length\tTotal weight\n";
	for(int ThisMemberType=0;ThisMemberType<=MaxLastMemberType;ThisMemberType++)
	{
		if(MemberTypePropsDef[ThisMemberType]==1)
		{
			double WeightThisType=TotalLengthofMembers[ThisMemberType]*rogA[ThisMemberType];
			SumTotalWeights+=WeightThisType;
			Results<<ThisMemberType<<"\t\t"<<TotalLengthofMembers[ThisMemberType]<<"\t\t"<<WeightThisType<<"\n";
		}
	}
	
	Results<<"\nSum of all member weights = "<<SumTotalWeights<<"\n\n";
	
	Results<<"x, y, z components of total factored load:\n\n"<<TotalFactoredLoad[0]<<"\t"<<TotalFactoredLoad[1]<<"\t"<<TotalFactoredLoad[2]<<"\n";
	Results<<"\nx, y, z components of total reaction:\n\n"<<TotalReaction[0]<<"\t"<<TotalReaction[1]<<"\t"<<TotalReaction[2]<<"\n";
	Results<<"\nDisplacements:\n\nNode\tx\ty\tz\n\n";
	for(int Node=0;Node<=LastNode;Node++)
	{
		if(NodeExists[Node]==1)
		{
			Results<<Node;
			for(int i=0;i<=2;i++)Results<<"\t"<<Coord[Node][i]-InitialCoord[Node][i];
			Results<<"\n";
		}
	}
	Results.close();
}