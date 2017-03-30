//Note on Macintosh run program once from 2DPinned.xcodeproj.
//This will create a folder 'build' and a folder 'Release' within 'build'.
//Put the data files in the folder Release and run again.
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "size.h"

using namespace std;

#define   MaxLastMember 199
#define   MaxLastNode    99

void Invert(int, double[][Max_n + 1], double[][Max_n + 1]);
void DrawStructure(void);
void dxfSetUp(void);
void dxfMember(int);
void dxfLoad(int, int);
void dxfMemberNo(int);
void dxfNodeNo(int);
void dxfNodeCircle(int);
void dxfFinishOff(void);
void printMDArray(int, double[][Max_n + 1]);

int End[2][MaxLastMember + 1], Fixed[Max_n + 1], LastMember, LastNode, LastDoF;
double x[2][MaxLastNode + 1], Load[Max_n + 1], delta[Max_n + 1], K[Max_n + 1][Max_n + 1], Flex[Max_n + 1][Max_n + 1], EA[MaxLastMember + 1], xDiff[2][MaxLastMember + 1], L[MaxLastMember + 1];
ifstream mems("members.txt"), nods("nodes.txt"), lds("loads.txt"), sup("supports.txt");
ofstream Picy("structure.dxf");
ofstream Tension("MemberTensions.txt");
int main(void)
{
    if (2 * MaxLastNode + 2 > Max_n + 1) {
	cout << "Problem with sizes\n";
	return 0;
    }
    mems >> LastMember;
    if (LastMember > MaxLastMember) {
	cout << "Too many members\n";
	return 0;
    }
    for (int tempint = 0; tempint <= LastMember; tempint++) {
	int Member;
	mems >> Member;
	mems >> End[0][Member] >> End[1][Member] >> EA[Member];
    }
    mems.close();
    nods >> LastNode;
    if (LastNode > MaxLastNode) {
	cout << "Too many Nodes\n";
	return 0;
    }
    for (int tempint = 0; tempint <= LastNode; tempint++) {
	int Node;
	nods >> Node;
	nods >> x[0][Node] >> x[1][Node];
    }
    nods.close();
    LastDoF = 2 * LastNode + 1;
    for (int ThisDoF = 0; ThisDoF <= LastDoF; ThisDoF++) {
	Load[ThisDoF] = 0.0;
	delta[ThisDoF] = 0.0;
	Fixed[ThisDoF] = 0;
	for (int ThatDoF = 0; ThatDoF <= LastDoF; ThatDoF++)
	    K[ThisDoF][ThatDoF] = 0.0;
    }
    for (;;) {
	int Node, Direction;
	lds >> Node;
	if (Node < 0)
	    break;
	lds >> Direction;
	lds >> Load[2 * Node + Direction];
    }
    for (;;) {
	int Node, Direction;
	sup >> Node;
	if (Node < 0)
	    break;
	sup >> Direction;
	Fixed[2 * Node + Direction] = 1;
	Load[2 * Node + Direction] = 0.0;
    }
    for (int Member = 0; Member <= LastMember; Member++) {
	double Lsq = 0.0;
	for (int ThisCoord = 0; ThisCoord <= 1; ThisCoord++) {
	    xDiff[ThisCoord][Member] = x[ThisCoord][End[1][Member]] - x[ThisCoord][End[0][Member]];
	    Lsq += xDiff[ThisCoord][Member] * xDiff[ThisCoord][Member];
	}
	L[Member] = sqrt(Lsq);
	for (int ThisCoord = 0; ThisCoord <= 1; ThisCoord++) {
	    for (int ThatCoord = 0; ThatCoord <= 1; ThatCoord++) {
		double ThisStiffness = (EA[Member] / (L[Member] * L[Member] * L[Member])) * xDiff[ThisCoord][Member] * xDiff[ThatCoord][Member];
		for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
		    for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {
			int ThisDoF = 2 * End[ThisEnd][Member] + ThisCoord;
			int ThatDoF = 2 * End[ThatEnd][Member] + ThatCoord;
			if (ThisEnd == ThatEnd)
			    K[ThisDoF][ThatDoF] += ThisStiffness;
			else
			    K[ThisDoF][ThatDoF] -= ThisStiffness;
		    }
		}
	    }
	}
    }
    for (int ThisDoF = 0; ThisDoF <= LastDoF; ThisDoF++) {
	if (Fixed[ThisDoF] == 1) {
	    K[ThisDoF][ThisDoF] = 1.0;
	    for (int ThatDoF = 0; ThatDoF <= LastDoF; ThatDoF++) {
		if (ThatDoF != ThisDoF) {
		    K[ThatDoF][ThisDoF] = 0.0;
		    K[ThisDoF][ThatDoF] = 0.0;
		}
	    }
	}
    }
    Invert(LastDoF, K, Flex);
    cout << "Global stifness matrix\n";
    printMDArray(LastDoF,K);
    cout << "Inverted matrix\n";
    printMDArray(LastDoF,Flex);
    for (int ThisDoF = 0; ThisDoF <= LastDoF; ThisDoF++) {
    for (int ThatDoF = 0; ThatDoF <= LastDoF; ThatDoF++){
	    delta[ThisDoF] += Flex[ThisDoF][ThatDoF] * Load[ThatDoF];
    }
    }
    for (int Member = 0; Member <= LastMember; Member++) {
	double MemberTension = (EA[Member] / (L[Member] * L[Member]))
	    * (xDiff[0][Member] * (delta[2 * End[1][Member] + 0] - delta[2 * End[0][Member] + 0])
	       + xDiff[1][Member] * (delta[2 * End[1][Member] + 1] - delta[2 * End[0][Member] + 1]));
	Tension << Member << "  " << MemberTension << "\n";
    }
    Tension.close();
    DrawStructure();
    cout << "dxf file written, end of program.\n";
    system("pause");
    return 0;
}

// END OF MAIN PROCEDURE
// -----------------------------------------------------------------------------------------------------------------------------


// SUB PROCEDURES
// -----------------------------------------------------------------------------------------------------------------------------

void printMDArray(int Max, double MDArray[][Max_n + 1])
{
    cout.precision(1);
    for (int i = 0; i <= Max; i++) {
    cout << fixed << i;
	for (int j = 0; j <= Max; j++) {
        cout << "\t";
//        cout << j << ":";
        cout << fixed << MDArray[i][j];
    }
	cout << "\n";
    }
}

void dxfSetUp(void)
{
    Picy << "0\n" << "SECTION\n" << "2\n" << "ENTITIES\n";
}
void DrawStructure(void)
{
    dxfSetUp();
    for (int Member = 0; Member <= LastMember; Member += 1)
	dxfMember(Member);
    for (int Node = 0; Node <= LastNode; Node += 1)
	dxfNodeCircle(Node);
    for (int Member = 0; Member <= LastMember; Member += 1)
	dxfMemberNo(Member);
    for (int Node = 0; Node <= LastNode; Node += 1)
	dxfNodeNo(Node);
    for (int Node = 0; Node <= LastNode; Node += 1) {
	for (int Direction = 0; Direction <= 1; Direction++) {
	    if (Load[2 * Node + Direction] != 0.0)
		dxfLoad(Node, Direction);
	}
    }
    dxfFinishOff();
}
void dxfMember(int ThisMember)
{
    Picy << "0\nLINE\n8\nUndeflected\n" << "10\n" << x[0][End[0][ThisMember]] << "\n" << "20\n" << x[1][End[0][ThisMember]] << "\n" << "11\n" << x[0][End[1][ThisMember]] << "\n" << "21\n" << x[1][End[1][ThisMember]] << "\n" << "62\n0\n";
    Picy << "0\nLINE\n8\nDeflected\n" << "10\n" << x[0][End[0][ThisMember]] + delta[2 * End[0][ThisMember] + 0] << "\n" << "20\n" << x[1][End[0][ThisMember]] + delta[2 * End[0][ThisMember] + 1] << "\n" << "11\n" << x[0][End[1][ThisMember]] + delta[2 * End[1][ThisMember] + 0] << "\n" << "21\n" << x[1][End[1][ThisMember]] + delta[2 * End[1][ThisMember] + 1] << "\n" << "62\n1\n";
}
void dxfNodeCircle(int ThisNode)
{
    double CircleRadius = 0.02;
    Picy << "0\nCIRCLE\n8\nUndeflected\n" << "10\n" << x[0][ThisNode] << "\n" << "20\n" << x[1][ThisNode] << "\n" << "40\n" << CircleRadius << "\n" << "62\n0\n";
    Picy << "0\nCIRCLE\n8\nDeflected\n" << "10\n" << x[0][ThisNode] + delta[2 * ThisNode + 0] << "\n" << "20\n" << x[1][ThisNode] + delta[2 * ThisNode + 1] << "\n" << "40\n" << CircleRadius << "\n" << "62\n1\n";
}
void dxfLoad(int ThisNode, int ThisDirection)
{
    Picy << "0\nLINE\n8\nInformation\n" << "10\n" << x[0][ThisNode] << "\n" << "20\n" << x[1][ThisNode] << "\n";
    if (ThisDirection == 0) {
	Picy << "11\n" << x[0][ThisNode] + Load[2 * ThisNode + ThisDirection] << "\n" << "21\n" << x[1][ThisNode] << "\n";
    } else {
	Picy << "11\n" << x[0][ThisNode] << "\n" << "21\n" << x[1][ThisNode] + Load[2 * ThisNode + ThisDirection] << "\n";
    }
    Picy << "62\n1\n";
}
void dxfMemberNo(int ThisMember)
{
    double TextSize = 0.05;
    Picy << "0\nTEXT\n8\nInformation\n" << "10\n" << (x[0][End[0][ThisMember]] + x[0][End[1][ThisMember]]) / 2.0 << "\n" << "20\n" << (x[1][End[0][ThisMember]] + x[1][End[1][ThisMember]]) / 2.0 << "\n" << "40\n" << TextSize << "\n" << "62\n0\n" << "1\n" << ThisMember << "\n";
}
void dxfNodeNo(int ThisNode)
{
    double TextSize = 0.05;
    Picy << "0\nTEXT\n8\nInformation\n" << "10\n" << x[0][ThisNode] + TextSize << "\n" << "20\n" << x[1][ThisNode] + TextSize << "\n" << "40\n" << TextSize << "\n" << "62\n0\n" << "1\n" << ThisNode << "\n";
}
void dxfFinishOff(void)
{
    Picy << "0\n" << "ENDSEC\n" << "0\n" << "EOF\n";
    Picy.close();
}
