//Written by Chris J K Williams, University of Bath, UK

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include "Graphics.h"

// LastOrdinaryMemberType has to be 3, 4 or 5
#define LastOrdinaryMemberType 3

#define DefinedMembersPerLevel  6 + 3 * LastOrdinaryMemberType
#define LastLevel 12
#define LastNode 6 * (LastLevel + 1) + 2 - 3
#define LastMember (6 + 3 * LastOrdinaryMemberType) * (LastLevel - 1) + 15 + 6 + 2 - 3 - 3 + 3
#define RightOrLeftHanded - 1.0

void InitialValues(void);
void StartingGeometry(void);
void ForcesFromMembers(void);
double CalculateMemberLength(int);
double NodalRadius(int);
int Permute(int);
int LastLevelPermute(int);
void RotateNode(int);
void WriteDXF(void);
void StrutProperties(int);
void TieProperties(int);
void MemberSetUp(int);
double RelativeLevelFromTop(int);

double Initial_x[LastNode + 1][3],x[LastNode + 1][3],Force[LastNode + 1][3],Velocity[LastNode + 1][3],Load[LastNode + 1][3],
NodalArea[LastNode + 1],Mass[LastNode + 1];

int NodeType[LastNode + 1],RotationalSymmetryType[LastNode + 1];

double InitialLength[LastMember + 1],Diameter[LastMember + 1],
EA[LastMember + 1],EI[LastMember + 1],RoA[LastMember + 1],EulerLoad[LastMember + 1],MemberStrengthInTension[LastMember + 1],Tension[LastMember + 1];

int End[LastMember + 1][2],MemberType[LastMember + 1],MemberLevel[LastMember + 1],MemberBuckled[LastMember + 1],MemberBroken[LastMember + 1];

double a[LastLevel + 1],c[LastLevel + 1],Height[LastLevel + 1];

int MembersPerLevel,PrestressMember;

double PointPosition[2],delta_t,root3,Prestress,MaxBucklingOverYield,UnstressFactor,gravity,WindPressureTimesDragCoeff,TotalNodalArea,CurrentShrinkPerLevel;

char NumericalValue[101];
char *MyText;

int main(int argc, char *  argv[])
{
	//Note: Use consistent system of units e.g. kilogramme - metre - second in which case forces are in newtons
	
	MembersPerLevel = DefinedMembersPerLevel;
	
	PI = 4.0 * atan(1.0);
	root3 = sqrt(3.0);
	
	Prestress = 1000.0;
	
	MaxBucklingOverYield = 0.7;
	
	UnstressFactor = 1.0;
	
	gravity = 9.81;
	WindPressureTimesDragCoeff = 1200.0;//In newtons per square metre
	
	delta_t = 0.0001;
	
	InitialValues();
	StartingGeometry();
	WriteDXF_Now = 0;
	
	TotalNodalArea = 0.0;
	for(int Node = 0;Node <= LastNode; Node ++)TotalNodalArea += NodalArea[Node];
	
	glutInit(&argc,argv);
	MyGraphics(1.0,1.0,1.0);
	return 0;
}

static void Draw(void)
{
	if(RestartWhenReady == 1)StartingGeometry();
	
	if(InitialStateWhenReady == 1)
	{
		InitialValues();
		StartingGeometry();
	}
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	glLoadIdentity();
	
	glPointSize(10.0);
	
	for(int Button = 0;Button <= LastButton; Button ++)
	{
		if(OverThisButton[Button] == 1)glColor4f(1.0,0.0,0.0,0.5);else glColor4f(0.0,0.0,0.0,0.2);
		glBegin(GL_POINTS);
		PointPosition[0] = ButtonX;
		PointPosition[1] = ButtonY[Button];
		glVertex2dv(PointPosition);
		glEnd();
		
		if(Button == 0)MyText="Zoom";
		if(Button == 1)MyText="Relative ring cable length";
		if(Button == 2)MyText="Relative prestress";
		if(Button == 3)MyText="Taper (press Restart to take effect)";
		if(Button == 4)MyText="Relative size (press Restart to take effect)";
		if(Button == 5)MyText="Unstress";
		if(Button == 6)MyText="Load";
		if(Button == 7)MyText="Restart";
		if(Button == 8)MyText="Reset";
		if(Button == 9)MyText="Write DXF - overwrites any existing file";
		
		glColor4f(0.0,0.0,0.0,0.5);
		glRasterPos2d(ButtonX + 0.02,ButtonY[Button] - 0.01);
		
		int StringLength=strlen(MyText);
		for(int MyChar=0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
		if(Button > 0 && Button < 5)
		{
			if(Button == 1)sprintf(NumericalValue,"%8.2f",RelativeRingLength);
			else
			{
				if(Button == 2)sprintf(NumericalValue,"%8.2f",RelativePrestress);
				else 
				{
				if(Button == 3)sprintf(NumericalValue,"%8.2f",ShrinkPerLevel);
				else sprintf(NumericalValue,"%8.2f",SizeMultiplier);
				}
			}
			int LengthOfString=strlen(NumericalValue);
			for (int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
				glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
		}
	}
	
	glColor4f(0.0,0.0,0.0,0.5);
	{
		glRasterPos2d(ButtonX + 0.02,0.45);
		MyText="Height at top =";
		int StringLength=strlen(MyText);
		for(int MyChar=0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
		sprintf(NumericalValue,"%8.2f",x[LastNode - 3][2]);
		int LengthOfString=strlen(NumericalValue);
		for (int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	}
	
	{
		glRasterPos2d(ButtonX + 0.02,0.35);
		MyText="Height at support =";
		int StringLength=strlen(MyText);
		for(int MyChar=0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
		sprintf(NumericalValue,"%8.2f",x[1][2]);
		int LengthOfString=strlen(NumericalValue);
		for (int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	}
	
	{
		glRasterPos2d(ButtonX + 0.02,0.25);
		MyText="Height at bottom of lowest struts - NOT attached to a support =";
		int StringLength=strlen(MyText);
		for(int MyChar=0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
		sprintf(NumericalValue,"%8.2f",x[0][2]);
		int LengthOfString=strlen(NumericalValue);
		for (int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	}
	
	{
		glRasterPos2d(ButtonX + 0.02,0.15);
		MyText="Total frontal area =";
		int StringLength=strlen(MyText);
		for(int MyChar=0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
		sprintf(NumericalValue,"%8.2f",TotalNodalArea);
		int LengthOfString=strlen(NumericalValue);
		for (int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	}
	
	glPopMatrix();
	
	for(int Member = 0;Member <= LastMember;Member ++)
	{
		if(MemberType[Member] == 0)
		{
			glLineWidth(2.0);
		glColor4f(0.0,0.0,1.0,1.0);}
		else
		{
			glLineWidth(1.0);
			glColor4f(0.0,0.0,0.0,1.0);
		}
		if(FormfindUnstressLoad == 2)//That is apply load
		{
			glBegin(GL_LINE_STRIP);
			glVertex3dv(Initial_x[End[Member][0]]);
			glVertex3dv(Initial_x[End[Member][1]]);
			glEnd();
		}
		if(MemberBuckled[Member] != 1 && MemberBroken[Member] != 1)
		{
			if(MemberType[Member] == 0)glColor4f(0.0,0.0,1.0,1.0);
			else glColor4f(0.0,0.0,0.0,0.5);
		}
		else glColor4f(1.0,0.0,0.0,1.0);
		if(FormfindUnstressLoad == 1)//That is unstress
		{
			if(abs(Tension[Member]) > 1.0e-6 * Prestress)glColor4f(1.0,0.0,0.0,1.0);else glColor4f(0.0,0.0,1.0,1.0);
			if(Member == PrestressMember)
			{
				glColor4f(1.0,0.0,1.0,1.0);
				glLineWidth(3.0);
			}
		}
		glBegin(GL_LINE_STRIP);
		glVertex3dv(x[End[Member][0]]);
		glVertex3dv(x[End[Member][1]]);
		glEnd();
	}
	
	glPointSize(3.0);
	glBegin(GL_POINTS);
	for(int Node = 0;Node <= LastNode; Node ++)
	{
		if(NodeType[Node] == 0)glColor4f(0.0,0.0,0.0,0.5);
		else glColor4f(1.0,0.0,0.0,1.0);
		if(FormfindUnstressLoad == 2)//That is apply load
			glVertex3dv(Initial_x[Node]);
		glVertex3dv(x[Node]);
	}
	glEnd();
	
	glutSwapBuffers();
	
	if(WriteDXF_Now == 1)
	{
		WriteDXF();
		ofstream Julia("Information.txt");
		Julia<<"Last node = "<<LastNode<<"\n";
		Julia<<"Last member = "<<LastMember<<"\n";
		for(int whichMemberType = 0;whichMemberType <= 100; whichMemberType++)
		{
			int FoundOne = 0;
			for(int Member = 0;Member <= LastMember;Member ++)
			{
				if(MemberType[Member] == whichMemberType)FoundOne = 1;
			}
			if(FoundOne == 1)
			{
				Julia<<"\nMember type "<<whichMemberType<<"\n\n";
				for(int Member = 0;Member <= LastMember;Member ++)
				{
					if(MemberType[Member] == whichMemberType)
					{
						double BucklingOverYield = EulerLoad[Member] / MemberStrengthInTension[Member];
						Julia<<"\tMember "<<Member<<"\n\t\tInitial length "<<InitialLength[Member]<<"m Diameter "<<1000.0 * Diameter[Member]<<"mm\n";
						Julia<<"\t\tMember strength in tension "<<MemberStrengthInTension[Member] / 1000.0<<"kN Euler load "<<EulerLoad[Member] / 1000.0<<"kN Ratio "<<EulerLoad[Member] / MemberStrengthInTension[Member]<<"\n";
						if(BucklingOverYield > MaxBucklingOverYield)
						{
							cout<<"\nMember "<<Member<<" is not slender\n";cout << "Press ENTER to continue... ";
							cin.get();
						}
					}
				}
			}
		}
		for(int Node = 0;Node <= 1; Node ++)
		{
			Julia<<"Base radius type "<<Node<<" = "<<NodalRadius(Node)<<"\n";
		}
		Julia<<"Height of structure = "<<x[LastNode - 3][2]<<"m\n";//3 extra nodes at base
		Julia.close();
		WriteDXF_Now = 0;
	}
	
	for(int cycle = 0;cycle <= 1000; cycle ++)
	{
		ForcesFromMembers();
		
		for(int Node = 0;Node <= LastNode; Node ++)
		{
			if(NodeType[Node] == 0)
			{
				if(FormfindUnstressLoad == 0  && (RotationalSymmetryType[Node] != 0 && RotationalSymmetryType[Node] != 1))//That is formfinding
					RotateNode(Node);
				else
				{
					for(int xyz = 0;xyz <= 2; xyz ++)
					{
						double CarryOver = 0.999;
						if(FormfindUnstressLoad == 2)//That is apply load
						{
							Force[Node][xyz] += Load[Node][xyz];
							CarryOver = 0.9999;
						}
						
						Velocity[Node][xyz] = CarryOver * Velocity[Node][xyz] + Force[Node][xyz] * delta_t / Mass[Node];
						x[Node][xyz] += Velocity[Node][xyz] * delta_t;
					}
				}
			}
			if(FormfindUnstressLoad == 0)//That is formfinding
			{
				for(int xyz = 0;xyz <= 2; xyz ++)Initial_x[Node][xyz] = x[Node][xyz];
			}
		}
		
		if(FormfindUnstressLoad == 0)//That is formfinding
		{
			double BaseLevel = x[0][2];
			for(int Node = 0;Node <= LastNode; Node ++)x[Node][2] -= BaseLevel;
			
			double RadiusNode0 = NodalRadius(0);
			double RadiusNode1 = NodalRadius(1);
			for(int ONLYxy = 0; ONLYxy <=2; ONLYxy ++)
				x[1][ONLYxy] *= (1.0 + 0.001 * (RadiusNode0 - RadiusNode1) / RadiusNode0);
			RotateNode(3);
			RotateNode(5);
			if(NodeType[1] != 1 || NodeType[3] != 1 || NodeType[5] != 1)cout<<"Moving wrong node\n";
			
			double RadiusLevel0 = NodalRadius(0);
			double RadiusLevel1 = NodalRadius(6);
			double RadiusTopMinus0 = NodalRadius(LastNode - 3);//3 extra nodes at base, 3 nodes on top
			double RadiusTopMinus1 = (NodalRadius(LastNode - 6) + NodalRadius(LastNode - 7)) / 2.0;//3 extra nodes at base, 3 nodes on top
			//double RadiusTopMinus1 = NodalRadius(LastNode - 7);//3 extra nodes at base, 3 nodes on top
			for(int Member = 0;Member <= LastMember;Member ++)
			{
				if(MemberType[Member] == LastOrdinaryMemberType + 1) InitialLength[Member] *= (1.0 + 0.0001 * (RadiusLevel1    / CurrentShrinkPerLevel - RadiusLevel0   ) / RadiusLevel0);
				if(MemberType[Member] == LastOrdinaryMemberType + 2) InitialLength[Member] *= (1.0 + 0.0001 * (RadiusTopMinus1 * CurrentShrinkPerLevel - RadiusTopMinus0) / RadiusTopMinus0);
			}
			
			for(int whichNode = 0;whichNode <= 2;whichNode ++)
			{
				for(int xyz = 0;xyz <= 2; xyz ++)x[LastNode - 2 + whichNode][xyz] = x[2 * whichNode][xyz];
				x[LastNode - 2 + whichNode][2] -= 0.5;
				NodeType[LastNode - 2 + whichNode] = 1;
			}
		}
	}
}

double NodalRadius(int Node)
{
	return sqrt(x[Node][0] * x[Node][0] + x[Node][1] * x[Node][1]);
}

void ForcesFromMembers(void)
{
	for(int Node = 0;Node <= LastNode; Node ++)for(int xyz = 0;xyz <= 2; xyz ++)Force[Node][xyz] = 0;
	
	for(int Member = 0;Member <= LastMember;Member ++)
	{
		double Length = CalculateMemberLength(Member);
		
		if(
		   FormfindUnstressLoad != 0 //That is NOT formfinding
		   || MemberType[Member] <= 1
		   || (MemberType[Member] >= LastOrdinaryMemberType + 1 && MemberType[Member] <= LastOrdinaryMemberType + 2))
		{
			double ThisInitialLength = InitialLength[Member];
			if(MemberType[Member] == 1) ThisInitialLength *= RelativeRingLength;
			Tension[Member] = EA[Member] * (Length - ThisInitialLength) / ThisInitialLength;
		}
		else
		{
			double PrestressAtThisLevel = Prestress * RelativeLevelFromTop(Member) * RelativeLevelFromTop(Member);
			if(MemberType[Member] == 2)
			{
				Tension[Member] = PrestressAtThisLevel;
				if(MemberLevel[Member] == LastLevel -1)Tension[Member] *= 1.6;
			}
			else
			{
				Tension[Member] = RelativePrestress * PrestressAtThisLevel;
				if(MemberLevel[Member] == LastLevel -1)Tension[Member] *= 0.2;
			}
			
			InitialLength[Member] = Length / (Tension[Member] / EA[Member] + 1.0);
		}
		
		if(Tension[Member] < - EulerLoad[Member])
		{
			Tension[Member] = - EulerLoad[Member];
			MemberBuckled[Member] = 1;
		}
		else MemberBuckled[Member] = 0;
		
		if(Tension[Member] > MemberStrengthInTension[Member])
		{
			Tension[Member] = 0.0;
			MemberBroken[Member] = 1;
		}
		else MemberBroken[Member] = 0;
		
		if(Member == PrestressMember && FormfindUnstressLoad == 1)
		{
			UnstressFactor *= 0.999;
			Tension[Member] *= UnstressFactor;
		}
		
		double TensionCoefficient = Tension[Member] / Length;
		
		for(int xyz = 0;xyz <= 2; xyz ++)
		{
			Force[End[Member][0]][xyz] += TensionCoefficient * (x[End[Member][1]][xyz] - x[End[Member][0]][xyz]);
			Force[End[Member][1]][xyz] += TensionCoefficient * (x[End[Member][0]][xyz] - x[End[Member][1]][xyz]);
		}
	}
}

double CalculateMemberLength(int Member)
{
	double LengthSq = 0.0;
	for(int xyz = 0;xyz <= 2; xyz ++)LengthSq += (x[End[Member][1]][xyz] - x[End[Member][0]][xyz]) * (x[End[Member][1]][xyz] - x[End[Member][0]][xyz]);
	return sqrt(LengthSq);
}

int Permute(int MyNumber)
{
	int MyAnswer = MyNumber;
	while (MyAnswer > 5)
	{
		MyAnswer -= 6;
	}
	while (MyAnswer < 0)
	{
		MyAnswer += 6;
	}
	return MyAnswer;
}

int LastLevelPermute(int MyNumber)
{
	int MyAnswer = MyNumber;
	while (MyAnswer > 2)
	{
		MyAnswer -= 3;
	}
	while (MyAnswer < 0)
	{
		MyAnswer += 3;
	}
	return MyAnswer;
}

void RotateNode(int Node)
{
	if(RotationalSymmetryType[Node] >= 2)
	{
		if(RotationalSymmetryType[Node] == 2 || RotationalSymmetryType[Node] == 3 || RotationalSymmetryType[Node] == 6)//6 is at top only
		{
			int otherNode = Node - 2;
			if(RotationalSymmetryType[Node] == 6)otherNode = Node - 1;//6 is at top only
			x[Node][0] = ( - x[otherNode][0] - RightOrLeftHanded * root3 * x[otherNode][1])/ 2.0;
			x[Node][1] = ( - x[otherNode][1] + RightOrLeftHanded * root3 * x[otherNode][0])/ 2.0;
			x[Node][2] =     x[otherNode][2];
		}
		else
		{
			int otherNode = Node - 4;
			if(RotationalSymmetryType[Node] == 7)otherNode = Node - 2;//7 is at top only
			x[Node][0] = ( - x[otherNode][0] + RightOrLeftHanded * root3 * x[otherNode][1])/ 2.0;
			x[Node][1] = ( - x[otherNode][1] - RightOrLeftHanded * root3 * x[otherNode][0])/ 2.0;
			x[Node][2] =     x[otherNode][2];
		}
	}
}

void WriteDXF(void)
{
	ofstream Madeleine("TensegrityTower.dxf");
	
	Madeleine<<"0\nSECTION\n2\nENTITIES\n";
	
	for(int Member = 0;Member <= LastMember;Member ++)
	{
		double memberlength,Xx,Xy,Xz,Yx,Yy,Yz,Zx,Zy,Zz,thislength;
		
		memberlength=sqrt(
						  (x[End[Member][1]][0]-x[End[Member][0]][0])*
						  (x[End[Member][1]][0]-x[End[Member][0]][0])+
						  (x[End[Member][1]][1]-x[End[Member][0]][1])*
						  (x[End[Member][1]][1]-x[End[Member][0]][1])+
						  (x[End[Member][1]][2]-x[End[Member][0]][2])*
						  (x[End[Member][1]][2]-x[End[Member][0]][2]));
		Zx=(x[End[Member][1]][0]-x[End[Member][0]][0])/memberlength;
		Zy=(x[End[Member][1]][1]-x[End[Member][0]][1])/memberlength;
		Zz=(x[End[Member][1]][2]-x[End[Member][0]][2])/memberlength;
		if(fabs(Zx)<1.0/64.0&&fabs(Zy)<1.0/64.0)
		{
			Xx=Zz;Xy=0.0;Xz=-Zx;
		}
		else
		{
			Xx=-Zy;Xy=+Zx;Xz=0.0;
		}
		thislength=sqrt(Xx*Xx+Xy*Xy+Xz*Xz);
		Xx=Xx/thislength;Xy=Xy/thislength;Xz=Xz/thislength;
		Yx=Zy*Xz-Zz*Xy;
		Yy=Zz*Xx-Zx*Xz;
		Yz=Zx*Xy-Zy*Xx;
		Madeleine<<"0\nCIRCLE\n8\n"<<MemberType[Member]<<"\n";
		Madeleine<<"39\n"<<memberlength<<"\n";
		Madeleine<<"10\n"<<x[End[Member][0]][0]*Xx+x[End[Member][0]][1]*Xy+x[End[Member][0]][2]*Xz<<"\n";
		Madeleine<<"20\n"<<x[End[Member][0]][0]*Yx+x[End[Member][0]][1]*Yy+x[End[Member][0]][2]*Yz<<"\n";
		Madeleine<<"30\n"<<x[End[Member][0]][0]*Zx+x[End[Member][0]][1]*Zy+x[End[Member][0]][2]*Zz<<"\n";
		Madeleine<<"40\n"<<Diameter[Member] / 2.0<<"\n";
		Madeleine<<"210\n"<<Zx<<"\n";
		Madeleine<<"220\n"<<Zy<<"\n";
		Madeleine<<"230\n"<<Zz<<"\n";
		Madeleine<<"62\n"<<MemberType[Member]<<"\n";
	}
	
	double SupportRadius = 0.05;
	for(int Node = 0;Node <= LastNode; Node ++)
	{
		int m = 18;
		int n = 9;
		if(NodeType[Node] != 0)
		{
			for(int i = 0; i <= m - 1; i ++)
			{
				for(int j = 0; j <= n - 1; j ++)
				{
					Madeleine<<"0\n3DFACE\n8\nSupports\n";
					int istep = 0;
					int jstep = 0;
					Madeleine<<"10\n"<<x[Node][0] + SupportRadius * cos(PI * double(2 * (j + jstep) - n) / double(2 * n)) * cos(PI * double(2 * (i + istep)) / double(m))<<"\n";
					Madeleine<<"20\n"<<x[Node][1] + SupportRadius * cos(PI * double(2 * (j + jstep) - n) / double(2 * n)) * sin(PI * double(2 * (i + istep)) / double(m))<<"\n";
					Madeleine<<"30\n"<<x[Node][2] + SupportRadius * sin(PI * double(2 * (j + jstep) - n) / double(2 * n))<<"\n";
					
					istep = 1;
					Madeleine<<"11\n"<<x[Node][0] + SupportRadius * cos(PI * double(2 * (j + jstep) - n) / double(2 * n)) * cos(PI * double(2 * (i + istep)) / double(m))<<"\n";
					Madeleine<<"21\n"<<x[Node][1] + SupportRadius * cos(PI * double(2 * (j + jstep) - n) / double(2 * n)) * sin(PI * double(2 * (i + istep)) / double(m))<<"\n";
					Madeleine<<"31\n"<<x[Node][2] + SupportRadius * sin(PI * double(2 * (j + jstep) - n) / double(2 * n))<<"\n";
					
					jstep = 1;
					Madeleine<<"12\n"<<x[Node][0] + SupportRadius * cos(PI * double(2 * (j + jstep) - n) / double(2 * n)) * cos(PI * double(2 * (i + istep)) / double(m))<<"\n";
					Madeleine<<"22\n"<<x[Node][1] + SupportRadius * cos(PI * double(2 * (j + jstep) - n) / double(2 * n)) * sin(PI * double(2 * (i + istep)) / double(m))<<"\n";
					Madeleine<<"32\n"<<x[Node][2] + SupportRadius * sin(PI * double(2 * (j + jstep) - n) / double(2 * n))<<"\n";
					
					istep = 0;
					Madeleine<<"13\n"<<x[Node][0] + SupportRadius * cos(PI * double(2 * (j + jstep) - n) / double(2 * n)) * cos(PI * double(2 * (i + istep)) / double(m))<<"\n";
					Madeleine<<"23\n"<<x[Node][1] + SupportRadius * cos(PI * double(2 * (j + jstep) - n) / double(2 * n)) * sin(PI * double(2 * (i + istep)) / double(m))<<"\n";
					Madeleine<<"33\n"<<x[Node][2] + SupportRadius * sin(PI * double(2 * (j + jstep) - n) / double(2 * n))<<"\n";
					
					Madeleine<<"62\n1\n";
				}
			}
		}
	}
	
	Madeleine<<"0\nENDSEC\n0\nEOF\n";
	Madeleine.close();
	cout<<"DXF file written.\n";
}

void StrutProperties(int Member)
{
	InitialLength[Member] = CalculateMemberLength(Member);
	
	double EStrut = 210.0e9;
	double SigmaMaxStrut = 350.0e6;
	double SlendernessControl = 2.0;
	double SlendernessRatio = 120.0 * (SlendernessControl - RelativeLevelFromTop(Member)) / (SlendernessControl - 1.0);
	double DiaStrut = 2.0 * sqrt(2.0) * InitialLength[Member] / SlendernessRatio;
	double StrutWallThickness = 0.003 * RelativeLevelFromTop(Member);
	double StrutInnerDiameter = DiaStrut - 2.0 * StrutWallThickness;
	double AreaStrut = PI * (DiaStrut * DiaStrut - StrutInnerDiameter * StrutInnerDiameter) / 4.0;
	double IStrut    = PI * (DiaStrut * DiaStrut * DiaStrut * DiaStrut - StrutInnerDiameter * StrutInnerDiameter * StrutInnerDiameter * StrutInnerDiameter) / 64.0;
	double DensitySrut = 7850.0;
	
	Diameter[Member] = DiaStrut;
	EA[Member] = EStrut * AreaStrut;
	EI[Member] = EStrut * IStrut;
	MemberStrengthInTension[Member] = SigmaMaxStrut * AreaStrut;
	RoA[Member] = DensitySrut * AreaStrut;
	
	MemberSetUp(Member);
}

void TieProperties(int Member)
{
	InitialLength[Member] = CalculateMemberLength(Member);
	
	double ETie = 210.0e9;
	double SigmaMaxTie = 1000.0e6;
	//double DiaTie = 0.012 * RelativeLevelFromTop(Member);
	double DiaTie = 0.007 * RelativeLevelFromTop(Member);
	double AreaTie   = PI *  DiaTie   * DiaTie / 4.0;
	double DensityTie = 7850.0;
	
	Diameter[Member] = DiaTie;
	EA[Member] = ETie * AreaTie;
	EI[Member] = 0.0;
	MemberStrengthInTension[Member] = SigmaMaxTie * AreaTie;
	RoA[Member] = DensityTie * AreaTie;
	
	MemberSetUp(Member);
}

void MemberSetUp(int Member)
{
	MemberBuckled[Member] = 0;
	MemberBroken[Member] = 0;
	
	Mass[End[Member][0]] += RoA[Member] * InitialLength[Member] / 2.0;
	Mass[End[Member][1]] += RoA[Member] * InitialLength[Member] / 2.0;
	
	NodalArea[End[Member][0]] += Diameter[Member] * InitialLength[Member] / 2.0;
	NodalArea[End[Member][1]] += Diameter[Member] * InitialLength[Member] / 2.0;
	
	EulerLoad[Member] = PI * PI *  EI[Member] / (InitialLength[Member] * InitialLength[Member]);
}

void StartingGeometry(void)
{
	RestartWhenReady = 0;
	InitialStateWhenReady = 0;
	FormfindUnstressLoad = 0;
	
	CurrentShrinkPerLevel = ShrinkPerLevel;
	
	double a_basic = 1.0 * SizeMultiplier;
	
	double heightratio = 1.45;
	
	double aminus1 = a_basic / ShrinkPerLevel;
	
	a[0] = aminus1 * ShrinkPerLevel;
	Height[0] = 0.0;
	
	for(int Level = 0; Level <= LastLevel; Level ++)
	{
		if(Level != 0)
		{
			a[Level] = a[Level - 1] * ShrinkPerLevel;
			Height[Level] = Height[Level - 1] + c[Level - 1];
		}
		c[Level] = a[Level] * heightratio;
	}
	
	int LevelEven = 0;
	for(int Level = 0; Level <= LastLevel; Level ++)
	{
		if(LevelEven == 1)LevelEven = 0;else LevelEven = 1;
		for(int Rotation = 0; Rotation <= 5; Rotation ++)
		{
			if(Level == LastLevel && Rotation == 3)break;
			
			double angle = double(Rotation) * PI / 3.0;
			
			if(Level == LastLevel) angle = 2.0 * angle - double(LevelEven) * PI / 3.0;
			
			int Node = 6 * Level + Rotation;
			
			Mass[Node] = 0.0;
			NodalArea[Node] = 0.0;
			
			if(Level == 0 && (Rotation == 1 || Rotation == 3 || Rotation == 5))NodeType[Node] = 1;else NodeType[Node] = 0;
			
			RotationalSymmetryType[Node] = Rotation;
			
			if(Level == LastLevel && RotationalSymmetryType[Node] != 0) RotationalSymmetryType[Node] += 5;
			
			angle *= RightOrLeftHanded;
			
			if(Rotation == 0 || Rotation == 1)
			{
				x[Node][0] = a[Level] * cos(angle);
				x[Node][1] = a[Level] * sin(angle);
				x[Node][2] = Height[Level];
				if(Level == LastLevel)MaxHalfDimension = 0.8 * x[Node][2] / 2.0;
			}
			else RotateNode(Node);
		}
	}
	
	for(int Level = 0; Level <= LastLevel; Level ++)
	{
		for(int Rotation = 0; Rotation <= 5; Rotation ++)
		{
			if(Level == LastLevel && Rotation == 3)break;
			
			int Member = MembersPerLevel * Level + Rotation;
			if(Level == LastLevel)Member -= 3 * (LastOrdinaryMemberType - 2);
			End[Member][0] = 6 * Level + Rotation;
			
			if(Level != LastLevel)End[Member][1] = 6 * Level + Permute(Rotation + 1);
			else End[Member][1] = 6 * Level + LastLevelPermute(Rotation + 1);
			
			MemberLevel[Member] = Level;TieProperties(Member);
			
			if(Level != 0 && Level != LastLevel)MemberType[Member] = 1;
			else
			{
				if(Level == 0)MemberType[Member] = LastOrdinaryMemberType + 1;
				else MemberType[Member] = LastOrdinaryMemberType + 2;
			}
		}
	}
	
	int MemberRotate = 1;
	int WhichWay = - 1;
	PrestressMember = - 1;
	for(int Level = 0; Level <= LastLevel - 1; Level ++)
	{
		MemberRotate ++;if(MemberRotate > 1)MemberRotate -= 2;
		WhichWay = - WhichWay;
		for(int DoubleRotation = 0; DoubleRotation <= 2; DoubleRotation ++)
		{
			int Member = MembersPerLevel * Level + 6 + DoubleRotation;
			End[Member][0] = 6 * Level + Permute(MemberRotate + 2 * DoubleRotation);
			
			if(Level != LastLevel - 1)End[Member][1] = 6 * (Level + 1) + Permute(MemberRotate + 2 * WhichWay + 2 * DoubleRotation);
			else End[Member][1] = 6 * (Level + 1) + LastLevelPermute(MemberRotate + 1 * WhichWay + 1 * DoubleRotation);
			
			MemberLevel[Member] = Level;StrutProperties(Member);
			
			MemberType[Member] = 0;
		}
		
		for(int DoubleRotation = 0; DoubleRotation <= 2; DoubleRotation ++)
		{
			int Member = MembersPerLevel * Level + 9 + DoubleRotation;
			End[Member][0] = 6 * Level + Permute(MemberRotate + 2 * DoubleRotation);
			
			if(Level != LastLevel - 1)End[Member][1] = 6 * (Level + 1) + Permute(MemberRotate + 2 * DoubleRotation);
			else End[Member][1] = 6 * (Level + 1) + LastLevelPermute(MemberRotate + 1 * DoubleRotation);
			
			MemberLevel[Member] = Level;TieProperties(Member);
			
			MemberType[Member] = 2;
		}
		
		if(Level <= LastLevel - 2)
		{
			for(int DoubleRotation = 0; DoubleRotation <= 2; DoubleRotation ++)
			{
				int Member = MembersPerLevel * Level + 12 + DoubleRotation;
				if(PrestressMember == - 1)PrestressMember = Member;
				End[Member][0] = 6 * Level + Permute(MemberRotate + 2 * DoubleRotation + 1);
				End[Member][1] = 6 * (Level + 1) + Permute(MemberRotate + 2 * DoubleRotation + 1);
				
				MemberLevel[Member] = Level;TieProperties(Member);
				
				MemberType[Member] = 3;
			}
			
			if(LastOrdinaryMemberType >= 4)
			{
				for(int DoubleRotation = 0; DoubleRotation <= 2; DoubleRotation ++)
				{
					int Member = MembersPerLevel * Level + 15 + DoubleRotation;
					End[Member][0] = 6 * Level + Permute(MemberRotate + 2 * DoubleRotation - WhichWay);
					End[Member][1] = 6 * (Level + 1) + Permute(MemberRotate + 2 * DoubleRotation - 2 * WhichWay);
					
					MemberLevel[Member] = Level;TieProperties(Member);
					
					MemberType[Member] = 4;
				}
			}
			if(LastOrdinaryMemberType >= 5)
			{
				for(int DoubleRotation = 0; DoubleRotation <= 2; DoubleRotation ++)
				{
					int Member = MembersPerLevel * Level + 18 + DoubleRotation;
					End[Member][0] = 6 * Level + Permute(MemberRotate + 2 * DoubleRotation);
					End[Member][1] = 6 * (Level + 1) + Permute(MemberRotate + 2 * DoubleRotation - WhichWay);
					
					MemberLevel[Member] = Level;TieProperties(Member);
					
					MemberType[Member] = 5;
				}
			}
		}
	}
	
	for(int whichNode = 0;whichNode <= 2;whichNode ++)
	{
		int Member = LastMember - 5 + whichNode;
		End[Member][0] = LastNode - whichNode - 3;
		int whichtopcables = 1;
		if(whichtopcables == 0)
		{
		if(whichNode == 0)End[Member][1] = LastNode - 7;
		 if(whichNode == 1)End[Member][1] = LastNode - 9;
		 if(whichNode == 2)End[Member][1] = LastNode - 11;
		 }
		 else
		 {
		if(whichNode == 0)End[Member][1] = LastNode - 9;
		if(whichNode == 1)End[Member][1] = LastNode - 11;
		if(whichNode == 2)End[Member][1] = LastNode - 7;
		}
		
		MemberLevel[Member] = LastLevel - 1;TieProperties(Member);
		
		MemberType[Member] = 3;
	}
	
	for(int whichNode = 0;whichNode <= 2;whichNode ++)
	{
		int Member = LastMember - 2 + whichNode;
		End[Member][0] = LastNode - 2 + whichNode;
		End[Member][1] = 2 * whichNode;
		
		int Level = 0;
		MemberLevel[Member] = Level;TieProperties(Member);
		
		MemberType[Member] = LastOrdinaryMemberType + 3;
	}
	
	for(int Node = 0;Node <= LastNode; Node ++)
	{
		for(int xyz = 0;xyz <= 2; xyz ++)Load[Node][xyz] = 0.0;
		Load[Node][2] = - gravity * Mass[Node];
		Load[Node][1] = WindPressureTimesDragCoeff * NodalArea[Node];
		
		for(int xyz = 0;xyz <= 2; xyz ++)Initial_x[Node][xyz] = x[Node][xyz];
	}
}

double RelativeLevelFromTop(int Member)
{
	//int TopAdd = 1;
	int TopAdd = 3;
	return double(LastLevel + TopAdd - MemberLevel[Member]) / double(LastLevel + TopAdd);
}

void InitialValues(void)
{
	ShrinkPerLevel = 0.96;
	RelativeRingLength = 1.0;
	RelativePrestress = 1.0;
	SizeMultiplier = 1.0;
}
