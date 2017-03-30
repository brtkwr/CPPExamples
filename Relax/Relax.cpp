//Written by Chris J K Williams, University of Bath, UK

#include "data.h"
#include "Graphics.h"

double myFunction(double);
double NonLinearElastic(int,int,double,double);

void aToA(void);
void MemIntForMo(int,double&,double&,double&,double&,double&,double&);
void MemberContributions(void);
void NodeMovement(void);
void ZeroStiffness(int);

int main(int argc, char* argv[])
{
	PI=4.0*atan(1.0);
	if(ReadBasicData()==0)
	{cout<<"Press 'Return' to close\n";char Pause;Pause=getchar();return 0;}
	MaxHalfDimension=0.0;
	for(int i=0;i<=2;i++)
	{
		if(MaxHalfDimension<fabs(MaxCoord[i]-MinCoord[i])/2.0)MaxHalfDimension=fabs(MaxCoord[i]-MinCoord[i])/2.0;
	}
	LoadScale=0.1*MaxHalfDimension/TypicalLoad;
	
	zRotation=0.0;
	HorizontalRotation=0.0;
	zRotationIncrement=45.0;HorizontalRotationIncrement=-180.0*acos(1.0/sqrt(3.0))/PI;
	
	xTrans=0.0;yTrans=0.0;zTrans=0.0;
	
	for(int i=0;i<=2;i++)AverageCoord[i]=(MaxCoord[i]+MinCoord[i])/2.0;
	
	PictureCycle=-1;//Negative value stops picture saving
	LastCalculationCycle=5;
	FirstCycle=1;
		
		glutInit(&argc,argv);
		SetUpGraphics(1.0,1.0,1.0);
		cout<<"Press 'Return' to close\n";char Pause;Pause=getchar();return 0;
}

static void Draw(void)
{
	double CoordinatesToPlot[3],OtherCoordinatesToPlot[3];
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	for(CalculationCycle=0;CalculationCycle<=LastCalculationCycle;CalculationCycle++)
	{
		aToA();
		MemberContributions();
		NodeMovement();
		FirstCycle=0;
	}
	
	glLineWidth(2.0);
	for(int Member=0;Member<=LastMember;Member++)
	{
		if(MemberExists[Member]==1)
		{
			if(NonLinear==0||ShowPeak==0)
			{
			if(ShowTensionOnlyElements==1&&TensionOnlyMember[MemberType[Member]]==1)glColor4f(0.0,0.0,1.0,1.0);
			else
			{
				if(MemberType[Member]==MemberTypeToShow)glColor4f(1.0,0.0,0.0,1.0);
				else glColor4f(0.0,0.0,0.0,1.0);
			}
			}
			else
			{
			if(MemberReachedPeak[Member]==1)glColor4f(1.0,0.0,0.0,1.0);
			else glColor4f(0.0,0.0,0.0,0.3);
			}
			glBegin(GL_LINE_STRIP);
			for(int i=0;i<=2;i++)CoordinatesToPlot[i]=Coord[End[0][Member]][i]-AverageCoord[i];
			glVertex3dv(CoordinatesToPlot);
			for(int i=0;i<=2;i++)CoordinatesToPlot[i]=Coord[End[1][Member]][i]-AverageCoord[i];
			glVertex3dv(CoordinatesToPlot);
			glEnd();
			
			if(ShowMemberLocalAxis==1&&TensionOnlyMember[MemberType[Member]]==0)
			{
				double MemberLocalAxisScale=MaxHalfDimension/100.0;
				glColor4f(0.0,0.0,1.0,0.5);
				glBegin(GL_LINE_STRIP);
				for(int i=0;i<=2;i++)
					CoordinatesToPlot[i]=Coord[End[0][Member]][i]+MemberLocalAxisScale*Rotatedy[0][Member][i]-AverageCoord[i];
				for(int i=0;i<=2;i++)
					OtherCoordinatesToPlot[i]=Coord[End[0][Member]][i]-MemberLocalAxisScale*Rotatedy[0][Member][i]-AverageCoord[i];
				glVertex3dv(CoordinatesToPlot);glVertex3dv(OtherCoordinatesToPlot);
				glEnd();
				glPointSize(2.0);
				glBegin(GL_POINTS);
				glVertex3dv(CoordinatesToPlot);glVertex3dv(OtherCoordinatesToPlot);
				glEnd();
				
				glBegin(GL_LINE_STRIP);
				for(int i=0;i<=2;i++)
					CoordinatesToPlot[i]=Coord[End[1][Member]][i]-MemberLocalAxisScale*Rotatedy[1][Member][i]-AverageCoord[i];
				for(int i=0;i<=2;i++)
					OtherCoordinatesToPlot[i]=Coord[End[1][Member]][i]+MemberLocalAxisScale*Rotatedy[1][Member][i]-AverageCoord[i];
				glVertex3dv(CoordinatesToPlot);glVertex3dv(OtherCoordinatesToPlot);
				glEnd();
				glBegin(GL_POINTS);
				glVertex3dv(CoordinatesToPlot);glVertex3dv(OtherCoordinatesToPlot);
				glEnd();
			}
			if(ShowInitial==1)
			{
			glColor4f(1.0,0.5,0.5,0.8);
			glBegin(GL_LINE_STRIP);
			for(int i=0;i<=2;i++)CoordinatesToPlot[i]=InitialCoord[End[0][Member]][i]-AverageCoord[i];
			glVertex3dv(CoordinatesToPlot);
			for(int i=0;i<=2;i++)CoordinatesToPlot[i]=InitialCoord[End[1][Member]][i]-AverageCoord[i];
			glVertex3dv(CoordinatesToPlot);
			glEnd();
			}
		}
	}
	
	for(int Node=0;Node<=LastNode;Node++)
	{
		if(NodeExists[Node]==1)
		{
			glColor4f(0.0,0.0,0.0,1.0);
			if(NodeDispType[Node]==0)glPointSize(3.0);
			else glPointSize(5.0);
			glBegin(GL_POINTS);
			for(int i=0;i<=2;i++)CoordinatesToPlot[i]=Coord[Node][i]-AverageCoord[i];
			glVertex3dv(CoordinatesToPlot);
			glEnd();
			if(NodeDispType[Node]!=0)
			{
				glPointSize(3.0);
				
				if(NodeDispType[Node]==1)glColor4f(1.0,0.0,0.0,1.0);
				if(NodeDispType[Node]==2)glColor4f(0.0,1.0,0.0,1.0);
				if(NodeDispType[Node]==4)glColor4f(0.0,0.0,1.0,1.0);
				if(NodeDispType[Node]>4)glColor4f(1.0,1.0,1.0,1.0);
				
				glBegin(GL_POINTS);
				glVertex3dv(CoordinatesToPlot);
				glEnd();
			}
		if(ShowLoads==1)
		{
			glColor4f(0.0,1.0,0.0,0.5);
			glBegin(GL_LINE_STRIP);
			for(int i=0;i<=2;i++)
				CoordinatesToPlot[i]=Coord[Node][i]-AverageCoord[i];
			for(int i=0;i<=2;i++)
				OtherCoordinatesToPlot[i]=CoordinatesToPlot[i]+LoadScale*AppliedSafetyFact*AppliedNodeLoad[Node][i];
			glVertex3dv(CoordinatesToPlot);glVertex3dv(OtherCoordinatesToPlot);
			glEnd();
			glPointSize(2.0);
			glBegin(GL_POINTS);
			glVertex3dv(OtherCoordinatesToPlot);
			glEnd();
		}
		if(ShowInitial==1&&NodeDispType[Node]==0)
			{
			glPointSize(3.0);
			glColor4f(1.0,0.5,0.5,0.8);
			glBegin(GL_POINTS);
			for(int i=0;i<=2;i++)CoordinatesToPlot[i]=InitialCoord[Node][i]-AverageCoord[i];
			glVertex3dv(CoordinatesToPlot);
			glEnd();
			}
		}
	}
	
	Controls();
	
	glutSwapBuffers();
	
	if(PictureCycle>=0)
	{
		if(PictureCycle==0)SavePicture();
		PictureCycle++;
		if(PictureCycle>100)PictureCycle=0;
	}
}

void aToA(void)
{
	for(int Node=0;Node<=LastNode;Node++)
	{
		if(NodeExists[Node]==1)
		{
			double asquared=0.0;
			for(int i=0;i<=2;i++)asquared+=a[Node][i]*a[Node][i];
			double tempdouble=1.0-asquared;
			for(int i=0;i<=2;i++)A[Node][i][i]=tempdouble;
			
			A[Node][0][1]=-2.0*a[Node][2];
			A[Node][1][2]=-2.0*a[Node][0];
			A[Node][2][0]=-2.0*a[Node][1];
			
			A[Node][1][0]=-A[Node][0][1];
			A[Node][2][1]=-A[Node][1][2];
			A[Node][0][2]=-A[Node][2][0];
			
			tempdouble=1.0+asquared;
			for(int i=0;i<=2;i++)
			{
				for(int j=0;j<=2;j++)
				{
					A[Node][i][j]+=2.0*a[Node][i]*a[Node][j];
					A[Node][i][j]=A[Node][i][j]/tempdouble;
				}
			}
		}
	}
}

void MemIntForMo(int ThisMember,double& ThisT,double& ThisMx0,double& ThisMy0,double& ThisMx1,double& ThisMy1,double& ThisMphi)
{
	if(OrigLength[ThisMember]==0.0)cout<<"Member "<<ThisMember<<" has zero ooriginal length\n";
	for(int i=0;i<=2;i++)
	{
		p[i]=Coord[End[1][ThisMember]][i]-Coord[End[0][ThisMember]][i];
		Rotatedx[0][ThisMember][i]=0.0;
		Rotatedy[0][ThisMember][i]=0.0;
		Rotatedx[1][ThisMember][i]=0.0;
		Rotatedy[1][ThisMember][i]=0.0;
		for(int j=0;j<=2;j++)
		{
			Rotatedx[0][ThisMember][i]+=A[End[0][ThisMember]][i][j]*x[0][ThisMember][j];
			Rotatedy[0][ThisMember][i]+=A[End[0][ThisMember]][i][j]*y[0][ThisMember][j];
			Rotatedx[1][ThisMember][i]+=A[End[1][ThisMember]][i][j]*x[1][ThisMember][j];
			Rotatedy[1][ThisMember][i]+=A[End[1][ThisMember]][i][j]*y[1][ThisMember][j];
		}
	}
	
	double Thetax1=0.0;
	double Thetay1=0.0;
	double Thetax2=0.0;
	double Thetay2=0.0;
	for(int i=0;i<=2;i++)
	{
		Thetax1+=p[i]*Rotatedy[0][ThisMember][i];
		Thetay1-=p[i]*Rotatedx[0][ThisMember][i];
		Thetax2+=p[i]*Rotatedy[1][ThisMember][i];
		Thetay2-=p[i]*Rotatedx[1][ThisMember][i];
	}
	Thetax1=Thetax1/OrigLength[ThisMember];
	Thetay1=Thetay1/OrigLength[ThisMember];
	Thetax2=Thetax2/OrigLength[ThisMember];
	Thetay2=Thetay2/OrigLength[ThisMember];
	
	double esuba=(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]-OrigLength[ThisMember]*OrigLength[ThisMember])
		/(2.0*OrigLength[ThisMember]);
	
	double esubb=(4.0*(Thetax1*Thetax1+Thetay1*Thetay1+Thetax2*Thetax2+Thetay2*Thetay2)
				  -2.0*(Thetax1*Thetax2+Thetay1*Thetay2))*OrigLength[ThisMember]/60.0;
	
	double phi=0.0;
	for(int i=0;i<=2;i++)phi+=Rotatedx[0][ThisMember][i]*Rotatedy[1][ThisMember][i]-Rotatedx[1][ThisMember][i]*Rotatedy[0][ThisMember][i];
	phi=phi/2.0;
	
	ThisT=EA[MemberType[ThisMember]]*(esuba+esubb)/OrigLength[ThisMember];
	if(TensionOnlyMember[MemberType[ThisMember]]==1&&ThisT<0.0)ThisT=0.0;
	
	ThisMx0=(EIxx[MemberType[ThisMember]]/OrigLength[ThisMember])*(4.0*Thetax1+2.0*Thetax2);
	ThisMy0=(EIyy[MemberType[ThisMember]]/OrigLength[ThisMember])*(4.0*Thetay1+2.0*Thetay2);
	ThisMx1=(EIxx[MemberType[ThisMember]]/OrigLength[ThisMember])*(4.0*Thetax2+2.0*Thetax1);
	ThisMy1=(EIyy[MemberType[ThisMember]]/OrigLength[ThisMember])*(4.0*Thetay2+2.0*Thetay1);
	
	double tempdouble=ThisT*OrigLength[ThisMember]/30.0;
	ThisMx0+=tempdouble*(4.0*Thetax1-Thetax2);
	ThisMy0+=tempdouble*(4.0*Thetay1-Thetay2);
	ThisMx1+=tempdouble*(4.0*Thetax2-Thetax1);
	ThisMy1+=tempdouble*(4.0*Thetay2-Thetay1);
	
	ThisMphi=GJ[MemberType[ThisMember]]*phi/OrigLength[ThisMember];
	
	if(NonLinear!=0)
	{
		MemberReachedPeak[ThisMember]=0;
		
		ThisT=NonLinearElastic(ThisMember,MemberType[ThisMember],ThisT,AxialPeak[MemberType[ThisMember]]);
		
		ThisMx0=NonLinearElastic(ThisMember,MemberType[ThisMember],ThisMx0,MPeakx[MemberType[ThisMember]]);
		ThisMy0=NonLinearElastic(ThisMember,MemberType[ThisMember],ThisMy0,MPeaky[MemberType[ThisMember]]);
		
		ThisMx1=NonLinearElastic(ThisMember,MemberType[ThisMember],ThisMx1,MPeakx[MemberType[ThisMember]]);
		ThisMy1=NonLinearElastic(ThisMember,MemberType[ThisMember],ThisMy1,MPeaky[MemberType[ThisMember]]);
		
		ThisMphi=
			  NonLinearElastic(ThisMember,MemberType[ThisMember],ThisMphi,MPeakTorsion[MemberType[ThisMember]]);
	}
}

void MemberContributions(void)
{
	double CoordinatesToPlot[3],OtherCoordinatesToPlot[3];
	
	for(int i=0;i<=2;i++)TotalFactoredLoad[i]=0.0;
	
	for(int Node=0;Node<=LastNode;Node++)
	{
		if(NodeExists[Node]==1)
		{
			for(int i=0;i<=2;i++)
			{
				double ThisLoad=AppliedSafetyFact*AppliedNodeLoad[Node][i];
				NodeForce[Node][i]=ThisLoad;
				TotalFactoredLoad[i]+=ThisLoad;
				NodeMom[Node][i]=AppliedSafetyFact*AppliedNodeMom[Node][i];
				for(int j=0;j<=2;j++)Stiff[Node][i][j]=0.0;
			}
			RotStiff[Node]=0.0;
		}
	}
	
	for(int Member=0;Member<=LastMember;Member++)
	{
		if(MemberExists[Member]==1)
		{
			if(End[0][Member]!=End[1][Member])
			{
				double T,Mx0,My0,Mx1,My1,Mphi;
				
				double MemberWeight=OwnWtSafetyFact*rogA[MemberType[Member]]*OrigLength[Member];
				NodeForce[End[0][Member]][2]-=MemberWeight/2.0;
				NodeForce[End[1][Member]][2]-=MemberWeight/2.0;
				TotalFactoredLoad[2]-=MemberWeight;
				
				if(ShowLoads==1&&CalculationCycle==LastCalculationCycle)
				{
					glColor4f(0.0,1.0,0.0,0.5);
					glBegin(GL_LINE_STRIP);
					for(int i=0;i<=2;i++)
						CoordinatesToPlot[i]=(Coord[End[0][Member]][i]+Coord[End[1][Member]][i])/2.0-AverageCoord[i];
					for(int i=0;i<=2;i++)
						OtherCoordinatesToPlot[i]=CoordinatesToPlot[i];
					OtherCoordinatesToPlot[2]-=LoadScale*MemberWeight;
					glVertex3dv(CoordinatesToPlot);glVertex3dv(OtherCoordinatesToPlot);
					glEnd();
					glPointSize(2.0);
					glBegin(GL_POINTS);
					glVertex3dv(OtherCoordinatesToPlot);
					glEnd();
				}
				
				MemIntForMo(Member,T,Mx0,My0,Mx1,My1,Mphi);
				
				for(int i=0;i<=2;i++)
				{
					double tempdouble=(T*p[i]+Mx0*Rotatedy[0][Member][i]+Mx1*Rotatedy[1][Member][i]-My0*Rotatedx[0][Member][i]-My1*Rotatedx[1][Member][i])/OrigLength[Member];
					NodeForce[End[0][Member]][i]+=tempdouble;
					NodeForce[End[1][Member]][i]-=tempdouble;
					
					int j,k;
					
					if(i==0){j=1;k=2;}
					if(i==1){j=2;k=0;}
					if(i==2){j=0;k=1;}
					
					tempdouble=Mphi*(Rotatedx[0][Member][j]*Rotatedy[1][Member][k]-Rotatedy[0][Member][j]*Rotatedx[1][Member][k]
								   -(Rotatedx[0][Member][k]*Rotatedy[1][Member][j]-Rotatedy[0][Member][k]*Rotatedx[1][Member][j]))/2.0;
					
					NodeMom[End[0][Member]][i]-=(Mx0*(p[k]*Rotatedy[0][Member][j]-p[j]*Rotatedy[0][Member][k])
												-My0*(p[k]*Rotatedx[0][Member][j]-p[j]*Rotatedx[0][Member][k]))/OrigLength[Member]
						+tempdouble;
					NodeMom[End[1][Member]][i]-=(Mx1*(p[k]*Rotatedy[1][Member][j]-p[j]*Rotatedy[1][Member][k])
												-My1*(p[k]*Rotatedx[1][Member][j]-p[j]*Rotatedx[1][Member][k]))/OrigLength[Member]
						-tempdouble;
				}
				
				if(T>0.0)
				{
					double tempdouble=T/OrigLength[Member];
					for(int i=0;i<=2;i++)
					{
						Stiff[End[0][Member]][i][i]+=tempdouble;
						Stiff[End[1][Member]][i][i]+=tempdouble;
					}
					
					tempdouble=T*OrigLength[Member]/7.5;
					RotStiff[End[0][Member]]+=tempdouble;
					RotStiff[End[1][Member]]+=tempdouble;
				}
				
				double tempdouble=1.0/(OrigLength[Member]*OrigLength[Member]*OrigLength[Member]);
				for(int i=0;i<=2;i++)
				{
					if(AxStiffType==1)
					{
						Stiff[End[0][Member]][i][i]+=tempdouble*EA[MemberType[Member]];
						Stiff[End[1][Member]][i][i]+=tempdouble*EA[MemberType[Member]];
					}
					for(int j=0;j<=2;j++)
					{
						if(AxStiffType==1)
						{
							Stiff[End[0][Member]][i][j]+=tempdouble*(-EA[MemberType[Member]]*(Rotatedx[0][Member][i]*Rotatedx[0][Member][j]
																							 +Rotatedy[0][Member][i]*Rotatedy[0][Member][j])
																	 +12.0*(EIxx[MemberType[Member]]*Rotatedy[0][Member][i]*Rotatedy[0][Member][j]
																		   +EIyy[MemberType[Member]]*Rotatedx[0][Member][i]*Rotatedx[0][Member][j]));
							Stiff[End[1][Member]][i][j]+=tempdouble*(-EA[MemberType[Member]]*(Rotatedx[1][Member][i]*Rotatedx[1][Member][j]
																							 +Rotatedy[1][Member][i]*Rotatedy[1][Member][j])
																	 +12.0*(EIxx[MemberType[Member]]*Rotatedy[1][Member][i]*Rotatedy[1][Member][j]
																		   +EIyy[MemberType[Member]]*Rotatedx[1][Member][i]*Rotatedx[1][Member][j]));
						}
						else
						{
							Stiff[End[0][Member]][i][j]+=tempdouble*(EA[MemberType[Member]]*p[i]*p[j]
																	 +12.0*(EIxx[MemberType[Member]]*Rotatedy[0][Member][i]*Rotatedy[0][Member][j]
																		   +EIyy[MemberType[Member]]*Rotatedx[0][Member][i]*Rotatedx[0][Member][j]));
							Stiff[End[1][Member]][i][j]+=tempdouble*(EA[MemberType[Member]]*p[i]*p[j]
																	 +12.0*(EIxx[MemberType[Member]]*Rotatedy[1][Member][i]*Rotatedy[1][Member][j]
																		   +EIyy[MemberType[Member]]*Rotatedx[1][Member][i]*Rotatedx[1][Member][j]));
						}
					}
				}
				
				if(EIxx[MemberType[Member]]>EIyy[MemberType[Member]])
					tempdouble=4.0*EIxx[MemberType[Member]];
				else tempdouble=4.0*EIyy[MemberType[Member]];
				if(GJ[MemberType[Member]]>tempdouble)tempdouble=GJ[MemberType[Member]];
				tempdouble=tempdouble/OrigLength[Member];
				
				RotStiff[End[0][Member]]+=tempdouble;
				RotStiff[End[1][Member]]+=tempdouble;
			}
		}
	}
}

void NodeMovement(void)
{
SumForce=0.0;
for(int i=0;i<=2;i++)TotalReaction[i]=0.0;
	for(int Node=0;Node<=LastNode;Node++)
	{
		if(NodeExists[Node]==1)
		{
			if(NodeDispType[Node]!=0)
			{
				int tempint=NodeDispType[Node];
				if(tempint>=4)
				{
					Stiff[Node][0][2]=0.0;Stiff[Node][2][0]=0.0;
					Stiff[Node][1][2]=0.0;Stiff[Node][2][1]=0.0;
					TotalReaction[2]-=NodeForce[Node][2];
					NodeForce[Node][2]=0.0;tempint-=4;
				}
				if(tempint>=2)
				{
					Stiff[Node][2][1]=0.0;Stiff[Node][1][2]=0.0;
					Stiff[Node][0][1]=0.0;Stiff[Node][1][0]=0.0;
					TotalReaction[1]-=NodeForce[Node][1];
					NodeForce[Node][1]=0.0;tempint-=2;
				}
				if(tempint>=1)
				{
					Stiff[Node][1][0]=0.0;Stiff[Node][0][1]=0.0;
					Stiff[Node][2][0]=0.0;Stiff[Node][0][2]=0.0;
					TotalReaction[0]-=NodeForce[Node][0];
					NodeForce[Node][0]=0.0;
				}
			}
			
			double ForceMagSq=0.0;
			for(int i=0;i<=2;i++)ForceMagSq+=NodeForce[Node][i]*NodeForce[Node][i];
			SumForce+=sqrt(ForceMagSq);
			
			if(NodeRotType[Node]!=0)
			{
				int tempint=NodeRotType[Node];
				if(tempint>=4){NodeMom[Node][2]=0.0;tempint-=4;}
				if(tempint>=2){NodeMom[Node][1]=0.0;tempint-=2;}
				if(tempint>=1){NodeMom[Node][0]=0.0;}
			}
			
			if(NodalMatrixControl==0)
			{
				double Determinant=Stiff[Node][0][0]*Stiff[Node][1][1]*Stiff[Node][2][2]
				+2.0*Stiff[Node][0][1]*Stiff[Node][1][2]*Stiff[Node][2][0]
				-Stiff[Node][0][0]*Stiff[Node][1][2]*Stiff[Node][2][1]
				-Stiff[Node][1][1]*Stiff[Node][2][0]*Stiff[Node][0][2]
				-Stiff[Node][2][2]*Stiff[Node][0][1]*Stiff[Node][1][0];
				
				if(Determinant!=0.0)
				{
				for(int i=0;i<=2;i++)
				{
					int j,k;
					
					if(i==0){j=1;k=2;}
					if(i==1){j=2;k=0;}
					if(i==2){j=0;k=1;}
					Displacement[Node][i]=
						((Stiff[Node][1][j]*Stiff[Node][2][k]-Stiff[Node][1][k]*Stiff[Node][2][j])*NodeForce[Node][0]
						+(Stiff[Node][2][j]*Stiff[Node][0][k]-Stiff[Node][2][k]*Stiff[Node][0][j])*NodeForce[Node][1]
						+(Stiff[Node][0][j]*Stiff[Node][1][k]-Stiff[Node][0][k]*Stiff[Node][1][j])*NodeForce[Node][2])
						/(Determinant*DispMassMult)
						+DispCarryOver*Displacement[Node][i];
				}
				}
				else ZeroStiffness(Node);
			}
			else
			{
				double SumSq=0.0;
				for(int i=0;i<=2;i++)
				{
					for(int j=0;j<=2;j++)SumSq+=Stiff[Node][i][j]*Stiff[Node][i][j];
				}
				double RootSumSq=sqrt(SumSq);
				if(RootSumSq!=0.0)
				{
				for(int i=0;i<=2;i++)Displacement[Node][i]=NodeForce[Node][i]/(RootSumSq*DispMassMult)+DispCarryOver*Displacement[Node][i];
				}
				else ZeroStiffness(Node);
			}
			
			for(int i=0;i<=2;i++)Coord[Node][i]+=Displacement[Node][i];
			
			if(RotStiff[Node]!=0.0)//This allows cable structures
			{
				double tempdouble=RotInertiaMult*RotStiff[Node];
				for(int i=0;i<=2;i++)Rotation[Node][i]=NodeMom[Node][i]/tempdouble+RotCarryOver*Rotation[Node][i];
			}
			else
			{
				for(int i=0;i<=2;i++)Rotation[Node][i]=0.0;
			}
			double tempdouble=1.0;
			for(int i=0;i<=2;i++)tempdouble-=a[Node][i]*Rotation[Node][i]/2.0;
			for(int i=0;i<=2;i++)
			{
				int j,k;
				
				if(i==0){j=1;k=2;}
				if(i==1){j=2;k=0;}
				if(i==2){j=0;k=1;}
				a[Node][i]=
					(a[Node][i]+(Rotation[Node][i]+Rotation[Node][j]*a[Node][k]-Rotation[Node][k]*a[Node][j])/2.0)/tempdouble;
			}
		}
	}
}

double myFunction(int ThisMemberType,double argument)
{
	double absargument=fabs(argument);
	double doubleValue=absargument/FunctionIncrement;
	int FunctionInteger=int(doubleValue);
	double weighting=doubleValue-double(FunctionInteger);
	if(weighting<0.0||weighting>1.0)cout<<"\n\nWeighting problem!\n\n";
	double answer;
	if(FunctionInteger<=FunctionLastInteger-1)
	answer=(1.0-weighting)*FunctionValue[ThisMemberType][FunctionInteger]+weighting*FunctionValue[ThisMemberType][FunctionInteger+1];
	else answer=1.0;
	if(argument<0.0)answer=-answer;
	return answer;
}

double NonLinearElastic(int WhichMember,int WhichMemberType,double ElasticValue,double PeakValue)
{
	double Proportion;
	if(MemberTypePropsDef[WhichMemberType]!=1)cout<<"\n\nNon-elastic member type not defined!\n\n";
	if(PeakValue!=0.0)
	{
	Proportion=myFunction(WhichMemberType,ElasticValue/PeakValue);
	if(fabs(Proportion)>0.95)MemberReachedPeak[WhichMember]=1;
	}
	else Proportion=0.0;
	return Proportion*PeakValue;
}

void ZeroStiffness(int ThisNode)
{
		double CoordinatesToPlot[3];
		glColor4f(1.0,0.0,0.0,1.0);glPointSize(10.0);
		glBegin(GL_POINTS);
		for(int i=0;i<=2;i++)CoordinatesToPlot[i]=Coord[ThisNode][i]-AverageCoord[i];
		glVertex3dv(CoordinatesToPlot);glEnd();
		if(FirstCycle==1)
		{
		cout<<"\nNode "<<ThisNode<<" has zero stiffness. Stiffness matrix:\n";
		for(int i=0;i<=2;i++)
		{
		cout<<"\t";
		for(int j=0;j<=2;j++)cout<<"\t"<<Stiff[ThisNode][i][j];
		cout<<"\n";
		}
		for(int Member=0;Member<=LastMember;Member++)
		{
		if(MemberExists[Member]==1)
		{
		if(ThisNode==End[0][Member]||ThisNode==End[1][Member])cout<<"     It is connected to member "<<Member<<" of type "<<MemberType[Member];
		if(ThisNode==End[0][Member])cout<<" at end 0";
		if(ThisNode==End[1][Member])cout<<" at end 1";
		if(ThisNode==End[0][Member]||ThisNode==End[1][Member])cout<<"\n";
		}
		}
		}
}
