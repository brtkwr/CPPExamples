//Written by Chris J K Williams, University of Bath, UK

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include <stdio.h>
#include <string.h>
#include "Graphics.h"

#define   mMax 2250
#define   nMax 1275
#define   LastStreakParticle  1000000
#define   LastCylinderPoint 360
void ContourLines(void);
void Contour(double[3],double[3],double[3]);
int SetUp(void);
double PsiBelow(int,int);
double PsiAbove(int,int);
int jp1(int);
int jm1(int);

int NodeOnCylinder[mMax+1][nMax+1],StreakParticleInUse[LastStreakParticle+1];
int m,n,Lastdotcycle,StreakParticlesPerCycle,LastStreakParticleSoFar,
CylinderEdgeNodes,Cycle;
int DrawingSaveInterval=500;//Never saves if negative
double Psi[mMax+1][nMax+1],
delSqPsi[mMax+1][nMax+1],PsiDot[mMax+1][nMax+1],
BasicError[mMax+1][nMax+1],Error[mMax+1][nMax+1],PsiDotIncrement[mMax+1][nMax+1],
xStreakParticle[LastStreakParticle+1],yStreakParticle[LastStreakParticle+1],StreakParticleVorticity[LastStreakParticle+1],
PsiStreakParticle[LastStreakParticle+1],
xCylinderEdgePoint[LastCylinderPoint+1],yCylinderEdgePoint[LastCylinderPoint+1],
xTri[3],yTri[3],psiTri[3],
DeltaXorY,DeltaXorYSq,WindSpeed,DeltaT,ReynoldsNumber,MuOverRo,SumSqError,ContourInterval,
CylinderDiameter,CylinderRadius,CylinderRadiusSq,InkRadius,InkRadiusSq,xCylinderCentre,yCylinderCentre,CylinderPsi,
CylinderPsiDot,CylinderPsiDotIncrement,PsiChangeAcrossChannel,
xShift,yShift,psiScale;

int main(int argc, char* argv[])
{
	PI=4.0*atan(1.0);
	if(SetUp()!=0)return 0;
	int DataControl;
	cout<<"Type 0 to continue with new data, or any other number to read from file.\n";
	cin>>DataControl;
	if(DataControl!=0)Readdata();
	SumSqError=0.0;
	
	glutInit(&argc,argv);
	SetUpGraphics(1.0,1.0,1.0);
	
	cout<<"Finished.\n";
	return 0;
}
int SetUp(void)
{
	FileNumber=-1;
	Scale=0.5;//1.5
	
	WindSpeed=200.0;
	ReynoldsNumber=1.0*WindSpeed*Scale;
	
	/*From SPH program:
	ZoneSize=5.0;
	LastZone[0]=int(300*Scale)-1;
	LastZone[1]=int(170*Scale)-1;
	Boundary[xyz]=0.5*(LastZone[xyz]+1)*ZoneSize;
	CylinderRadius=0.075*Boundary[0];*/
	
	float BasicNumber=2.0;
	
	m=int(BasicNumber*300.0*Scale);
	n=int(BasicNumber*170.0*Scale);
	if(m>mMax||n>nMax){cout<<"\nScale too big, m = "<<m<<",n = "<<n<<"\n";return 1;}
	
	DeltaXorY=5.0/double(BasicNumber);
	
	psiScale=1.0/WindSpeed;
	
	CylinderRadius=0.075*float(m)*DeltaXorY/2.0;
	CylinderDiameter=2.0*CylinderRadius;
	CylinderRadiusSq=CylinderRadius*CylinderRadius;
	InkRadius=(1.0+0.0001/Scale)*CylinderDiameter/2.0;
	InkRadiusSq=InkRadius*InkRadius;
	StreakParticlesPerCycle=int(200.0*Scale);
	
	MuOverRo=CylinderDiameter*WindSpeed/ReynoldsNumber;
	xShift=double(m)*DeltaXorY/2.0;
	yShift=double(n)*DeltaXorY/2.0;
	
	DeltaT=0.05*DeltaXorY/WindSpeed;
	
	ofstream Julia;
	Julia.open("Information.txt");
	Julia<<"ReynoldsNumber's number = "<<ReynoldsNumber<<"\n";
	Julia.close();
	
	Lastdotcycle=20;
	
	LastStreakParticleSoFar=-1;
	
	DeltaXorYSq=DeltaXorY*DeltaXorY;
	ContourInterval=1000.0*DeltaXorY;
	xCylinderCentre=0.25*double(m)*DeltaXorY;
	yCylinderCentre=0.5*DeltaXorY*double(n);
	
	for(int CylinderPoint=0;CylinderPoint<=LastCylinderPoint;CylinderPoint++)
	{
		xCylinderEdgePoint[CylinderPoint]=-xShift+xCylinderCentre+CylinderRadius*cos(2.0*PI*double(CylinderPoint)/double(LastCylinderPoint));
		yCylinderEdgePoint[CylinderPoint]=-yShift+yCylinderCentre+CylinderRadius*sin(2.0*PI*double(CylinderPoint)/double(LastCylinderPoint));
	}
	
	PsiChangeAcrossChannel=WindSpeed*DeltaXorY*float(n);
	CylinderPsi=0.001*PsiChangeAcrossChannel;CylinderPsiDot=0.0;CylinderPsiDotIncrement=0.0;
	for(int i=0;i<=m;i++)
	{
		for(int j=0;j<=n;j++)
		{
			NodeOnCylinder[i][j]=0;
			double zReal=DeltaXorY*double(i)-xCylinderCentre;
			double zImag=DeltaXorY*(double(j)-0.5*double(n));
			Psi[i][j]=WindSpeed*zImag;
			
			double eta=CylinderRadius/sqrt(zReal*zReal+zImag*zImag);
			double Correction=CylinderPsi-WindSpeed*zImag*eta*eta;
			
			Correction*=
				(1.0+zReal/xCylinderCentre)/(1.0+eta*zReal/xCylinderCentre)*
				(1.0-zReal/(double(m)*DeltaXorY-xCylinderCentre))/(1.0-eta*zReal/(double(m)*DeltaXorY-xCylinderCentre))*
				(1.0+zImag/(0.5*double(n)*DeltaXorY))/(1.0+eta*zImag/(0.5*double(n)*DeltaXorY))*
				(1.0-zImag/(0.5*double(n)*DeltaXorY))/(1.0-eta*zImag/(0.5*double(n)*DeltaXorY));
			Psi[i][j]+=Correction;
			
			delSqPsi[i][j]=0.0;
			PsiDot[i][j]=0.0;
			PsiDotIncrement[i][j]=0.0;
		}
	}
	for(int i=0;i<=m;i++)
	{
		for(int j=0;j<=n;j++)
		{
			if(sqrt(
					(DeltaXorY*i-xCylinderCentre)*(DeltaXorY*i-xCylinderCentre)+
					(DeltaXorY*j-yCylinderCentre)*(DeltaXorY*j-yCylinderCentre))<CylinderRadius)
			{
				NodeOnCylinder[i][j]=1;Psi[i][j]=CylinderPsi;
			}
		}
	}
	CylinderEdgeNodes=0;
	for(int i=1;i<=m-1;i++)
	{
		for(int j=1;j<=n-1;j++)
		{
			if(NodeOnCylinder[i][j]==1&&
			   (NodeOnCylinder[i-1][j]==0||NodeOnCylinder[i+1][j]==0||NodeOnCylinder[i][j-1]==0||NodeOnCylinder[i][j+1]==0))
			{
				NodeOnCylinder[i][j]=2;
				CylinderEdgeNodes++;
			}
		}
	}
	for(int StreakParticle=0;StreakParticle<=LastStreakParticle;StreakParticle++)StreakParticleInUse[StreakParticle]=0;
	return 0;
}

void ContourLines(void)
{
	for(int i=0;i<=m-1;i++)
	{
		for(int j=0;j<=n-1;j++)
		{
			xTri[2]=DeltaXorY*(double(i)+0.5);yTri[2]=DeltaXorY*(double(j)+0.5);
			psiTri[2]=(Psi[i+0][j+0]+Psi[i+0][j+1]+Psi[i+1][j+0]+Psi[i+1][j+1])/4.0;
			
			xTri[0]=DeltaXorY*double(i+0);xTri[1]=DeltaXorY*double(i+1);
			yTri[0]=DeltaXorY*double(j+0);yTri[1]=DeltaXorY*double(j+0);
			psiTri[0]=Psi[i+0][j+0];psiTri[1]=Psi[i+1][j+0];
			Contour(xTri,yTri,psiTri);
			
			xTri[0]=DeltaXorY*double(i+1);xTri[1]=DeltaXorY*double(i+1);
			yTri[0]=DeltaXorY*double(j+0);yTri[1]=DeltaXorY*double(j+1);
			psiTri[0]=Psi[i+1][j+0];psiTri[1]=Psi[i+1][j+1];
			Contour(xTri,yTri,psiTri);
			
			xTri[0]=DeltaXorY*double(i+1);xTri[1]=DeltaXorY*double(i+0);
			yTri[0]=DeltaXorY*double(j+1);yTri[1]=DeltaXorY*double(j+1);
			psiTri[0]=Psi[i+1][j+1];psiTri[1]=Psi[i+0][j+1];
			Contour(xTri,yTri,psiTri);
			
			xTri[0]=DeltaXorY*double(i+0);xTri[1]=DeltaXorY*double(i+0);
			yTri[0]=DeltaXorY*double(j+1);yTri[1]=DeltaXorY*double(j+0);
			psiTri[0]=Psi[i+0][j+1];psiTri[1]=Psi[i+0][j+0];
			Contour(xTri,yTri,psiTri);
		}
	}
}

void Contour(double xTriangle[3],double yTriangle[3],double psiTriangle[3])
{
	int CornerPlus1,CornerPlus2,intpsiStart,intpsiStop;
	double psiStart,psiStop,thisPsi;
	for(int Corner=0;Corner<=2;Corner++)
	{
		CornerPlus1=Corner+1;if(Corner+1>2)CornerPlus1-=3;
		CornerPlus2=Corner+2;if(Corner+2>2)CornerPlus2-=3;
		if((psiTriangle[CornerPlus1]-psiTriangle[Corner])*(psiTriangle[CornerPlus2]-psiTriangle[Corner])>
		   1.0e-6*ContourInterval*ContourInterval);
		{
			if(psiTriangle[CornerPlus1]-psiTriangle[Corner]>0.0)
			{
				psiStart=psiTriangle[Corner];
				if(psiTriangle[CornerPlus1]>psiTriangle[CornerPlus2])psiStop=psiTriangle[CornerPlus2];
				else psiStop=psiTriangle[CornerPlus1];
			}
			else
			{
				psiStop=psiTriangle[Corner];
				if(psiTriangle[CornerPlus1]<psiTriangle[CornerPlus2])psiStart=psiTriangle[CornerPlus2];
				else psiStart=psiTriangle[CornerPlus1];
			}
		}
		intpsiStart=int(psiStart/ContourInterval);if(psiStart<0.0)intpsiStart-=1;
		intpsiStop=int(psiStop/ContourInterval);  if(psiStop<0.0)intpsiStop-=1;
		if(intpsiStart<intpsiStop)
		{
			for(int intpsi=intpsiStart+1;intpsi<=intpsiStop;intpsi++)
			{
				double rv[3];
				glBegin(GL_LINE_STRIP);
				thisPsi=intpsi*ContourInterval;
				rv[0]=-xShift+xTriangle[Corner]+(xTriangle[CornerPlus1]-xTriangle[Corner])*(thisPsi-psiTriangle[Corner])/(psiTriangle[CornerPlus1]-psiTriangle[Corner]);
				rv[1]=-yShift+yTriangle[Corner]+(yTriangle[CornerPlus1]-yTriangle[Corner])*(thisPsi-psiTriangle[Corner])/(psiTriangle[CornerPlus1]-psiTriangle[Corner]);
				rv[2]=psiScale*thisPsi;
				glVertex3dv(rv);
				rv[0]=-xShift+xTriangle[Corner]+(xTriangle[CornerPlus2]-xTriangle[Corner])*(thisPsi-psiTriangle[Corner])/(psiTriangle[CornerPlus2]-psiTriangle[Corner]);
				rv[1]=-yShift+yTriangle[Corner]+(yTriangle[CornerPlus2]-yTriangle[Corner])*(thisPsi-psiTriangle[Corner])/(psiTriangle[CornerPlus2]-psiTriangle[Corner]);
				glVertex3dv(rv);glEnd();
			}
		}
	}
}

static void Draw(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	int NumberOfStreakParticlesInUseMinus1=-1;
	for(int StreakParticle=0;StreakParticle<=LastStreakParticleSoFar;StreakParticle++)
	{
		if(StreakParticleInUse[StreakParticle]==1)
		{
			if(
			   (xStreakParticle[StreakParticle]-xCylinderCentre)*
			   (xStreakParticle[StreakParticle]-xCylinderCentre)+
			   (yStreakParticle[StreakParticle]-yCylinderCentre)*
			   (yStreakParticle[StreakParticle]-yCylinderCentre)
			   <InkRadiusSq)StreakParticleInUse[StreakParticle]=0;
			else NumberOfStreakParticlesInUseMinus1++;
		}
	}
	
	for(int NewStreakParticle=1;NewStreakParticle<=StreakParticlesPerCycle;NewStreakParticle++)
	{
		int StreakParticle=-1;
		for(;;)
		{
			StreakParticle++;
			if(StreakParticle>LastStreakParticleSoFar)LastStreakParticleSoFar=StreakParticle;
			if(LastStreakParticleSoFar>LastStreakParticle)LastStreakParticleSoFar=LastStreakParticle;
			if(StreakParticle>LastStreakParticle)break;
			int FoundOne=0;
			if(StreakParticleInUse[StreakParticle]==0)
			{
				FoundOne=1;
				StreakParticleInUse[StreakParticle]=1;
				double Radius=CylinderRadius+(InkRadius-CylinderRadius)*double(rand())/double(RAND_MAX);
				double angle=2.0*PI*double(rand())/double(RAND_MAX);
				xStreakParticle[StreakParticle]=xCylinderCentre+Radius*cos(angle);
				yStreakParticle[StreakParticle]=yCylinderCentre+Radius*sin(angle);
				StreakParticleVorticity[StreakParticle]=0.0;
				PsiStreakParticle[StreakParticle]=0.0;
			}
			if(FoundOne==1)break;
		}
	}
	
	double rv[3];
	
	glLineWidth(3.0);
	glColor4d(1.0,0.0,0.0,0.5);
	glBegin(GL_LINE_STRIP);
	rv[0]=-0.53*DeltaXorY*double(m);rv[1]=-yShift;rv[2]=-psiScale*WindSpeed*DeltaXorY*double(n)/2.0;glVertex3dv(rv);
	rv[1]=-yShift+DeltaXorY*double(n)*double(NumberOfStreakParticlesInUseMinus1)/double(LastStreakParticle);
	rv[2]=psiScale*WindSpeed*DeltaXorY*double(n)*(double(NumberOfStreakParticlesInUseMinus1)/double(LastStreakParticle)-0.5);
	glVertex3dv(rv);
	glEnd();
	glPointSize(4.0);
	glBegin(GL_POINTS);
	glColor4d(0.0,0.0,0.0,1.0);
	glVertex3dv(rv);
	glEnd();
	glColor4d(0.0,0.0,1.0,0.5);
	glBegin(GL_LINE_STRIP);
	glVertex3dv(rv);
	rv[1]=+yShift;rv[2]=psiScale*WindSpeed*DeltaXorY*double(n)/2.0;glVertex3dv(rv);
	glEnd();
	
	glLineWidth(2.0);
	glColor4d(0.0,0.0,0.0,1.0);
	glBegin(GL_LINE_STRIP);
	rv[0]=-xShift;rv[1]=-yShift;rv[2]=psiScale*Psi[0][0];glVertex3dv(rv);
	rv[0]=DeltaXorY*double(m)-xShift;rv[2]=psiScale*Psi[m][0];glVertex3dv(rv);
	rv[1]=DeltaXorY*double(n)-yShift;rv[2]=psiScale*Psi[m][n];glVertex3dv(rv);
	rv[0]=-xShift;rv[2]=psiScale*Psi[0][n];glVertex3dv(rv);
	rv[1]=-yShift;rv[2]=psiScale*Psi[0][0];glVertex3dv(rv);
	glEnd();
	
	if(ShowGrid==1)
	{
		glLineWidth(1.0);
		glColor4d(0.0,0.0,1.0,0.5);
		for(int j=1;j<=n-1;j++)
		{
		glBegin(GL_LINE_STRIP);
		for(int i=0;i<=m;i++)
		{
		rv[0]=-xShift+DeltaXorY*double(i);rv[1]=-yShift+double(j)*DeltaXorY;rv[2]=psiScale*Psi[i][j];glVertex3dv(rv);
		}
		glEnd();
		}
		for(int i=1;i<=m-1;i++)
		{
			glBegin(GL_LINE_STRIP);
			for(int j=0;j<=n;j++)
			{
			rv[0]=-xShift+double(i)*DeltaXorY;rv[1]=-yShift+double(j)*DeltaXorY;rv[2]=psiScale*Psi[i][j];glVertex3dv(rv);
			}
						glEnd();
		}
	}
	
	glLineWidth(2.0);
	glColor4d(0.0,0.0,0.0,1.0);
	glBegin(GL_LINE_STRIP);
	for(int CylinderPoint=0;CylinderPoint<=LastCylinderPoint;CylinderPoint++)
	{
		rv[0]=xCylinderEdgePoint[CylinderPoint];
		rv[1]=yCylinderEdgePoint[CylinderPoint];
		rv[2]=psiScale*CylinderPsi;
		glVertex3dv(rv);
	}
	glEnd();
	
	if(ShowContours==1)
	{
		glLineWidth(1.0);
		glColor4d(0.0,0.0,0.0,0.5);
		ContourLines();
	}
	
	if(ShowStreak==1)
	{
		double ThisPointSize=0.4/Scale;
		if(ThisPointSize<1.0)ThisPointSize=1.0;
		glPointSize(ThisPointSize);
		glBegin(GL_POINTS);
		for(int StreakParticle=0;StreakParticle<=LastStreakParticleSoFar;StreakParticle++)
		{
			if(StreakParticleInUse[StreakParticle]==1)
			{
				double ColourControl=10.0*StreakParticleVorticity[StreakParticle]/(WindSpeed/DeltaXorY);
				double ColourShift=-1.0;
				glColor4d((1.0+tanh(ColourShift+ColourControl))/2.0,0.0,(1.0+tanh(ColourShift-ColourControl))/2.0,1.0);
				rv[0]=xStreakParticle[StreakParticle]-xShift;
				rv[1]=yStreakParticle[StreakParticle]-yShift;
				rv[2]=psiScale*PsiStreakParticle[StreakParticle];
				glVertex3dv(rv);
			}
		}
		glEnd();
	}
	
	char NumericalValue[101];
	char *MyText;
	glColor4f(0.0,0.0,0.0,0.5);
	glRasterPos2f(-550.0*Scale,-450.0*Scale);
	MyText="ReynoldsNumber number  ";
	int LengthOfString=strlen(MyText);
	for(int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
	sprintf(NumericalValue,"%i",int(ReynoldsNumber));
	LengthOfString=strlen(NumericalValue);
	for(int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	MyText="         Cycle  ";
	LengthOfString=strlen(MyText);
	for(int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
	sprintf(NumericalValue,"%i",Cycle);
	LengthOfString=strlen(NumericalValue);
	for(int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	MyText="         Non-dimensional time ";
	LengthOfString=strlen(MyText);
	for(int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,MyText[MyCharacter]);
	sprintf(NumericalValue,"%8.2f",DeltaT*Cycle*WindSpeed/CylinderDiameter);
	LengthOfString=strlen(NumericalValue);
	for(int MyCharacter=0;MyCharacter<LengthOfString;MyCharacter++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,NumericalValue[MyCharacter]);
	glutSwapBuffers();
	
	for(int i=1;i<=m-1;i++)
	{
		for(int j=0;j<=n;j++)
		{
			delSqPsi[i][j]=(Psi[i-1][j]+Psi[i+1][j]+PsiBelow(i,j)+PsiAbove(i,j)-4.0*Psi[i][j])/DeltaXorYSq;
		}
	}
	for(int i=2;i<=m-2;i++)
	{
		for(int j=0;j<=n;j++)
		{
			BasicError[i][j]=MuOverRo*(delSqPsi[i-1][j]+delSqPsi[i+1][j]+delSqPsi[i][jm1(j)]+delSqPsi[i][jp1(j)]-4.0*delSqPsi[i][j]);
			BasicError[i][j]-=(PsiAbove(i,j)-PsiBelow(i,j))*(delSqPsi[i+1][j]   -delSqPsi[i-1][j]   )/4.0;
			BasicError[i][j]+=(Psi[i+1][j]  -Psi[i-1][j]  )*(delSqPsi[i][jp1(j)]-delSqPsi[i][jm1(j)])/4.0;
		}
	}
	for(int dotcycle=0;dotcycle<=Lastdotcycle;dotcycle++)
	{
		double CylinderError=0.0;
		for(int i=2;i<=m-2;i++)
		{
			for(int j=0;j<=n;j++)
			{
				Error[i][j]=BasicError[i][j]-(PsiDot[i-1][j]+PsiDot[i+1][j]+PsiDot[i][jm1(j)]+PsiDot[i][jp1(j)]-4.0*PsiDot[i][j]);
				if(NodeOnCylinder[i][j]!=0)CylinderError+=Error[i][j];
			}
		}
		CylinderPsiDotIncrement=0.9*CylinderPsiDotIncrement-0.5*CylinderError/(2.0*double(CylinderEdgeNodes));
		CylinderPsiDot+=CylinderPsiDotIncrement;
		for(int i=2;i<=m-2;i++)
		{
			for(int j=0;j<=n;j++)
			{
				if(NodeOnCylinder[i][j]==0)
				{
					PsiDotIncrement[i][j]=0.9*PsiDotIncrement[i][j]-0.5*Error[i][j]/4.0;
					PsiDot[i][j]+=PsiDotIncrement[i][j];
				}
				else
				{
					PsiDotIncrement[i][j]=CylinderPsiDotIncrement;
					PsiDot[i][j]=CylinderPsiDot;
				}
			}
		}
	}
	CylinderPsi+=DeltaT*(CylinderPsiDot+CylinderPsiDotIncrement/2.0);
	for(int i=2;i<=m-2;i++)
	{
		for(int j=0;j<=n;j++)//Periodic boundaries top and bottom
		//for(int j=2;j<=n-2;j++)//Fixed position and slope top and bottom
		{
			if(NodeOnCylinder[i][j]==0)Psi[i][j]+=DeltaT*(PsiDot[i][j]+PsiDotIncrement[i][j]/2.0);
			else Psi[i][j]=CylinderPsi;
		}
		Psi[i][0]=(Psi[i][0]+Psi[i][n])/2.0-PsiChangeAcrossChannel/2.0;
		Psi[i][n]=Psi[i][0]+PsiChangeAcrossChannel;
	}
	//for(int j=0;j<=n;j++)Psi[m-1][j]=2.0*Psi[m-2][j]-Psi[m-3][j];
	//for(int j=0;j<=n;j++)Psi[m-0][j]=2.0*Psi[m-1][j]-Psi[m-2][j];
	for(int j=0;j<=n;j++)Psi[m-1][j]=(Psi[m-2][j]+Psi[m][j])/2.0;
	
	for(int StreakParticle=0;StreakParticle<=LastStreakParticleSoFar;StreakParticle++)
	{
		if(StreakParticleInUse[StreakParticle]==1)
		{
			double reali=xStreakParticle[StreakParticle]/DeltaXorY;int i=int(reali);if(i<1)i=1;if(i>m-2)i=m-2;double myu=reali-i;
			double realj=yStreakParticle[StreakParticle]/DeltaXorY;int j=int(realj);if(j<1)j=1;if(j>n-2)j=n-2;double myv=realj-j;
			xStreakParticle[StreakParticle]+=(
											  (Psi[i+0][j+1]-Psi[i+0][j-1])*(1.0-myu)*(1.0-myv)+
											  (Psi[i+1][j+1]-Psi[i+1][j-1])*myu*(1.0-myv)+
											  (Psi[i+0][j+2]-Psi[i+0][j-0])*(1.0-myu)*myv+
											  (Psi[i+1][j+2]-Psi[i+1][j-0])*myu*myv)*DeltaT/(2.0*DeltaXorY);
			yStreakParticle[StreakParticle]-=(
											  (Psi[i+1][j+0]-Psi[i-1][j+0])*(1.0-myu)*(1.0-myv)+
											  (Psi[i+2][j+0]-Psi[i+0][j+0])*myu*(1.0-myv)+
											  (Psi[i+1][j+1]-Psi[i-1][j+1])*(1.0-myu)*myv+
											  (Psi[i+2][j+1]-Psi[i+0][j+1])*myu*myv)*DeltaT/(2.0*DeltaXorY);
			StreakParticleVorticity[StreakParticle]=(
													 delSqPsi[i+0][j+0]*(1.0-myu)*(1.0-myv)+
													 delSqPsi[i+1][j+0]*myu*(1.0-myv)+
													 delSqPsi[i+0][j+1]*(1.0-myu)*myv+
													 delSqPsi[i+1][j+1]*myu*myv)/2.0;
			PsiStreakParticle[StreakParticle]=(
													 Psi[i+0][j+0]*(1.0-myu)*(1.0-myv)+
													 Psi[i+1][j+0]*myu*(1.0-myv)+
													 Psi[i+0][j+1]*(1.0-myu)*myv+
													 Psi[i+1][j+1]*myu*myv)/2.0;
			if(
			   xStreakParticle[StreakParticle]<0.0||xStreakParticle[StreakParticle]>DeltaXorY*double(m)||
			   yStreakParticle[StreakParticle]<0.0||yStreakParticle[StreakParticle]>DeltaXorY*double(n))StreakParticleInUse[StreakParticle]=0;
		}
	}
	if(DrawingSaveInterval>=0&&Cycle==DrawingSaveInterval*int(Cycle/DrawingSaveInterval))SavePicture();
	Cycle++;
}

double PsiBelow(int i,int j)
{
	if(j>0)
		return Psi[i][j-1];
	else
		return Psi[i][n-1]-PsiChangeAcrossChannel;
}

double PsiAbove(int i,int j)
{
	if(j<n)
		return Psi[i][j+1];
	else 
		return Psi[i][0+1]+PsiChangeAcrossChannel;
}

int jp1(int j)
{
	if(j<n)return j+1;else return 1;
}

int jm1(int j)
{
	if(j>0)return j-1;else return n-1;
}

void Writedata(void)
{
	ofstream Julia;
	Julia.open("Data.txt");
	for(int i=0;i<=m;i++)
	{
		for(int j=0;j<=n;j++)
			Julia<<NodeOnCylinder[i][j]<<" "<<Psi[i][j]<<" "<<delSqPsi[i][j]<<" "<<PsiDot[i][j]<<" "<<PsiDotIncrement[i][j]<<"\n";
	}
	for(int StreakParticle=0;StreakParticle<=LastStreakParticle;StreakParticle++)
	{
		Julia<<StreakParticleInUse[StreakParticle]<<"\n";
		if(StreakParticleInUse[StreakParticle]==1)Julia<<xStreakParticle[StreakParticle]<<" "<<yStreakParticle[StreakParticle]<<"\n";
	}
	Julia<<Cycle<<" "<<LastStreakParticleSoFar<<" "<<FileNumber<<"\n";
	Julia.close();
	Julia.open("StreamFunction.dxf");
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	Julia<<"0\nCIRCLE\n8\nCylinder\n";
	Julia<<"10\n"<<xCylinderCentre<<"\n";
	Julia<<"20\n"<<yCylinderCentre<<"\n";
	Julia<<"40\n"<<CylinderRadius<<"\n";
	for(int StreakParticle=0;StreakParticle<=LastStreakParticle;StreakParticle++)
	{
		if(StreakParticleInUse[StreakParticle]==1)
		{
			Julia<<"0\nCIRCLE\n8\nParticles\n";
			Julia<<"10\n"<<xStreakParticle[StreakParticle]<<"\n";
			Julia<<"20\n"<<yStreakParticle[StreakParticle]<<"\n";
			Julia<<"40\n"<<CylinderRadius/200.0<<"\n";
		}
	}
	Julia<<"0\nENDSEC\n0\nEOF\n\n";
	Julia.close();
}

void Readdata(void)
{
	ifstream Julia;
	Julia.open("Data.txt");
	for(int i=0;i<=m;i++)
	{
		for(int j=0;j<=n;j++)
			Julia>>NodeOnCylinder[i][j]>>Psi[i][j]>>delSqPsi[i][j]>>PsiDot[i][j]>>PsiDotIncrement[i][j];
	}
	for(int StreakParticle=0;StreakParticle<=LastStreakParticle;StreakParticle++)
	{
		Julia>>StreakParticleInUse[StreakParticle];
		if(StreakParticleInUse[StreakParticle]==1)Julia>>xStreakParticle[StreakParticle]>>yStreakParticle[StreakParticle];
	}
	Julia>>Cycle>>LastStreakParticleSoFar>>FileNumber;
	Julia.close();
}

