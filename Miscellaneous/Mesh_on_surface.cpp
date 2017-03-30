#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define m_Maximum 1201
#define NumEllipses 6

double x[m_Maximum + 1][m_Maximum + 1],y[m_Maximum + 1][m_Maximum + 1],z[m_Maximum + 1][m_Maximum + 1],
PI,L,Lsq,x1,y_1,z1,x2,y2,z2,x3,y3,z3,xvalue,yvalue,zvalue,zx,zy,
xsqui,ysqui,xmult,ymult,QuantitySq,Quantity,
dL1sqbydX,dL2sqbydX,dL1sqbydY,dL2sqbydY,determinant,error1,error2,
maxerrorLsq,stressratio,rootstressratio,
a[NumEllipses],b[NumEllipses],xcentre[NumEllipses],ycentre[NumEllipses],
angle23,angle34,temptheta1,temptheta2,xTemp,yTemp,existanceTemp;

int existance[m_Maximum + 1][m_Maximum + 1],m,repeat,step,vertmax,
numDiv,iinc,jinc;

void Function(void);
void NewPoint(int,int);
void PlanShape(void);
void HereOrNot(void);
void MyLine(int,int);

ofstream Julia("Mesh.dxf");

int main(void)
{
	PI=4.0*atan(1.0);
	stressratio=2.0;rootstressratio=sqrt(stressratio);
	
	step=10;
	L=1.0/(1.0*step);Lsq=L*L;
	m=120.0/L;if(m>m_Maximum){cout<<"Arrays too big\n";return 0;}
	vertmax=m/3;
	
	for(int i=0;i<=m;i++)
	{
		int j=m-i;
		xvalue=L*(i-j)/sqrt(1.0+stressratio);
		yvalue=0.0;
		Function();
        x[i][j]=xvalue;
        y[i][j]=yvalue;
        z[i][j]=zvalue;
		existance[i][j]=1;
	}
	
	maxerrorLsq=0.0;
	for(int toporbot=-1;toporbot<=1;toporbot+=2)
	{
		for(int vert=1;vert<=vertmax;vert++)
		{
			for(int hori=0;hori<=m-vert;hori++)
			{
				int i,j;
				if(toporbot==1)
				{
					i=hori+vert;
					j=m+vert-i;
				}
				else
				{
					i=m-hori-vert;
					j=m-vert-i;
				}
				
				if(i<0||i>m||j<0||j>m){cout<<"Out of range\n";return 0;}
				
				x1=x[i-toporbot][j];
				y_1=y[i-toporbot][j];
				z1=z[i-toporbot][j];
				x2=x[i][j-toporbot];
				y2=y[i][j-toporbot];
				z2=z[i][j-toporbot];
				x3=x[i-toporbot][j-toporbot];
				y3=y[i-toporbot][j-toporbot];
				z3=z[i-toporbot][j-toporbot];
				NewPoint(vert,toporbot);
				x[i][j]=xvalue;
				y[i][j]=yvalue;
				z[i][j]=zvalue;
			}
		}
	}
	
	cout<<"Maximum error (before final repeat) = "<<maxerrorLsq<<"\n";
	
	Julia<<"0\nSECTION\n2\nENTITIES\n";
	
	PlanShape();
	
	for(int i=0;i<=m;i+=step)
	{
		for(int jline=0;jline*step<=m-step;jline+=1)
		{
			if(i+jline*step+step<=m+vertmax&&i+jline*step>=m-vertmax)
			{
				int j=jline*step;iinc=0;jinc=step;
				if(existance[i][j]==existance[i+iinc][j+jinc])MyLine(i,j);
				else
				{
					for(int jbit=0;jbit<=step-1;jbit++)
					{
						j=jline*step+jbit;jinc=1;MyLine(i,j);
					}
				}
			}
		}
	}
	
	for(int j=0;j<=m;j+=step)
	{
		for(int iline=0;iline*step<=m-step;iline+=1)
		{
			if(iline*step+j+step<=m+vertmax&&iline*step+j>=m-vertmax)
			{
				int i=iline*step;iinc=step;jinc=0;
				if(existance[i][j]==existance[i+iinc][j+jinc])MyLine(i,j);
				else
				{
					for(int ibit=0;ibit<=step-1;ibit++)
					{
						int i=iline*step+ibit;iinc=1;MyLine(i,j);
					}
				}
			}
		}
	}
	
	Julia<<"0\nENDSEC\n0\nEOF\n";Julia.close();
	
	cout<<"Finished\n";
	return 0;
}

void MyLine(int i,int j)
{
	if(existance[i][j]*existance[i+iinc][j+jinc]!=0)
		Julia<<"0\nLINE\n8\nMesh\n";
	else
		Julia<<"0\nLINE\n8\nExess Mesh\n";
	Julia<<"10\n"<<x[i][j]<<"\n";
	Julia<<"20\n"<<y[i][j]<<"\n";
	Julia<<"30\n"<<z[i][j]<<"\n";
	Julia<<"11\n"<<x[i+iinc][j+jinc]<<"\n";
	Julia<<"21\n"<<y[i+iinc][j+jinc]<<"\n";
	Julia<<"31\n"<<z[i+iinc][j+jinc]<<"\n";
	if(existance[i][j]*existance[i+iinc][j+jinc]!=0)
		Julia<<"62\n0\n";
	else
		Julia<<"62\n1\n";
}

void NewPoint(int vert, int toporbot)
{
	xvalue=(x1+x2)/2.0;yvalue=(y_1+y2)/2.0;Function();
	xsqui=(x1-x2)/2.0;ysqui=(y_1-y2)/2.0;
	xmult=xsqui*(1.0+zx*zx)+ysqui*zx*zy;
	ymult=ysqui*(1.0+zy*zy)+xsqui*zx*zy;
	QuantitySq=(Lsq-(xsqui*xsqui*(1.0+zx*zx)+ysqui*ysqui*(1.0+zy*zy)+2.0*xsqui*ysqui*zx*zy))
		/(ymult*ymult*(1.0+zx*zx)+xmult*xmult*(1.0+zy*zy)-2.0*ymult*xmult*zx*zy);
	if(QuantitySq<0.0)cout<<"Problem point\n";
	else
	{
		Quantity=sqrt(QuantitySq);
		if(vert!=1){if(ymult*((x1+x2)/2.0-x3)-xmult*((y_1+y2)/2.0-y3)<0.0)Quantity=-Quantity;}
		else {if(toporbot*Quantity*xmult>0.0)Quantity=-Quantity;}
		xvalue=(x1+x2)/2.0+Quantity*ymult;
		yvalue=(y_1+y2)/2.0-Quantity*xmult;
		Function();
		
		dL1sqbydX=2.0*(xvalue-x1)+2.0*zx*((xvalue-x1)*zx+(yvalue-y_1)*zy);
		dL1sqbydY=2.0*(yvalue-y_1)+2.0*zy*((xvalue-x1)*zx+(yvalue-y_1)*zy);
		dL2sqbydX=2.0*(xvalue-x2)+2.0*zx*((xvalue-x2)*zx+(yvalue-y2)*zy);
		dL2sqbydY=2.0*(yvalue-y2)+2.0*zy*((xvalue-x2)*zx+(yvalue-y2)*zy);
		determinant=dL1sqbydX*dL2sqbydY-dL1sqbydY*dL2sqbydX;
		
		for(repeat=1;repeat<=5;repeat++)
		{
			error1=(xvalue-x1)*(xvalue-x1)+(yvalue-y_1)*(yvalue-y_1)+(zvalue-z1)*(zvalue-z1)-Lsq;
			error2=(xvalue-x2)*(xvalue-x2)+(yvalue-y2)*(yvalue-y2)+(zvalue-z2)*(zvalue-z2)-Lsq;
			
			xvalue-=(+error1*dL2sqbydY-error2*dL1sqbydY)/determinant;
			yvalue-=(-error1*dL2sqbydX+error2*dL1sqbydX)/determinant;
			Function();
		}
		if(maxerrorLsq<error1*error1)maxerrorLsq=error1*error1;
		if(maxerrorLsq<error2*error2)maxerrorLsq=error2*error2;
	}
}

void Function(void)
{
	double constant1,constant2,constant3,constant4,yfunct;
	
	constant1=1.0/1000.0;
	constant2=1.0/50.0;
	constant3=3.5;
	constant4=80.0;
	yfunct=yvalue-5.0;
	zvalue=-xvalue*xvalue*constant1-yfunct*yfunct*constant2+constant3*cos(2.0*PI*xvalue/constant4)*cosh(2.0*PI*yfunct/(rootstressratio*constant4));
    zx=   -2.0*xvalue*constant1     -constant3*(2.0*PI/constant4)*sin(2.0*PI*xvalue/constant4)*cosh(2.0*PI*yfunct/(rootstressratio*constant4));
    zy=   -2.0*yfunct*constant2
		+constant3*(2.0*PI/(rootstressratio*constant4))*cos(2.0*PI*xvalue/constant4)*sinh(2.0*PI*yfunct/(rootstressratio*constant4));
	//zvalue=0.0;zx=0.0;zy=0.0;
}

void PlanShape(void)
{
	numDiv=500;
	
	a[0]=35.0;b[0]=a[0]*rootstressratio;
	a[1]=a[0];b[1]=a[1]*rootstressratio;
	a[2]=a[0];b[2]=a[2]*rootstressratio;
	a[3]=3.75;b[3]=a[3]*rootstressratio;
	a[4]=a[0];b[4]=a[4]*rootstressratio;
	a[5]=a[0];b[5]=a[5]*rootstressratio;
	xcentre[0]=0.0;  ycentre[0]=68.5;
	xcentre[1]=35.0; ycentre[1]=63.5;
	xcentre[2]=71.0; ycentre[2]=41.5;
	angle23=57.0*PI/180.0;    xcentre[3]=xcentre[2]-(+a[2]+a[3])*cos(angle23);ycentre[3]=ycentre[2]-(+b[2]+b[3])*sin(angle23);
	angle34=109.0*PI/180.0;   xcentre[4]=xcentre[3]+(-a[3]+a[4])*cos(angle34);ycentre[4]=ycentre[3]+(-b[3]+b[4])*sin(angle34);
	xcentre[5]=xcentre[4]/3.0;ycentre[5]=ycentre[4];
	
	for(int ellipse=0;ellipse<=NumEllipses-1;ellipse++)
	{
		for(int Div=0;Div<=numDiv;Div++)
		{
			temptheta1=(2.0*PI*Div)/(1.0*(numDiv+1));
			temptheta2=(2.0*PI*(Div+1))/(1.0*(numDiv+1));
			Julia<<"0\nLINE\n8\nEllipses\n";
			Julia<<"10\n"<<xcentre[ellipse]+a[ellipse]*cos(temptheta1)<<"\n";
			Julia<<"20\n"<<ycentre[ellipse]+b[ellipse]*sin(temptheta1)<<"\n";
			Julia<<"11\n"<<xcentre[ellipse]+a[ellipse]*cos(temptheta2)<<"\n";
			Julia<<"21\n"<<ycentre[ellipse]+b[ellipse]*sin(temptheta2)<<"\n";
			Julia<<"62\n0\n";
			
			if(xcentre[ellipse]!=0.0)
			{
				Julia<<"0\nLINE\n8\nEllipses\n";
				Julia<<"10\n"<<-xcentre[ellipse]+a[ellipse]*cos(temptheta1)<<"\n";
				Julia<<"20\n"<<+ycentre[ellipse]+b[ellipse]*sin(temptheta1)<<"\n";
				Julia<<"11\n"<<-xcentre[ellipse]+a[ellipse]*cos(temptheta2)<<"\n";
				Julia<<"21\n"<<+ycentre[ellipse]+b[ellipse]*sin(temptheta2)<<"\n";
				Julia<<"62\n0\n";
			}
		}
	}
	
	for(int i=0;i<=m;i++)
	{
		for(int j=0;j<=m;j++)
		{
			if(i+j<=m+vertmax&&i+j>=m-vertmax)
			{
				xTemp=x[i][j];
				yTemp=y[i][j];
				HereOrNot();
				existance[i][j]=existanceTemp;
			}
		}
	}
}

void HereOrNot(void)
{
	existanceTemp=0;
	for(int mirror=-1;mirror<=1;mirror+=2)
	{
		for(int ellipse=4;ellipse<=NumEllipses-1;ellipse++)
			if((xTemp-mirror*xcentre[ellipse])*(xTemp-mirror*xcentre[ellipse])/(a[ellipse]*a[ellipse])
			   +(yTemp-ycentre[ellipse])*(yTemp-ycentre[ellipse])/(b[ellipse]*b[ellipse])<1.0)existanceTemp=1;
		
		if(((xTemp-mirror*xcentre[3])/(a[3]*cos(angle23))-(yTemp-ycentre[3])/(b[3]*sin(angle23)))*mirror>0.0
		   &&((xTemp-mirror*xcentre[3])/(a[3]*cos(angle34))-(yTemp-ycentre[3])/(b[3]*sin(angle34)))*mirror<0.0)existanceTemp=0;
	}
	
	int ellipse=3;
	
	for(int mirror=-1;mirror<=1;mirror+=2)
		if((xTemp-mirror*xcentre[ellipse])*(xTemp-mirror*xcentre[ellipse])/(a[ellipse]*a[ellipse])
		   +(yTemp-ycentre[ellipse])*(yTemp-ycentre[ellipse])/(b[ellipse]*b[ellipse])<1.0)existanceTemp=1;
	
	for(int mirror=-1;mirror<=1;mirror+=2)
	{
		for(int ellipse=0;ellipse<=2;ellipse++)
			if((xTemp-mirror*xcentre[ellipse])*(xTemp-mirror*xcentre[ellipse])/(a[ellipse]*a[ellipse])
			   +(yTemp-ycentre[ellipse])*(yTemp-ycentre[ellipse])/(b[ellipse]*b[ellipse])<1.0)existanceTemp=0;
	}
}