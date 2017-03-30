//Written by Dr Chris J K Williams, University of Bath, UK
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

void init_graphics(void);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
void Museum(void);
void MouseClick(int, int, int, int);

int argc;
char **argv;
bool LeftMouseDown=0,RightMouseDown=0;
float rv[2][2];	
GLUquadricObj *quadObj;

int DrawingDone,LastDrawingDone;

#define   surfacenum     40000
#define   linenum        20000
#define   maxnumrowsp1  801
#define   maxnocols     400
#define   delta         0.01
#define   centreheight 20.955
#define   edgeheight   19.71
#define   normalcalc    1



void setup(void);
void setupxy(void);
void dothebound(void);
void setbound(void);
void typicalnodexy(void);
void circumference(void);
void polar(void);
void calcxy(void);
void calczanalytic(void);
void findtheheight(void);
void findnormals(void);
void DrawTheGrid(void);
void drawagridline(void);

//WindowPtr   myWindow;
//Rect        theRect;
//PicHandle   myPicture;
///const       longZero=0;
//long        longCount,longWhere;
//short       globalRef;
//int         iforMac,picycounter,side,line,tempx,tempy,localsavepicts;
//SFReply     *reply;
//Point       wher;
//RGBColor    FillColour,FrameColour,LineColour;
//PolyHandle  Poly;

  /*  WindowRef 		fullScreenWindow;
    Ptr			oldState;
    CGrafPtr    drawingPort;
    Rect 	bounds;
    RGBColor	randomColor,FillColour;
    PolyHandle  Poly;
    
    short h,v,halfwidth,halfheight;*/


double PI,
       tempfloat,tempfloatread1,tempfloatread2,tempfloatread3,tempfloatread4,tempfloatread5,
       xline[linenum],yline[linenum],zline[linenum],
       x[surfacenum],y[surfacenum],z[surfacenum],
       lastx[surfacenum],lasty[surfacenum],lastz[surfacenum],
       xtoplot[surfacenum],ytoplot[surfacenum],
       dzbydx[surfacenum],dzbydy[surfacenum],
       xmovement[surfacenum],ymovement[surfacenum],zmovement[surfacenum],
       prevxmovement[surfacenum],prevymovement[surfacenum],prevzmovement[surfacenum],
       xcamera,ycamera,zcamera,//scale,
	   overallscale,
	   crot,srot,cameradistance,
       tx,ty,tz,
       a,perpdist[4],r,theta,
       ratio,
       xdistance,ydistance,grad,ConeAmount,
       xvalue,yvalue,zvalue,storezvalue[3][3],
       xnormal[surfacenum],ynormal[surfacenum],znormal[surfacenum],
       lambda,mu,
       previousKE,myKE,overrelax=0.0,carryover,radialx,radialy,radialz,
       Waagner,xshift,yshift;//,xshrink;
      
int   step,row,col,numrows,numcols,
      topcol[4],botcol[4],toprow[4],botrow[4],minbotrow,maxtoprow,
      boundaryrow[maxnumrowsp1],
      minnumrows,
      rowstodraw,colstodraw,colstodrawperside[4],noddy,Q,
      k,
      whichcol,tempint,
      row2,col2,row3,col3,
      phiorw,temprow,trow,tcol,
      cycle,maxnumxycycles,
      wherex,wherey,
      tempintread1,tempintread2,tempintread3,
      intxdraw[4],intydraw[4],GridsSoFar,
      end1,end2,diagonal,run,
      intermediate,intermediatecycle;

   // int halfwidth =  510;
   // int halfheight = 400;

short int linetype[linenum];

int main(int argc, char* argv[])
{

    PI=4.0*atan(1.0);
	glutInit(&argc,argv);
    init_graphics();
    DrawingDone=0;LastDrawingDone=50;
  
    glClear(GL_COLOR_BUFFER_BIT);
    glutKeyboardFunc(Key);
    glutIdleFunc(Draw);
    glutDisplayFunc(Draw);
    glutMainLoop();
      
    return 0;

}



void setup(void)
{
botcol[3]=2*(topcol[0]+topcol[1]+topcol[2]+topcol[3]-botcol[0]-botcol[1]-botcol[2]);

minnumrows=step*minnumrows;
for(k=0;k<=3;k+=1)
{
topcol[k]=step*topcol[k];
botcol[k]=step*botcol[k];
}

numcols=botcol[3];

toprow[0]=0;
for(k=0;k<=3;k+=1)
{
if(k!=0)toprow[k]=botrow[k-1]+topcol[k]-botcol[k-1];
botrow[k]=toprow[k]-(botcol[k]-topcol[k]);
}

minbotrow=botrow[0];
for(k=1;k<=3;k+=1){if(botrow[k]<minbotrow)minbotrow=botrow[k];}

for(k=0;k<=3;k+=1)
{
toprow[k]-=minbotrow;
botrow[k]-=minbotrow;
}

maxtoprow=toprow[0];
for(k=1;k<=3;k+=1){if(toprow[k]>maxtoprow)maxtoprow=toprow[k];}

numrows=maxtoprow+minnumrows;

for(k=0;k<=3;k+=1)
{
if(k==0){for(col=0;col<=topcol[0]-1;col+=1)boundaryrow[col]=botrow[3]+col;}
else
{for(col=botcol[k-1];col<=topcol[k]-1;col+=1)boundaryrow[col]=botrow[k-1]+col-botcol[k-1];}

for(col=topcol[k];col<=botcol[k]-1;col+=1)
{
boundaryrow[col]=toprow[k]-(col-topcol[k]);
}
if(k==0||k==2)
{
boundaryrow[topcol[k]-1]-=1;
boundaryrow[topcol[k]]-=2;
boundaryrow[topcol[k]+1]-=1;
}
if(k==1||k==3)
{
boundaryrow[topcol[k]-5]-=1;
boundaryrow[topcol[k]-4]-=2;
boundaryrow[topcol[k]-3]-=3;
boundaryrow[topcol[k]-2]-=4;
boundaryrow[topcol[k]-1]-=5;
boundaryrow[topcol[k]]-=6;
boundaryrow[topcol[k]+1]-=5;
boundaryrow[topcol[k]+2]-=4;
boundaryrow[topcol[k]+3]-=3;
boundaryrow[topcol[k]+4]-=2;
boundaryrow[topcol[k]+5]-=1;
}
}
boundaryrow[numcols]=botrow[3];
}

void setupxy(void)
{
for(col=0;col<=topcol[0];col+=1)y[Q*boundaryrow[col]+col]
=((1.0*col-1.0*topcol[0])/(1.0*topcol[0]))*perpdist[3];
for(col=topcol[0]+1;col<=botcol[0];col+=1)y[Q*boundaryrow[col]+col]
=((1.0*col-1.0*topcol[0])/(1.0*botcol[0]-1.0*topcol[0]))*perpdist[1];

for(col=botcol[0]+1;col<=topcol[1];col+=1)x[Q*boundaryrow[col]+col]
=((1.0*topcol[1]-1.0*col)/(1.0*topcol[1]-1.0*botcol[0]))*perpdist[0];
for(col=topcol[1]+1;col<=botcol[1];col+=1)x[Q*boundaryrow[col]+col]
=((1.0*topcol[1]-1.0*col)/(1.0*botcol[1]-1.0*topcol[1]))*perpdist[2];

for(col=botcol[1]+1;col<=topcol[2];col+=1)y[Q*boundaryrow[col]+col]
=((1.0*topcol[2]-1.0*col)/(1.0*topcol[2]-1.0*botcol[1]))*perpdist[1];
for(col=topcol[2]+1;col<=botcol[2];col+=1)y[Q*boundaryrow[col]+col]
=((1.0*topcol[2]-1.0*col)/(1.0*botcol[2]-1.0*topcol[2]))*perpdist[3];

for(col=botcol[2]+1;col<=topcol[3];col+=1)x[Q*boundaryrow[col]+col]
=((1.0*col-1.0*topcol[3])/(1.0*topcol[3]-1.0*botcol[2]))*perpdist[2];
for(col=topcol[3]+1;col<=botcol[3];col+=1)x[Q*boundaryrow[col]+col]
=((1.0*col-1.0*topcol[3])/(1.0*botcol[3]-1.0*topcol[3]))*perpdist[0];

for(col=0;col<=numcols;col+=1)
{
x[Q*numrows+col]=a*cos(2.0*PI*(col-topcol[0])/(1.0*numcols));
y[Q*numrows+col]=a*sin(2.0*PI*(col-topcol[0])/(1.0*numcols));
}

dothebound();

for(col=0;col<=numcols;col+=1)
{
for(row=boundaryrow[col]+1;row<=numrows-1;row+=1)
{
x[Q*row+col]=(x[Q*boundaryrow[col]+col]*(1.0*numrows-1.0*row)
             +x[Q*numrows+col]*(1.0*row-1.0*boundaryrow[col]))
             /(1.0*numrows-1.0*boundaryrow[col]);
y[Q*row+col]=(y[Q*boundaryrow[col]+col]*(1.0*numrows-1.0*row)
             +y[Q*numrows+col]*(1.0*row-1.0*boundaryrow[col]))
             /(1.0*numrows-1.0*boundaryrow[col]);
}
}
}

void dothebound(void)
{
for(k=0;k<=3;k+=1)
{
if(k==0){for(col=0;col<=topcol[0]-1;col+=1)setbound();}
else {for(col=botcol[k-1];col<=topcol[k]-1;col+=1)setbound();}

for(col=topcol[k];col<=botcol[k];col+=1)setbound();
}

for(row=0;row<=botrow[3];row+=1)
{
x[Q*row+numcols]=perpdist[0];
y[Q*row]=-perpdist[3];
}
}

void setbound(void)
{
for(row=0;row<=boundaryrow[col];row+=1)
{
if(k==0){x[Q*row+col]= perpdist[0];if(row!=boundaryrow[col])y[Q*row+col]=y[Q*boundaryrow[col]+col];}
if(k==1){y[Q*row+col]= perpdist[1];if(row!=boundaryrow[col])x[Q*row+col]=x[Q*boundaryrow[col]+col];}
if(k==2){x[Q*row+col]=-perpdist[2];if(row!=boundaryrow[col])y[Q*row+col]=y[Q*boundaryrow[col]+col];}
if(k==3){y[Q*row+col]=-perpdist[3];if(row!=boundaryrow[col])x[Q*row+col]=x[Q*boundaryrow[col]+col];}
}
}

void typicalnodexy(void)
{
if(col>=1)whichcol=col-1;else whichcol=numcols-1;

tx=x[Q*row+col];
ty=y[Q*row+col];
polar();

lambda=1.0-0.0004*(1.5*numrows-1.0*row)*(1.0+0.05*(1.0-cos(2.0*theta)));
//lambda=1.0-0.1*cos(8.0*theta);

mu=2.0-lambda;

tx=(x[Q*row+col+1]+x[Q*row+whichcol])/4.0-x[Q*row+col]/2.0;
ty=(y[Q*row+col+1]+y[Q*row+whichcol])/4.0-y[Q*row+col]/2.0;
tz=(z[Q*row+col+1]+z[Q*row+whichcol])/4.0-z[Q*row+col]/2.0;

tx+=(lambda*x[Q*(row+1)+col]+mu*x[Q*(row-1)+col])/4.0-x[Q*row+col]/2.0;
ty+=(lambda*y[Q*(row+1)+col]+mu*y[Q*(row-1)+col])/4.0-y[Q*row+col]/2.0;
tz+=(lambda*z[Q*(row+1)+col]+mu*z[Q*(row-1)+col])/4.0-z[Q*row+col]/2.0;

xmovement[Q*row+col]=overrelax*tx;
ymovement[Q*row+col]=overrelax*ty;
zmovement[Q*row+col]=overrelax*tz;

tempfloat=xmovement[Q*row+col]*xnormal[Q*row+col]
         +ymovement[Q*row+col]*ynormal[Q*row+col]
         +zmovement[Q*row+col]*znormal[Q*row+col];

xmovement[Q*row+col]-=xnormal[Q*row+col]*tempfloat;
ymovement[Q*row+col]-=ynormal[Q*row+col]*tempfloat;
zmovement[Q*row+col]-=znormal[Q*row+col]*tempfloat;
}

void calczanalytic(void)
{
for(row=0;row<=numrows;row+=1)
{
for(col=0;col<=numcols;col+=1)
{
for(wherex=0;wherex<=2;wherex+=1)
{
for(wherey=0;wherey<=2;wherey+=1)
{
if(wherex==1||wherey==1)
{
if(wherex==1&&wherey==1&&row<=boundaryrow[col])storezvalue[wherex][wherey]=edgeheight;
else
{
if(wherex==1&&wherey==1&&row==numrows)storezvalue[wherex][wherey]=centreheight;
else
{
xvalue=x[Q*row+col]+(1.0+1.0*wherex)*delta;
yvalue=y[Q*row+col]+(1.0+1.0*wherey)*delta;
findtheheight();
storezvalue[wherex][wherey]=zvalue;
}
}
}
}
}
z[Q*row+col]=storezvalue[1][1];
dzbydx[Q*row+col]=(storezvalue[2][1]-storezvalue[0][1])/(2.0*delta);
dzbydy[Q*row+col]=(storezvalue[1][2]-storezvalue[1][0])/(2.0*delta);
}
}

for(row=numrows;row>=0;row-=1)
{
for(col=0;col<=numcols;col+=1)
{
if((col==0||col==botcol[0]||col==botcol[1]||col==botcol[2]||col==botcol[3])
&&row<=boundaryrow[col])
{
dzbydx[Q*row+col]=dzbydx[Q*(row+1)+col];
dzbydy[Q*row+col]=dzbydy[Q*(row+1)+col];
}
if(((col>botcol[0]&&col<botcol[1])||(col>botcol[2]&&col<botcol[3]))
&&row<=boundaryrow[col])dzbydx[Q*row+col]=0.0;
if(((col>botcol[1]&&col<botcol[2])||(col>0&&col<botcol[0]))
&&row<=boundaryrow[col])dzbydy[Q*row+col]=0.0;
}
}
}

void calcxy(void)
{
//int pause;
double doublepause;

previousKE=0.0;
for(col=0;col<=numcols-1;col+=1)
{
for(row=boundaryrow[col]+1;row<=numrows-1;row+=1)
{
prevxmovement[Q*row+col]=0.0;
prevymovement[Q*row+col]=0.0;
prevzmovement[Q*row+col]=0.0;
}
}

//maxnumxycycles=200;
maxnumxycycles=1;//5;//5;

//overallscale=1.6;
/*overallscale=100/(perpdist[3]+perpdist[1]);
if(overallscale>100/(perpdist[0]+perpdist[2]))overallscale=200/(perpdist[0]+perpdist[2]);
overallscale=0.95*overallscale;*/

//overallscale=9.5;
//xshrink=0.9;
//xshrink=1.0;

//overrelax=1.0;//This works fine!!
//overrelax=1.0;//0.5;//This slows it down!
if(LeftMouseDown==1){overrelax=1.0;LeftMouseDown=0;RightMouseDown=0;}
if(RightMouseDown==1){overrelax=0.0;LeftMouseDown=0;RightMouseDown=0;
//setup();

//Q=numcols+1;

//for(run=0;run<=1;run+=1)
//{

  
setupxy();
}
carryover=0.99;//0.9999;

findnormals();
for(row=0;row<=numrows;row+=1)
{
x[Q*row+numcols]=x[Q*row];
y[Q*row+numcols]=y[Q*row];
z[Q*row+numcols]=z[Q*row];
}

//glColor3f(1.0,0.1*DrawingDone, 0.0);

//if(DrawingDone==LastDrawingDone)glColor3f(1.0,0.0,1.0);

doublepause=0.0;
//if(DrawingDone==1)for(pause=0;pause<=1000000000;pause++)doublepause=sqrt(fabs(cos(doublepause)))*sqrt(fabs(cos(doublepause)));

DrawTheGrid();

intermediate=10;
intermediatecycle=0;
cycle=-1;
for(;;)
{
cycle+=1;

if(intermediatecycle==0||cycle==maxnumxycycles)
{
for(col=0;col<=numcols;col+=1)
{
for(row=0;row<=numrows;row+=1)
{
lastx[Q*row+col]=x[Q*row+col];
lasty[Q*row+col]=y[Q*row+col];
lastz[Q*row+col]=z[Q*row+col];
}
}
}

//if(cycle==0||cycle==maxnumxycycles)while (!Button());
if(cycle==maxnumxycycles)break;

if(intermediatecycle==0)findnormals();

intermediatecycle+=1;
if(intermediatecycle>intermediate)intermediatecycle=0;

myKE=0.0;

for(col=0;col<=numcols-1;col+=1)
{
for(row=boundaryrow[col]+1;row<=numrows-1;row+=1)
{
typicalnodexy();
}
}

for(col=0;col<=numcols-1;col+=1)
{
for(row=boundaryrow[col]+1;row<=numrows-1;row+=1)
{
xmovement[Q*row+col]+=carryover*prevxmovement[Q*row+col];
ymovement[Q*row+col]+=carryover*prevymovement[Q*row+col];
zmovement[Q*row+col]+=carryover*prevzmovement[Q*row+col];

myKE+=xmovement[Q*row+col]*xmovement[Q*row+col]
   +ymovement[Q*row+col]*ymovement[Q*row+col]
   +zmovement[Q*row+col]*zmovement[Q*row+col];

x[Q*row+col]+=xmovement[Q*row+col];
y[Q*row+col]+=ymovement[Q*row+col];
z[Q*row+col]+=zmovement[Q*row+col];

if(col==topcol[0]||col==topcol[2])y[Q*row+col]=0.0;
if(col==topcol[1]||col==topcol[3])x[Q*row+col]=0.0;

prevxmovement[Q*row+col]=xmovement[Q*row+col];
prevymovement[Q*row+col]=ymovement[Q*row+col];
prevzmovement[Q*row+col]=zmovement[Q*row+col];
}
}

if(myKE<previousKE)
{
previousKE=0.0;
for(col=0;col<=numcols-1;col+=1)
{
for(row=boundaryrow[col]+1;row<=numrows-1;row+=1)
{
prevxmovement[Q*row+col]=0.0;
prevymovement[Q*row+col]=0.0;
prevzmovement[Q*row+col]=0.0;
}
}
}
else previousKE=myKE;

for(row=0;row<=numrows;row+=1)
{
x[Q*row+numcols]=x[Q*row];
y[Q*row+numcols]=y[Q*row];
z[Q*row+numcols]=z[Q*row];
}
}
}

void polar(void)
{
r=sqrt(tx*tx+ty*ty);
if(fabs(ty)>fabs(tx)){theta=acos(tx/r);if(ty<0.0)theta=-theta;}
else{theta=asin(ty/r);if(tx<0.0)theta=PI-theta;}
}

void findtheheight(void)
{
tx=xvalue;
ty=yvalue;
polar();

tempfloat=(1.0-xvalue/perpdist[0])*(1.0+xvalue/perpdist[2])
         *(1.0-yvalue/perpdist[1])*(1.0+yvalue/perpdist[3]);
         
zvalue=(centreheight-edgeheight)*tempfloat/
((1.0-(a/r)*(xvalue/perpdist[0]))*(1.0+(a/r)*(xvalue/perpdist[2]))
*(1.0-(a/r)*(yvalue/perpdist[1]))*(1.0+(a/r)*(yvalue/perpdist[3])))
+edgeheight;

grad=0.0;

xdistance=perpdist[0]-xvalue;ydistance=perpdist[1]-yvalue;
grad+=sqrt(xdistance*xdistance+ydistance*ydistance)/(xdistance*ydistance);

xdistance=perpdist[0]-xvalue;ydistance=perpdist[3]+yvalue;
grad+=sqrt(xdistance*xdistance+ydistance*ydistance)/(xdistance*ydistance);

xdistance=perpdist[2]+xvalue;ydistance=perpdist[3]+yvalue;
grad+=sqrt(xdistance*xdistance+ydistance*ydistance)/(xdistance*ydistance);

xdistance=perpdist[2]+xvalue;ydistance=perpdist[1]-yvalue;
grad+=sqrt(xdistance*xdistance+ydistance*ydistance)/(xdistance*ydistance);

ConeAmount=0.5;                          //This controls "foldiness" at corners
zvalue+=ConeAmount*
(
3.5*(1.0+cos(2.0*theta))/2.0
+3.0*(1.0-cos(2.0*theta))/2.0
+0.3*sin(theta)
)
*(1.0-(a/r))/grad;

zvalue+=(1.0-ConeAmount)*
(
(35.0+10.0*tempfloat)*(1.0+cos(2.0*theta))/2.0
//=1 when theta=0,180;=0 when theta=90,270
+24.0*(((1.0-cos(2.0*theta))/2.0)+sin(theta))/2.0
//=0 when theta=0,180, 270;=1 when theta=90
+(7.5+12.0*tempfloat)*(((1.0-cos(2.0*theta))/2.0)-sin(theta))/2.0
//=0 when theta=0,180,  90;=1 when theta=270
-1.6//Lowers whole thing
)
*(r/a-1.0)*tempfloat;

//This controls short span:
zvalue-=10.0*((1.0+cos(2.0*theta))/2.0)*(r/a-1.0)*tempfloat;

//This raises roof over North Portico:
zvalue+=10.0*pow(((((1.0-cos(2.0*theta))/2.0)+sin(theta))/2.0),2.0)
*(r/a-1.0)*tempfloat*(1.0-3.0*(r/a-1.0)*tempfloat);

//This raises roof over South Portico:
zvalue+=2.5*pow(((((1.0-cos(2.0*theta))/2.0)-sin(theta))/2.0),2.0)
*pow((r/a-1.0),3.0)*tempfloat;

//This controls corners:
//zvalue-=5.0*((1.0-cos(4.0*theta))/2.0)*(r/a-1.0)*tempfloat;

//This is Waagner Buro's four corners raise
Waagner=14.0;
zvalue+=1.05*((1.0-(a/r))/grad)*
(exp(Waagner*(-(1.0-xvalue/perpdist[0])))+
 exp(Waagner*(-(1.0+xvalue/perpdist[2]))))*
(exp(Waagner*(-(1.0-yvalue/perpdist[1])))+
 exp(Waagner*(-(1.0+yvalue/perpdist[3]))));
}

void findnormals(void)
{
if(normalcalc==1)
{
for(col=0;col<=numcols-1;col+=1)
{
for(row=0;row<=numrows;row+=1)
{
if(row<=boundaryrow[col])z[Q*row+col]=edgeheight;
else
{
if(row==numrows)z[Q*row+col]=centreheight;
else
{
xvalue=x[Q*row+col];
yvalue=y[Q*row+col];
findtheheight();
z[Q*row+col]=zvalue;
}
}
}
}
for(col=0;col<=numcols-1;col+=1)
{
if(col>=1)whichcol=col-1;else whichcol=numcols-1;
for(row=boundaryrow[col]+1;row<=numrows-1;row+=1)
{
xnormal[Q*row+col]=(y[Q*(row+1)+col]-y[Q*(row-1)+col])*(z[Q*row+col+1]-z[Q*row+whichcol])
                  -(z[Q*(row+1)+col]-z[Q*(row-1)+col])*(y[Q*row+col+1]-y[Q*row+whichcol]);
ynormal[Q*row+col]=(z[Q*(row+1)+col]-z[Q*(row-1)+col])*(x[Q*row+col+1]-x[Q*row+whichcol])
                  -(x[Q*(row+1)+col]-x[Q*(row-1)+col])*(z[Q*row+col+1]-z[Q*row+whichcol]);
znormal[Q*row+col]=(x[Q*(row+1)+col]-x[Q*(row-1)+col])*(y[Q*row+col+1]-y[Q*row+whichcol])
                  -(y[Q*(row+1)+col]-y[Q*(row-1)+col])*(x[Q*row+col+1]-x[Q*row+whichcol]);

tempfloat=xnormal[Q*row+col]*xnormal[Q*row+col]
         +ynormal[Q*row+col]*ynormal[Q*row+col]
         +znormal[Q*row+col]*znormal[Q*row+col];

tempfloat=sqrt(tempfloat);

xnormal[Q*row+col]=xnormal[Q*row+col]/tempfloat;
ynormal[Q*row+col]=ynormal[Q*row+col]/tempfloat;
znormal[Q*row+col]=znormal[Q*row+col]/tempfloat;
}
}
}

else
{
for(col=0;col<=numcols-1;col+=1)
{
for(row=boundaryrow[col]+1;row<=numrows-1;row+=1)
{
xnormal[Q*row+col]=0.0;
ynormal[Q*row+col]=0.0;
znormal[Q*row+col]=1.0;
}
}
}
}

void DrawTheGrid(void)
{
/*xcamera=40.0;
ycamera=60.0;
zcamera=40.0;

scale=4.3;
xshift=-80.0;
yshift=85.0;


cameradistance=sqrt(xcamera*xcamera+ycamera*ycamera);
crot=xcamera/cameradistance;
srot=ycamera/cameradistance;

for(row=0;row<=numrows;row+=1)
{
for(col=0;col<=numcols;col+=1)
{
tx=x[Q*row+col];
ty=y[Q*row+col];
xtoplot[Q*row+col]=+tx*crot+ty*srot;
ytoplot[Q*row+col]=-tx*srot+ty*crot;
tx=xtoplot[Q*row+col];
ty=ytoplot[Q*row+col];
xtoplot[Q*row+col]=ty*cameradistance/(cameradistance-tx);
ytoplot[Q*row+col]=(z[Q*row+col]-zcamera)*cameradistance/(cameradistance-tx)
                            +zcamera;
xtoplot[Q*row+col]=scale*xtoplot[Q*row+col];
ytoplot[Q*row+col]=scale*ytoplot[Q*row+col];
}
}*/

//scale=1.0;
xshift=(perpdist[3]-perpdist[1])/2.0;
yshift=0.0;

/*for(row=0;row<=numrows;row+=1)
{
for(col=0;col<=numcols;col+=1)
{
xtoplot[Q*row+col]=-scale*y[Q*row+col];
ytoplot[Q*row+col]=scale*x[Q*row+col];
}
}*/

//if(run==0)
{
glColor4f(0.0,0.0,1.0,1.0);glLineWidth(1.0);
//if(DrawingDone==0)glColor4f(0.0,0.0,1.0,0.5);
//if(DrawingDone==LastDrawingDone)glColor4f(1.0,0.0,0.0,1.0);
for(col=0;col<=numcols-1;col+=1)
{
for(row=0;row<=numrows-1;row+=1)
{
end1=Q*row+col;
end2=Q*(row+1)+col;
drawagridline();
}
}
for(row=0;row<=numrows;row+=1)
{
for(col=0;col<=numcols-1;col+=1)
{
end1=Q*row+col;
end2=Q*row+col+1;
drawagridline();
}
}
}
//else
{
glColor4f(0.0,0.0,0.0,1.0);glLineWidth(2.0);
//if(DrawingDone==0){glColor4f(1.0,0.0,0.0,1.0);glLineWidth(2.0);}
//if(DrawingDone==LastDrawingDone){glColor4f(1.0,0.0,0.0,1.0);glLineWidth(2.0);}
for(col=0;col<=numcols-step;col+=step)
{
if(
col!=topcol[0]&&
col!=topcol[1]&&
col!=topcol[2]&&
col!=topcol[3]&&
col!=topcol[1]-2*step&&col!=topcol[1]+2*step&&
col!=topcol[3]-2*step&&col!=topcol[3]+2*step
)
{
for(row=boundaryrow[col]+2*step;row<=numrows;row+=2*step)
{
end1=Q*(row-2*step)+col;
end2=Q*row+col;
drawagridline();
}
}
else
{
end1=Q*boundaryrow[col]+col;
end2=Q*(boundaryrow[col]+step)+col;
drawagridline();
for(row=boundaryrow[col]+3*step;row<=numrows;row+=2*step)
{
end1=Q*(row-2*step)+col;
end2=Q*row+col;
drawagridline();
}
}
}

diagonal=step;
for(col=0;col<=numcols-step;col+=step)
{
if(diagonal==step)diagonal=0;else diagonal=step;
for(row=diagonal;row<=numrows-step;row+=2*step)
{
if(row>=boundaryrow[col]&&row+step!=boundaryrow[col+step])
{
end1=Q*row+col;
if(col!=numcols-step)end2=Q*(row+step)+col+step;
else end2=Q*(row+step);
drawagridline();
}
}
}

diagonal=0;
for(col=step;col<=numcols;col+=step)
{
if(diagonal==step)diagonal=0;else diagonal=step;
for(row=diagonal;row<=numrows-step;row+=2*step)
{
if(row>=boundaryrow[col]&&row+step!=boundaryrow[col-step])
{
if(col!=numcols)end1=Q*row+col;
else end1=Q*row;
end2=Q*(row+step)+col-step;
drawagridline();
}
}
}

//Rectanglar in plan edge beam i.e. outside edge beam
col=0;
row=0;
end1=Q*row+col;
for(col=step;col<=numcols;col+=step)
{
end2=Q*row+col;
drawagridline();
end1=end2;
}

//Reading Room edge beam
row=numrows;
col=0;
end1=Q*row+col;
for(col=step;col<=numcols;col+=step)
{
end2=Q*row+col;
drawagridline();
end1=end2;
}
}
}




void drawagridline(void)
{
/*if(cycle!=0)
{
MainLineColour.red=int(65535.0*1.0);
MainLineColour.green=int(65535.0*0.1);
MainLineColour.blue=int(65535.0*0.1);
RGBForeColor(&MainLineColour);
intxdraw[0]=int(overallscale*(xshift+scale*lasty[end1]));
intydraw[0]=int(overallscale*(yshift-scale*lastx[end1]));
intxdraw[1]=int(overallscale*(xshift+scale*lasty[end2]));
intydraw[1]=int(overallscale*(yshift-scale*lastx[end2]));
     MoveTo(halfwidth+intxdraw[0],halfheight-intydraw[0]);
     LineTo(halfwidth+intxdraw[1],halfheight-intydraw[1]);
}

MainLineColour.red=int(65535.0*1.0);
MainLineColour.green=int(65535.0*1.0);
MainLineColour.blue=int(65535.0*1.0);
RGBForeColor(&MainLineColour);
intxdraw[0]=int(overallscale*(xshift+scale*y[end1]));
intydraw[0]=int(overallscale*(yshift-scale*x[end1]));
intxdraw[1]=int(overallscale*(xshift+scale*y[end2]));
intydraw[1]=int(overallscale*(yshift-scale*x[end2]));*/
    // MoveTo(halfwidth+intxdraw[0],halfheight-intydraw[0]);
    // LineTo(halfwidth+intxdraw[1],halfheight-intydraw[1]);

      rv[0][0] = overallscale*(xshift+y[end1]);//*xshrink;
      rv[0][1] = overallscale*(yshift-x[end1]);

      rv[1][0] = overallscale*(xshift+y[end2]);//*xshrink;
      rv[1][1] = overallscale*(yshift-x[end2]);

//cout<<rv[0][0]<<"  "<<rv[0][1]<<"  "<<rv[1][0]<<"  "<<rv[1][1]<<"\n";

#define LINE_STEP 4

glBegin(GL_LINE_STRIP);
	glVertex2fv(rv[0]);
	glVertex2fv(rv[1]);
glEnd();
    // MoveTo(halfwidth+intxdraw[0],halfheight-intydraw[0]);
   //  LineTo(halfwidth+intxdraw[1],halfheight-intydraw[1]);

/*if(run==0)
{
Chris<<"0\nLINE\n8\n0\n";
Chris<<"10\n"<<x[end1]<<"\n";
Chris<<"20\n"<<y[end1]<<"\n";
Chris<<"11\n"<<x[end2]<<"\n";
Chris<<"21\n"<<y[end2]<<"\n";
if(cycle==maxnumxycycles)Chris<<"62\n0\n";else Chris<<"62\n1\n";
}*/
}

  void init_graphics(void)
{
int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	
	cout<<"Your screen resolution is "<<screenWidth<<" x "<<screenHeight<<"\n";
	
	int WindowHalfWidth =int(double(screenWidth)/2.0)-5;
	int WindowHalfHeight=int(double(screenHeight)/2.0)-25;
	
	glutInitWindowSize(2*WindowHalfWidth, 2*WindowHalfHeight);
	
	float overallscalex=0.95*float(2*WindowHalfWidth)/(46.025+51.125);
	float overallscaley=0.95*float(2*WindowHalfHeight)/(36.625+36.625);
    if(overallscaley>overallscalex)overallscale=overallscalex;
	else overallscale=overallscaley;

//perpdist[0]=36.625;  perpdist[1]=46.025;  perpdist[2]=36.625;  perpdist[3]=51.125;

/*overallscale=100/(perpdist[3]+perpdist[1]);
if(overallscale>100/(perpdist[0]+perpdist[2]))overallscale=200/(perpdist[0]+perpdist[2]);
overallscale=0.95*overallscale;*/

 //   glutInitWindowSize(2*halfwidth, 2*halfheight);
   // if(Macintosh==1)    glutInit(&argc, argv);
   /* glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("British Museum");
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glViewport(0, 0, 2*halfwidth, 2*halfheight);
    gluOrtho2D(-halfwidth,halfwidth,-halfheight,halfheight);*/
	
	//InitialZoom=0.95*double(WindowHalfHeight)/MaxHalfDimension;
	//Zoom=InitialZoom;
	glutInitWindowPosition(0,0);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Relax");
	glClearColor(1.0, 1.0, 1.0, 0.0);
	//glOrtho(-WindowHalfWidth,WindowHalfWidth,-WindowHalfHeight,WindowHalfHeight,-10.0*WindowHalfWidth,10.0*WindowHalfWidth);
	gluOrtho2D(-WindowHalfWidth,WindowHalfWidth,-WindowHalfHeight,WindowHalfHeight);
	//glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
	glClear(GL_COLOR_BUFFER_BIT);
	glutKeyboardFunc(Key);
	glutIdleFunc(Draw);
	glutDisplayFunc(Draw);
	glutMouseFunc(MouseClick);
//	glutMotionFunc(MouseDrag);
//	glutPassiveMotionFunc(MouseMove);
	glEnable(GL_BLEND);
	glEnable(GL_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_TEXTURE);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
//	glScaled(Zoom,Zoom,Zoom);
	//Rotate();

	
}

static void Draw(void)
{



//if(DrawingDone<=LastDrawingDone)
{
//if(DrawingDone==0)
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    

Museum();

 /*   rv[0][0] = 40*(DrawingDone-5);
    rv[0][1] = 0.0;
    glColor3f(1.0,0.1*DrawingDone, 0.0);
    
    for(count=0;count<=10000;count++)
    {
      rv[1][0] = rv[0][0]+200.0*cos(count*2.0*PI/10000.0);
      rv[1][1] = rv[0][1]+200.0*sin(count*2.0*PI/10000.0);

#define LINE_STEP 4

glBegin(GL_LINE_STRIP);
	glVertex2fv(rv[0]);
	glVertex2fv(rv[1]);
glEnd();
}*/
glutSwapBuffers();
DrawingDone+=1;
}
}

static void Key(unsigned char key, int x, int y)
{
switch (key) {
case 'q':gluDeleteQuadric(quadObj);exit(0);
case '\033':gluDeleteQuadric(quadObj);exit(0);
}
}

void Museum(void)
{
//cout<<DrawingDone<<"\n";
if(DrawingDone==0)
{
run=1;

 
PI=4.0*atan(1.0);
step=2;

//Dim>>a;
//for(k=0;k<=3;k+=1)Dim>>perpdist[k];
//Dim.close();

a=22.245;
perpdist[0]=36.625;  perpdist[1]=46.025;  perpdist[2]=36.625;  perpdist[3]=51.125;

//Top>>minnumrows;
//for(k=0;k<=3;k+=1)Top>>topcol[k];
//for(k=0;k<=2;k+=1)Top>>botcol[k];
//Top.close();

minnumrows=18;
topcol[0]=17;    topcol[1]=47;    topcol[2]=77;    topcol[3]=107;
   botcol[0]=34;    botcol[1]=60;    botcol[2]=94;
   
setup();

Q=numcols+1;

//for(run=0;run<=1;run+=1)
//{

  
setupxy();
}
calcxy();
//}
//Chris<<"0\nENDSEC\n0\nEOF\n";Chris.close();
//while (!Button());
//    HideWindow(fullScreenWindow);
//    EndFullScreen(oldState,nil);

//    return 0;

}

void MouseClick(int button, int state, int x, int y) 
{
	if((button==GLUT_LEFT_BUTTON)&&(state==GLUT_DOWN))
	{
		LeftMouseDown=1;
		}
		if((button==GLUT_RIGHT_BUTTON)&&(state==GLUT_DOWN)) {
		RightMouseDown=1;
	}
		}