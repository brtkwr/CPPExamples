//Written by Chris J K Williams, University of Bath, UK
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include "Graphics.h"
#include "BritishMuseum.h"
string MyText;

int main(int argc, char * argv[])
{
	cout<<"m = "<<m<<"\n";
	cout<<"n = "<<n<<"\n";
	cout<<"The program will now run forever until 'q' or ESCAPE is pressed on the keyboard\n";
	cout<<"Left mouse button for pan\n";
	cout<<"Right mouse button (or CONTROL left button) for rotate\n";
	cout<<"Move spot to zoom\n";
	
	SetupData();
	
	glutInit(&argc,argv);
	
	MyGraphics(1.0,1.0,1.0);
	return 0;
}

static void Draw(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	ButtonControls();
	if(freshStart == 1)start = 0;
	
	CalculateAndDraw();
	
	glutSwapBuffers();
}

void ButtonControls(void)
{
	glPushMatrix();
	glLoadIdentity();
	
	glPointSize(4.0);
	for(int Button = 0; Button <= lastButton; Button ++)
	{
		glBegin(GL_POINTS);
		
		if(Button <= 3)
		{
			if(ButtonOn[Button] == 1)glColor4f(0.0,0.0,1.0,1.0);else glColor4f(0.0,0.0,0.0,1.0);
		}
		else {
			if(Button == 4){if(start == 0)glColor4f(0.0,1.0,0.0,1.0);else glColor4f(0.0,0.0,0.0,1.0);}
			if(Button == 5){if(start == 1)glColor4f(1.0,0.0,0.0,1.0);else glColor4f(0.0,0.0,0.0,1.0);}
			if(Button == 6)glColor4f(0.0,0.0,1.0,1.0);
		}
		
		PointPosition[0] = ButtonX[Button];
		PointPosition[1] = ButtonY[Button];
		glVertex2dv(PointPosition);
		glEnd();
		
		glColor4f(0.0,0.0,0.0,0.2);
		if(Button == 0)MyText = "zoom";
		if(Button == 1)MyText = "lambda";
		if(Button == 2)MyText = "alpha";
		if(Button == 3)MyText = "beta";
		if(Button == 4)MyText = "Start";
		if(Button == 5)MyText = "Stop";
		if(Button == 6)MyText = "Fresh start";
		
		glRasterPos2d(ButtonX[Button] + 0.02,ButtonY[Button] - 0.01);
		
		int StringLength = MyText.length();
		for(int MyChar = 0;MyChar<StringLength;MyChar++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,MyText[MyChar]);
	}
	
	glPopMatrix();
}