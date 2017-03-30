//Written by Chris J K Williams, University of Bath, UK

int LineExists(int,int,int,int);
void Axial(int,int,int,int);
void Fresh(void);
void Function(double,double);
void AddVertex(int,int,int,GLfloat,GLfloat,GLfloat,GLfloat);
void plotLinesOrPoints(int,int);
void ButtonControls(void);

#define  subdivision 2
#define  p (17 * subdivision)
#define  q (17 * subdivision)
#define  r_value (13 * subdivision)
#define  s (13 * subdivision)
#define  m (2 * (p + q + r_value + s)) - 1
#define  n (35 * subdivision)
#define  LastPossibleVertex 2 * m

GLfloat chris_vertices[3 * (LastPossibleVertex + 1)];
GLfloat chris_colours[4 * (LastPossibleVertex + 1)];

#define  CarryOver 0.99
#define  MovementFactor 0.5

#include "LinesOrPoints.h"
#include "Constants.h"

int  NodeType[m + 1][n + 1];

double x[m + 1][n + 1][3],CableNetForce[m + 1][n + 1][3],
velocity[m + 1][n + 1][3],Normal[m + 1][n + 1][3],
ReadingRoomPosition[m + 1][3],ReadingRoomHeight,
PointPosition[2],
CableNetTensionCoefficient,a,b,c,d,MotionStartControl,asq,
z,ThisNormal[3];

void CalculateAndDraw(void);
void SetupData(void);

void SetupData(void)
{
	CableNetTensionCoefficient = 0.1;
	
	a = 22.245;
	b = 36.625;
	c = 46.025;
	d = 51.125;
	ReadingRoomHeight = 20.955 - 19.71;
	
	asq = a * a;
	
	for(int WhichEighth = 0;WhichEighth <= 6; WhichEighth += 2)
	{
		int  i_start = 0;
		int i_peak = 0;
		int j_peak = 0;
		int i_stop = 0;
		int j_triangle = 0;
		
		if(WhichEighth == 0)
		{
			i_start = 0;
			i_peak = p;
			j_peak = p;
			i_stop = p + q - 1;
			j_triangle = 1 * subdivision;
		}
		if(WhichEighth == 2)
		{
			i_start = p + q;
			i_peak = p + q + r_value;
			j_peak = p - q + r_value;
			i_stop = p + q + 2 * r_value - 1;
			j_triangle = 3 * subdivision;
		}
		if(WhichEighth == 4)
		{
			i_start = p + q + 2 * r_value;
			i_peak = p + 2 * q + 2 * r_value;
			j_peak = p;
			i_stop = 2 * p + 2 * q + 2 * r_value - 1;
			j_triangle = 1 * subdivision;
		}
		if(WhichEighth == 6)
		{
			i_start = 2 * p + 2 * q + 2 * r_value;
			i_peak = 2 * p + 2 * q + 2 * r_value + s;
			j_peak = s;
			i_stop = m;
			j_triangle = 3 * subdivision;
		}
		
		for(int i = i_start;i <= i_stop; i ++)
		{
			for(int j = 0;j <= n;j ++)
			{
				if(j < n)NodeType[i][j] = 0;else NodeType[i][j] = 1;
				if(j == 0)
				{
					if(i == i_start)NodeType[i][j] = WhichEighth + 2;
					else NodeType[i][j] = WhichEighth + 3;
				}
				int CuttingPattern = 1;
				if(CuttingPattern == 1)
				{
					if(j <= j_peak - (i - i_peak) && j <=  j_peak + (i - i_peak) && j <= j_peak - j_triangle)
					{
						NodeType[i][j] = - 1;
						if(j == j_peak - (i - i_peak) || j ==  j_peak + (i - i_peak) || j == j_peak - j_triangle)
						{
							if(i == i_start)NodeType[i][j] = WhichEighth + 2;
							else NodeType[i][j] = WhichEighth + 3;
						}
					}
				}
			}
		}
	}
	
	double setupmultiplier = 8.0 * atan(1.0) / (double)(m + 1);
	
	for(int i = 0;i <= m;i ++)
	{
		double theta = (double) (i - p) * setupmultiplier;
		ReadingRoomPosition[i][0] = a * cos(theta);
		ReadingRoomPosition[i][1] = a * sin(theta);
		ReadingRoomPosition[i][2] = ReadingRoomHeight;
	}
	
	freshStart = 1;
}

void CalculateAndDraw(void)
{
	if(freshStart == 1)Fresh();
	
	glLineWidth(1.0);
	
	chris_vertices[ 0] = + b;chris_vertices[ 1] = - d;chris_vertices[ 2] = 0.0;
	chris_vertices[ 3] = + b;chris_vertices[ 4] = + c;chris_vertices[ 5] = 0.0;
	chris_vertices[ 6] = - b;chris_vertices[ 7] = + c;chris_vertices[ 8] = 0.0;
	chris_vertices[ 9] = - b;chris_vertices[10] = - d;chris_vertices[11] = 0.0;
	chris_vertices[12] = + b;chris_vertices[13] = - d;chris_vertices[14] = 0.0;
	
	for(int Corner = 0;Corner <= 4;Corner ++)
	{
		chris_colours[4 * Corner + 0] = 0.0;
		chris_colours[4 * Corner + 1] = 0.2;
		chris_colours[4 * Corner + 2] = 0.4;
		chris_colours[4 * Corner + 3] = 0.5;
	}
	
	plotLinesOrPoints(1,4);
	
	for(int itry = 0;itry <= m + 1;itry ++)
	{
		int i = itry;
		if(i > m)i -= m + 1;
		
		for(int xyz = 0; xyz <= 2;xyz ++)chris_vertices[3 * itry + xyz] = ReadingRoomPosition[i][xyz];
		
		chris_colours[4 * itry + 0] = 0.0;
		chris_colours[4 * itry + 1] = 0.2;
		chris_colours[4 * itry + 2] = 0.4;
		chris_colours[4 * itry + 3] = 0.5;
	}
	plotLinesOrPoints(1,m + 1);
	
	glLineWidth(1.0);
	
	int LineDrawStep = subdivision - 1;
	for(int j = 0;j <= n;j ++)
	{
		LineDrawStep ++;
		if(LineDrawStep >= subdivision)LineDrawStep = 0;
		
		int NowDrawing = 1;
		int myVertex = - 1;
		for(int itry = 0;itry <= m + start;itry ++)
		{
			int i = itry;
			if(i > m)i -= m + 1;
			if(NodeType[i][j] >= 0)
			{
				if(NowDrawing == 0)NowDrawing = 1;
				myVertex ++;
				if(LineDrawStep == 0)AddVertex(myVertex,i,j,1.0,0.0,0.0,0.5);
				else AddVertex(myVertex,i,j,1.0,0.0,0.0,0.25);
			}
			else
			{
				if(NowDrawing == 1)
				{
					plotLinesOrPoints(1,myVertex);
					myVertex = - 1;
					NowDrawing = 0;
				}
			}
			
		}
		plotLinesOrPoints(1,myVertex);
	}
	
	LineDrawStep = subdivision - 1;
	for(int i = 0;i <= m;i ++)
	{
		LineDrawStep ++;
		if(LineDrawStep >= subdivision)LineDrawStep = 0;
		
		int myVertex = - 1;
		for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] >= 0)
			{
				myVertex ++;
				if(LineDrawStep == 0)AddVertex(myVertex,i,j,1.0,0.0,0.0,0.5);
				else AddVertex(myVertex,i,j,1.0,0.0,0.0,0.25);
			}
		}
		plotLinesOrPoints(1,myVertex);
	}
	
	for(int i = subdivision;i <= m + 1 - subdivision;i += 2 * subdivision)
	{
		int myVertex = - 1;
		for(int j = n;j >= 0;j -= subdivision)
		{
			int Actual_i = i + (n - j);
			if(Actual_i > m && Actual_i - subdivision <= m && start == 0)
			{
				plotLinesOrPoints(1,myVertex);
				myVertex = - 1;
			}
			if(Actual_i > m)Actual_i -= m + 1;
			if(NodeType[Actual_i][j] >= 0)
			{
				myVertex ++;
				AddVertex(myVertex,Actual_i,j,0.0,0.0,0.0,0.5);
			}
		}
		plotLinesOrPoints(1,myVertex);
	}
	
	for(int i = subdivision;i <= m + 1 - subdivision;i += 2 * subdivision)
	{
		int myVertex = - 1;
		for(int j = n;j >= 0;j -= subdivision)
		{
			int Actual_i = i - (n - j);
			if(Actual_i < 0 && Actual_i + subdivision >= 0 && start == 0)
			{
				plotLinesOrPoints(1,myVertex);
				myVertex = - 1;
			}
			if(Actual_i < 0)Actual_i += m + 1;
			if(NodeType[Actual_i][j] >= 0)
			{
				myVertex ++;
				AddVertex(myVertex,Actual_i,j,0.0,0.0,0.0,0.5);
			}
		}
		plotLinesOrPoints(1,myVertex);
	}
	
	int LineStart = 0;
	for(int i = 0;i <= m;i += subdivision)
	{
		int myVertex = - 1;
		if(LineStart == subdivision)LineStart = 0;else LineStart = subdivision;
		for(int j = n - LineStart;j >= 0;j -= 2 * subdivision)
		{
			if(NodeType[i][j] >= 0)
			{
				myVertex ++;
				AddVertex(myVertex,i,j,0.0,0.0,0.0,0.5);
			}
			else
			{
				if(NodeType[i][j + subdivision] >= 0)
				{
					myVertex ++;
					AddVertex(myVertex,i,j + subdivision,0.0,0.0,0.0,0.5);
				}
			}
		}
		plotLinesOrPoints(1,myVertex);
	}
	
	int ShowNormals = 0;
	if(start == 1 && ShowNormals == 1)
	{
		glColor4f(0.0,0.0,1.0,1.0);
		
		for(int i = 0;i <= m;i += subdivision)
		{
			for(int j = n;j >= 0;j -= subdivision)
			{
				if(NodeType[i][j] == 0)
				{
					for(int xyz = 0; xyz <= 2;xyz ++)
					{
						chris_vertices[xyz] = x[i][j][xyz];
						chris_vertices[xyz + 3] = x[i][j][xyz] + 0.5 * Normal[i][j][xyz];
					}
					for(int NormalEnd = 0;NormalEnd <= 4;NormalEnd +=4)
					{
						chris_colours[NormalEnd + 0] = 0.0;
					chris_colours[NormalEnd + 1] = 0.2;
					chris_colours[NormalEnd + 2] = 1.0;
					chris_colours[NormalEnd + 3] = 0.5;
				}
					plotLinesOrPoints(1,1);
			}
		}
		}
	}
	
	glPointSize(2.0);
	int myVertex = - 1;
	for(int i = 0;i <= m;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			if(NodeType[i][j] > 0)
			{
				myVertex ++;
				AddVertex(myVertex,i,j,0.0,0.0,1.0,1.0);			}
		}
	}
	plotLinesOrPoints(0,myVertex);
	
	if(start == 1)
	{
		for(int CalculationLoop = 0;CalculationLoop <= 0;CalculationLoop++)
		{
			if(MotionStartControl < 1.0)MotionStartControl += 0.0001;
			
			for(int i = 0;i <= m;i ++)
			{
				for(int j = 0;j <= n;j ++)
				{
					for(int xyz = 0; xyz <= 2;xyz ++)CableNetForce[i][j][xyz] = 0.0;
				}
			}
			for(int i = 0;i <= m;i ++)
			{
				for(int j = 0;j <= n;j ++)
				{
					Axial(i,i + 1,j,j);
					if(j + 1 <= n)Axial(i,i,j,j + 1);
				}
			}
			for(int i = 0;i <= m;i ++)
			{
				for(int j = 0;j <= n - 1;j ++)
				{
					if(NodeType[i][j] == 0)
					{
						Function(x[i][j][0],x[i][j][1]);
						
						double ScalarProduct = 0.0;
						for(int xyz = 0;xyz <= 2;xyz ++)
						{
							Normal[i][j][xyz] = ThisNormal[xyz];
							ScalarProduct += ThisNormal[xyz] * CableNetForce[i][j][xyz];
						}
						
						for(int xyz = 0;xyz <= 2;xyz ++)CableNetForce[i][j][xyz] -= ScalarProduct * ThisNormal[xyz];
						
						for(int xyz = 0;xyz <= 2;xyz ++)
						{
							velocity[i][j][xyz] = CarryOver * velocity[i][j][xyz] + MovementFactor * CableNetForce[i][j][xyz] / (1.0 + lambda);
							x[i][j][xyz] += MotionStartControl * velocity[i][j][xyz];
						}
						
						double zcorrection = (x[i][j][2] - z) * ThisNormal[2];
						
						for(int xyz = 0;xyz <= 2;xyz ++)x[i][j][xyz] -= zcorrection * ThisNormal[xyz];
					}
				}
			}
			
			for(int i = 0;i <= m;i ++)
			{
				for(int xyz = 0; xyz <= 2;xyz ++)x[i][n][xyz] = ReadingRoomPosition[i][xyz];
				for(int j = 0;j <= n;j ++)
				{
					if(NodeType[i][j] == 2 || NodeType[i][j] == 3 || NodeType[i][j] == 4){x[i][j][0] = + b;x[i][j][2] = 0.0;}
					if(NodeType[i][j] == 4 || NodeType[i][j] == 5 || NodeType[i][j] == 6){x[i][j][1] = + c;x[i][j][2] = 0.0;}
					if(NodeType[i][j] == 6 || NodeType[i][j] == 7 || NodeType[i][j] == 8){x[i][j][0] = - b;x[i][j][2] = 0.0;}
					if(NodeType[i][j] == 8 || NodeType[i][j] == 9 || NodeType[i][j] == 2){x[i][j][1] = - d;x[i][j][2] = 0.0;}
					
					if(NodeType[i][j] == 3)
					{
						if(i <= p)x[i][j][1] = d * (double) (i - p) / (double) p;
						else x[i][j][1] = c * (double) (i - p) / (double) q;
					}
					
					if(NodeType[i][j] == 7)
					{
						if(i > p + 2 * q + 2 * r_value)x[i][j][1] = d * (double) (p + 2 * q + 2 * r_value - i) / (double) p;
						else x[i][j][1] = c * (double) (p + 2 * q + 2 * r_value - i) / (double) q;
					}
					
					if(NodeType[i][j] == 5)x[i][j][0] = b - 2.0 * b * (double) (i - (p + q)) / (double)(2 * r_value);
					if(NodeType[i][j] == 9)x[i][j][0] = b - 2.0 * b * (double) (m + 1 - i) / (double)(2 * r_value);
				}
			}
		}
	}
}

int LineExists(int i,int inext,int j,int jnext)
{
	if((NodeType[i][j] == 0 && NodeType[inext][jnext] >= 0) ||
	   (NodeType[i][j] >= 0 && NodeType[inext][jnext] == 0) ||
	   (NodeType[i][j] >= 0 && NodeType[inext][jnext] >= 0 && (j == 0 || j == n)))return 1;
	else return 0;
}

void Axial(int i,int inext_trial,int j,int jnext)
{
	double deltax[3];
	int inext = inext_trial;
	if(inext > m)inext = 0;
	if(LineExists(i,inext,j,jnext))
	{
		for(int xyz = 0;xyz <= 2;xyz ++)deltax[xyz] = x[inext][jnext][xyz] - x[i][j][xyz];
		
		double tensioncoefficient = CableNetTensionCoefficient;
		
		int CableDirection = 0;
		
		if(jnext == j){CableDirection = 1;tensioncoefficient *= lambda;}
		
		for(int xyz = 0;xyz <= 2;xyz ++)
		{
			double thisForce = tensioncoefficient * deltax[xyz];
			
			CableNetForce[i][j][xyz] += thisForce;
			CableNetForce[inext][jnext][xyz] -= thisForce;
		}
	}
}

void Fresh(void)
{
	for(int i = 0;i <= m;i ++)
	{
		for(int j = 0;j <= n;j ++)
		{
			x[i][j][1] = + 0.8 * (c + d) * (((double)(2 * i - m) / (double) (2 * m)) - 0.02);
			x[i][j][0] = - 0.8 * (c + d) * (((double)(2 * j - n) / (double) (2 * m)) - 0.3);
			x[i][j][2] = 0.0;
			
			for(int xyz = 0; xyz <= 2;xyz ++)velocity[i][j][xyz] = 0.0;
		}
	}
	
	lambda = 0.5;
	alpha = 0.8;
	beta = 2.0;
	
	MotionStartControl = 0.0;
	freshStart = 0;
}

void Function(double x, double y)
{
	double rsq = x * x + y * y;
	if(rsq < 0.5 * asq)rsq = 0.5 * asq;//This to avoid nodes geting stuck in the middle;
	double r = sqrt(rsq);
	double r_over_a = r / a;
	double z1,z2,z3,z1x,z2x,z3x,z1y,z2y,z3y;
	
	double T1 = b - x;
	double T2 = b + x;
	double T3 = c - y;
	double T4 = d + y;
	
	double B1 = b - x / r_over_a;
	double B2 = b + x / r_over_a;
	double B3 = c - y / r_over_a;
	double B4 = d + y / r_over_a;
	
	double Q1 = sqrt(1.0 / (T1 * T1) + 1.0 / (T3 * T3));
	double Q2 = sqrt(1.0 / (T2 * T2) + 1.0 / (T3 * T3));
	double Q3 = sqrt(1.0 / (T1 * T1) + 1.0 / (T4 * T4));
	double Q4 = sqrt(1.0 / (T2 * T2) + 1.0 / (T4 * T4));
	
	z1 = ReadingRoomHeight * T1 * T2 * T3 * T4 / (B1 * B2 * B3 * B4);
	
	z2 = (r - a) * T1 * T2 * T3 * T4 / (b * b * c * d);
	
	z3 =  (1.0 - a / r) / (Q1 + Q2 + Q3 + Q4);
	
	z1x = z1 * (- 1.0 / T1 + 1.0 / T2
				+ (a / r) * (1.0 - x * x / rsq) * (1.0 / B1 - 1.0 / B2)
				- (a / r) * (      x * y / rsq) * (1.0 / B3 - 1.0 / B4));
	z1y = z1 * (- 1.0 / T3 + 1.0 / T4
				+ (a / r) * (1.0 - y * y / rsq) * (1.0 / B3 - 1.0 / B4)
				- (a / r) * (      x * y / rsq) * (1.0 / B1 - 1.0 / B2));
	
	z2x = z2 * ((x / (r * a)) / (r_over_a - 1.0) - 1.0 / T1 + 1.0 / T2);
	z2y = z2 * ((y / (r * a)) / (r_over_a - 1.0) - 1.0 / T3 + 1.0 / T4);
	
	z3x = (z3 * z3 / (1.0 - a / r)) * ((a * x / (r * rsq)) / z3
									   - (1.0 / Q1 + 1.0 / Q3) / (T1 * T1 * T1)
									   + (1.0 / Q2 + 1.0 / Q4) / (T2 * T2 * T2));
	
	z3y = (z3 * z3 / (1.0 - a / r)) * ((a * y / (r * rsq)) / z3
									   - (1.0 / Q1 + 1.0 / Q2) / (T3 * T3 * T3)
									   + (1.0 / Q3 + 1.0 / Q4) / (T4 * T4 * T4));
	
	z = z1 + alpha * z2 + beta * z3;
	double zx = z1x + alpha * z2x + beta * z3x;
	double zy = z1y + alpha * z2y + beta * z3y;
	
	double NormalMagnitude = sqrt(1.0 + zx * zx + zy * zy);
	ThisNormal[0] = - zx / NormalMagnitude;
	ThisNormal[1] = - zy / NormalMagnitude;
	ThisNormal[2] =  1.0 / NormalMagnitude;
}

void AddVertex(int myVertex,int i,int j,GLfloat Red,GLfloat Green,GLfloat Blue,GLfloat alpha)
{
	if(myVertex <= LastPossibleVertex)
	{
		for(int xyz = 0; xyz <= 2;xyz ++)
			chris_vertices[3 * myVertex + xyz] = (GLfloat)x[i][j][xyz];
		chris_colours[4 * myVertex + 0] = Red;
		chris_colours[4 * myVertex + 1] = Green;
		chris_colours[4 * myVertex + 2] = Blue;
		chris_colours[4 * myVertex + 3] = alpha;
	}
}
