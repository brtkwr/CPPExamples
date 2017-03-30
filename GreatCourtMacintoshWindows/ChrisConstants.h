//Written by Chris J K Williams, University of Bath, UK

#define iPhone 0

#if (iPhone == 1)

//Use the following for iPhone

#define ChrisScale				4.0
#define xDirection 2
#define yDirection 0
#define zDirection 1
#define  subdivision 1

int freshStart = 0;
int start =1;
double lambda,alpha,beta;
int LastCalculationLoop = 5;
int ShowNormals = 0;
int ShowMathematicalGrid = 0;

GLfloat StructureRed   = 0.0;
GLfloat StructureGreen = 0.0;
GLfloat StructureBlue  = 0.0;
GLfloat StructureAlpha = 0.2;

GLfloat GridRed   = 1.0;
GLfloat GridGreen = 0.0;
GLfloat GridBlue  = 0.0;
GLfloat GridAlpha = 0.1;
GLfloat GridEmpahsisAlpha = 0.05;

double MotionStartControlIncrement = 0.001;

#else

//Use the following for Macintosh (or Windows)

 #define xDirection 0
 #define yDirection 1
 #define zDirection 2
 #define  subdivision 2
 
 int LastCalculationLoop = 5;
 int ShowNormals = 0;
 int ShowMathematicalGrid = 1;
 
 GLfloat StructureRed   = 0.0;
 GLfloat StructureGreen = 0.0;
 GLfloat StructureBlue  = 0.0;
 GLfloat StructureAlpha = 0.4;
 
 GLfloat GridRed   = 1.0;
 GLfloat GridGreen = 0.0;
 GLfloat GridBlue  = 0.0;
 GLfloat GridAlpha = 0.5;
 GLfloat GridEmpahsisAlpha = 0.2;
 
 
 double MotionStartControlIncrement = 0.0001;

#endif