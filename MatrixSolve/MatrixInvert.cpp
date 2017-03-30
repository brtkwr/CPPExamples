#include <cmath>
#include "size.h"

void Invert(int n,double C[][nplus1],double B[][nplus1])
{
double A[nplus1][nplus1],divide,mult,tempdouble;
int i,j,k,biggest;

for(i=0;i<=n;i++)
{
for(j=0;j<=n;j++)A[i][j]=C[i][j];//This is so that C[i][j] is returned unchanged
}

for(i=0;i<=n;i++)
{
for(j=0;j<=n;j++)B[i][j]=0.0;
B[i][i]=1.0;
}

for(i=0;i<=n;i++)
{
if(i!=n)
{
biggest=i;
for(k=i+1;k<=n;k++){if(fabs(A[k][i])>fabs(A[biggest][i]))biggest=k;}
for(j=0;j<=n;j++)
{
tempdouble=A[i][j];A[i][j]=A[biggest][j];A[biggest][j]=tempdouble;
tempdouble=B[i][j];B[i][j]=B[biggest][j];B[biggest][j]=tempdouble;
}
}

divide=A[i][i];
for(j=0;j<=n;j++)
{
A[i][j]=A[i][j]/divide;
B[i][j]=B[i][j]/divide;
}

if(i!=n)
{
   for(k=i+1;k<=n;k++)
   {
   mult=A[k][i];
   for(j=0;j<=n;j++)
   {
   A[k][j]-=mult*A[i][j];
   B[k][j]-=mult*B[i][j];
   }
   }
}
}

for(i=n;i>=1;i-=1)
{
   for(k=i-1;k>=0;k-=1)
   {
   mult=A[k][i];
   for(j=0;j<=n;j++)
   {
   A[k][j]-=mult*A[i][j];
   B[k][j]-=mult*B[i][j];
   }
   }
}
}