#include "size.h"

void Mult(int n,double A[][nplus1],double B[][nplus1],double C[][nplus1])
{
int i,j,k;

for(i=0;i<=n;i++)
{
for(j=0;j<=n;j++)
{
C[i][j]=0.0;
for(k=0;k<=n;k++)C[i][j]+=A[i][k]*B[k][j];
}
}
}