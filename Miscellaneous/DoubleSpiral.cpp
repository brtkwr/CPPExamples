#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
#define m 2000
#define n 25
int i,j;
double x[m+1],y[m+1],PI,a,phi,psi,
theta,costheta,sintheta,alpha,beta,B;
ofstream Julia("Spiral.dxf");
int main(void){
PI=4.0*atan(1.0);a=1000.0;
theta=60.0*PI/180.0;
costheta=cos(theta);
sintheta=sin(theta);
Julia<<"0\nSECTION\n2\nENTITIES\n";
for(j=0;j<=n;j++){cout<<n-j<<"\n";
beta=(PI*(2*j-n)*costheta)/(1.0*(n+1));
for(i=0;i<=m;i++){
alpha=30.0*(1.0*i-0.5*m)/(1.0*m);
phi=alpha*costheta-beta*sintheta;
psi=alpha*sintheta+beta*costheta;
B=cosh(phi)+cos(psi);
x[i]=a*sinh(phi)/B;y[i]=a*sin(psi)/B;}
for(i=0;i<=m-1;i++){
Julia<<"0\nLINE\n8\nSpirals\n";
Julia<<"10\n"<<x[i]<<"\n";
Julia<<"20\n"<<y[i]<<"\n";
Julia<<"11\n"<<x[i+1]<<"\n";
Julia<<"21\n"<<y[i+1]<<"\n";}}
Julia<<"0\nENDSEC\n0\nEOF\n";
Julia.close();
cout<<"DXF file written.\n";
return 0;}