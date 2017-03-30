#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
using namespace std;
#define W 3600
#define H 3000
short int point[W][H];
int iteration,a,b,i;
double r,x;char s;
ofstream tga;
int main(void){
tga.open("LogisticMap.tga",
ofstream::binary);
s=0;tga.put(s);s=0;tga.put(s);
s=2;tga.put(s);
i=0;tga.put(i%256);tga.put(i/256);
i=0;tga.put(i%256);tga.put(i/256);
s=0;tga.put(s);
i=0;tga.put(i%256);tga.put(i/256);
i=0;tga.put(i%256);tga.put(i/256);
tga.put(W%256);tga.put(W/256);
tga.put(H%256);tga.put(H/256);
s=24;tga.put(s);s=0x20;tga.put(s);
for(a=0;a<=W-1;a++){
if(100*int(a/100)==a)cout<<a<<"\n";
r=2.9+(1.1*a)/(1.0*(W-1));
for(b=0;b<=H-1;b++)point[a][b]=0;
x=0.5;iteration=0;b=0;
for(;;){
iteration++;x=r*x*(1.0-x);
if(iteration>100000){
b=int(x*(H-1)+0.5);
if(b<0)b=0;if(b>H-1)b=H-1;
point[a][b]++;}
if(point[a][b]>=2000)break;}}
for(b=H-1;b>=0;b-=1){
for(a=0;a<=W-1;a++){
i=int(255.0*(1.0-
tanh(1.0*point[a][b]/255.0)));
tga.put(i);tga.put(i);tga.put(i);}}
tga.close();cout<<"File written\n";
return 0;}