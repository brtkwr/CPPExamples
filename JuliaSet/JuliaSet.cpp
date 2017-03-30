#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
using namespace std;
int W,H,i,j,Colour,Number,MaxNumber;
double cReal,cImag,x,y,xSq,ySq;char s;
ofstream tga;
int main(void){
cout<<"Real and imaginary parts of c (for example: -0.16 -0.651)\n";
cin>>cReal>>cImag;
cout<<"Height of image in pixels\n";cin>>H;W=int((4.0*H)/3.0);
tga.open("JuliaSet.tga",ofstream::binary);
s=0;tga.put(s);s=0;tga.put(s);s=2;tga.put(s);
i=0;tga.put(i%256);tga.put(i/256);
i=0;tga.put(i%256);tga.put(i/256);
s=0;tga.put(s);
i=0;tga.put(i%256);tga.put(i/256);
i=0;tga.put(i%256);tga.put(i/256);
tga.put(W%256);tga.put(W/256);
tga.put(H%256);tga.put(H/256);
s=24;tga.put(s);s=0x20;tga.put(s);
MaxNumber=0;
for(j=H-1;j>=0;j-=1){if(100*int(j/100)==j)cout<<j<<"\n";
for(i=0;i<=W-1;i++){Number=0;
x=2.0*(2.0*i-1.0*(W-1))/(1.0*(W-1));
y=2.0*(2.0*j-1.0*(H-1))/(1.0*(W-1));
do{xSq=x*x;ySq=y*y;y=2.0*x*y+cImag;x=xSq-ySq+cReal;Number++;
if(MaxNumber<Number)MaxNumber=Number;}while(xSq+ySq<4.0&&Number<10000);
Colour=int(255/cosh((1.0*Number)/500.0));tga.put(Colour);
Colour=int(255/cosh((1.0*Number)/250.0));tga.put(Colour);tga.put(Colour);}}
tga.close();cout<<"Maximum number = "<<MaxNumber<<"\nFile written\n";
return 0;}