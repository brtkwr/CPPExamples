#include "ChrisNet.h"
int main(void)
{
ReadNet();
ofstream Madeleine("Net.dxf");
Madeleine<<"0\nSECTION\n2\nENTITIES\n";
for(element=0;element<=lastelement;element++)
{
Madeleine<<"0\n3DLINE\n8\nNet\n";
Madeleine<<"10\n"<<x[end[element][0]][0]<<"\n";
Madeleine<<"20\n"<<x[end[element][0]][1]<<"\n";
Madeleine<<"30\n"<<x[end[element][0]][2]<<"\n";
Madeleine<<"11\n"<<x[end[element][1]][0]<<"\n";
Madeleine<<"21\n"<<x[end[element][1]][1]<<"\n";
Madeleine<<"31\n"<<x[end[element][1]][2]<<"\n";
}
Madeleine<<"0\nENDSEC\n0\nEOF\n";
Madeleine.close();cout<<"Finished\n";return 0;
}
