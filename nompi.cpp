
#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;
int main (int argc, char *argv[])
{
  ofstream fout;
  fout.open("1.txt", ios_base::app);
  fout<<argc<<endl;
  for (int i=0;i<argc;i++)
  fout<<argv[i]<<endl; 
  fout<<"___";
  fout.close();
  return 0;
}