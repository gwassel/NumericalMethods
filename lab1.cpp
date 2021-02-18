#include <iostream>
#include "equation.cpp"

int main(){
RightPart rp;
int systemNum = 2;
DifferentialEquation eq1(rp, systemNum);
eq1.solveWithRungeKutta();

system("pause");
return 0;
}