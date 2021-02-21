#include <iostream>
#include "equation.cpp"

int main()
{
     RightPart rp;
     int systemNum = 2;
     DifferentialEquation eq1(rp, systemNum);
     //eq1.solveWithRungeKutta();
     eq1.solveWithPredictorCorrector();
     
     cout << "end" << endl;
     return 0;
}