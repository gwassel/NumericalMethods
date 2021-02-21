#include "equation.cpp"

int main()
{
     RightPart rp;
     int systemNum = 5;
     DifferentialEquation eq1(rp, systemNum);
     //eq1.solveWithRungeKutta();
     //eq1.solveWithAdams();
     eq1.solveWithPredictorCorrector();
     
     cout << "end" << endl;
     return 0;
}