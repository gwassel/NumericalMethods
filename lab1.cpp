#include "equation.cpp"

int main()
{
     RightPart rp;
     int systemNum = 3;
     DifferentialEquation eq1(rp, systemNum);
     eq1.solveWithRungeKutta();
     //eq1.solveWithAdams();
     //eq1.solveWithPredictorCorrector();
     //eq1.solveWithRungeKutta2Order();

     cout << "end" << endl;
     return 0;
}