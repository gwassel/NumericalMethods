#include "equation.cpp"

int main()
{
     RightPart rp;
     int systemNum = 2;
     DifferentialEquationSolver eq1(rp, systemNum);
     //eq1.solveWithRungeKutta();
     //eq1.solveWithAdams();
     //eq1.solveWithPredictorCorrector();
     eq1.solveWithRungeKutta2Order();
     return 0;
}