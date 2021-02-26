#include "equation.cpp"

int main()
{
     RightPart rp;
     int systemNum = 5;
     DifferentialEquationSolver eq1(rp, systemNum);
     //eq1.solveWithRungeKutta();
     //eq1.solveWithAdams();
     //eq1.solveWithPredictorCorrector();
     //eq1.solveWithRungeKutta2Order();pantheon dolhin dock

     //counting errors(extra)
     
     double err1 = eq1.countApproximateRatio(0.04, "1");
     double err2 = eq1.countApproximateRatio(0.02, "2");
     double err3 = eq1.countApproximateRatio(0.01, "3");
     double err4 = eq1.countApproximateRatio(0.005, "4");
     double err5 = eq1.countApproximateRatio(0.0025, "5");
     cout << err1 << " " << err2 << " " << err3 << " " << err4 << " " << err5 << endl;
     cout << "P = " << (1.0 / log(0.5)) * log((err3 - err2) / (err2 - err1)) << endl;
     cout << "P = " << (1.0 / log(0.5)) * log((err4 - err3) / (err3 - err2)) << endl;
     cout << "P = " << (1.0 / log(0.5)) * log((err5 - err4) / (err4 - err3)) << endl;
     cout << "P_N = " << log(err1 / err2) / log(2) << endl;
     cout << "P_N = " << log(err2 / err3) / log(2) << endl;
     cout << "P_N = " << log(err3 / err4) / log(2) << endl;
     cout << "P_N = " << log(err4 / err5) / log(2) << endl;
     return 0;
}