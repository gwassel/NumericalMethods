#include "equation.cpp"

void readConfig(DifferentialEquationSolver &solver, RightPart &rp)
{
     const std::string FILE_NAME_CFG = "configs.json";
     std::ifstream fileJsonInput;
     fileJsonInput.open(FILE_NAME_CFG);

     nlohmann::json objJson;
     fileJsonInput >> objJson;

     fileJsonInput.close();

     solver.setEps(objJson["eps"]);
     solver.setMethod(objJson["method"]);
     solver.setRightPart(rp, objJson["system"]);
     solver.setNumberOfPoints(objJson["points_number"]);
}

int main()
{
     RightPart rp;
     DifferentialEquationSolver solver;
     readConfig(solver, rp);
     solver.solve();
     //solver.compareApproximateRatioAndAccuracyRatio();
     return 0;
}
