#include <iostream>
#include <fstream>
#include "RightPart.h"

using namespace std;
class DifferentialEquation{
private:
void (*f)(double, double*, double*);
double* initialConditions;
int equationsCount;
double first;
double last;
int pointsNumber = 11;
double* tOutPoints;
double** xOutMatrix;

public:
DifferentialEquation(RightPart rp, int systemNum);
int solveWithRungeKutta();
void outputFile();
};