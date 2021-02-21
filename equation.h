#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "RightPart.h"

class vec{
public:
double* p;
vec(int size){
 p = new double[size];
}
};

using namespace std;
class DifferentialEquation
{
private:
    void (*f)(double, double *, double *);
    double *initialConditions;
    int equationsCount;
    double first;
    double last;
    int pointsNumber = 10;
    double tauStart;
    vector<double> tOutPoints;
    vector<vec> xOutMatrix;
    double eps = 1e-6;

    double calculateError(double *y1, double* y2);
    void stepWithRungeKutta(double* k1,double* k2,double* k3,double* k4,double *varX, double *tmpX, double tau, double& t);

public:
    DifferentialEquation(RightPart rp, int systemNum);
    int solveWithRungeKutta(bool flag = false, int count = 4);
    int solveWithAdams(int methodOrder = 4, int numberOfPoints = 1000);
    void outputFile();
    void printResult();
};