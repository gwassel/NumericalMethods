#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "RightPart.h"

class MyVector
{
    public:
    double *array;
    MyVector(int size)
    {
        array = new double[size];
    }
};

using namespace std;
class DifferentialEquationSolver
{
private:
    double last;
    double first;
    double tauStart;
    double eps = 1e-6;
    int equationsCount;
    int pointsNumber = 1000;
    vector<double> tOutPoints;
    double *initialConditions;
    vector<MyVector> xOutMatrix;
    void (*f)(double, double *, double *);

    double calculateError(double *y1, double *y2);
    void assignInitialValues(double*& varX1,double*& varX2);
    void stepWithRungeKutta2Order(double *k1, double *k2, double *varX, double *tmpX, double tau, double &t);
    void deleteMemory(double*& k1, double*& k2, double*& k3, double*& k4, double*& varX1, double*& varX2, double*& tmpX);
    void initializeMemory(double*& k1, double*& k2, double*& k3, double*& k4, double*& varX1, double*& varX2, double*& tmpX);
    void stepWithRungeKutta(double *k1, double *k2, double *k3, double *k4, double *varX, double *tmpX, double tau, double &t);

public:
    DifferentialEquationSolver(RightPart rp, int systemNum);
    int solveWithRungeKutta(bool flag = false, int count = 4);
    int solveWithRungeKutta2Order();
    int solveWithAdams(int methodOrder = 4, int numberOfPoints = 1000);
    int solveWithPredictorCorrector(int methodOrder = 4, int numberOfPoints = 1000);
    void outputFile();
    void printResult();
};