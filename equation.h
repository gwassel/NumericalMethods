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
    double eps = 1e-7;
    int equationsCount;
    vector<double> tOutPoints;
    double *initialConditions;
    vector<MyVector> xOutMatrix;
    void (*f)(double, double *, double *);

    double calculateErrorNorm(double *y1, double *y2, int p);
    void assignInitialValues(double *&varX1, double *&varX2);
    void stepWithRungeKutta2Order(double *k1, double *k2, double *varX, double *tmpX, double tau, double &t);
    void deleteMemory(double *&k1, double *&k2, double *&k3, double *&k4, double *&varX1, double *&varX2, double *&tmpX);
    void initializeMemory(double *&k1, double *&k2, double *&k3, double *&k4, double *&varX1, double *&varX2, double *&tmpX);
    void stepWithRungeKutta(double *k1, double *k2, double *k3, double *k4, double *varX, double *tmpX, double tau, double &t);

public:
    void outputFile();
    void printResult();
    double countApproximateRatio(double tau, string fname);
    DifferentialEquationSolver(RightPart rp, int systemNum);
    int solveWithRungeKutta2Order(int pointsNumber = 10000);
    int solveWithAdams(int methodOrder = 4, int numberOfPoints = 1000);
    int solveWithPredictorCorrector(int methodOrder = 4, int numberOfPoints = 1000);
    int solveWithRungeKutta(bool flag = false, int count = 4, int pointsNumber = 10000);
};