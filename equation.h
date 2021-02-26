#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "RightPart.h"
#include "json.hpp"

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
    int method;
    double last;
    double first;
    double eps = 1e-7;
    int equationsCount;
    vector<double> tOutPoints;
    double *initialConditions;
    vector<MyVector> xOutMatrix;
    void (*f)(double, double *, double *);

    double countApproximateRatio(double tau, string fname);
    int solveWithRungeKutta2Order(int pointsNumber = 10000);
    double calculateErrorNorm(double *y1, double *y2, int p);
    void assignInitialValues(double *&varX1, double *&varX2);
    int solveWithAdams(int methodOrder = 4, int numberOfPoints = 1000);
    int solveWithPredictorCorrector(int methodOrder = 4, int numberOfPoints = 1000);
    int solveWithRungeKutta(bool flag = false, int count = 4, int pointsNumber = 10000);
    void stepWithRungeKutta2Order(double *k1, double *k2, double *varX, double *tmpX, double tau, double &t);
    void deleteMemory(double *&k1, double *&k2, double *&k3, double *&k4, double *&varX1, double *&varX2, double *&tmpX);
    void stepWithRungeKutta(double *k1, double *k2, double *k3, double *k4, double *varX, double *tmpX, double tau, double &t);
    void initializeMemory(double *&k1, double *&k2, double *&k3, double *&k4, double *&varX1, double *&varX2, double *&tmpX,
                          double *&saveX1);

public:
    int solve();
    double getEps();
    void outputFile();
    string getMethod();
    void printResult();
    void setEps(double userEps);
    DifferentialEquationSolver();
    ~DifferentialEquationSolver();
    void setMethod(string userMethod = "RK4");
    void setRightPart(RightPart rp, int systemNum);
    void compareApproximateRatioAndAccuracyRatio();
    DifferentialEquationSolver(RightPart rp, int systemNum);
};