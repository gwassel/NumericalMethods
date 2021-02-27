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
    double eps;
    int equationsCount;
    int numberOfPoints;
    vector<double> tauOnStep;
    double *initialConditions;
    vector<double> tOutPoints;
    vector<MyVector> xOutMatrix;
    void (*f)(double, double *, double *);

    void solveWithExplicitEuler();
    int solveWithRungeKutta2Order();
    int solveWithAdams(int methodOrder = 4);
    int solveWithPredictorCorrector(int methodOrder = 4);
    double countApproximateRatio(double tau, string fname);
    double calculateErrorNorm(double *y1, double *y2, int p);
    void assignInitialValues(double *&varX1, double *&varX2);
    int solveWithRungeKutta4Order(bool flag = false, int count = 4);
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
    DifferentialEquationSolver();
    ~DifferentialEquationSolver();
    void setEps(double userEps = 1e-7);
    void setMethod(string userMethod = "RK4");
    void setRightPart(RightPart rp, int systemNum);
    void compareApproximateRatioAndAccuracyRatio();
    void setNumberOfPoints(int userNumberOfPoints = 1000);
    DifferentialEquationSolver(RightPart rp, int systemNum);
};