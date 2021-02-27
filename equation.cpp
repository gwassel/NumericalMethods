#include "equation.h"
int DifferentialEquationSolver::solveWithRungeKutta4Order(bool flag, int count)
{
    double *k1;
    double *k2;
    double *k3;
    double *k4;
    double maxTau;

    double *tmpX;
    double *varX1;
    double *varX2;
    double *saveX1;

    initializeMemory(k1, k2, k3, k4, varX1, varX2, tmpX, saveX1);

    double T;
    double t = first;
    int readyPoints = 1;
    double tau = (last - first) / (numberOfPoints - 1);
    maxTau = tau;

    assignInitialValues(varX1, varX2);
    tauOnStep.push_back(tau);

    while (t < last)
    {
        T = t;
        for (int i = 0; i < equationsCount; i++)
            saveX1[i] = varX1[i];

        while (true)
        {
            for (int i = 0; i < equationsCount; i++)
                varX1[i] = varX2[i] = saveX1[i];

            t = T;

            stepWithRungeKutta(k1, k2, k3, k4, varX1, tmpX, tau, t);

            t = T;

            stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, tau / 2.0, t);
            stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, tau / 2.0, t);

            if (calculateErrorNorm(varX1, varX2, 4) < eps)
            {
                readyPoints++;
                tOutPoints.push_back(t);
                tauOnStep.push_back(tau);
                xOutMatrix.push_back(MyVector(equationsCount));
                for (int i = 0; i < equationsCount; i++)
                    xOutMatrix[xOutMatrix.size() - 1].array[i] = varX1[i];
                t = T + tau;

                if (readyPoints == count && flag == true)
                {
                    deleteMemory(k1, k2, k3, k4, varX1, varX2, tmpX);
                    return 0;
                }

                if (calculateErrorNorm(varX1, varX2, 4) < eps)
                {
                    tau *= 2;
                    if (tau > maxTau)
                        maxTau = tau;
                }
                break;
            }
            else
            {
                t = T;
                tau /= 2.0;
            }
            if (tau < 1e-20)
            {
                deleteMemory(k1, k2, k3, k4, varX1, varX2, tmpX);
                cout << "fatal error" << endl;
                return 4000;
            }
        };
    }
    deleteMemory(k1, k2, k3, k4, varX1, varX2, tmpX);
    outputFile();
    cout << "Max tau needed to achieve eps = " << eps << ", " << maxTau << endl;
    return 0;
}

void DifferentialEquationSolver::assignInitialValues(double *&varX1, double *&varX2)
{
    tauOnStep.clear();
    tOutPoints.clear();
    xOutMatrix.clear();
    tOutPoints.push_back(first);
    xOutMatrix.push_back(MyVector(equationsCount));

    for (int i = 0; i < equationsCount; i++)
    {
        xOutMatrix[xOutMatrix.size() - 1].array[i] = initialConditions[i];
        varX1[i] = varX2[i] = initialConditions[i];
    }
}

void DifferentialEquationSolver::deleteMemory(double *&k1, double *&k2, double *&k3, double *&k4, double *&varX1, double *&varX2, double *&tmpX)
{
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] varX1;
    delete[] varX2;
    delete[] tmpX;
}

void DifferentialEquationSolver::initializeMemory(double *&k1, double *&k2, double *&k3, double *&k4, double *&varX1,
                                                  double *&varX2, double *&tmpX, double *&saveX1)
{
    varX1 = new double[equationsCount];
    varX2 = new double[equationsCount];
    tmpX = new double[equationsCount];
    k1 = new double[equationsCount];
    k2 = new double[equationsCount];
    k3 = new double[equationsCount];
    k4 = new double[equationsCount];
    saveX1 = new double[equationsCount];
}

DifferentialEquationSolver::DifferentialEquationSolver(RightPart rp, int systemNum)
{
    setRightPart(rp, systemNum);
}

DifferentialEquationSolver::DifferentialEquationSolver()
{
    method = 0;
    eps = 0.0;
    first = 0.0;
    last = 0.0;
    equationsCount = 0;
    initialConditions = nullptr;
}

DifferentialEquationSolver::~DifferentialEquationSolver()
{
    xOutMatrix.clear();
    tOutPoints.clear();
}

void DifferentialEquationSolver::setRightPart(RightPart rp, int systemNum)
{
    f = rp.getF(systemNum, equationsCount, initialConditions, first, last);
}

void DifferentialEquationSolver::printResult()
{
    cout << "Array of independent variable values:" << endl;
    for (int j = 0; j < tOutPoints.size(); j++)
        cout << tOutPoints[j] << " ";
    cout << endl;

    cout << "Array of function values:" << endl;
    for (int i = 0; i < xOutMatrix.size(); i++)
    {
        for (int j = 0; j < equationsCount; j++)
            cout << xOutMatrix[i].array[j] << " ";
        cout << endl;
    }
    cout << endl;
}

double DifferentialEquationSolver::calculateErrorNorm(double *y1, double *y2, int p)
{
    int s, u;
    switch (p)
    {
    case 2:
        s = 3;
        break;
    case 4:
        s = 15;
        break;
    }

    double max = fabs(y2[0] - y1[0]);
    u = 0;
    for (int i = 1; i < equationsCount; i++)
    {
        if (max < fabs((y2[i] - y1[i])))
        {
            max = fabs((y2[i] - y1[i]));
            u = i;
        }
    }
    return max / s;
}

void DifferentialEquationSolver::stepWithRungeKutta(double *k1, double *k2, double *k3, double *k4, double *varX, double *tmpX,
                                                    double tau, double &t)
{
    double ONE_SIXTH = (1.0 / 6.0);
    double HALF = 0.5;
    //k1
    f(t, varX, k1);
    for (int j = 0; j < equationsCount; j++)
        tmpX[j] = varX[j] + tau * HALF * k1[j];

    //k2
    t += HALF * tau;
    f(t, tmpX, k2);
    for (int j = 0; j < equationsCount; j++)
        tmpX[j] = varX[j] + tau * HALF * k2[j];

    //k3
    f(t, tmpX, k3);
    for (int j = 0; j < equationsCount; j++)
        tmpX[j] = varX[j] + tau * k3[j];

    //k4
    t += HALF * tau;
    f(t, tmpX, k4);
    for (int j = 0; j < equationsCount; j++)
        varX[j] += ONE_SIXTH * tau * (k1[j] + k2[j] + k2[j] + k3[j] + k3[j] + k4[j]);
}

int DifferentialEquationSolver::solveWithAdams(int methodOrder)
{
    solveWithRungeKutta4Order(true);

    double t = tOutPoints[tOutPoints.size() - 1];
    double step = fabs(last - t) / numberOfPoints;
    double coef = step / 24.0;

    double C0 = 55.0;
    double C1 = -59.0;
    double C2 = 37.0;
    double C3 = -9.0;

    double **tmpf = new double *[methodOrder];
    for (int i = 0; i < methodOrder; i++)
        tmpf[i] = new double[equationsCount];

    for (int i = 0; i < numberOfPoints; i++)
    {
        t = t + step;
        tOutPoints.push_back(t);
        xOutMatrix.push_back(MyVector(equationsCount));

        for (int i1 = 1; i1 <= methodOrder; i1++)
            f(tOutPoints[tOutPoints.size() - 1 - i1] /*t - step * i1*/, xOutMatrix[xOutMatrix.size() - i1 - 1].array, tmpf[i1 - 1]);

        for (int j = 0; j < equationsCount; j++)
            xOutMatrix[xOutMatrix.size() - 1].array[j] = xOutMatrix[xOutMatrix.size() - 2].array[j] + coef * (C0 * tmpf[0][j] + C1 * tmpf[1][j] + C2 * tmpf[2][j] + C3 * tmpf[3][j]);
    };
    outputFile();
    return 0;
}

void DifferentialEquationSolver::outputFile()
{
    cout << "Outputting file" << endl;
    ofstream fOut;
    string fileName;

    for (int i = 0; i < equationsCount; i++)
    {
        fileName = "res" + to_string(i) + ".txt";
        fOut.open(fileName, ios::out);

        for (int j = 0; j < xOutMatrix.size(); j++)
            fOut << xOutMatrix[j].array[i] << endl;
        fOut.close();
    }

    fileName = "t.txt";
    fOut.open(fileName, ios::out);

    for (int j = 0; j < tOutPoints.size(); j++)
        fOut << tOutPoints[j] << endl;

    fOut.close();

    fileName = "tau.txt";
    fOut.open(fileName, ios::out);

    for (int j = 0; j < tauOnStep.size(); j++)
        fOut << tauOnStep[j] << endl;

    fOut.close();
}

int DifferentialEquationSolver::solveWithPredictorCorrector(int methodOrder)
{
    solveWithRungeKutta4Order(true);

    double t = tOutPoints[tOutPoints.size() - 1];
    double step = fabs(last - t) / numberOfPoints;

    double C0 = step / 24.0;
    double C1 = 55.0;
    double C2 = -59.0;
    double C3 = 37.0;
    double C4 = -9.0;

    double Q1 = 9.0;
    double Q2 = 19.0;
    double Q3 = -5.0;

    double *tmpX = new double[equationsCount];

    double **tmpf = new double *[methodOrder];
    for (int i = 0; i < methodOrder; i++)
        tmpf[i] = new double[equationsCount];

    for (int i = 0; i < numberOfPoints; i++)
    {
        t += step;
        tOutPoints.push_back(t);
        xOutMatrix.push_back(MyVector(equationsCount));

        for (int i1 = 1; i1 <= methodOrder; i1++)
            f(tOutPoints[tOutPoints.size() - 1 - i1] /*t - step * i1*/, xOutMatrix[xOutMatrix.size() - i1 - 1].array, tmpf[i1 - 1]);

        for (int j = 0; j < equationsCount; j++)
            xOutMatrix[xOutMatrix.size() - 1].array[j] = xOutMatrix[xOutMatrix.size() - 2].array[j] + C0 * (C1 * tmpf[0][j] + C2 * tmpf[1][j] + C3 * tmpf[2][j] + C4 * tmpf[3][j]);

        f(tOutPoints[tOutPoints.size() - 1], xOutMatrix[xOutMatrix.size() - 1].array, tmpX);

        for (int j = 0; j < equationsCount; j++)
            xOutMatrix[xOutMatrix.size() - 1].array[j] = xOutMatrix[xOutMatrix.size() - 2].array[j] + C0 * (Q1 * tmpX[j] + Q2 * tmpf[0][j] + Q3 * tmpf[1][j] + tmpf[2][j]);
    };
    outputFile();

    return 0;
}

int DifferentialEquationSolver::solveWithRungeKutta2Order()
{
    double *varX1 = new double[equationsCount];
    double *saveX1 = new double[equationsCount];
    double *varX2 = new double[equationsCount];
    double *tmpX = new double[equationsCount];

    double *k1 = new double[equationsCount];
    double *k2 = new double[equationsCount];

    double T;
    double maxTau;
    double t = first;
    double tau = (last - first) / (numberOfPoints - 1);
    maxTau = tau;

    tOutPoints.clear();
    xOutMatrix.clear();
    tOutPoints.push_back(first);
    xOutMatrix.push_back(MyVector(equationsCount));

    for (int i = 0; i < equationsCount; i++)
    {
        xOutMatrix[xOutMatrix.size() - 1].array[i] = initialConditions[i];
        varX1[i] = varX2[i] = initialConditions[i];
    }

    while (t < last)
    {
        T = t;
        for (int i = 0; i < equationsCount; i++)
            saveX1[i] = varX1[i];
        while (true)
        {
            for (int i = 0; i < equationsCount; i++)
            {
                varX1[i] = saveX1[i];
                varX2[i] = varX1[i];
            }

            t = T;

            stepWithRungeKutta2Order(k1, k2, varX1, tmpX, tau, t);

            t = T;

            stepWithRungeKutta2Order(k1, k2, varX2, tmpX, tau / 2.0, t);
            stepWithRungeKutta2Order(k1, k2, varX2, tmpX, tau / 2.0, t);

            if (calculateErrorNorm(varX1, varX2, 2) < eps)
            {
                tOutPoints.push_back(t);
                xOutMatrix.push_back(MyVector(equationsCount));
                for (int i = 0; i < equationsCount; i++)
                    xOutMatrix[xOutMatrix.size() - 1].array[i] = varX1[i];
                t = T + tau;

                if (calculateErrorNorm(varX1, varX2, 2) < eps)
                {
                    tau *= 2;
                    if (maxTau < tau)
                        maxTau = tau;
                }

                break;
            }
            else
            {
                t = T;
                tau /= 2.0;
            }
            if (tau < 1e-20)
            {
                cout << "fatal error" << endl;
                return 4000;
            }
        };
    }
    outputFile();
    cout << "Max tau needed to achieve eps = " << eps << ", " << maxTau << endl;
    return 0;
}

void DifferentialEquationSolver::stepWithRungeKutta2Order(double *k1, double *k2, double *varX, double *tmpX, double tau, double &t)
{
    //k1
    f(t, varX, k1);
    for (int j = 0; j < equationsCount; j++)
        tmpX[j] = varX[j] + k1[j] * tau * 0.5;

    //k2
    t += 0.5 * tau;
    f(t, tmpX, k2);

    for (int j = 0; j < equationsCount; j++)
        varX[j] += tau * k2[j];
}

double DifferentialEquationSolver::countApproximateRatio(double mtau, string fname)
{
    double *k1;
    double *k2;
    double *k3;
    double *k4;

    double *tmpX;
    double *varX1;
    double *varX2;
    double *saveX1; //not needed

    initializeMemory(k1, k2, k3, k4, varX1, varX2, tmpX, saveX1);

    double T;
    double maxErr;
    double t = first;

    assignInitialValues(varX1, varX2);

    //one init step
    T = t;
    for (int i = 0; i < equationsCount; i++)
        varX2[i] = varX1[i];

    t = T;

    stepWithRungeKutta(k1, k2, k3, k4, varX1, tmpX, mtau, t);

    t = T;

    stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, mtau / 2.0, t);
    stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, mtau / 2.0, t);

    t = T + mtau;
    maxErr = calculateErrorNorm(varX1, varX2, 4);

    //all other steps
    while (t < last)
    {
        T = t;
        for (int i = 0; i < equationsCount; i++)
            varX2[i] = varX1[i];

        t = T;

        stepWithRungeKutta(k1, k2, k3, k4, varX1, tmpX, mtau, t);

        t = T;

        stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, mtau / 2.0, t);
        stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, mtau / 2.0, t);

        if (maxErr < calculateErrorNorm(varX1, varX2, 4))
            maxErr = calculateErrorNorm(varX1, varX2, 4);
        t = T + mtau;
    }
    deleteMemory(k1, k2, k3, k4, varX1, varX2, tmpX);

    return maxErr;
}

void DifferentialEquationSolver::setEps(double userEps)
{
    eps = userEps;
}

double DifferentialEquationSolver::getEps()
{
    return eps;
}

void DifferentialEquationSolver::setMethod(string userMethod)
{
    if (userMethod == "RK2")
        method = 1;
    if (userMethod == "RK4")
        method = 2;
    if (userMethod == "A")
        method = 3;
    if (userMethod == "PC")
        method = 4;
    if (userMethod == "EE")
        method = 5;
    if (userMethod == "EB")
        method = 6;
}

void DifferentialEquationSolver::setNumberOfPoints(int userNumberOfPoints)
{
    numberOfPoints = userNumberOfPoints;
}

string DifferentialEquationSolver::getMethod()
{
    switch (method)
    {
    case 1:
        return "Runge Kutta 2 order";
        break;
    case 2:
        return "Runge Kutta 4 order";
        break;
    case 3:
        return "Adams";
        break;
    case 4:
        return "Predictor-corrector";
        break;
    case 5:
        return "Euler explicit";
        break;
    case 6:
        return "Euler backword";
        break;
    }
}

int DifferentialEquationSolver::solve()
{
    switch (method)
    {
    case 1:
        cout << "RK2" << endl;
        solveWithRungeKutta2Order();
        break;
    case 2:
        cout << "RK4" << endl;
        solveWithRungeKutta4Order();
        break;
    case 3:
        cout << "A" << endl;
        solveWithAdams();
        break;
    case 4:
        cout << "PD" << endl;
        solveWithPredictorCorrector();
        break;
    case 5:
        cout << "EE" << endl;
        solveWithExplicitEuler();
        break;
    case 6:
        //backword euler
        break;
    }
    return 0;
}

void DifferentialEquationSolver::compareApproximateRatioAndAccuracyRatio()
{
    double err1 = countApproximateRatio(0.04, "1");
    double err2 = countApproximateRatio(0.02, "2");
    double err3 = countApproximateRatio(0.01, "3");
    double err4 = countApproximateRatio(0.005, "4");
    double err5 = countApproximateRatio(0.0025, "5");
    cout << "Showing maximum errors in system with different tau" << endl;
    cout << err1 << " " << err2 << " " << err3 << " " << err4 << " " << err5 << endl;
    cout << "Showing approximate ratio" << endl;
    cout << "P = " << (1.0 / log(0.5)) * log((err3 - err2) / (err2 - err1)) << endl;
    cout << "P = " << (1.0 / log(0.5)) * log((err4 - err3) / (err3 - err2)) << endl;
    cout << "P = " << (1.0 / log(0.5)) * log((err5 - err4) / (err4 - err3)) << endl;
    cout << "Showing accuracy ratio" << endl;
    cout << "P_N = " << log(err1 / err2) / log(2) << endl;
    cout << "P_N = " << log(err2 / err3) / log(2) << endl;
    cout << "P_N = " << log(err3 / err4) / log(2) << endl;
    cout << "P_N = " << log(err4 / err5) / log(2) << endl;
}

void DifferentialEquationSolver::solveWithExplicitEuler()
{
    double t = first;
    double tau = (last - first) / (numberOfPoints - 1);

    xOutMatrix.clear();
    tOutPoints.clear();

    double *varX = new double[equationsCount];
    double *tmpX = new double[equationsCount];

    for (int i = 0; i < equationsCount; i++)
        varX[i] = initialConditions[i];

    xOutMatrix.push_back(MyVector(equationsCount));
    tOutPoints.push_back(t);
    for (int i = 0; i < equationsCount; i++)
        xOutMatrix[xOutMatrix.size() - 1] = varX[i] + tau * tmpX[i];

    while (t < last)
    {
        xOutMatrix.push_back(MyVector(equationsCount));
        tOutPoints.push_back(t);

        f(t, varX, tmpX);

        for (int i = 0; i < equationsCount; i++)
            xOutMatrix[xOutMatrix.size() - 1].array[i] = varX[i] + tau * tmpX[i];

        for (int i = 0; i < equationsCount; i++)
            varX[i] = xOutMatrix[xOutMatrix.size() - 1].array[i];

        t += tau;
    }

    outputFile();
}