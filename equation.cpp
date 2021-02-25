#include "equation.h"
int DifferentialEquationSolver::solveWithRungeKutta(bool flag, int count)
{
    double *k1;
    double *k2;
    double *k3;
    double *k4;

    double *tmpX;
    double *varX1;
    double *saveX1 = new double[equationsCount];
    double *varX2;

    initializeMemory(k1, k2, k3, k4, varX1, varX2, tmpX);

    double T;
    double norm;
    double badTau;
    double t = first;
    int readyPoints = 1;
    double tau = tauStart;

    assignInitialValues(varX1, varX2);

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

            stepWithRungeKutta(k1, k2, k3, k4, varX1, tmpX, tau, t);

            t = T;

            stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, tau / 2.0, t);
            stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, tau / 2.0, t);

            //cout << calculateErrorNorm(varX1, varX2, 4) << " " << tau << endl;
            //cout << varX1[0] << " " << varX1[1] << " " << varX1[2] << " " << varX1[3] << endl;
            //cout << varX2[0] << " " << varX2[1] << " " << varX2[2] << " " << varX2[3] << endl;
            //cin.get();
            //cout << tau<< endl;

            if (calculateErrorNorm(varX1, varX2, 4) < eps)
            {
                readyPoints++;
                tOutPoints.push_back(t);
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
                    tau *= 2;

                break;
            }
            else
            {
                t = T;
                tau /= 2.0;
                //break;
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
    return 0;
}

void DifferentialEquationSolver::assignInitialValues(double *&varX1, double *&varX2)
{
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

void DifferentialEquationSolver::initializeMemory(double *&k1, double *&k2, double *&k3, double *&k4, double *&varX1, double *&varX2, double *&tmpX)
{
    varX1 = new double[equationsCount];
    varX2 = new double[equationsCount];
    tmpX = new double[equationsCount];
    k1 = new double[equationsCount];
    k2 = new double[equationsCount];
    k3 = new double[equationsCount];
    k4 = new double[equationsCount];
}

DifferentialEquationSolver::DifferentialEquationSolver(RightPart rp, int systemNum)
{
    f = rp.getF(systemNum, equationsCount, initialConditions, first, last);
    tauStart = (last - first) / (pointsNumber - 1);
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
    /*double sum = 0.0;
    for (int i = 0; i < equationsCount; i++)
        sum += (y2[i] - y1[i]) * ((y2[i] - y1[i]) / (225.0));
    return sqrt(sum);*/

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
    //cout << "max:" << max << endl;
    u = 0;
    for (int i = 1; i < equationsCount; i++)
    {
        //cout << fabs((y2[i] - y1[i])) << endl;
        if (max < fabs((y2[i] - y1[i])))
        {
            max = fabs((y2[i] - y1[i]));
            u = i;
        }
    }
    return max / s;

    /*double sum = 0.0;
    for (int i = 0; i < equationsCount; i++)
        sum += fabs(y2[i] - y1[i]) / s;
    return sum;*/
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

int DifferentialEquationSolver::solveWithAdams(int methodOrder, int numberOfPointsAdams)
{
    solveWithRungeKutta(true);

    double t = tOutPoints[tOutPoints.size() - 1];
    double step = fabs(last - t) / numberOfPointsAdams;
    double coef = step / 24.0;

    double C0 = 55.0;
    double C1 = -59.0;
    double C2 = 37.0;
    double C3 = -9.0;

    double **tmpf = new double *[methodOrder];
    for (int i = 0; i < methodOrder; i++)
        tmpf[i] = new double[equationsCount];

    for (int i = 0; i < numberOfPointsAdams; i++)
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

    for (int j = 0; j < xOutMatrix.size(); j++)
        fOut << tOutPoints[j] << endl;

    fOut.close();
}

int DifferentialEquationSolver::solveWithPredictorCorrector(int methodOrder, int numberOfPoints)
{
    solveWithRungeKutta(true);

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
    double t = first;
    double tau = tauStart;

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

            //cout << calculateErrorNorm(varX1, varX2, 2) << " " << t << " " << tau << endl;
            //cout << varX1[0] << " " << varX1[1] << varX1[2] << " " << varX1[3] << endl;
            //cout << varX2[0] << " " << varX2[1] << varX2[2] << " " << varX2[3] << endl;
            //cin.get();

            if (calculateErrorNorm(varX1, varX2, 2) < eps)
            {
                tOutPoints.push_back(t);
                xOutMatrix.push_back(MyVector(equationsCount));
                for (int i = 0; i < equationsCount; i++)
                    xOutMatrix[xOutMatrix.size() - 1].array[i] = varX1[i];
                t = T + tau;

                if (calculateErrorNorm(varX1, varX2, 2) < eps)
                    tau *= 2;

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
