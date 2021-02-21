#include "equation.h"
int DifferentialEquation::solveWithRungeKutta(bool flag, int count)
{
    //delete memory
    double *varX1 = new double[equationsCount];
    double *varX2 = new double[equationsCount];
    double *tmpX = new double[equationsCount];

    double *k1 = new double[equationsCount];
    double *k2 = new double[equationsCount];
    double *k3 = new double[equationsCount];
    double *k4 = new double[equationsCount];

    double T;
    double norm;
    double t = first;
    double tau = tauStart;
    int readyPoints = 1;

    tOutPoints.clear();
    xOutMatrix.clear();
    tOutPoints.push_back(first);
    xOutMatrix.push_back(vec(equationsCount));

    for (int i = 0; i < equationsCount; i++)
    {
        xOutMatrix[xOutMatrix.size() - 1].p[i] = initialConditions[i];
        varX1[i] = varX2[i] = initialConditions[i];
    }

    while (t < last)
    {
        T = t;
        norm = -1.0;
        while (true)
        {
            t = T;

            stepWithRungeKutta(k1, k2, k3, k4, varX1, tmpX, tau, t);

            t = T;

            stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, tau / 2.0, t);
            stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, tau / 2.0, t);

            if (calculateError(varX1, varX2) < eps || fabs(calculateError(varX1, varX2) - norm) < 1e-8)
            {
                readyPoints++;
                tOutPoints.push_back(t);
                xOutMatrix.push_back(vec(equationsCount));
                for (int i = 0; i < equationsCount; i++)
                    xOutMatrix[xOutMatrix.size() - 1].p[i] = varX1[i];
                t = T + tau;

                if (readyPoints == count && flag == true)
                {
                    stepAfterRK = tau;
                    return 0;
                }

                if (calculateError(varX1, varX2) < eps / 100.0)
                    tau *= 2;

                break;
            }
            else
            {
                readyPoints = 1;
                xOutMatrix.clear();
                tOutPoints.clear();
                tOutPoints.push_back(first);
                xOutMatrix.push_back(vec(equationsCount));

                for (int i = 0; i < equationsCount; i++)
                {
                    xOutMatrix[xOutMatrix.size() - 1].p[i] = initialConditions[i];
                    varX1[i] = varX2[i] = initialConditions[i];
                }

                tauStart /= 10.0;
                tau = tauStart;
                t = first;
                break;
            }
            norm = calculateError(varX1, varX2);
        };

        if (tauStart < 1e-20)
        {
            return 4000;
        }
    }
    outputFile();
    return 0;
}

DifferentialEquation::DifferentialEquation(RightPart rp, int systemNum)
{
    f = rp.getF(systemNum, equationsCount, initialConditions, first, last);
    tauStart = (last - first) / (pointsNumber - 1);
}

void DifferentialEquation::printResult()
{
    cout << "T:" << endl;
    for (int j = 0; j < tOutPoints.size(); j++)
        cout << tOutPoints[j] << " ";
    cout << endl;

    cout << "X:" << endl;
    for (int i = 0; i < xOutMatrix.size(); i++)
    {
        for (int j = 0; j < equationsCount; j++)
        {
            cout << xOutMatrix[i].p[j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

double DifferentialEquation::calculateError(double *y1, double *y2)
{
    double sum = 0.0;
    for (int i = 0; i < equationsCount; i++)
        sum += (y2[i] - y1[i]) * ((y2[i] - y1[i]) / (225.0));
    return sqrt(sum);
}

void DifferentialEquation::stepWithRungeKutta(double *k1, double *k2, double *k3, double *k4, double *varX, double *tmpX,
                                              double tau, double &t)
{
    double oneSixth = (1.0 / 6.0);
    //k1
    f(t, varX, k1);
    for (int j = 0; j < equationsCount; j++)
        tmpX[j] = varX[j] + k1[j] * tau * 0.5;

    //k2
    t += 0.5 * tau;
    f(t, tmpX, k2);
    for (int j = 0; j < equationsCount; j++)
        tmpX[j] = varX[j] + k2[j] * tau * 0.5;

    //k3
    f(t, tmpX, k3);
    for (int j = 0; j < equationsCount; j++)
        tmpX[j] = varX[j] + k3[j] * tau;

    //k4
    t += 0.5 * tau;
    f(t, tmpX, k4);
    for (int j = 0; j < equationsCount; j++)
    {
        varX[j] += oneSixth * tau * (k1[j] + k2[j] + k2[j] + k3[j] + k3[j] + k4[j]);
    }
}

int DifferentialEquation::solveWithAdams(int methodOrder, int numberOfPoints)
{
    solveWithRungeKutta(true);

    double t = tOutPoints[tOutPoints.size() - 1];
    double step = fabs(last - t) / numberOfPoints;
    //double step = stepAfterRK;
    //int pnum = (int)(fabs(last - t) / step);
    double coef = step / 24.0;

    double **tmpf = new double *[methodOrder];
    for (int i = 0; i < methodOrder; i++)
    {
        tmpf[i] = new double[equationsCount];
    }

    for (int i = 0; i < numberOfPoints; i++)
    {
        t = t + step;
        tOutPoints.push_back(t);
        xOutMatrix.push_back(vec(equationsCount));

        for (int i1 = 1; i1 <= methodOrder; i1++)
            f(tOutPoints[tOutPoints.size() - 1 - i1] /*t - step * i1*/, xOutMatrix[xOutMatrix.size() - i1 - 1].p, tmpf[i1 - 1]);

        for (int j = 0; j < equationsCount; j++)
            xOutMatrix[xOutMatrix.size() - 1].p[j] = xOutMatrix[xOutMatrix.size() - 2].p[j] + coef * (55.0 * tmpf[0][j] - 59.0 * tmpf[1][j] + 37.0 * tmpf[2][j] - 9.0 * tmpf[3][j]);
    };
    outputFile();
    return 0;
}

void DifferentialEquation::outputFile()
{
    ofstream fOut;
    string fileName;
    for (int i = 0; i < equationsCount; i++)
    {
        fileName = "res" + to_string(i) + ".txt";
        fOut.open(fileName, ios::out);
        for (int j = 0; j < xOutMatrix.size(); j++)
        {
            fOut << xOutMatrix[j].p[i] << endl;
        }
        fOut.close();
    }
    fileName = "t.txt";
    fOut.open(fileName, ios::out);
    for (int j = 0; j < xOutMatrix.size(); j++)
    {
        fOut << tOutPoints[j] << endl;
    }
    fOut.close();
}

int DifferentialEquation::solveWithPredictorCorrector(int methodOrder, int numberOfPoints)
{
    solveWithRungeKutta(true);

    double t = tOutPoints[tOutPoints.size() - 1];
    double step = fabs(last - t) / numberOfPoints;
    double coef = step / 24.0;

    double *tmpX = new double[equationsCount];

    double **tmpf = new double *[methodOrder];
    for (int i = 0; i < methodOrder; i++)
    {
        tmpf[i] = new double[equationsCount];
    }

    for (int i = 0; i < numberOfPoints; i++)
    {
        t = t + step;
        tOutPoints.push_back(t);
        xOutMatrix.push_back(vec(equationsCount));

        for (int i1 = 1; i1 <= methodOrder; i1++)
            f(tOutPoints[tOutPoints.size() - 1 - i1] /*t - step * i1*/, xOutMatrix[xOutMatrix.size() - i1 - 1].p, tmpf[i1 - 1]);

        for (int j = 0; j < equationsCount; j++)
            xOutMatrix[xOutMatrix.size() - 1].p[j] = xOutMatrix[xOutMatrix.size() - 2].p[j] + coef * (55.0 * tmpf[0][j] - 59.0 * tmpf[1][j] + 37.0 * tmpf[2][j] - 9.0 * tmpf[3][j]);

        f(tOutPoints[tOutPoints.size() - 1], xOutMatrix[xOutMatrix.size() - 1].p, tmpX);

        for (int j = 0; j < equationsCount; j++)
            xOutMatrix[xOutMatrix.size() - 1].p[j] = xOutMatrix[xOutMatrix.size() - 2].p[j] + coef * (9.0 * tmpX[j] + 19.0 * tmpf[0][j] - 5.0 * tmpf[1][j] + tmpf[2][j]);
    };
    outputFile();

    return 0;
}