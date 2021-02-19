#include "equation.h"

int DifferentialEquation::solveWithRungeKutta()
{
    double *varX1 = new double[equationsCount];
    double *varX2 = new double[equationsCount];
    double *tmpX = new double[equationsCount];

    double *k1 = new double[equationsCount];
    double *k2 = new double[equationsCount];
    double *k3 = new double[equationsCount];
    double *k4 = new double[equationsCount];

    double tau = (last - first) / (pointsNumber - 1);
    //double TAU =;

    double t = first;
    double T;
    double norm;
    int index;

    int trys;

    tOutPoints.push_back(first);
    xOutMatrix.push_back(vec(equationsCount));

    for (int i = 0; i < equationsCount; i++)
    {
        xOutMatrix[xOutMatrix.size() - 1].p[i] = initialConditions[i];
        varX1[i] = initialConditions[i];
        varX2[i] = initialConditions[i];
    }

    while (t < last)
    {
        //printResult();
        //tau = TAU;
        T = t;
        norm = -1.0;
        index = 0;
        trys = 0;
        while (true)
        {
            trys++;
            index++;
            //cout << "Iteration: " << index++ << " began with t = " << t << " tau = " << tau << endl;
            t = T;

            stepWithRungeKutta(k1, k2, k3, k4, varX1, tmpX, tau, t);

            t = T;

            stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, tau / 2.0, t);
            stepWithRungeKutta(k1, k2, k3, k4, varX2, tmpX, tau / 2.0, t);

            /*cout << "Err: " << index << " " << calculateError(varX1, varX2) << endl;
            for (int i = 0; i < equationsCount; i++)
            {
                cout << varX1[i] << " ";
            }
            cout << endl;

            for (int i = 0; i < equationsCount; i++)
            {
                cout << varX2[i] << " ";
            }
            cout << endl;

            cout << "F: " << varX1[3] - varX2[3] << endl;*/

            if (calculateError(varX1, varX2) < 1e-6 /*|| fabs(calculateError(varX1, varX2) - norm) < 1e-8*/)
            {
                //cout << "Added t = " << t << " tau = " << tau << endl;
                //cout << "Added " << t << endl;
                tOutPoints.push_back(t);
                xOutMatrix.push_back(vec(equationsCount));
                for (int i = 0; i < equationsCount; i++)
                    xOutMatrix[xOutMatrix.size() - 1].p[i] = varX1[i];
                tau *= 2;
                break;
            }
            else
            {
                //cout << "Divided" << endl;
                tau /= 2.0;
            }
            //cout << "Tau = " << tau << endl;
            norm = calculateError(varX1, varX2);
        };
        cout << endl;
        t = T + tau;
    }
    printResult();
    return 0;
}

DifferentialEquation::DifferentialEquation(RightPart rp, int systemNum)
{
    f = rp.getF(systemNum, equationsCount, initialConditions, first, last);
    //xOutMatrix = new double *[equationsCount];
    //for (int i = 0; i < pointsNumber; i++)
    //    xOutMatrix[i] = new double[pointsNumber];

    //for (int j = 0; j < pointsNumber; j++)
    //    tOutPoints[j] = 0.0;

    //for (int i = 0; i < equationsCount; i++)
    //{
    //   for (int j = 0; j < pointsNumber; j++)
    // {
    //   xOutMatrix[i][j] = 0.0;
    //}
    //}
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
    {
        sum += (y2[i] - y1[i]) * ((y2[i] - y1[i]) / (225.0));
    }

    /*double max = fabs(y1[0] - y2[0]);
    for (int i = 1; i < equationsCount; i++)
    {
        if (max < fabs(y1[i] - y2[i]))
        {
            max = fabs(y1[i] - y2[i]);
        }
    }*/

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