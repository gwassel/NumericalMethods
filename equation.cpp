#include "equation.h"

int DifferentialEquation::solveWithRungeKutta(){
double* varX = new double[equationsCount];
double* tmpX = new double[equationsCount];

double* k1 = new double[equationsCount];
double* k2 = new double[equationsCount];
double* k3 = new double[equationsCount];
double* k4 = new double[equationsCount];

double tau = (last - first) / (pointsNumber - 1);
double oneSixth = (1.0/6.0);
double t = first;
tOutPoints[0] = first;
int iOut = 1;
for(int i = 0; i < equationsCount; i++){
    xOutMatrix[i][0] = initialConditions[i];
    varX[i] = initialConditions[i];
}

for(int i = 1; i < pointsNumber; i++){
//k1
f(t, varX, k1);
for(int j = 0; j < equationsCount; j++)
    tmpX[j] = varX[j] + k1[j] * tau * 0.5;

//k2
t += 0.5 * tau;
f(t, tmpX, k2);
for(int j = 0; j < equationsCount; j++)
    tmpX[j] = varX[j] + k2[j] * tau * 0.5;

//k3
f(t, tmpX, k3);
for(int j = 0; j < equationsCount; j++)
    tmpX[j] = varX[j] + k3[j] * tau;

//k4
t += 0.5 * tau;
f(t, tmpX, k4);
for(int j = 0; j < equationsCount; j++)
    varX[j] += oneSixth * tau * (k1[j] + k2[j] + k2[j] + k3[j] + k3[j] + k4[j]);

tOutPoints[iOut] = t;
for(int j = 0; j < equationsCount; j++)
    xOutMatrix[j][iOut] = varX[j];
iOut++;

cout << "X:"<< i <<endl;
for(int i1 = 0; i1 < equationsCount; i1++){
    for(int j = 0; j < pointsNumber; j++){
        cout << xOutMatrix[i1][j] << " ";
    }
    cout << endl;
}

}

cout << "T:" <<endl;
for(int j = 0; j < pointsNumber; j++)
        cout << tOutPoints[j] << " ";
cout << endl;

cout << "X:" <<endl;
for(int i = 0; i < equationsCount; i++){
    for(int j = 0; j < pointsNumber; j++){
        cout << xOutMatrix[i][j] << " ";
    }
    cout << endl;
}
return 0;
}


DifferentialEquation::DifferentialEquation(RightPart rp, int systemNum){
f = rp.getF(systemNum, equationsCount, initialConditions, first, last);
tOutPoints = new double[pointsNumber];
xOutMatrix = new double*[equationsCount];
for(int i = 0; i < pointsNumber;i++)
    xOutMatrix[i] = new double[pointsNumber];

for(int j = 0; j < pointsNumber; j++)
        tOutPoints[j] = 0.0;

for(int i = 0; i < equationsCount; i++){
    for(int j = 0; j < pointsNumber; j++){
        xOutMatrix[i][j] = 0.0;
    }
}

}