//исправить на возврат указателя на функцию
typedef void (*func)(double, double*, double*);
class RightPart{

public:
func getF(int systemNum, int& equationsCount, double*& initConds, double& first, double& last){
switch (systemNum){
case 1:
    equationsCount = 2;
    initialConditions = new double[equationsCount];
    for(int i = 0; i < equationsCount; i++)
        initialConditions[i] = 0.0;
        initConds = initialConditions;
    first = 0.0;
    last = 10.0;
    return f1;
case 2:
    equationsCount = 4;
    initialConditions = new double[equationsCount];
    for(int i = 0; i < equationsCount; i++)
        initialConditions[i] = 0.0;
    initConds = initialConditions;
    first = 0.0;
    last = 10.0;
    return f2;
}
}

private:
double* initialConditions;
static void f2(double t, double* x, double* res){
    res[0] = 2.0*t;
    res[1] = 3.0*(x[0]+t*t)/2.0;
    res[2] = 4.0*(x[1]+t*x[0]+t*t*t)/3.0;
    res[3] = 5.0*(x[2]+t*x[1]+t*t*x[0]+t*t*t*t)/4.0;
}

static void f1(double t, double* x, double* res){
    res[0] = 2.0 * t;
    res[1] = 3.0 * x[0];
};
};