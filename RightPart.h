typedef void (*func)(double, double *, double *);
class RightPart
{
public:
    func getF(int systemNum, int &equationsCount, double *&initConds, double &first, double &last)
    {
        switch (systemNum)
        {
        case 1:
            equationsCount = 2;
            initialConditions = new double[equationsCount];
            for (int i = 0; i < equationsCount; i++)
                initialConditions[i] = 0.0;
            initConds = initialConditions;
            first = 0.0;
            last = 10.0;
            return f1;
        case 2:
            equationsCount = 4;
            initialConditions = new double[equationsCount];
            for (int i = 0; i < equationsCount; i++)
                initialConditions[i] = 0.0;
            initConds = initialConditions;
            first = 0.0;
            last = 10.0;
            return f2;
        case 3:
            equationsCount = 2;
            initialConditions = new double[equationsCount];
            for (int i = 0; i < equationsCount; i++)
                initialConditions[i] = 0.0;
            initConds = initialConditions;
            first = 0.0;
            last = 1.0;
            return f3;
        case 4:
            equationsCount = 2;
            initialConditions = new double[equationsCount];
            initialConditions[0] = 1.0;
            initialConditions[1] = 0.1;
            initConds = initialConditions;
            first = 0.0;
            last = 50.0;
            return f4;
        case 5:
            equationsCount = 2;
            initialConditions = new double[equationsCount];
            initialConditions[0] = 1.0;
            initialConditions[1] = 0.0;
            initConds = initialConditions;
            first = 0.0;
            last = 1.0;
            return f5;
        }
    }

private:
    double *initialConditions;
    static void f1(double t, double *x, double *res)
    {
        res[0] = 2.0 * t;
        res[1] = 3.0 * x[0];
    };

    static void f2(double t, double *x, double *res)
    {
        res[0] = 2.0 * t;
        res[1] = 3.0 * (x[0] + t * t) / 2.0;
        res[2] = 4.0 * (x[1] + t * x[0] + t * t * t) / 3.0;
        res[3] = 5.0 * (x[2] + t * x[1] + t * t * x[0] + t * t * t * t) / 4.0;
    }

    static void f3(double t, double *x, double *res)
    {
        res[0] = 2 * x[0] + x[1] * x[1] - 1;
        res[1] = 6 * x[0] - x[1] * x[1] + 1;
    }

    static void f4(double t, double *x, double *res)
    {
        double H = (x[0] * x[0] + x[1] * x[1]) * (x[0] * x[0] + x[1] * x[1]) - 0.02 * (x[0] * x[0] - x[1] * x[1]);
        double H_x = -0.04 * x[0] + 4.0 * x[0] * (x[0] * x[0] + x[1] * x[1]);
        double H_y = 0.04 * x[1] + 4.0 * x[1] * (x[0] * x[0] + x[1] * x[1]);
        res[0] = H_y - 0.1 * H * H_x + 0.15 * x[1] * sin(0.25 * t);
        res[1] = -H_x - 0.1 * H * H_y + 0.15 * x[1] * sin(0.25 * t);
    }

    static void f5(double t, double *x, double *res)
    {
        double k = 20.0;
        double m = 0.3;
        res[0] = x[1];
        res[1] = -(k / m) * x[0];
    }
};