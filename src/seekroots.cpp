#include "seekroots.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

SeekRoots::SeekRoots()
{
    
}

SeekRoots::~SeekRoots()
{
}

SeekRoots::SeekRoots(const vector<double> &coeffs)
{
    coefficients = coeffs;
    degree = coeffs.size() - 1;
}

double SeekRoots::Function(double x)
{
    double result = 0;
    for (size_t i = 0; i <= degree; i++)
    {
        result += coefficients[degree-i] * pow(x, i);
    }
    return result;
}

double SeekRoots::DerivedFunction(double x)
{
    double result = 0;
    for (size_t i = 1; i <= degree; i++)
    {
        result += i * coefficients[degree-i] * pow(x, i - 1);
    }
    return result;
}

void SeekRoots::PrintFunction()
{
    cout << "方程为：";
    for (size_t i = 0; i < degree; i++)
    {
        if (i == 0)
        {
            cout << coefficients[i] << "x^" << degree-i;
        }
        else
        {
            if (coefficients[i] > 0)
            {
                cout << " + " << coefficients[i] << "x^" << degree-i;
            }
            else
            {
                cout << " - " << -coefficients[i] << "x^" << degree-i;
            }
        }
    }
    if (coefficients[degree] > 0)
    {
        cout << " + " << coefficients[degree] << " = 0" << endl;
    }
    else
    {
        cout << " - " << -coefficients[degree] << " = 0" << endl;
    }
}



void SeekRoots::FindRoots()
{
    if (coefficients[0] < 0) // 如果首项系数为负，则将多项式的所有系数取反
    {
        for (int i = 0; i < degree; i++)
        {
            coefficients[i] = -coefficients[i];
        }
    }
    if (coefficients[degree] == 0)
    {
        vector<double> newcoeff = coefficients;
        newcoeff.pop_back();
        SeekRoots newfunction(newcoeff);
        newfunction.FindRoots();
        roots = newfunction.roots;
        roots.push_back(0);
        return;
    }
    if (degree >= 5)
    {
        NewtonAllRoots();
    }
    else if (degree == 4)
    {
        QuarticEquation();
    }
    else if (degree == 3)
    {
        CubicEquation();
    }
    else if (degree == 2)
    {
        QuadraticEquation();
    }
    else if (degree == 1)
    {
        LinearEquation();
    }
    else
    {
        cout << "error" << endl;
    }
}

void SeekRoots::FindRoots(const vector<double> &coeffs)
{
    coefficients = coeffs;
    degree = coeffs.size() - 1;
    if (coefficients[0] < 0) // 如果首项系数为负，则将多项式的所有系数取反
    {
        for (int i = 0; i < degree; i++)
        {
            coefficients[i] = -coefficients[i];
        }
    }
    if (coefficients[degree] == 0)
    {
        vector<double> newcoeff = coefficients;
        newcoeff.pop_back();
        SeekRoots newfunction(newcoeff);
        newfunction.FindRoots();
        roots = newfunction.roots;
        roots.push_back(0);
        return;
    }
    if (degree >= 5)
    {
        NewtonAllRoots();
    }
    else if (degree == 4)
    {
        QuarticEquation();
    }
    else if (degree == 3)
    {
        CubicEquation();
    }
    else if (degree == 2)
    {
        QuadraticEquation();
    }
    else if (degree == 1)
    {
        LinearEquation();
    }
    else
    {
        cout << "error" << endl;
    }
}

void SeekRoots::PrintRoots()
{
    if (roots_imaginary.size() == 0)
    {
        cout << "方程的根为：";
        for (const auto &root : roots)
        {
            if (abs(root - int(root)) < 1e-6)
            {
                cout << int(root) << "  ";
            }
            else
            {
                cout << fixed << setprecision(3);
                cout << root << "  ";
            }
        }
        cout << endl;
    }
    else if (roots.size() == 0)
    {
        cout << fixed << setprecision(4);
        cout << "方程的复根为：" << endl;
        for (size_t i = 0; i < roots_real.size(); i++)
        {
            if (roots_imaginary[i] > 0)
            {
                cout << setw(5) << roots_real[i] << " + " << setw(5) << roots_imaginary[i] << " i" << endl;
            }
            else
            {
                cout << setw(5) << roots_real[i] << " - " << setw(5) << -roots_imaginary[i] << " i" << endl;
            }
        }
    }
    else
    {
        cout << "方程的实根为：";
        for (const auto &root : roots)
        {
            if (abs(root - int(root)) < 1e-6)
            {
                cout << int(root) << "  ";
            }
            else
            {
                cout << fixed << setprecision(3);
                cout << root << "  ";
            }
        }
        cout << endl;
        cout << fixed << setprecision(4);
        cout << "方程的复根为：" << endl;
        for (size_t i = 0; i < roots_real.size(); i++)
        {
            if (roots_imaginary[i] > 0)
            {
                cout << setw(5) << roots_real[i] << " + " << setw(5) << roots_imaginary[i] << " i" << endl;
            }
            else
            {
                cout << setw(5) << roots_real[i] << " - " << setw(5) << -roots_imaginary[i] << " i" << endl;
            }
        }
    }
}

void SeekRoots::LinearEquation()
{
    roots.push_back(-coefficients[1]/coefficients[0]);
}

void SeekRoots::QuadraticEquation()
{
    double a = coefficients[0];
    double b = coefficients[1];
    double c = coefficients[2];
    double delta = b * b - 4 * a * c;

    if (delta > 0)
    {
        roots.push_back((-b - sqrt(delta)) / (2 * a));
        roots.push_back((-b + sqrt(delta)) / (2 * a));
    }
    else if (delta == 0)
    {
        roots.push_back(-b / (2 * a));
        roots.push_back(-b / (2 * a)); // 重根
    }
    else
    {
        roots_real.push_back(-b / (2 * a));
        roots_imaginary.push_back(sqrt(-delta) / (2 * a));
        roots_real.push_back(-b / (2 * a));
        roots_imaginary.push_back(-sqrt(-delta) / (2 * a));
    }
}

void SeekRoots::CubicEquation()
{
    // 卡尔丹公式
    // 三次方程 x^3 + bx^2 + cx + d = 0
    double a = coefficients[0];
    double b = coefficients[1];
    double c = coefficients[2];
    double d = coefficients[3];

    // 标准化
    if (a != 1)
    {
        b /= a;
        c /= a;
        d /= a;
    }

    // 化为无二次项的形式
    double p = 0.0;
    double q = 0.0;
    if (b == 0)
    {
        p = c;
        q = d;
    }
    else
    {
        p = (3 * c - b * b) / (3 * a * a);
        q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);
    }
    double delta = q * q / 4 + p * p * p / 27; // 判别式

    if (delta == 0)
    {
        if (p == 0 && q == 0) // 三重零根
        {
            roots.push_back(0);
            roots.push_back(0);
            roots.push_back(0);
        }
        else // 三实根
        {
            double u = cbrt(-q / 2);
            roots.push_back(2 * u);
            roots.push_back(-u);
            roots.push_back(-u);
        }
    }
    else if (delta > 0) // 一实两复根
    {
        double u = cbrt(-q / 2 + sqrt(delta));
        double v = cbrt(-q / 2 - sqrt(delta));
        roots.push_back(u + v); // 实根
        roots_real.push_back(-u / 2 - v / 2);
        roots_imaginary.push_back(sqrt(3) * (u - v) / 2); // 复根2
        roots_real.push_back(-u / 2 - v / 2);
        roots_imaginary.push_back(-sqrt(3) * (u - v) / 2); // 复根3
    }
    else // 三实根
    {
        double r = sqrt(-p / 3);
        double phi = acos(-q / (2 * r * r * r));
        roots.push_back(2 * r * cos(phi / 3)); // 实根1
        roots.push_back(2 * r * cos((phi + 2 * M_PI) / 3)); // 实根2
        roots.push_back(2 * r * cos((phi + 4 * M_PI) / 3)); // 实根3
    }
}

void SeekRoots::QuarticEquation()
{
    // 天衍公式
    double a = coefficients[0];
    double b = coefficients[1];
    double c = coefficients[2];
    double d = coefficients[3];
    double e = coefficients[4];

    double delta_D = 3 * b * b - 8 * a * c;
    double delta_E = -b * b * b + 4 * a * b * c - 8 * a * a * d;
    double delta_F = 3 * b * b * b * b + 16 * a * a * c * c - 16 * a * b * b * c + 16 * a * a * b * d - 64 * a * a * a * e;

    double delta_A = delta_D * delta_D - 3 * delta_F;
    double delta_B = delta_D * delta_F - 9 * delta_E * delta_E;
    double delta_C = delta_F * delta_F - 3 * delta_D * delta_E * delta_E;

    double delta = delta_B * delta_B - 4 * delta_A * delta_C;

    if (delta_D == 0 && delta_E == 0 && delta_F == 0) // 四重根
    {
        double x = -b / (4 * a);
        roots.push_back(x);
        roots.push_back(x);
        roots.push_back(x);
        roots.push_back(x);
    }
    else if (delta_A == 0 && delta_B == 0 && delta_C == 0 && delta_D * delta_E * delta_F != 0) // 四实根三重根
    {
        roots.push_back((-b * delta_D + 9 * delta_E) / (4 * a * delta_D));
        double x = (-b * delta_D - 3 * delta_E) / (4 * a * delta_D);
        roots.push_back(x);
        roots.push_back(x);
        roots.push_back(x); // 三重根
    }
    else if (delta_E == 0 && delta_F == 0 && delta_D != 0) // 两对二重根
    {

        double x_real = -b / (4 * a);
        if (delta_D > 0) // 两对实二重根
        {
            double x_tail = sqrt(delta_D) / (4 * a);
            roots.push_back(x_real+x_tail);
            roots.push_back(x_real+x_tail);
            roots.push_back(x_real-x_tail);
            roots.push_back(x_real-x_tail);
        }
        else // 两对复二重根
        {
            double x_imaginary = sqrt(-delta_D) / (4 * a);
            roots_real.push_back(x_real);
            roots_imaginary.push_back(x_imaginary);
            roots_real.push_back(x_real);
            roots_imaginary.push_back(x_imaginary);
            roots_real.push_back(x_real);
            roots_imaginary.push_back(-x_imaginary);
            roots_real.push_back(x_real);
            roots_imaginary.push_back(-x_imaginary);
        }
    }
    else if (delta_A * delta_B * delta_C != 0 && delta == 0) // 一对二重实根
    {
        double x_real = (-b + 2 * delta_A * delta_E / delta_B) / (4 * a);
        if (delta_A * delta_B > 0) // 其余两根为不等实根
        {
            double x_tail = sqrt(2 * delta_B /delta_A) / (4 * a);
            roots.push_back(x_real+x_tail);
            roots.push_back(x_real-x_tail);
        }
        else // 其余两根为共轭复根
        {
            double x_imaginary = sqrt(-2 * delta_B /delta_A) / (4 * a);
            roots_real.push_back(x_real);
            roots_imaginary.push_back(x_imaginary);
            roots_real.push_back(x_real);
            roots_imaginary.push_back(-x_imaginary);
        }
        double x3 = (-b - 2 * delta_A * delta_E / delta_B) / (4 * a);
        roots.push_back(x3);
        roots.push_back(x3);
    }
    else if (delta > 0) // 两个不等实根一队共轭复根
    {
        double z_real = delta_A * delta_D - 3 * delta_B / 2;
        double z1 = z_real + 3 * sqrt(delta) / 2;
        double z2 = z_real - 3 * sqrt(delta) / 2;
        double temp = cbrt(z1) + cbrt(z2);
        double z = delta_D * delta_D - delta_D * temp + temp * temp - 3 * delta_A;
        double sgn_E = (delta_E > 0) ? 1 : -1;
        double x_real = (-b + sgn_E * sqrt((delta_D + temp) / 3)) / (4 * a);
        double x_tail = sqrt((-2 * delta_D + temp + 2 * sqrt(z)) / 3) / (4 * a);
        roots.push_back(x_real + x_tail);
        roots.push_back(x_real - x_tail);
        roots_real.push_back(x_real);
        roots_imaginary.push_back(x_tail); // 复根1
        roots_real.push_back(x_real);
        roots_imaginary.push_back(-x_tail); // 复根2
    }
    else if (delta < 0)
    {
        double x_real = -b/(4 * a);
        double x_tail1 = sqrt(delta_D + 2 * sqrt(delta_F)) / (4 * a);
        double x_tail2 = sqrt(delta_D - 2 * sqrt(delta_F)) / (4 * a);
        if (delta_E == 0 && delta_D > 0 && delta_F > 0)
        {
            roots.push_back(x_real + x_tail1);
            roots.push_back(x_real - x_tail1);
            roots.push_back(x_real + x_tail2);
            roots.push_back(x_real - x_tail2);
        }
        else if (delta_E == 0 && delta_D < 0 && delta_F > 0)
        {
            roots_real.push_back(x_real);
            roots_imaginary.push_back(x_tail1);
            roots_real.push_back(x_real);
            roots_imaginary.push_back(-x_tail1);
            roots_real.push_back(x_real);
            roots_imaginary.push_back(x_tail2);
            roots_real.push_back(x_real);
            roots_imaginary.push_back(-x_tail2);
        }
        else if (delta_E == 0 && delta_F < 0)
        {
            double x_tail = sqrt(2 * delta_D + 2 * sqrt(delta_A - delta_F)) / (8 * a);
            double x_imaginary = sqrt(-2 * delta_D + 2 * sqrt(delta_A - delta_F)) / (8 * a);
            roots_real.push_back(x_real+x_tail);
            roots_imaginary.push_back(x_imaginary);
            roots_real.push_back(x_real+x_tail);
            roots_imaginary.push_back(-x_imaginary);
            roots_real.push_back(x_real-x_tail);
            roots_imaginary.push_back(x_imaginary);
            roots_real.push_back(x_real-x_tail);
            roots_imaginary.push_back(-x_imaginary);
        }
        else if (delta_E != 0)
        {
            double phi = acos((3 * delta_B - 2 * delta_A * delta_D) / (2 * delta_A * sqrt(delta_A)));
            double y1 = (delta_D - 2 * sqrt(delta_A) * cos(phi / 3)) / 3;
            double y2 = (delta_D + sqrt(delta_A) * (cos(phi / 3) + sqrt(3) * sin(phi / 3))) / 3;
            double y3 = (delta_D + sqrt(delta_A) * (cos(phi / 3) - sqrt(3) * sin(phi / 3))) / 3;
            double sgn_E = (delta_E > 0) ? 1 : -1;
            double y1_sqrt = sqrt(abs(y1));
            double y2_sqrt = sqrt(abs(y2));
            double y3_sqrt = sqrt(abs(y3));
            if (delta_D > 0 && delta_F > 0)
            {
                roots.push_back(x_real + ((sgn_E * y1_sqrt) + y2_sqrt + y3_sqrt) / (4 * a));
                roots.push_back(x_real + ((sgn_E * y1_sqrt) - y2_sqrt - y3_sqrt) / (4 * a));
                roots.push_back(x_real - ((sgn_E * y1_sqrt) + y2_sqrt - y3_sqrt) / (4 * a));
                roots.push_back(x_real - ((sgn_E * y1_sqrt) - y2_sqrt + y3_sqrt) / (4 * a));
            }
            else
            {
                roots_real.push_back(x_real + y2_sqrt / (4 * a));
                roots_imaginary.push_back((sgn_E * y1_sqrt - y3_sqrt) / (4 * a));
                roots_real.push_back(x_real + y2_sqrt / (4 * a));
                roots_imaginary.push_back(-(sgn_E * y1_sqrt - y3_sqrt) / (4 * a));
                roots_real.push_back(x_real - y2_sqrt / (4 * a));
                roots_imaginary.push_back((sgn_E * y1_sqrt + y3_sqrt) / (4 * a));
                roots_real.push_back(x_real - y2_sqrt / (4 * a));
                roots_imaginary.push_back(-(sgn_E * y1_sqrt + y3_sqrt) / (4 * a));
            }
        }
    }
}

double SeekRoots::NewtonMethod(double x0)
{
    double takeEqualBoundaries = 1e-4;
    if (coefficients[degree] > pow(10,degree))
    {
        takeEqualBoundaries = 1e-2;
    }
    double eps = 1e-6; // 精度
    double x = x0;
    double fx = Function(x);
    double dfx = DerivedFunction(x); // 计算导数
    for (int i = 0; i < 1000; i++)
    {
        if (abs(fx) < eps)
        {
            if (abs(dfx) < takeEqualBoundaries)
            {
                for (double temp_x = x - 0.1; temp_x < x + 0.1; temp_x += takeEqualBoundaries/100)
                {
                    if (abs(fx) > abs(Function(temp_x)))
                    {
                        x = temp_x;
                        fx = Function(x);
                    }
                    else
                    {
                        temp_x += eps;
                    }
                }
            }
            break;
        }
        else if (abs(dfx) < takeEqualBoundaries) // 导数为0，选择新的初始值迭代
        {
            while (abs(DerivedFunction(x)) < takeEqualBoundaries)
            {
                if (x < x0)
                {
                    x -= takeEqualBoundaries * 10;
                }
                else
                {
                    x += takeEqualBoundaries * 10;
                }
            }
            return NewtonMethod(x);
        }
        x -= fx / dfx; // 牛顿迭代公式
        fx = Function(x);
        dfx = (Function(x + eps) - Function(x)) / eps; // 计算导数
    }

    if (abs(fx) > eps)
    {
        cout << "牛顿法求根失败" << endl;
        return 0;
    }
    else
    {
        return x;
    }
}

void SeekRoots::NewtonAllRoots()
{
    double takeEqualBoundaries = 1e-4;
    int multiple = 10;
    if (coefficients[degree] > pow(10,degree))
    {
        multiple = (int)(coefficients[degree] / pow(10,degree+1));
        takeEqualBoundaries /= multiple;
    }

    vector<double> coeffs;
    for (int i = 0; i < degree; i++)
    {
        coeffs.push_back(coefficients[i] * (degree - i));
    }
    SeekRoots derivedFunction(coeffs);

    double zero_point = NewtonMethod(0); // 先求出一个根
    roots.push_back(zero_point);

    if (roots.size() == 0)
    {
        cout << "牛顿法求根失败" << endl;
        return;
    }
    double a = zero_point - multiple; // 左边界
    double fa = Function(a);
    double dfa = DerivedFunction(a); // 计算导数
    double b = zero_point + multiple; // 右边界
    double fb = Function(b);
    double dfb = DerivedFunction(b); // 计算导数

    while (true)
    {
        if (degree %2 == 0 && fa > 0 && dfa < 0 && dfa > DerivedFunction(a-1))
        {
            double derivedroot = derivedFunction.NewtonMethod(a);
            if (derivedroot >= a)
            {
                a = NewtonMethod(a);
                if (abs(a - zero_point) > takeEqualBoundaries)
                {
                    roots.push_back(a);
                }
                a = max(a, derivedroot);
                break;
            }
        }
        else if (degree %2 == 1 && fa < 0 && dfa > 0 && dfa < DerivedFunction(a-1))
        {
            double derivedroot = derivedFunction.NewtonMethod(a);
            if (derivedroot >= a)
            {
                a = NewtonMethod(a);
                if (abs(a - zero_point) > takeEqualBoundaries)
                {
                    roots.push_back(a);
                }
                a = max(a, derivedroot);
                break;
            }
        }
        else
        {
            a -= multiple;
            fa = Function(a);
            dfa = DerivedFunction(a); // 计算导数
        }
    }
    while (true)
    {
        if (fb > 0 && dfb > 0 && dfb < DerivedFunction(b+1))
        {
            double derivedroot = derivedFunction.NewtonMethod(b);
            if (derivedroot <= b)
            {
                b = NewtonMethod(b);
                if (abs(b - zero_point) > takeEqualBoundaries)
                {
                    roots.push_back(b);
                }
                b = min(b, derivedroot);
                break;
            }
        }
        else
        {
            b += multiple;
            fb = Function(b);
            dfb = DerivedFunction(b); // 计算导数
        }
    }

    vector<double> edge;
    for (double x = a; x < b; x += 0.005 * multiple)
    {
        edge.push_back(x);
    }

    // 调整顺序,每10个x移动到前面,从而尽可能减少遍历次数
    for (int i = 1; i < (int)(edge.size()/multiple); i += 1)
    {
        double temp = edge[i * multiple];
        edge.erase(edge.begin() + i * multiple);
        edge.insert(edge.begin() + i, temp);
    }

    for (const double &x : edge)
    {
        double root = NewtonMethod(x);
        if (root == 0 && abs(Function(x)) > takeEqualBoundaries)
        {
            continue;
        }

        bool isUnique = true;
        for (double existroot : roots)
        {
            if (abs(existroot - root) < takeEqualBoundaries)
            {
                isUnique = false;
                break;
            }
        }
        if (isUnique)
        {
            if (abs(root - (int)root) < takeEqualBoundaries)
            {
                roots.push_back((int)root);
            }
            else
            {
                roots.push_back(root);
            }
            if (roots.size() == degree)
            {
                break;
            }
        }
    }

    sort(roots.begin(), roots.end());
}