#ifndef SEEKROOTS_H
#define SEEKROOTS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <algorithm>

using namespace std;

class SeekRoots
{
protected:
    int degree;
    vector<double> coefficients;
    vector<double> roots;
    vector<double> roots_real;
    vector<double> roots_imaginary;

public:
    SeekRoots();
    ~SeekRoots();
    SeekRoots(const vector<double> &coeffs);
    double Function(double x); // 计算函数值
    double DerivedFunction(double x); // 计算导数值

    void LinearEquation(); // 一次
    void QuadraticEquation(); // 二次
    void CubicEquation(); // 三次
    void QuarticEquation(); // 四次

    double NewtonMethod(double x0); // 牛顿法求根
    void NewtonAllRoots(); // 牛顿法求所有根

    void FindRoots();
    void FindRoots(const vector<double> &coeffs); // 求根
    void PrintFunction(); // 打印函数
    void PrintRoots();
};

#endif