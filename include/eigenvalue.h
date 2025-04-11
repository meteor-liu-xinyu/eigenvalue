#ifndef EIGENVALUE_H
#define EIGENVALUE_H

#include "seekroots.h"

#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

class SeekRoots;

class Eigenvalue : public SeekRoots
{
private:
    int n;
    vector<double> coefficient; //系数
    vector<vector<float>> matrix;//矩阵的存储
public:
    Eigenvalue();
    Eigenvalue(const vector<vector<float>> &input);
    ~Eigenvalue();
    int MasterChildFormula(int k, const vector<int>& row_col_num); // 求主子式
    void Eigenpolynomials(); // 求特征多项式系数

    int Check0(vector<vector<float>> &temp, int p, int q, int &exchange_num);
    vector<vector<float>> TransSteppedType(const vector<vector<float>> &input);//矩阵的阶梯形化简
    double Det(const vector<vector<float>> &input); // 求行列式

    void SetMatrix();
    void Run();
    void printFeaturePolynomials();
    void printEigenvalues(); // 打印特征值

};






#endif