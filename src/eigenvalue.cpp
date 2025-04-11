#include "eigenvalue.h"

Eigenvalue::Eigenvalue()
{
}

Eigenvalue::Eigenvalue(const vector<vector<float>> &input)
{
    n = input.size();
    matrix = input;
}

Eigenvalue::~Eigenvalue()
{
}

int Eigenvalue::Check0(vector<vector<float>> &temp, int p, int q, int &exchange_num)
{
    int n = temp.size();
    for (int i = p; i < n; i++)
    {
        if (temp[i][q]!=0)//判断第i行的第q个元素是否为0
        {
            if (i==p)//如果第p行的第q个元素不为0，则直接返回1
            {
                return 1;
            }
            else//如果第p行的第q个元素为0，则继续向下寻找，直至第i行的第q个元素不为0
            {
                //将第i行元素与第p行元素交换
                vector<float> exchange = temp[i];
                temp[i] = temp[p];
                temp[p] = exchange;
                exchange_num++;
                return 1;
            }
        }
    }
    return 0;
}

vector<vector<float>> Eigenvalue::TransSteppedType(const vector<vector<float>> &input)
{
    int n = input.size();
    int m = input[0].size();
    vector<vector<float>> temp = input;
    int exchange_num = 0;
    int min = n < m ? n : m;
    int i = 0,ii = 0;
    for (i = 0; i < min && ii < m; i++,ii++)
    {
        if (Check0(temp, i, ii, exchange_num) == 1)//保证第i行第i个元素不为0
        {
            //做初等行变换将后面的行第i个元素变为0
            for (int j = i + 1; j < n; j++)
            {
                float multiple = temp[j][ii] / temp[i][ii];
                for (int k = 0; k < m; k++)
                {
                    temp[j][k] = temp[j][k] - temp[i][k] * multiple;
                }
            }
        }
        else//否则第i行起每一行第i列全为0
        {
            i--;
        }
        for (int s = 0; s < n; s++)
        {
            for (int t = 0; t < m; t++)
            {
                if (temp[s][t]>-0.0000001 && temp[s][t]<0.0000001)
                {
                    temp[s][t] = 0;
                }
            }
        }
    }

    if (exchange_num % 2 == 1)//如果交换行数为奇数，则在最后一个非零行乘-1
    {
        for (int i = n-1; i >= 0; i--)
        {
            if (temp[i][m-1] != 0)
            {
                for (int j = 0; j < m; j++)
                {
                    temp[i][j] = -temp[i][j];
                }
                break;
            }
        }
    }

    return temp;
}

double Eigenvalue::Det(const vector<vector<float>> &input)
{
    int n = input.size();
    int m = input[0].size();
    if (n != m)
    {
        // 矩阵不是方阵，无法计算行列式
        return -1;
    }
    double det = 1;
    vector<vector<float>> temp = TransSteppedType(input);
    for (int i = 0; i < n; i++)
    {
        det *= temp[i][i];
    }
    if (det>-0.0000001 && det<0.0000001)
    {
        det = 0;
    }
    return det;
}

int Eigenvalue::MasterChildFormula(int k, const vector<int>& row_col_num)
{
    vector<vector<float>> temp = matrix;
    for (int i = 0; i < k; i++)
    {
        // 删除第row_col_num[i]行
        temp.erase(temp.begin() + row_col_num[i]);
    }
    for (int i = 0; i < k; i++)
    {
        // 删除第row_col_num[i]列
        for (int j = 0; j < n-k; j++)
        {
            temp[j].erase(temp[j].begin() + row_col_num[i]);
        }
    }

    return Det(temp);
}

void Eigenvalue::Eigenpolynomials()
{
    vector<int> row_col_num;
    coefficient.push_back(1); // n次项系数

    for (int k = n-1; k >= 0; k--)
    {
        row_col_num.clear();
        vector<int> combination(k, 0);
        for (int i = 0; i < k; i++)
        {
            combination[i] = i;
        }

        float result = 0.0;
        while (true)
        {
            // 从0到n-1从小到大生成所有可能的k个数组合
            row_col_num = combination;
            int term = MasterChildFormula(k, row_col_num);
            result += term;

            int i;
            for (i = k - 1; i >= 0; i--)
            {
                if (combination[i] != i + n - k)
                {
                    break;
                }
            }

            if (i < 0)
            {
                break;
            }

            combination[i]++;
            for (int j = i + 1; j < k; j++)
            {
                combination[j] = combination[j - 1] + 1;
            }
        }

        if (k % 2 == 0)
        {
            coefficient.push_back(result);
        }
        else
        {
            coefficient.push_back(-result);
        }
    }
}

void Eigenvalue::SetMatrix()
{
    cout << "请输入矩阵的阶数n:";
    cin >> n;
    matrix.resize(n, vector<float>(n, 0));
    cout << "请输入矩阵的元素:"<< endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> matrix[i][j];
        }
    }
    cout << endl;
}

void Eigenvalue::printFeaturePolynomials()
{
    string lambda = "λ";
    cout << "特征多项式为：" << lambda << "^" << n;
    for (int i = 1; i < coefficient.size(); i++)
    {
        double abstemp = coefficient[i];
        if (coefficient[i] < 0)
        {
            cout << " - ";
            abstemp = -abstemp;
        }
        else
        {
            cout << " + ";
        }
        // 如果为整数，则不显示小数点
        if (abstemp == (int)abstemp)
        {
            cout << (int)abstemp;
        }
        else
        {
            cout << fixed << setprecision(3) << abstemp;
        }
        cout << " ";
        if (i == coefficient.size() - 1)
        {
            break;
        }
        cout << lambda << "^" << n - i;
    }
    cout << endl;
}

void Eigenvalue::printEigenvalues()
{
    if (n >= 5 && coefficient[n - 1] > pow (10, n))
    {
        cout << "粗略的";
    }
    cout << "特征值为：";
    if (roots_imaginary.size() == 0) // 只有实根
    {
        for (const auto &root : roots)
        {
            if (root == (int)root)
            {
                cout << (int)root << "  ";
            }
            else
            {
                cout << fixed << setprecision(3) << root << "  ";
            }
        }
        cout << endl;
    }
    else if (roots.size() == 0) // 只有复根
    {
        cout << endl;
        cout << fixed << setprecision(4);
        for (int i = 0; i < roots_real.size(); i++)
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
    else // 既有实根又有复根
    {
        for (const auto &root : roots)
        {
            if (root == (int)root)
            {
                cout << (int)root << "  ";
            }
            else
            {
                cout << fixed << setprecision(3) << root << "  ";
            }
        }
        cout << endl;
        cout << fixed << setprecision(4);
        for (int i = 0; i < roots_real.size(); i++)
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

void Eigenvalue::Run()
{
    Eigenpolynomials();
    printFeaturePolynomials();
    FindRoots(coefficient);
    printEigenvalues();
}