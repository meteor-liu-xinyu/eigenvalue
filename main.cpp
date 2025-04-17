#include "eigenvalue.h"
#include <windows.h>
#include <conio.h> // 用于 _getch() 函数

int main()
{
    // 设置控制台输入和输出编码为 UTF-8
    SetConsoleCP(65001); // 设置控制台输入代码页为 UTF-8
    SetConsoleOutputCP(65001); // 设置控制台输出代码页为 UTF-8

    while (true)
    {
        Eigenvalue eigenvalue;
        eigenvalue.SetMatrix();
        eigenvalue.Run();
        cout << "是否继续？(y/n):  ";
        char c;
        cin >> c;
        if (c != 'y' && c != 'Y')
        {
            break;
        }
        cout << endl;
    }

    // Eigenvalue test(
    //     {{1,1,1,1},
    //     {1,1,-1,-1},
    //     {1,-1,1,-1},
    //     {1,-1,-1,1}}
    // );
    // test.Run();

    cout << "按任意键退出...";
    _getch();
    return 0;
}