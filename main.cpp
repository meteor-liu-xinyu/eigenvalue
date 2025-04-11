#include "eigenvalue.h"
#include<windows.h>

int main()
{
    SetConsoleOutputCP(65001);

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

    return 0;
}