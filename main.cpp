#include <iostream>
#include <ostream>
#include"Matrix2d.h"
#include"matmath.h"
using namespace std;

int main()
{
    Matrix2d<float> mat0(4, 3, 1);
    mat0 << 2, 2, 3,
            4, 5, 6,
            4, 5, 6,
            4, 5, 6;
    cout << mat0.reshape(2,6).Tr();

    Matrix2d<float> mat1(4, 3, 5.6);
    Matrix2d<float> mat2(4, 3, 7.77777), mat3;
    mat3 = mat1 + sqrt(mat2) * (mat0);
    cout <<mat3;

    mat3 = hypot(mat1, mat2);
    cout <<mat3;

    mat3 = abs(mat2);
    cout <<mat3;
    cout<<rand();
    return 0;
}
