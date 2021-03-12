#include <iostream>
#include <ostream>
#include"Matrix2d.h"
#include"matmath.h"
using namespace std;

int main()
{
    Matrix2d<float> mat1(4, 3, 5.6);
    Matrix2d<float> mat2(4, 3, -7.77777), mat3;
    mat3 = mat1 + sqrt(mat2) * sin(mat1);
    cout <<mat3;
    mat3 = hypot(mat1, mat2);
    cout <<mat3;
    mat3 = abs(mat2);
    cout <<mat3;
    cout<<rand();
    return 0;
}
