#include <iostream>
#include <ostream>
#include"Matrix2d.h"
#include"matmath.h"
#include"linalg.h"
using namespace std;

int main()
{
    Matrix2d<float> mat0(4, 3, 0);
    mat0 << 0, 2, 3,
            0, 5, 7,
            0, 1, 6,
            0, 5, 0;
    Matrix2d<float> mat00 = Matrix2d<float>::eye(3);
    //double q = mat_triu(mat0);
    //cout<<det(mat0);
    cout<<mat0;
    cout<<Matrix2d<float>::vstack(mat0, mat00);

/*
    cout << mat0.reshape(2,6).Tr();
    cout << mat0.block(1,0,3,3);

    Matrix2d<float> mat1(4, 3, 5.6);
    Matrix2d<float> mat2(4, 3, 7.77777), mat3;
    mat3 = mat1 + sqrt(mat2) * (mat0);
    cout <<mat3;

    mat3 = hypot(mat1, mat2);
    cout <<mat3;

    mat3 = abs(mat2);
    cout <<mat3;
*/
    return 0;
}
