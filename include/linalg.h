#ifndef LINALG_H
#define LINALG_H

#include "Matrix2d.h"
/********* 矩阵操作 **********/
template<class T> Matrix2d<T> matmul(Matrix2d<T> &m1,  Matrix2d<T> &m2) {
    assertm(m1.shape[1] == m2.shape[0], "Matrix Shape: matmul");
    Matrix2d<T> ans(m1.shape[0], m2.shape[1], 0);
    for(int i=0; i<m1.shape[0]; ++i) {
        for(int j=0; j<m2.shape[1]; ++j){
            for(int k=0; k<m2.shape[0]; ++k) {
                ans[i][j] += (m1[i][k] + m2[k][j]);
            }
        }
    }
    return ans;
}

/* 两行相加: line(j) = line(j) + s * line(i)  */

template<typename T> void mat_lines_add(Matrix2d<T> &m, int i, int j, double s) {
    for(int k=0; k<m.shape[1]; ++k) {
        m[j][k] += s*m[i][k];
    }
}

/*   交换两行   */
template<typename T> void mat_lines_swap(Matrix2d<T> &m, int i, int j) {
    swap(m[j], m[i]);
}

/*  上三角阵 upper triangular matrix */
template<typename T> int mat_triu(Matrix2d<T> &m) {
    assertm(m.shape[0] == m.shape[1], "Matrix Shape");
    int q = 1;
    for(int j=0; j<m.shape[1]-1; ++j) {  //水平方向
        if(m[j][j] == 0) {  //尝试交换行，保证对角线没有0
            int k=j+1;
            while(k < m.shape[0] && m[k][j] == 0) {k++;}  //寻找非零元素
            if(k < m.shape[0]) {mat_lines_swap(m, k, j); q = -q;} else {continue;} //找不到非零元素，整一列都是0
        }

        for(int i=j+1; i < m.shape[0]; ++i) {  //垂直方向
            if(m[i][j] == 0) continue;
            double s = m[i][j]/m[j][j];
            mat_lines_add(m, j, i, -s);
        }

    }
    return q;

}
/*  det  */
// 先三角阵，然后求det
template<typename T> double det(Matrix2d<T> m) {
    assertm(m.shape[0] == m.shape[1], "Matrix Shape");
    int q = mat_triu(m);
    for(int i=0; i<m.shape[0]; ++i) {
        q *= m[i][i];
    }
    return q;
}

/*  逆 */
//高斯消元法


#endif // LINALG_H
