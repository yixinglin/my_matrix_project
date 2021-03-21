#ifndef MATRIX2D_H
#define MATRIX2D_H

#include<iostream>
#include<vector>
#include<cmath>
#include<assert.h>
#include<ostream>
#include<istream>

#define assertm(exp, msg) assert(((void)msg, exp))

template <class T> class Matrix2d
{
    public:
        /*  构造函数  */
         Matrix2d(){};
         Matrix2d(const Matrix2d<T> &other);
         Matrix2d(int r, int c, T val);
         Matrix2d(std::vector<std::vector<T>> m);
        virtual ~Matrix2d() {}

    public:
        /*  运算符重载  */
        std::vector<T> &operator[] (int i) {return mat[i];};
        /* 一元运算符 */
        void unary_operaton(T (*func)(T));
        Matrix2d<T> &operator++(int n);
        Matrix2d<T> &operator+=(T n);
        Matrix2d<T> &operator-=(T n);
        Matrix2d<T> &operator*=(T n);
        Matrix2d<T> &operator/=(T n);
        Matrix2d<T> &operator+=(const Matrix2d<T> &other);
        Matrix2d<T> &operator-=(const Matrix2d<T> &other);
        Matrix2d<T> &operator*=(const Matrix2d<T> &other);
        Matrix2d<T> &operator/=(const Matrix2d<T> &other);

        /*  %,   */

        /* 二元运算符 */
        Matrix2d<T> binary_operation(const T n, T (*func)(T, T));
        Matrix2d<T> binary_operation(const Matrix2d<T> &other, T (*func)(T, T));
        Matrix2d<T> operator+(T n);
        Matrix2d<T> operator-(T n);
        Matrix2d<T> operator*(T n);
        Matrix2d<T> operator/(T n);
        Matrix2d<T> operator+(const Matrix2d<T> &other);
        Matrix2d<T> operator-(const Matrix2d<T> &other);
        Matrix2d<T> operator*(const Matrix2d<T> &other);
        Matrix2d<T> operator/(const Matrix2d<T> &other);

        // 矩阵赋值，e.g. matrix_obj << 1, 2, 3, 4, 5, 6;
        Matrix2d<T> &operator<<(const T &n);
        Matrix2d<T> &operator,(const T &n);
    public:
        int shape[2];
        void print_matrix();
        template <class A> friend std::ostream &operator<<(std::ostream &output, const Matrix2d<A> &m);
        template <class A> friend std::istream &operator>>(std::istream &input, const Matrix2d<A> &m);


        Matrix2d<T> reshape(int row, int col);
        Matrix2d<T> Tr(); // transpose
        Matrix2d<T> block(int i, int j, int p, int q); // block operation. Block of size (p,q), starting at (i,j)
        T at(int i, int j) const {return mat[i][j];}

    public:
        static Matrix2d<T> eye(int i);
        static Matrix2d<T> hstack(Matrix2d<T> &m, Matrix2d<T> &n);
        static Matrix2d<T> vstack(Matrix2d<T> &m, Matrix2d<T> &n);

    private:
        std::vector<std::vector<T>> mat;
        int input_ptr = 0;
    private:
        void __fill(int r, int c, T val);
};

/* 类模板要写在头文件里面 */
/*********Matrix Construction**********/
template <class T> void Matrix2d<T>::__fill(int r, int c, T val) {
    for(int i=0; i < r; i++) {
        std::vector<T> row(c, val);
        mat.push_back(row);
    }
    shape[0] = r;
    shape[1] = c;
}

template <class T> Matrix2d<T>::Matrix2d(int r, int c, T val) {
    std::cout << "Construct Matrix2d<T>::__fill(int r, int c, T val)"<<std::endl;
    __fill(r, c, val);
}

template <class T> Matrix2d<T>::Matrix2d(const Matrix2d<T> &other) {
    std::cout << "Construct Matrix2d<T>::Matrix2d(const Matrix2d<T> &other)"<<std::endl;
    this->mat = other.mat;  // vector<vector<T>> 是深复制
    shape[0] = other.mat.size();
    shape[1] = other.mat[0].size();
}

template <class T> Matrix2d<T>::Matrix2d(std::vector<std::vector<T>> m) {
    this->mat = m;  // vector 是深复制
    shape[0] = this->mat.size();
    shape[1] = this->mat[0].size();
}

/* Matrix assignment  */
template <class T> Matrix2d<T> &Matrix2d<T>::operator<<(const T &n) {
    mat[0][0] = n;
    input_ptr = 1;
    return *this;
}

template <class T> Matrix2d<T> &Matrix2d<T>::operator,(const T &n) {
     assert(input_ptr < shape[1]*shape[0]);
     int r = input_ptr / shape[1];
     int c = input_ptr % shape[1];
     mat[r][c] = n;
     input_ptr += 1;
     return *this;
}

/********* Display matrix **********/
template<class A> std::ostream &operator<<(std::ostream &output, const Matrix2d<A> &m) {
    output <<std::endl << "[";
    for(int i=0; i < m.shape[0]; i++) {
        output <<"[";
        for(int j=0; j < m.shape[1]; j++) {
            output << m.mat[i][j];
            if(j < m.shape[1]-1) output << ", ";
        }
        output <<"]";
        if(i < m.shape[0]-1) output<<std::endl;

    }
    output <<"]"<<std::endl;
    return output;
}

template <class A> std::istream &operator>>(std::istream &input, const Matrix2d<A> &m) {

    for(int i=0; i<m.shape[0]; ++i) {
        for(int j=0; j<m.shape[1]; ++j) {
            input >> (A)m.mat[i][j];
        }

    }
    return input;
}
/********* Basic matrix operations **********/
/*  重载运算符  */
/*
template<class T> std::vector<T> &Matrix2d<T>::operator[] (int i){
    return mat[i];
}
*/
/* 重载一元运算符 */
template<class T> void Matrix2d<T>::unary_operaton(T (*func)(T)) {
    //typename std::vector<std::vector<T>>::iterator i;  //利用迭代器提高效率
    //typename std::vector<T>::iterator j;    // i and j are pointers
    for(auto i=mat.begin(); i != mat.end(); ++i) {
        for(auto j = i->begin(); j != i->end(); ++j) {
            *j = func(*j);
        }
    }
}

template<class T> Matrix2d<T> &Matrix2d<T>::operator++(int n) {
    unary_operaton([](T a) -> T {return a+1;});
    return *this;  // lambda expression
}

/*2次复制*/
template<class T> Matrix2d<T> &Matrix2d<T>::operator+=(T n) {
    Matrix2d<T> N(shape[0], shape[1], n);
    *this = this->binary_operation(N, [](T a, T b) -> T {return a+b;});  // 复制
    return *this;
}

template<class T> Matrix2d<T> &Matrix2d<T>::operator-=(T n) {
    Matrix2d<T> N(shape[0], shape[1], n);
    *this = this->binary_operation(N, [](T a, T b) -> T {return a-b;});  // 复制
    return *this;
}

template<class T> Matrix2d<T> &Matrix2d<T>::operator*=(T n) {
    Matrix2d<T> N(shape[0], shape[1], n);
    *this = this->binary_operation(N, [](T a, T b) -> T {return a*b;});  // 复制
    return *this;
}

template<class T> Matrix2d<T> &Matrix2d<T>::operator/=(T n) {
    Matrix2d<T> N(shape[0], shape[1], n);
    *this = this->binary_operation(N, [](T a, T b) -> T {return a/b;});  // 复制
    return *this;
}
/*一次复制*/
template<class T> Matrix2d<T> &Matrix2d<T>::operator+=(const Matrix2d<T> &other) {
    // this 地址不改变
    *this = this->binary_operation(other, [](T a, T b) -> T {return a+b;});  // 复制
    return *this;
}

template<class T> Matrix2d<T> &Matrix2d<T>::operator-=(const Matrix2d<T> &other) {
    *this = this->binary_operation(other, [](T a, T b) -> T {return a-b;});  // 复制
    return *this;
}

template<class T> Matrix2d<T> &Matrix2d<T>::operator*=(const Matrix2d<T> &other) {
    *this = this->binary_operation(other, [](T a, T b) -> T {return a*b;});  // 复制
    return *this;
}

template<class T> Matrix2d<T> &Matrix2d<T>::operator/=(const Matrix2d<T> &other) {
    *this = this->binary_operation(other, [](T a, T b) -> T {return a/b;});  // 复制
    return *this;
}

/* 双元运算符 */
template<class T> Matrix2d<T> Matrix2d<T>::binary_operation(const Matrix2d<T> &other, T (*func)(T, T)) {
    assertm (other.shape[0] == this->shape[0],
             "Shapes should be the same");
    assertm (other.shape[1] == this->shape[1],
             "Shapes should be the same");

    Matrix2d<T> ans(other);
    auto i = this->mat.begin();  //利用迭代器提高效率
    auto ans_i = ans.mat.begin();
    typename std::vector<T>::iterator ans_j, j;   // i and j are pointers

    for(; i != mat.end(); ++i, ++ans_i) {
        j = i->begin();
        ans_j = ans_i->begin();
        for(; j != i->end(); ++j, ++ans_j) {
            *ans_j = func(*j, *ans_j);
        }
    }
    return ans;
}

template<class T> Matrix2d<T> Matrix2d<T>::binary_operation(const T n, T (*func)(T, T)) {
    Matrix2d<T> other(shape[0], shape[1], n);
    return this->binary_operation(other, func);
}

/*  一次复制 */
template<class T> Matrix2d<T> Matrix2d<T>::operator+(const Matrix2d<T> &other) {
    return binary_operation(other,
                [](T a, T b) -> T {return a+b;});  // lambda expression
}

template<class T> Matrix2d<T> Matrix2d<T>::operator-(const Matrix2d<T> &other) {
    return binary_operation(other,
                [](T a, T b) -> T {return a-b;});  // lambda expression
}

template<class T> Matrix2d<T> Matrix2d<T>::operator*(const Matrix2d<T> &other) {
    return binary_operation(other,
                [](T a, T b) -> T {return a*b;});  // lambda expression
}

template<class T> Matrix2d<T> Matrix2d<T>::operator/(const Matrix2d<T> &other) {
    return binary_operation(other,
                [](T a, T b) -> T {return a/b;});  // lambda expression
}

template<class T> Matrix2d<T> Matrix2d<T>::operator+(T n) {
    return binary_operation(n,
                [](T a, T b) -> T {return a+b;});  // lambda expression
}

template<class T> Matrix2d<T> Matrix2d<T>::operator-(T n) {
    Matrix2d<T> N(shape[0], shape[1], n);
    return binary_operation(N,
                [](T a, T b) -> T {return a-b;});  // lambda expression
}

template<class T> Matrix2d<T> Matrix2d<T>::operator*(T n) {
    return binary_operation(n,
                [](T a, T b) -> T {return a*b;});  // lambda expression
}

template<class T> Matrix2d<T> Matrix2d<T>::operator/(T n) {
    return binary_operation(n,
                [](T a, T b) -> T {return a/b;});  // lambda expression
}

// transpose
template<class T> Matrix2d<T> Matrix2d<T>::Tr() {
    Matrix2d<T> ans(shape[1], shape[0], 0);
    for(int i=0; i<shape[0]; ++i) {
        for(int j=0; j<shape[1]; ++j) {
            ans[j][i] = mat[i][j];
        }
    }
    return ans;
}

// reshape
template<class T> Matrix2d<T> Matrix2d<T>::reshape(int row, int col) {
    Matrix2d<T> ans(row, col, 0);
    assert(row*col == shape[0]*shape[1]);
    int ans_ptr = 0;
    for(auto i=mat.begin(); i != mat.end(); ++i) {
        for(auto j = i->begin(); j != i->end(); ++j, ++ans_ptr) {
            int r = ans_ptr / col;
            int c = ans_ptr % col;
            ans[r][c] = *j;
        }
    }
    return ans;
}

// block operation,
// Block of size (p,q), starting at (i,j)
template<class T> Matrix2d<T> Matrix2d<T>::block(int i, int j, int p, int q) {
    assertm(i >= 0 && i+p <= shape[0],
            "i exceeded scope");
    assertm(j >= 0 && j+q <= shape[1],
            "j exceeded scope");
    Matrix2d<T> ans(p, q, 0);
    int ans_ptr = 0;
    for(int ii = i; ii<i+p; ++ii) {
        for(int jj = j; jj<j+q; ++jj, ++ans_ptr) {
            int r = ans_ptr / q;
            int c = ans_ptr % q;
            ans[r][c] = mat[ii][jj];
        }

    }
    return ans;
}

template<class T>  Matrix2d<T> Matrix2d<T>::eye(int i) {
    Matrix2d<T> ans(i, i, 0);
    while(i--) {
        ans[i][i] = 1;
    }
    return ans;
}

template<class T>  Matrix2d<T> Matrix2d<T>::hstack(Matrix2d<T> &m, Matrix2d<T> &n) {
    assertm(m.shape[0] == n.shape[0], "hstack");
    Matrix2d<T> ans(m.shape[0], m.shape[1] + n.shape[1], 0);
    for(int i=0; i<ans.shape[0]; ++i) {
        for(int j=0; j<m.shape[1]; ++j) {
            ans[i][j] = m[i][j];
        }

        for(int j=0; j<n.shape[1]; ++j) {
            ans[i][j+m.shape[1]] = n[i][j];
        }

    }
    return ans;
}

template<class T>  Matrix2d<T> Matrix2d<T>::vstack(Matrix2d<T> &m, Matrix2d<T> &n) {
    assertm(m.shape[1] == n.shape[1], "vstack");
    Matrix2d<T> ans(m.shape[0] + n.shape[0], m.shape[1], 0);
    for(int i=0; i<ans.shape[1]; ++i) {
        for(int j=0; j<m.shape[0]; ++j) {
            ans[j][i] = m[j][i];
        }

        for(int j=0; j<n.shape[0]; ++j) {
            ans[j+m.shape[0]][i] = n[j][i];
        }

    }
    return ans;
}

/*
template<class T>
Matrix2d<T> matmul(const Matrix2d<T> &m1,  const Matrix2d<T> &m2) {
    assertm(m1.shape[1] == m2.shape[0], "Matrix Shape: matmul");

    Matrix2d<T> ans(m1.shape[0], m2.shape[1], 0);
    for(int i=0; i<m1.shape[0]; ++i) {
        for(int j=0; j<m2.shape[1]; ++j){
            for(int k=0; k<m2.shape[0]; ++k) {
                ans[i][j] += (m1[i][k] * m2[k][j]);

            }
        }
    }
    return ans;
}
*/

#endif // MATRIX2D_H



