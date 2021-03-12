#ifndef MATRIX2D_H
#define MATRIX2D_H

#include<iostream>
#include<vector>
#include<cmath>
#include<assert.h>
#include<ostream>

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
        std::vector<T> &operator[](int i);
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

    public:
        int shape[2];
        void print_matrix();
        template <class A> friend std::ostream &operator<<(std::ostream &output, const Matrix2d<A> &m);

    private:
        std::vector<std::vector<T>> mat;

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

/********* Basic matrix operations **********/
/*  重载运算符  */
template<class T> std::vector<T> &Matrix2d<T>::operator[](int i) {
    return mat[i];
}

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
    assert (other.shape[0] == this->shape[0]);
    assert (other.shape[1] == this->shape[1]);

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


#endif // MATRIX2D_H


