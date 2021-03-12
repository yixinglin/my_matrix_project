#ifndef MATMATH_H
#define MATMATH_H
#include"Matrix2d.h"
#include<cmath>

/**********Math Functions****************/
/*一元运算*/
template<class T> Matrix2d<T> sin(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::sin(a);});
    return ans;
}

template<class T> Matrix2d<T> cos(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::cos(a);});
    return ans;
}

template<class T> Matrix2d<T> tan(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::tan(a);});
    return ans;
}

template<class T> Matrix2d<T> asin(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::asin(a);});
    return ans;
}

template<class T> Matrix2d<T> acos(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::acos(a);});
    return ans;
}

template<class T> Matrix2d<T> atan(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::atan(a);});
    return ans;
}

template<class T> Matrix2d<T> sqrt(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::sqrt(a);});
    return ans;
}

template<class T> Matrix2d<T> abs(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::abs(a);});
    return ans;
}

template<class T> Matrix2d<T> exp(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::exp(a);});
    return ans;
}

template<class T> Matrix2d<T> log(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::log(a);});
    return ans;
}

template<class T> Matrix2d<T> ceil(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::ceil(a);});
    return ans;
}

template<class T> Matrix2d<T> floor(const Matrix2d<T> &m) {
    Matrix2d<T> ans(m);
    ans.unary_operaton([](T a) -> T {return std::floor(a);});
    return ans;
}

/*二元运算*/
// 2次 复制
template<class T> Matrix2d<T> pow(Matrix2d<T> &base, double power) {
    return base.binary_operation(power,
            [](T a, T b) -> T {return std::pow(a, b);});
}

//三角形斜边
template<class T> Matrix2d<T> hypot(Matrix2d<T> &x, Matrix2d<T> &y) {
    return x.binary_operation(y,
            [](T a, T b) -> T {return std::hypot(a, b);});
}

template<class T> Matrix2d<T> atan2(Matrix2d<T> &y, Matrix2d<T> &x) {
    return y.binary_operation(x,
            [](T a, T b) -> T {return std::hypot(a, b);});
}


class matmath
{
    public:
        matmath();
        virtual ~matmath();

    protected:

    private:
};

#endif // MATMATH_H
