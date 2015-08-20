#ifndef UTILS_MATH_HPP
#define UTILS_MATH_HPP

namespace utils{

template<class T>
T factorial(int n){
    return (n == 1 || n == 0) ? 1 : factorial<T>(n - 1) * n;
}
template<>
mpfr::mpreal factorial(int n){
    return mpfr::fac_ui(n);
}

template<class T>
T const_pi(){
    return M_PI;
}
template<>
mpfr::mpreal const_pi(){
    return mpfr::const_pi();
}

template<class T>
int sgn(const T& val){
    return (T(0) < val) - (val < T(0));
}
template<>
int sgn(const mpfr::mpreal& val){
    return mpfr::sgn(val);
}

template<class T>
T abs(const T& val){
    return std::abs(val);
}
template<>
mpfr::mpreal abs(const mpfr::mpreal& val){
    return mpfr::abs(val);
}

}

#endif
