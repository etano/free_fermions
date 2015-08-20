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
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<class T>
T running_abs(T sum_so_far, T x){
    return sum_so_far + sqrt(x*x);
}

template<class T>
T opt_abs(T val, bool opt){
    if(opt)
        return sqrt(val*val);
    else
        return val;
}

}

#endif
