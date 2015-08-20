#ifndef UTILS_RANDOM_HPP
#define UTILS_RANDOM_HPP

namespace utils{

template<typename RealType>
struct rng
{
    // Constructor
    rng(int seed)
    { mpfr::random(seed); }

    // Random Functions

    // Generate a random number between 0 and 1
    // return a uniform number in [0,1].
    RealType unif_rand() { return mpfr::random(); }

    // Generate a random number in a real interval.
    // param a one end point of the interval
    // param b the other end of the interval
    // return a inform rand numberin [a,b].
    template <typename T>
    RealType unif_rand(const T a, const T b) { return (b-a)*unif_rand() + a; }

    // Generate a random integer between 1 and a given value.
    // param n the largest value return a uniform random value in [1,...,n]
    RealType unif_rand (const int n) { return floor(unif_rand()*n+1); }

    // Generate a normal distribution random number
    RealType norm_rand() { return sqrt(-2*log(mpfr::random()))*cos(2*utils::const_pi<RealType>()*mpfr::random()); }

    // Generate a normal distribution random number of mean m and variance s
    template <typename T>
    RealType norm_rand(const T m, const T s) { return norm_rand()*s + m; }
};

} // namespace

#endif
