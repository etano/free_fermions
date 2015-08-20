#include "main.hpp"

int main(int argc, char** argv)
{
     // Parse arguments
    utils::parameters params(argc,argv);

    // Required arguments
    int N = params.pop<int>("N");
    int D = params.pop<int>("D");
    bool polarized = params.pop<bool>("polarized");

    // Optional arguments
    bool use_egas_units = params.pop<bool>("units") or false;
    int n_max = params.pop<int>("n_max") or 10;
    const int digits = params.pop<int>("precision") or 8;

    // Set precision
    using mpfr::mpreal;
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    std::cout.precision(digits);    // Show all the digits

    // The rest of the required arguments
    mpreal L, T, lambda;
    if(use_egas_units){
        mpreal rs = params.pop<mpreal>("rs");
        mpreal theta = params.pop<mpreal>("theta");
        L = egas_units::calc_L(N,D,rs);
        T = egas_units::calc_T(N,D,rs,theta,polarized);
        lambda = egas_units::calc_lambda(D,rs);
    }else{
        L = params.pop<mpreal>("L");
        T = params.pop<mpreal>("T");
        lambda = params.pop<mpreal>("lambda");
    }

    // Create free fermion system
    FreeFermions<mpreal> ff(N,D,L,T,lambda,n_max);

    // Calculate partition functions
    mpreal ZB = ff.calc_Z(0);
    mpreal ZF = ff.calc_Z(1);

    // Write data to file
    write("PpB.dat", ff.calc_Pps(0,ZB), digits);
    write("PpF.dat", ff.calc_Pps(1,ZF), digits);
    write("EpB.dat", ff.calc_Eps(0,ZB), digits);
    write("EpF.dat", ff.calc_Eps(1,ZF), digits);
    write("sgn.dat", ff.calc_sign(ZB), digits);
    //write("Pk.dat", ff.calc_Pk());
    write("EB.dat", ff.calc_dZdB(0)/ZB, digits);
    write("EF.dat", -ff.calc_dZdB(1)/ZF, digits);

    return 0;
}

