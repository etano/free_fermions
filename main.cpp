#include "free_fermions.h"

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
    using RealType = mpfr::mpreal;
    RealType::set_default_prec(mpfr::digits2bits(digits));
    std::cout.precision(digits);    // Show all the digits

    // The rest of the required arguments
    RealType L, T, lambda;
    if(use_egas_units){
        RealType rs = params.pop<RealType>("rs");
        RealType theta = params.pop<RealType>("theta");
        L = egas_units::calc_L(N,D,rs);
        T = egas_units::calc_T(N,D,rs,theta,polarized);
        lambda = egas_units::calc_lambda(D,rs);
    }else{
        L = params.pop<RealType>("L");
        T = params.pop<RealType>("T");
        lambda = params.pop<RealType>("lambda");
    }

    // Create free fermion system
    std::cout << "creating free fermion system..." << std::endl;
    free_fermions::FreeFermions<RealType> ff(N,D,L,T,lambda,n_max);

    // Write data to file
    std::cout << "calculating quantities..." << std::endl;
    utils::write("PpEpB.dat", utils::zip(ff.calc_Pps(0),ff.calc_Eps(0)), digits);
    utils::write("PpEpF.dat", utils::zip(ff.calc_Pps(1),ff.calc_Eps(1)), digits);
    utils::write("PpB.dat", ff.calc_Pps(0), digits);
    utils::write("PpF.dat", ff.calc_Pps(1), digits);
    utils::write("EpB.dat", ff.calc_Eps(0), digits);
    utils::write("EpF.dat", ff.calc_Eps(1), digits);
    utils::write("EB.dat", ff.calc_E(0), digits);
    utils::write("EF.dat", ff.calc_E(1), digits);
    utils::write("sgn.dat", ff.calc_sign(), digits);
    utils::write("Pk.dat", ff.calc_Pks(), digits);
    utils::write("Pl.dat", ff.calc_Pls(), digits);
    utils::write("Plm1.dat", ff.calc_Plm1s(), digits);

    return 0;
}
