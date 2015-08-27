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
    bool add_noise = params.pop<bool>("add_noise") or false;

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

    // Calculate things
    std::cout << "calculating quantities..." << std::endl;
    auto PpB(ff.calc_Pps(0)); auto PpF(ff.calc_Pps(1));
    auto EpB(ff.calc_Eps(0)); auto EpF(ff.calc_Eps(1));
    auto EB(ff.calc_E(0)); auto EF(ff.calc_E(1));
    auto sign(ff.calc_sign()); auto Pk(ff.calc_Pks());
    auto Pl(ff.calc_Pls()); auto Plm1(ff.calc_Plm1s());
    auto El(ff.calc_Els());

    // Write data to file
    std::cout << "writing to file..." << std::endl;
    utils::write("PpEpB.dat", utils::zip(PpB,EpB), digits);
    utils::write("PpEpF.dat", utils::zip(PpF,EpF), digits);
    utils::write("PpB.dat", PpB, digits);
    utils::write("PpF.dat", PpF, digits);
    utils::write("EpB.dat", EpB, digits);
    utils::write("EpF.dat", EpF, digits);
    utils::write("EB.dat", EB, digits);
    utils::write("EF.dat", EF, digits);
    utils::write("sign.dat", sign, digits);
    utils::write("Pk.dat", Pk, digits);
    utils::write("Pl.dat", Pl, digits);
    utils::write("Plm1.dat", Plm1, digits);
    utils::write("El.dat", El, digits);
    if(add_noise){
        auto noise_ps(PpB);
        for(unsigned p=0; p<noise_ps.size(); p++)
            noise_ps[p] = 1e-100;
        auto noise_ks(Pk);
        for(unsigned k=0; k<noise_ks.size(); k++)
            noise_ks[k] = 1e-100;
        utils::write("PpB.dat", utils::zip(PpB,noise_ps), digits);
        utils::write("EpB.dat", utils::zip(EpB,noise_ps), digits);
        utils::write("Pk.dat", utils::zip(Pk,noise_ks), digits);
    }

    return 0;
}
