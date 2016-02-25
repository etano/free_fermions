#include "free_fermions.h"

int main(int argc, char** argv) {
    // Parse arguments
    utils::parameters params(argc, argv);

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
    std::cout.precision(digits);  // Show all the digits

    // The rest of the required arguments
    RealType L, T, T_F, lambda, rs;
    if (use_egas_units) {
        rs = params.pop<RealType>("rs");
        RealType theta = params.pop<RealType>("theta");
        L = egas_units::calc_L(N, D, rs);
        T = egas_units::calc_T(D, rs, theta, polarized);
        T_F = egas_units::calc_T_F(D, rs, polarized);
        lambda = egas_units::calc_lambda(D, rs);
        std::cout << "L : " << L << std::endl;
        std::cout << "T : " << T << std::endl;
        std::cout << "T_F : " << T_F << std::endl;
        std::cout << "lambda : " << lambda << std::endl;
    } else {
        L = params.pop<RealType>("L");
        T = params.pop<RealType>("T");
        lambda = params.pop<RealType>("lambda");
        T_F = egas_units::calc_T_F(D, N, L, polarized);
        rs = egas_units::calc_rs(D, N, L);
        std::cout << "T_F : " << T_F << std::endl;
        std::cout << "rs : " << rs << std::endl;
    }

    // Adjust if unpolarized
    if (!polarized)
        N /= 2;

    // Create free fermion system
    std::cout << "creating free fermion system..." << std::endl;
    std::vector<unsigned> used_sectors;
    free_fermions::FreeFermions<RealType> ff(N, D, L, T, T_F, rs, lambda, n_max, true, used_sectors);

    // Calculate things
    std::cout << "calculating quantities..." << std::endl;
    auto PpB(ff.calc_Pps(0));
    auto PpF(ff.calc_Pps(1));
    auto EpB(ff.calc_Eps(0));
    auto EpF(ff.calc_Eps(1));
    auto EB(ff.calc_E(0));
    auto EF(ff.calc_E(1));
    auto sign(ff.calc_sign());
    auto Pk(ff.calc_Pks());
    auto Pl(ff.calc_Pls());
    auto Plm1(ff.calc_Plm1s());
    auto El(ff.calc_Els());
    auto EB_thermo(ff.calc_E_thermo(0));
    auto EF_thermo(ff.calc_E_thermo(1));

    // Adjust if unpolarized
    if (!polarized) {
        N *= 2;
        EB *= 2.;
        EF *= 2.;
    }

    // Write data to file
    std::cout << "writing to file..." << std::endl;
    utils::write("PpEpB.dat", utils::zip(PpB, EpB), digits);
    utils::write("PpEpF.dat", utils::zip(PpF, EpF), digits);
    utils::write("PpB.dat", PpB, digits);
    utils::write("PpF.dat", PpF, digits);
    utils::write("EpB.dat", EpB, digits);
    utils::write("EpF.dat", EpF, digits);
    utils::write("EB.dat", EB, digits);
    utils::write("EBN.dat", EB / N, digits);
    utils::write("EB_thermo.dat", EB_thermo, digits);
    utils::write("EF.dat", EF, digits);
    utils::write("EFN.dat", EF / N, digits);
    utils::write("EF_thermo.dat", EF_thermo, digits);
    utils::write("sign.dat", sign, digits);
    utils::write("Pk.dat", Pk, digits);
    utils::write("Pl.dat", Pl, digits);
    utils::write("Plm1.dat", Plm1, digits);
    utils::write("El.dat", El, digits);
    if (add_noise) {
        auto noise_ps(PpB);
        for (unsigned p = 0; p < noise_ps.size(); p++)
            noise_ps[p] = 1e-8;
        auto noise_ks(Pk);
        for (unsigned k = 0; k < noise_ks.size(); k++)
            noise_ks[k] = 1e-8;
        utils::write("PpB.dat", utils::zip(PpB, noise_ps), digits);
        utils::write("EpB.dat", utils::zip(EpB, noise_ps), digits);
        utils::write("Pk.dat", utils::zip(Pk, noise_ks), digits);
    }

    return 0;
}
