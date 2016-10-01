#ifndef FREE_FERMIONS_HPP
#define FREE_FERMIONS_HPP

#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

namespace free_fermions {

/// Parameters for finding the Fermi density
struct fermi_params {
    double theta;
};

/// Function whose root is the Fermi density
double fermi_root(double x, void* params) {
    struct fermi_params* p = (struct fermi_params*)params;
    return gsl_sf_gamma(3./2.) * gsl_sf_fermi_dirac_half(x) - (2./3.) * pow(p->theta, -3./2.);
}

/// Class that builds a free fermion system and computes quantities of it.
template <typename RealType>
class FreeFermions {
    using VecRealType = std::vector<RealType>;
    using VecIntType = std::vector<unsigned>;
    using VecVecRealType = std::vector<std::vector<RealType>>;
    using VecVecIntType = std::vector<std::vector<unsigned>>;

   public:
    /// Constructor computes \Beta and V, and generates permutation sectors, k vectors, Z0s, dZ0dBs, Zps, dZpdBs, ZB, ZF, EB, and EF.
    FreeFermions(const unsigned N, const unsigned D, const RealType L, const RealType T, const RealType T_F, const RealType rs, const RealType lambda, const unsigned n_max, const bool generate_all_quantities, VecIntType& used_sectors)
        : N_(N), D_(D), n_max_(n_max), L_(L), T_(T), T_F_(T_F), rs_(rs), lambda_(lambda), beta_(1 / T), V_(pow(L, D)), used_sectors_(used_sectors), use_all_sectors_(false) {
        std::cout << "...generating permutation sectors..." << std::endl;
        generate_sectors();
        std::cout << "...generating k vectors..." << std::endl;
        std::vector<std::pair<RealType, unsigned>> k2s(generate_k2s(n_max_));
        std::cout << "...generating Z0s..." << std::endl;
        generate_Z0s(k2s);
        std::cout << "...generating dZ0dBs..." << std::endl;
        generate_dZ0dBs(k2s);
        std::cout << "...generating constants..." << std::endl;
        generate_constants();
        if (generate_all_quantities) {
            std::cout << "...generating Zps..." << std::endl;
            generate_Zps();
            std::cout << "...generating dZpdBs..." << std::endl;
            generate_dZpdBs();
        }
    }

    /// Calculate sign of sector given cycles
    static int calc_sector_sign(const std::vector<unsigned>& c) {
        unsigned n_odd = 0;
        for (unsigned l = 0; l < c.size(); l++)
            if (l % 2)
                n_odd += c[l];
        return (!(n_odd % 2)) ? 1 : -1;
    }

    /// Calculate energy given Eps, Pps and cs
    static VecRealType calc_energy_from_Eps_Pps(const VecRealType& Eps, const VecRealType& Pps, const FreeFermions& ff) {
        RealType num(0), dem(0);
        auto& cs = ff.get_cs();
        for (unsigned p = 0; p < Pps.size(); p++) {
            int sign = calc_sector_sign(cs[p]);
            num += sign * Pps[p] * Eps[p];
            dem += sign * Pps[p];
        }
        return VecRealType{num / dem};
    }

    /// Calculate Pps from Pls
    static VecRealType calc_Pps_from_Pls(const VecRealType& Pls, const FreeFermions& ff) {
        // Calculate constants
        auto N = ff.get_N();
        VecVecRealType Ppls(N);
        for (unsigned l = 0; l < N; l++) {
            VecRealType Ppls_l(N + 1);
            for (unsigned i = 0; i <= N; i++)
                Ppls_l[i] = ff.get_inverse_factorial()[i] * pow(Pls[l] / (l + 1), i);
            Ppls[l] = Ppls_l;
        }

        // Calculate Pps
        auto& cs = ff.get_cs();
        VecRealType Pps(cs.size());
        RealType Pp_tot(0);
        for (unsigned p = 0; p < cs.size(); p++) {
            RealType tot(1);
            for (unsigned l = 0; l < N; l++)
                tot *= Ppls[l][cs[p][l]];
            Pps[p] = tot;
            Pp_tot += tot;
        }
        for (unsigned p = 0; p < Pps.size(); p++)
            Pps[p] /= Pp_tot;
        return std::move(Pps);
    }

    /// Calculate Grad Pps from Pls
    static VecVecRealType calc_grad_Pps_from_Pls(const VecRealType& Pls, const VecRealType& Pps, const FreeFermions& ff) {
        // Calculate grad_Pps
        VecVecRealType grad_Pps(Pps.size());
        auto& cs = ff.get_cs();
        for (unsigned p = 0; p < Pps.size(); p++) {
            VecRealType grad_Pps_p(Pls.size());
            for (unsigned l = 0; l < Pls.size(); l++)
                grad_Pps_p[l] = (cs[p][l] * Pps[p] / Pls[l]) * (1 - Pps[p]);
            grad_Pps[p] = grad_Pps_p;
        }
        return std::move(grad_Pps);
    }

    /// Calculate average sign from Pps
    static VecRealType calc_sign_from_Pps(const VecRealType& Pps, const FreeFermions& ff) {
        RealType tot(0);
        auto& cs = ff.get_cs();
        for (unsigned p = 0; p < Pps.size(); p++)
            tot += calc_sector_sign(cs[p]) * Pps[p];
        return VecRealType{tot};
    }

    /// Calculate average sign from Pls
    static VecRealType calc_sign_from_Pls(const VecRealType& Pls, const FreeFermions& ff) {
        return calc_sign_from_Pps(calc_Pps_from_Pls(Pls, ff), ff);
    }

    /// Calculate cycle probabilities
    static VecRealType calc_Pks_from_Pps(const VecRealType& Pps, const FreeFermions& ff) {
        auto N = ff.get_N();
        VecRealType num(N, 0), dem(N, 0);
        auto& cs = ff.get_cs();
        for (unsigned p = 0; p < Pps.size(); p++) {
            for (unsigned l = 0; l < N; l++) {
                num[l] += Pps[p] * cs[p][l];
                dem[l] += Pps[p];
            }
        }
        VecRealType Pks(N);
        RealType tot(0);
        for (unsigned l = 0; l < N; l++) {
            Pks[l] = num[l] / dem[l];
            tot += Pks[l];
        }
        for (unsigned l = 0; l < N; l++)
            Pks[l] /= tot;
        return Pks;
    }

    /// Calculate Pks from Pls
    static VecRealType calc_Pks_from_Pls(const VecRealType& Pls, const FreeFermions& ff) {
        return calc_Pks_from_Pps(calc_Pps_from_Pls(Pls, ff), ff);
    }

    /// Calculate Grad Pks from Pls
    static VecVecRealType calc_grad_Pks_from_Pls(const VecRealType& Pls, const VecRealType& Pks, const FreeFermions& ff) {
        // Calculate Pps
        const auto Pps = calc_Pps_from_Pls(Pls, ff);

        // Calculate (dPk/dPl)/Pk
        VecVecRealType num(Pks.size());
        VecRealType dem(Pks.size());
        for (unsigned k = 0; k < Pks.size(); k++) {
            VecRealType num_k(Pls.size());
            for (unsigned l = 0; l < Pls.size(); l++)
                num_k[l] = 0;
            num[k] = num_k;
            dem[k] = 0;
        }
        auto& cs = ff.get_cs();
        for (unsigned p = 0; p < Pps.size(); p++) {
            for (unsigned k = 0; k < Pks.size(); k++) {
                const auto Pps_p_cs_p_k = Pps[p] * cs[p][k];
                for (unsigned l = 0; l < Pls.size(); l++)
                    num[k][l] += cs[p][l] * Pps_p_cs_p_k / Pls[l];
                dem[k] += Pps_p_cs_p_k;
            }
        }

        // Calculate grad_Pks
        VecVecRealType grad_Pks(Pks.size());
        for (unsigned k = 0; k < Pks.size(); k++) {
            const auto Pk_1_m_Pk_i_dem = Pks[k] * (1 - Pks[k]) / dem[k];
            VecRealType grad_Pks_k(Pls.size());
            for (unsigned l = 0; l < Pls.size(); l++)
                grad_Pks_k[l] = num[k][l] * Pk_1_m_Pk_i_dem;
            grad_Pks[k] = grad_Pks_k;
        }
        return std::move(grad_Pks);
    }

    /// Calculate Eps from Els
    static VecRealType calc_Eps_from_Els(const VecRealType& Els, const FreeFermions& ff) {
        const auto& cs = ff.get_cs();
        VecRealType Eps(cs.size());
        for (unsigned p = 0; p < cs.size(); p++) {
            RealType tot(0);
            for (unsigned l = 0; l < Els.size(); l++)
                tot += cs[p][l] * Els[l] * (l + 1);
            Eps[p] = tot;
        }
        return std::move(Eps);
    }

    /// Calculate Grad Eps from Els
    static VecVecRealType calc_grad_Eps_from_Els(const VecRealType& Els, const VecRealType& Eps, const FreeFermions& ff) {
        const auto& cs = ff.get_cs();
        VecVecRealType grad_Eps(Eps.size());
        for (unsigned p = 0; p < Eps.size(); p++) {
            VecRealType grad_Eps_p(Els.size());
            for (unsigned l = 0; l < Els.size(); l++)
                grad_Eps_p[l] = cs[p][l] * (l + 1);
            grad_Eps[p] = grad_Eps_p;
        }
        return std::move(grad_Eps);
    }

    /// Calculate E from Els
    static VecRealType calc_energy_from_Els_Pps(const VecRealType& Els, const VecRealType& Pps, const FreeFermions& ff) {
        return calc_energy_from_Eps_Pps(calc_Eps_from_Els(Els, ff), Pps, ff);
    }

    /// Calculate Pls from exponential model with parameters ps
    static VecRealType calc_Pls_from_Pl_exp_ps(const VecRealType& ps, const FreeFermions& ff) {
        VecRealType Pls(ff.get_N());
        for (unsigned l = 0; l < Pls.size(); l++) {
            RealType tot(0);
            for (unsigned i = 1; i < ps.size(); i += 2)
                tot += ps[i] * pow(l + 1, ps[i + 1]);
            Pls[l] = 1 + ps[0] * exp(-tot);
        }
        return std::move(Pls);
    }

    // the function is
    // f(l) = d + c*exp(-a*(l+1)**b)
    //      = d + c*exp(-a*exp(ln(l+1)*b))
    // the derivative w.r.t. d
    // df(l)/dd = 1
    // the derivative w.r.t. c
    // df(l)/dc = exp(-a*(l+1)**b)
    //          = (f(l) - d)/c
    // making the derivative w.r.t. a
    // df(l)/da = -c*(l+1)**b * exp(-a*(l+1)**b)
    //          = -c*(l+1)**b * df(l)/dc
    // and the derivative w.r.t. b
    // df(l)/db = -c*a*ln(l+1)*exp(ln(l+1)*b) * exp(-a*(l+1)**b)
    //          = -c*a*ln(l+1)*(l+1)**b * df(l)/dc
    //          = a*ln(l+1) * df(l)/da

    /// Calculate Grad Pls from exponential model with parameters ps
    static VecVecRealType calc_grad_Pls_from_Pl_exp_ps(const VecRealType& ps, const VecRealType& Pls, const FreeFermions& ff) {
        VecVecRealType grad_Pls_ps(Pls.size());
        for (unsigned l = 0; l < Pls.size(); l++) {
            VecRealType grad_Pls_ps_l(ps.size());
            grad_Pls_ps_l[0] = (Pls[l] - 1) / ps[0];
            for (unsigned i = 1; i < ps.size(); i += 2) {
                grad_Pls_ps_l[i] = -ps[0] * pow(RealType(l + 1), ps[i + 1]) * (Pls[l] - 1) / ps[0];
                grad_Pls_ps_l[i + 1] = -ps[0] * ps[i] * log(RealType(l + 1)) * pow(RealType(l + 1), ps[i + 1]) * (Pls[l] - 1) / ps[0];
            }
            grad_Pls_ps[l] = grad_Pls_ps_l;
        }
        return std::move(grad_Pls_ps);
    }

    /// Calculate Pps from exponential model of Pls with parameters ps
    static VecRealType calc_Pps_from_Pl_exp_ps(const VecRealType& ps, const FreeFermions& ff) {
        return calc_Pps_from_Pls(calc_Pls_from_Pl_exp_ps(ps, ff), ff);
    }

    /// Calculate Pps from exponential model of Pls with parameters ps
    static VecVecRealType calc_grad_Pps_from_Pl_exp_ps(const VecRealType& ps, const VecRealType& Pps, const FreeFermions& ff) {
        // Calculate Pls
        auto Pls = calc_Pls_from_Pl_exp_ps(ps, ff);
        // Calculate grad_Pps_Pls
        auto grad_Pps_Pls = calc_grad_Pps_from_Pls(Pls, Pps, ff);
        // Calculate grad_Pls_ps
        auto grad_Pls_ps = calc_grad_Pls_from_Pl_exp_ps(ps, Pls, ff);
        // Calculate grad_Pps_ps
        VecVecRealType grad_Pps_ps(Pps.size());
        for (unsigned p = 0; p < Pps.size(); p++) {
            VecRealType grad_Pps_ps_p(ps.size());
            for (unsigned i = 0; i < ps.size(); i++)
                grad_Pps_ps_p[i] = 0;
            grad_Pps_ps[p] = grad_Pps_ps_p;
        }
        for (unsigned p = 0; p < Pps.size(); p++)
            for (unsigned i = 0; i < ps.size(); i++)
                for (unsigned l = 0; l < Pls.size(); l++)
                    grad_Pps_ps[p][i] += grad_Pls_ps[l][i] * grad_Pps_Pls[p][l];
        return std::move(grad_Pps_ps);
    }

    /// Calculate Pp
    VecRealType calc_Pps(const bool fermi) {
        if (Zps_.size() != 0) {
            VecRealType Pps(Zps_);
            for (unsigned p = 0; p < Pps.size(); p++)
                Pps[p] = fermi ? Pps[p] / ZF_ : utils::abs(Pps[p]) / ZB_;
            return std::move(Pps);
        } else {
            std::cerr << "ERROR: Zps not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate Ep
    VecRealType calc_Eps(const bool fermi) {
        if ((Zps_.size() != 0) && (dZpdBs_.size() != 0)) {
            VecRealType Eps(dZpdBs_);
            for (unsigned p = 0; p < Eps.size(); p++)
                Eps[p] = fermi ? Eps[p] / Zps_[p] : utils::abs(Eps[p] / Zps_[p]);
            return std::move(Eps);
        } else {
            std::cerr << "ERROR: Zps or dZpdBs not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate cycle probabilities
    VecRealType calc_Pks() {
        if (Zps_.size() != 0) {
            return calc_Pks_from_Pps(calc_Pps(0), *this);
        } else {
            std::cerr << "ERROR: Zps not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate Pls
    VecRealType calc_Pls() {
        if (Z0s_.size() != 0) {
            return Z0s_;
        } else {
            std::cerr << "ERROR: Z0s not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate Plm1s
    VecRealType calc_Plm1s() {
        if (Z0s_.size() != 0) {
            VecRealType Plm1s(Z0s_);
            for (unsigned l = 0; l < Plm1s.size(); l++)
                Plm1s[l] -= 1;
            return std::move(Plm1s);
        } else {
            std::cerr << "ERROR: Z0s not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate Els
    VecRealType calc_Els() {
        if ((Z0s_.size() != 0) && (dZ0dBs_.size() != 0)) {
            VecRealType Els(Z0s_.size());
            for (unsigned l = 0; l < Els.size(); l++)
                Els[l] = -dZ0dBs_[l] / (Z0s_[l] * (l + 1));
            return std::move(Els);
        } else {
            std::cerr << "ERROR: Z0s and/or dZ0/dBs not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate average sign
    RealType calc_sign() {
        if (Zps_.size() != 0) {
            return calc_sign_from_Pps(calc_Pps(0), *this)[0];
        } else {
            std::cerr << "ERROR: Zps not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate total energy
    RealType calc_E(const bool fermi) {
        if ((Zps_.size() != 0) && (dZpdBs_.size() != 0)) {
            return fermi ? -dZFdB_ / ZF_ : dZBdB_ / ZB_;
        } else {
            std::cerr << "ERROR: Zps or dZpdBs not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate total energy in the thermodynamic limit
    RealType calc_E_thermo(const bool fermi) {
        if (D_ != 3) {
            std::cerr << "WARNING: Thermodynamic limit data only for 3D systems!" << std::endl;
            return 0.;
        }
        if (!fermi) {
            std::cerr << "WARNING: Thermodynamic limit data only for fermions!" << std::endl;
            return 0.;
        }
        int status;
        int iter = 0, max_iter = 1000;
        const gsl_root_fsolver_type* T;
        gsl_root_fsolver* s;
        double n = 0;
        double x_lo = -10.0, x_hi = 1000.0;
        gsl_function F;
        struct fermi_params params = {double(T_ / T_F_)};
        F.function = &fermi_root;
        F.params = &params;
        T = gsl_root_fsolver_bisection;
        s = gsl_root_fsolver_alloc(T);
        gsl_root_fsolver_set(s, &F, x_lo, x_hi);
        do {
            iter++;
            status = gsl_root_fsolver_iterate(s);
            n = gsl_root_fsolver_root(s);
            x_lo = gsl_root_fsolver_x_lower(s);
            x_hi = gsl_root_fsolver_x_upper(s);
            status = gsl_root_test_interval(x_lo, x_hi, 0, 1.e-8);
        } while (status == GSL_CONTINUE && iter < max_iter);
        gsl_root_fsolver_free(s);
        return pow(rs_, 3) * gsl_sf_gamma(5. / 2.) * gsl_sf_fermi_dirac_3half(n) / (3 * utils::const_pi<RealType>() * pow(0.5 * beta_, 5. / 2.));
    }

    /// Get number of particles
    const unsigned get_N() const { return N_; }
    /// Get number of permutation sectors
    const unsigned get_n_sectors() const { return n_sectors_; }
    /// Get number of used permutation sectors
    const unsigned get_n_used_sectors() const { return used_sectors_.size(); }
    /// Get map of permutation sectors
    const VecVecIntType& get_cs() const { return use_all_sectors_ ? all_cs_ : used_cs_; }
    /// Get inverse factorial
    const VecRealType& get_inverse_factorial() const { return inverse_factorial_; }
    /// Get use_all_sectors flag
    const bool get_use_all_sectors() const { return use_all_sectors_; }
    /// Set use_all_sectors flag
    void set_use_all_sectors(const bool use_all_sectors) { use_all_sectors_ = use_all_sectors; }

   private:
    const unsigned N_;               ///< number of particles
    const unsigned D_;               ///< dimension of system
    const unsigned n_max_;           ///< number of periodic images included
    unsigned n_sectors_;             ///< number of permutation sectors
    const RealType L_;               ///< box side length
    const RealType T_;               ///< temperature of system
    const RealType T_F_;             ///< Fermi temperature of system
    const RealType rs_;              ///< Wigner-Seitz radius
    const RealType lambda_;          ///< \hbar^2/2m
    const RealType beta_;            ///< 1/(k_B T)
    const RealType V_;               ///< volume of system
    RealType ZB_;                    ///< bosonic many-particle partition function
    RealType ZF_;                    ///< fermionic many-particle partition function
    RealType dZBdB_;                 ///< bosonic many-particle partition function
    RealType dZFdB_;                 ///< fermionic many-particle partition function
    VecRealType Z0s_;                ///< map of single particle partition functions
    VecRealType dZ0dBs_;             ///< map of beta derivative of single particle partition functions
    VecRealType Zps_;                ///< map of sector partition functions
    VecRealType dZpdBs_;             ///< map of beta derivative of sector partition functions
    VecRealType inverse_factorial_;  ///< map of inverse factorial
    VecVecIntType cs_;               ///< map of permutation sectors
    VecIntType used_sectors_;        ///< which sectors to use

    VecVecIntType all_cs_; ///< map of all permutation sectors
    VecVecIntType used_cs_; ///< map of used permutation sectors
    bool use_all_sectors_; ///< whether or not to use all sectors

    /// Initialize permutation sectors for a given species
    void generate_sectors(const unsigned sectors_max = 0) {
        std::vector<int> a(N_, 0);
        int k = 1;
        int y = N_ - 1;
        std::vector<std::vector<int>> perms;
        while (k != 0 && (sectors_max > perms.size() || !sectors_max)) {
            int x = a[k - 1] + 1;
            k -= 1;
            while (2 * x <= y) {
                a[k] = x;
                y -= x;
                k += 1;
            }
            int l = k + 1;
            while (x <= y && (sectors_max > perms.size() || !sectors_max)) {
                a[k] = x;
                a[l] = y;
                std::vector<int> b(a.begin(), a.begin() + k + 2);
                perms.push_back(b);
                x += 1;
                y -= 1;
            }
            a[k] = x + y;
            y = x + y - 1;
            std::vector<int> c(a.begin(), a.begin() + k + 1);
            perms.push_back(c);
        }
        n_sectors_ = perms.size();
        unsigned p(0);
        VecIntType included_sectors;
        for (const auto& perm : perms) {
            if ((std::find(used_sectors_.begin(), used_sectors_.end(), p) != used_sectors_.end()) or (used_sectors_.empty())) {
                std::vector<unsigned> c(N_, 0);
                for (const auto& l : perm)
                    c[l - 1]++;
                cs_.push_back(c);
                used_cs_.push_back(c);
                included_sectors.push_back(p);
                all_cs_.push_back(c);
            }
            p++;
        }
        if (used_sectors_.empty())
            used_sectors_ = included_sectors;
    }

    /// Generates vector of k^2
    std::vector<std::pair<RealType, unsigned>> generate_k2s(const int n_max) {
        std::vector<unsigned> all_n2s;
        for (int n_0 = -n_max; n_0 <= n_max; n_0++) {
            if (D_ == 1)
                all_n2s.push_back(n_0 * n_0);
            else {
                for (int n_1 = -n_max; n_1 <= n_max; n_1++) {
                    if (D_ == 2)
                        all_n2s.push_back(n_0 * n_0 + n_1 * n_1);
                    else {
                        for (int n_2 = -n_max; n_2 <= n_max; n_2++)
                            all_n2s.push_back(n_0 * n_0 + n_1 * n_1 + n_2 * n_2);
                    }
                }
            }
        }
        std::unordered_map<unsigned, unsigned> n2s;
        for (const auto& n2 : all_n2s) {
            auto iterator = n2s.find(n2);
            if (iterator == n2s.end())
                n2s[n2] = 1;
            else
                n2s[n2] += 1;
        }
        RealType k_box_2 = pow(2 * utils::const_pi<RealType>() / L_, 2);
        std::vector<std::pair<RealType, unsigned>> k2s;
        for (const auto& n2 : n2s)
            k2s.push_back(std::make_pair(k_box_2 * n2.first, n2.second));
        return std::move(k2s);
    }

    /// Initialize Z0s
    void generate_Z0s(const std::vector<std::pair<RealType, unsigned>>& k2s) {
        Z0s_.resize(N_);
        for (unsigned i = 0; i < N_; i++) {
            RealType i_plus_1_beta_lambda = (i + 1) * beta_ * lambda_;
            RealType tot(0);
            for (const auto& k2 : k2s)
                tot += k2.second * exp(-i_plus_1_beta_lambda * k2.first);
            Z0s_[i] = tot;
        }
    }

    /// Initialize dZ0/d\Betas
    void generate_dZ0dBs(const std::vector<std::pair<RealType, unsigned>>& k2s) {
        dZ0dBs_.resize(N_);
        for (unsigned i = 0; i < N_; i++) {
            RealType i_plus_1_lambda = (i + 1) * lambda_;
            RealType i_plus_1_beta_lambda = i_plus_1_lambda * beta_;
            RealType tot(0);
            for (const auto& k2 : k2s)
                tot += k2.second * (-i_plus_1_lambda * k2.first * exp(-i_plus_1_beta_lambda * k2.first));
            dZ0dBs_[i] = tot;
        }
    }

    /// Generate constants
    void generate_constants() {
        inverse_factorial_.resize(N_ + 1);
        for (unsigned i = 0; i <= N_; i++)
            inverse_factorial_[i] = 1 / utils::factorial<RealType>(i);
    }

    /// Initialize Zps
    void generate_Zps() {
        // Generate Zpws
        VecVecRealType Zpws(N_ + 1);
        for (unsigned i = 0; i <= N_; i++) {
            std::vector<RealType> Zp(N_, 0);
            for (unsigned j = 0; j < N_; j++)
                Zp[j] = inverse_factorial_[i] * pow(pow(-1, j) * Z0s_[j] / (j + 1), i);
            Zpws[i] = Zp;
        }

        // Generate Zps
        ZB_ = 0;
        ZF_ = 0;
        Zps_.resize(cs_.size());
        for (unsigned p = 0; p < Zps_.size(); p++) {
            RealType tot(1);
            for (unsigned l = 0; l < N_; l++)
                tot *= Zpws[cs_[p][l]][l];
            Zps_[p] = tot;
            ZB_ += utils::abs(tot);
            ZF_ += tot;
        }
    }

    /// Initialize dZpdBs
    void generate_dZpdBs() {
        dZBdB_ = 0;
        dZFdB_ = 0;
        dZpdBs_.resize(cs_.size());
        for (unsigned p = 0; p < dZpdBs_.size(); p++) {
            RealType tot(0);
            for (unsigned l = 0; l < N_; l++)
                tot += cs_[p][l] * dZ0dBs_[l] / Z0s_[l];
            tot *= Zps_[p];
            dZpdBs_[p] = tot;
            dZBdB_ += utils::abs(tot);
            dZFdB_ += tot;
        }
    }
};
}

#endif
