#ifndef FREE_FERMIONS_HPP
#define FREE_FERMIONS_HPP


/// Class that builds a free fermion system and computes quantities of it.
template<typename RealType>
class FreeFermions{
public:
    /// Constructor computes \Beta and V, and generates k vectors
    FreeFermions(unsigned N, unsigned D, RealType L, RealType T, RealType lambda, unsigned n_max)
        : N_(N), D_(D), L_(L), T_(T), lambda_(lambda)
    {
        beta_ = 1/T;
        V_ = pow(L,D);

        std::cout << "generating permutation sectors..." << std::endl;
        generate_sectors();
        std::cout << "generating k vectors..." << std::endl;
        std::vector<RealType> k2s(generate_k2s(n_max));
        std::cout << "generating Z0s..." << std::endl;
        generate_Z0s(k2s);
        std::cout << "generating dZ0dBs..." << std::endl;
        generate_dZ0dBs(k2s);
        std::cout << "generating Zps..." << std::endl;
        generate_Zps();
        std::cout << "generating dZpdBs..." << std::endl;
        generate_dZpdBs();
    }

    /// Calculate Z
    RealType calc_Z(bool fermi){
        RealType tot(0);
        for(const auto& Zp : Zps_)
            tot += utils::opt_abs(Zp,!fermi);
        return tot;
    }

    /// Calculate dZ/d\Beta
    RealType calc_dZdB(bool fermi){
        RealType tot(0);
        for(const auto& dZpdB : dZpdBs_)
            tot += utils::opt_abs(dZpdB,!fermi);
        return tot;
    }

    /// Calculate Pp
    std::vector<RealType> calc_Pps(bool fermi, RealType Z){
        std::vector<RealType> Pps(Zps_);
        for(auto& Pp : Pps)
            Pp = utils::opt_abs(Pp,!fermi)/Z;
        return Pps;
    }

    /// Calculate Ep
    std::vector<RealType> calc_Eps(bool fermi, RealType Z){
        std::vector<RealType> Eps(dZpdBs_);
        for(auto& Ep : Eps)
            Ep = utils::opt_abs(Ep,!fermi)/Z;
        return Eps;
    }

    /// Calculate average sign
    RealType calc_sign(RealType Z){
        RealType tot(0);
        for(const auto& Zp : Zps_)
            tot += utils::sgn(Zp)*utils::opt_abs(Zp,1)/Z;
        return tot;
    }

    /// Calculate cycle probabilities
    std::vector<RealType> calc_Pks(RealType Z){
        std::vector<RealType> num(N_,0), dem(N_,0);
        for(unsigned p=0; p<n_sectors_; p++){
            for(unsigned l=0; l<N_; l++){
                RealType PpB = utils::opt_abs(Zps_[p],1)/Z;
                num[l] += PpB * cs_[p][l];
                dem[l] += PpB;
            }
        }
        std::vector<RealType> Pks(N_);
        RealType tot(0);
        for(unsigned l=0; l<N_; l++){
            Pks[l] = num[l]/dem[l];
            tot += Pks[l];
        }
        for(auto& Pk : Pks)
            Pk /= tot;
        return Pks;
    }

private:
    unsigned N_; ///< number of particles
    unsigned D_; ///< dimension of system
    unsigned n_sectors_; ///< number of permutation sectors
    RealType L_; ///< box side length
    RealType T_; ///< temperature of system
    RealType lambda_; ///< \hbar^2/2m
    RealType beta_; ///< 1/(k_B T)
    RealType V_; ///< volume of system
    std::vector<RealType> Z0s_; ///< vector of single particle partition functions
    std::vector<RealType> dZ0dBs_; ///< vector of beta derivative of single particle partition functions
    std::vector<RealType> Zps_; ///< vector of sector partition functions
    std::vector<RealType> dZpdBs_; ///< vector of beta derivative of sector partition functions
    std::vector<std::vector<unsigned>> cs_; ///< vector of all permutation sectors

    /// Initialize permutation sectors for a given species
    void generate_sectors(const unsigned sectors_max=0)
    {
        std::vector<int> a(N_,0);
        int k = 1;
        int y = N_-1;
        std::vector<std::vector<int>> perms;
        while (k != 0 && (sectors_max > perms.size() || !sectors_max)) {
            int x = a[k-1] + 1;
            k -= 1;
            while (2*x <= y) {
                a[k] = x;
                y -= x;
                k += 1;
            }
            int l = k+1;
            while (x <= y && (sectors_max > perms.size() || !sectors_max)) {
                a[k] = x;
                a[l] = y;
                std::vector<int> b(a.begin(),a.begin()+k+2);
                perms.push_back(b);
                x += 1;
                y -= 1;
            }
            a[k] = x+y;
            y = x+y-1;
            std::vector<int> c(a.begin(),a.begin()+k+1);
            perms.push_back(c);
        }
        cs_.clear();
        n_sectors_ = perms.size();
        for(const auto& p : perms){
            std::vector<unsigned> c(N_,0);
            for(const auto& l : p)
                c[l-1]++;
            cs_.push_back(c);
        }
    }

    /// Generates vector of k^2
    std::vector<RealType> generate_k2s(int n_max){
        std::vector<std::vector<RealType>> ks; ///< vector of k vectors
        for(int n_0=-n_max; n_0<=n_max; n_0++){
            if(D_==1)
                ks.push_back({(RealType)n_0});
            else{
                for(int n_1=-n_max; n_1<=n_max; n_1++){
                    if(D_==2)
                        ks.push_back({(RealType)n_0,(RealType)n_1});
                    else{
                        for(int n_2=-n_max; n_2<=n_max; n_2++)
                            ks.push_back({(RealType)n_0,(RealType)n_1,(RealType)n_2});
                    }
                }
            }
        }
        RealType k_box = 2*utils::const_pi<RealType>()/L_;
        for(auto& k : ks)
            for(auto& k_x : k)
                k_x *= k_box;
        std::vector<RealType> k2s;
        for(auto& k : ks){
            RealType k2(0);
            for(auto& k_x : k)
                k2 += k_x*k_x;
            k2s.push_back(k2);
        }
        return std::move(k2s);
    }

    /// Initialize Z0s
    void generate_Z0s(std::vector<RealType>& k2s){
        Z0s_.resize(N_);
        for(unsigned i=0; i<N_; i++){
            RealType i_plus_1_beta_lambda = (i+1)*beta_*lambda_;
            RealType tot(0);
            for(const auto& k2 : k2s)
                tot += exp(-i_plus_1_beta_lambda*k2);
            Z0s_[i] = tot;
        }
    }

    /// Initialize dZ0/d\Betas
    void generate_dZ0dBs(std::vector<RealType>& k2s){
        dZ0dBs_.resize(N_);
        for(unsigned i=0; i<N_; i++){
            RealType i_plus_1_lambda = (i+1)*lambda_;
            RealType i_plus_1_beta_lambda = i_plus_1_lambda*beta_;
            RealType tot(0);
            for(const auto& k2 : k2s)
                tot += -i_plus_1_lambda*k2*exp(-i_plus_1_beta_lambda*k2);
            dZ0dBs_[i] = tot;
        }
    }

    /// Initialize Zps
    void generate_Zps(){
        // Generate Zpws
        std::vector<std::vector<RealType>> Zpws;
        for(unsigned i=0; i<=N_; i++){
            std::vector<RealType> Zp(N_,0);
            RealType i_factorial_i = 1/utils::factorial<RealType>(i);
            for(unsigned j=0; j<N_; j++)
                Zp[j] = i_factorial_i*pow(pow(-1,j)*Z0s_[j]/(j+1),i);
            Zpws.push_back(Zp);
        }

        // Generate Zps
        Zps_.resize(n_sectors_);
        for(unsigned p=0; p<n_sectors_; p++){
            RealType tot(1);
            for(unsigned l=0; l<N_; l++)
                tot *= Zpws[cs_[p][l]][l];
            Zps_[p] = tot;
        }
    }

    /// Initialize dZpdBs
    void generate_dZpdBs(){
        dZpdBs_.resize(n_sectors_);
        for(unsigned p=0; p<n_sectors_; p++){
            RealType tot(0);
            for(unsigned l=0; l<N_; l++)
                tot += cs_[p][l]*dZ0dBs_[l]/Z0s_[l];
            tot *= Zps_[p];
            dZpdBs_[p] = tot;
        }
    }


};

#endif
