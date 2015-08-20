#ifndef FREE_FERMIONS_HPP
#define FREE_FERMIONS_HPP

namespace free_fermions{

/// Class that builds a free fermion system and computes quantities of it.
template<typename RealType>
class FreeFermions{
public:
    /// Constructor computes \Beta and V, and generates permutation sectors, k vectors, Z0s, dZ0dBs, Zps, dZpdBs, ZB, ZF, EB, and EF.
    FreeFermions(const unsigned N, const unsigned D, const RealType L, const RealType T, const RealType lambda, const unsigned n_max, const bool generate_all=true)
        : N_(N), D_(D), L_(L), T_(T), lambda_(lambda), beta_(1/T), V_(pow(L,D))
    {
        std::cout << "...generating permutation sectors..." << std::endl;
        generate_sectors();
        if (generate_all){
            std::cout << "...generating k vectors..." << std::endl;
            std::vector<std::pair<RealType,unsigned>> k2s(generate_k2s(n_max));
            std::cout << "...generating Z0s..." << std::endl;
            generate_Z0s(k2s);
            std::cout << "...generating dZ0dBs..." << std::endl;
            generate_dZ0dBs(k2s);
            std::cout << "...generating Zps..." << std::endl;
            generate_Zps();
            std::cout << "...generating dZpdBs..." << std::endl;
            generate_dZpdBs();
        }
    }

    /// Calculate Pp
    std::vector<RealType> calc_Pps(const bool fermi){
        if(Zps_.size() != 0){
            std::vector<RealType> Pps(Zps_);
            for(unsigned p=0; p<n_sectors_; p++)
                Pps[p] = fermi ? Pps[p]/ZF_ : utils::abs(Pps[p])/ZB_;
            return Pps;
        }else{
            std::cerr << "ERROR: Zps not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate Ep
    std::vector<RealType> calc_Eps(const bool fermi){
        if((Zps_.size()!=0)&&(dZpdBs_.size()!=0)){
            std::vector<RealType> Eps(dZpdBs_);
            for(unsigned p=0; p<n_sectors_; p++)
                Eps[p] = fermi ? Eps[p]/Zps_[p] : utils::abs(Eps[p]/Zps_[p]);
            return Eps;
        }else{
            std::cerr << "ERROR: Zps or dZpdBs not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate total energy
    RealType calc_E(const bool fermi){
        if((Zps_.size()!=0)&&(dZpdBs_.size()!=0)){
            return fermi ? -dZFdB_/ZF_ : dZBdB_/ZB_;
        }else{
            std::cerr << "ERROR: Zps or dZpdBs not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate energy given Eps, Pps and cs
    static std::vector<RealType> calc_energy_from_Eps_Pps(const std::vector<RealType>& Eps, const std::vector<RealType>& Pps, const std::vector<std::vector<unsigned>>& cs){
        RealType num(0), dem(0);
        for(unsigned p=0; p<Pps.size(); p++){
            int sign = calc_sector_sign(cs[p]);
            num += sign*Pps[p]*Eps[p];
            dem += sign*Pps[p];
        }
        return std::vector<RealType>{num/dem};
    }

    /// Calculate sign of sector given cycles
    static int calc_sector_sign(const std::vector<unsigned>& c){
        unsigned n_even = 0;
        for(const auto& l : c)
            if((l+1)%2==0)
                n_even++;
        return (n_even%2==0) ? 1 : -1;
    }

    /// Calculate average sign
    static std::vector<RealType> calc_sign_from_Pps(const std::vector<RealType>& Pps, const std::vector<std::vector<unsigned>>& cs){
        RealType tot(0);
        for(unsigned p=0; p<Pps.size(); p++)
            tot += calc_sector_sign(cs[p])*Pps[p];
        return std::vector<RealType>{tot};
    }

    /// Calculate average sign
    RealType calc_sign(){
        if(Zps_.size() != 0){
            return calc_sign_from_Pps(calc_Pps(0), cs_)[0];
        }else{
            std::cerr << "ERROR: Zps not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate cycle probabilities
    static std::vector<RealType> calc_Pks_from_Pps(const std::vector<RealType>& Pps, const std::vector<std::vector<unsigned>>& cs){
        unsigned N = cs[0].size();
        std::vector<RealType> num(N,0), dem(N,0);
        for(unsigned p=0; p<Pps.size(); p++){
            for(unsigned l=0; l<N; l++){
                num[l] += Pps[p] * cs[p][l];
                dem[l] += Pps[p];
            }
        }
        std::vector<RealType> Pks(N);
        RealType tot(0);
        for(unsigned l=0; l<N; l++){
            Pks[l] = num[l]/dem[l];
            tot += Pks[l];
        }
        for(auto& Pk : Pks) // TODO: Not sure if normalization is needed
            Pk /= tot;
        return Pks;
    }

    /// Calculate cycle probabilities
    std::vector<RealType> calc_Pks(){
        if(Zps_.size() != 0){
            return calc_Pks_from_Pps(calc_Pps(0),cs_);
        }else{
            std::cerr << "ERROR: Zps not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate Pls
    std::vector<RealType> calc_Pls(){
        if(Z0s_.size() != 0){
            return Z0s_;
        }else{
            std::cerr << "ERROR: Z0s not yet generated!" << std::endl;
            exit(1);
        }
    }

    /// Calculate Plm1s
    std::vector<RealType> calc_Plm1s(){
        if(Z0s_.size() != 0){
            std::vector<RealType> Plm1s(Z0s_);
            for (auto& Pl : Plm1s)
                Pl -= 1;
            return Plm1s;
        }else{
            std::cerr << "ERROR: Z0s not yet generated!" << std::endl;
            exit(1);
        }
    }

    const unsigned get_N() { return N_; }
    const unsigned get_n_sectors() { return n_sectors_; }
    const std::vector<std::vector<unsigned>>& get_cs() { return cs_; }

private:
    const unsigned N_; ///< number of particles
    const unsigned D_; ///< dimension of system
    unsigned n_sectors_; ///< number of permutation sectors
    const RealType L_; ///< box side length
    const RealType T_; ///< temperature of system
    const RealType lambda_; ///< \hbar^2/2m
    const RealType beta_; ///< 1/(k_B T)
    const RealType V_; ///< volume of system
    RealType ZB_; ///< bosonic many-particle partition function
    RealType ZF_; ///< fermionic many-particle partition function
    RealType dZBdB_; ///< bosonic many-particle partition function
    RealType dZFdB_; ///< fermionic many-particle partition function
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
    std::vector<std::pair<RealType,unsigned>> generate_k2s(const int n_max){
        std::vector<unsigned> all_n2s;
        for(int n_0=-n_max; n_0<=n_max; n_0++){
            if(D_==1)
                all_n2s.push_back(n_0*n_0);
            else{
                for(int n_1=-n_max; n_1<=n_max; n_1++){
                    if(D_==2)
                        all_n2s.push_back(n_0*n_0 + n_1*n_1);
                    else{
                        for(int n_2=-n_max; n_2<=n_max; n_2++)
                            all_n2s.push_back(n_0*n_0 + n_1*n_1 + n_2*n_2);
                    }
                }
            }
        }
        std::unordered_map<unsigned,unsigned> n2s;
        for(const auto& n2 : all_n2s){
            auto iterator = n2s.find(n2);
            if(iterator == n2s.end()) n2s[n2] = 1;
            else n2s[n2] += 1;
        }
        RealType k_box_2 = pow(2*utils::const_pi<RealType>()/L_,2);
        std::vector<std::pair<RealType,unsigned>> k2s;
        for(const auto& n2 : n2s)
            k2s.push_back(std::make_pair(k_box_2*n2.first,n2.second));
        return std::move(k2s);
    }

    /// Initialize Z0s
    void generate_Z0s(const std::vector<std::pair<RealType,unsigned>>& k2s){
        Z0s_.resize(N_);
        for(unsigned i=0; i<N_; i++){
            RealType i_plus_1_beta_lambda = (i+1)*beta_*lambda_;
            RealType tot(0);
            for(const auto& k2 : k2s)
                tot += k2.second*exp(-i_plus_1_beta_lambda*k2.first);
            Z0s_[i] = tot;
        }
    }

    /// Initialize dZ0/d\Betas
    void generate_dZ0dBs(const std::vector<std::pair<RealType,unsigned>>& k2s){
        dZ0dBs_.resize(N_);
        for(unsigned i=0; i<N_; i++){
            RealType i_plus_1_lambda = (i+1)*lambda_;
            RealType i_plus_1_beta_lambda = i_plus_1_lambda*beta_;
            RealType tot(0);
            for(const auto& k2 : k2s)
                tot += k2.second*(-i_plus_1_lambda*k2.first*exp(-i_plus_1_beta_lambda*k2.first));
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
        ZB_ = 0; ZF_ = 0;
        Zps_.resize(n_sectors_);
        for(unsigned p=0; p<n_sectors_; p++){
            RealType tot(1);
            for(unsigned l=0; l<N_; l++)
                tot *= Zpws[cs_[p][l]][l];
            Zps_[p] = tot;
            ZB_ += utils::abs(tot);
            ZF_ += tot;
        }
    }

    /// Initialize dZpdBs
    void generate_dZpdBs(){
        dZBdB_ = 0; dZFdB_ = 0;
        dZpdBs_.resize(n_sectors_);
        for(unsigned p=0; p<n_sectors_; p++){
            RealType tot(0);
            for(unsigned l=0; l<N_; l++)
                tot += cs_[p][l]*dZ0dBs_[l]/Z0s_[l];
            tot *= Zps_[p];
            dZpdBs_[p] = tot;
            dZBdB_ += utils::abs(tot);
            dZFdB_ += tot;
        }
    }


};

}

#endif
