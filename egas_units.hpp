#ifndef EGAS_UNITS_HPP
#define EGAS_UNITS_HPP

namespace egas_units{

// rs
template<typename RealType>
RealType calc_rs(int D, int N, RealType L){
    RealType rs;
    if(D==3){
        rs = pow(N*3/(4*L*L*L*utils::const_pi<RealType>()), 1./3.);
    }else if(D==2){
        rs = sqrt(N/(L*L*utils::const_pi<RealType>()));
    }else{
        std::cerr << "ERROR: D must equal 2 or 3!" << std::endl;
        exit(1);
    }
    return rs;
}

// Wigner seitz radius in Bohr radii
template<typename RealType>
RealType calc_a(RealType rs){
    return rs;
}

// Fermi temperature in Hartree
template<typename RealType>
RealType calc_T_F(int D, RealType rs, bool polarized){
    RealType a = calc_a(rs);
    RealType T_F;
    if(D==3){
        if(polarized)
            T_F = pow(9*utils::const_pi<RealType>()/2,RealType(2)/RealType(3)) / (2*a*a);
        else
            T_F = pow(9*utils::const_pi<RealType>()/4,RealType(2)/RealType(3)) / (2*a*a);
    }else if(D==2){
        if(polarized)
            T_F = 4/(a*a);
        else
            T_F = 1/(a*a);
    }else{
        std::cerr << "ERROR: D must equal 2 or 3!" << std::endl;
        exit(1);
    }
    return T_F;
}

// Fermi temperature in Hartree
template<typename RealType>
RealType calc_T_F(int D, int N, RealType L, bool polarized){
    RealType rs = calc_rs(D,N,L);
    return calc_T_F(D, rs, polarized);
}

// Temperature in Hartree
template<typename RealType>
RealType calc_T(int D, RealType rs, RealType theta, bool polarized){
    return theta*calc_T_F(D,rs,polarized);
}

// Length of box in Bohr radii
template<typename RealType>
RealType calc_L(int N, int D, RealType rs){
    RealType a = calc_a(rs);
    RealType L;
    if(D==3)
        L = a * pow(4*N*utils::const_pi<RealType>()/3,RealType(1)/RealType(3));
    else if(D==2)
        L = a * sqrt(utils::const_pi<RealType>()*N);
    else{
        std::cerr << "ERROR: D must equal 2 or 3!" << std::endl;
        exit(1);
    }
    return L;
}

// \hbar^2/2m
template<typename RealType>
RealType calc_lambda(int D, RealType rs){
    return RealType(1)/RealType(2);
}

} // namespace

#endif
