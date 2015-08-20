#ifndef EGAS_UNITS_HPP
#define EGAS_UNITS_HPP

namespace egas_units{

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
            T_F = pow(9*mpfr::const_pi()/2,RealType(2)/RealType(3)) / (2*a*a);
        else
            T_F = pow(9*mpfr::const_pi()/4,RealType(2)/RealType(3)) / (2*a*a);
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

// Temperature in Hartree
template<typename RealType>
RealType calc_T(int N, int D, RealType rs, RealType theta, bool polarized){
    return theta*calc_T_F(D,rs,polarized);
}

// Length of box in Bohr radii
template<typename RealType>
RealType calc_L(int N, int D, RealType rs){
    RealType a = calc_a(rs);
    RealType L;
    if(D==3)
        L = a * pow(4*N*mpfr::const_pi()/3,RealType(1)/RealType(3));
    else if(D==2)
        L = a * sqrt(mpfr::const_pi()*N);
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
