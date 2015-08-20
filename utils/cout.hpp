#ifndef UTILS_COUT_HPP
#define UTILS_COUT_HPP

template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vs){
    for(unsigned i=0; i<vs.size(); i++)
        out << i << " " << vs[i] << std::endl;
    return out;
}

/// Compare the values of two vectors
template<class T>
struct CompareVec
{
    bool operator() (const std::vector<T> &a, const std::vector<T> &b) {
        for (uint32_t i=0; i<a.size(); i++)
            if (a[i] != b[i])
                return (a[i] > b[i]);
        return (a[0] > b[0]);
    }
};

template<class T, class U>
std::ostream& operator<<(std::ostream& out, const std::map<T,U,CompareVec<U>>& m){
    for(const auto& kvp : m)
        out << kvp.first << " : " << kvp.second << std::endl;
    return out;
}


template<class T>
void write(const std::string& f, const T val, const unsigned digits=10){
    std::ofstream f_stream;
    f_stream.precision(digits);
    f_stream.open(f);
    f_stream << val;
    f_stream.close();
}

#endif
