#ifndef UTILS_OSTREAM_HPP
#define UTILS_OSTREAM_HPP

template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vs){
    for(unsigned i=0; i<vs.size(); i++)
        out << i << " " << vs[i] << std::endl;
    return out;
}

template<class T, class U>
std::ostream& operator<<(std::ostream& out, const std::vector<std::pair<T,U>>& vs){
    for(unsigned i=0; i<vs.size(); i++)
        out << i << " " << vs[i].first << " " << vs[i].second << std::endl;
    return out;
}

namespace utils{

template<class T>
void write(const std::string& f, const T val, const unsigned digits=10){
    std::ofstream f_stream;
    f_stream.precision(digits);
    f_stream.open(f);
    f_stream << val;
    f_stream.close();
}

}

#endif
