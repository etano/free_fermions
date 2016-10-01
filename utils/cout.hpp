#ifndef UTILS_COUT_HPP
#define UTILS_COUT_HPP

template <class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vs) {
    for (unsigned i = 0; i < vs.size(); i++)
        out << i << " " << vs[i] << std::endl;
    return out;
}

template <class T, class U>
std::ostream& operator<<(std::ostream& out, const std::vector<std::pair<T, U>>& vs) {
    for (const auto& v : vs)
        out << v.first << " : " << v.second << std::endl;
    return out;
}

namespace utils {

template <class T>
void write(const std::string& f, const T val, const unsigned digits = 10) {
    std::ofstream f_stream;
    f_stream.precision(digits);
    f_stream.open(f);
    f_stream << val << "\n";
    f_stream.close();
}
}

#endif
