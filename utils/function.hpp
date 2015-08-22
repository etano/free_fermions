#ifndef UTILS_FUNCTION_HPP
#define UTILS_FUNCTION_HPP

namespace utils{

template<class T, class U, class V>
std::map<T,std::pair<U,V>> zip(const std::map<T,U>& a, const std::map<T,V>& b){
    assert(a.size() == b.size());
    std::map<T,std::pair<U,V>> c;
    for(const auto& e : a)
        c[e.first] = std::make_pair(a.at(e.first),b.at(e.first));
    return std::move(c);
}

}

#endif
