#ifndef UTILS_FUNCTION_HPP
#define UTILS_FUNCTION_HPP

namespace utils{

template<class T, class U, class V>
void unzip(const std::unordered_map<T,U>& ab, std::unordered_map<T,V>& b, std::unordered_map<T,V>& c){
    for(const auto& e : ab){
        b[e.first] = e.second.first;
        c[e.first] = e.second.second;
    }
}

template<class T, class U, class V>
std::unordered_map<T,std::pair<U,V>> zip(const std::unordered_map<T,U>& a, const std::unordered_map<T,V>& b){
    assert(a.size() == b.size());
    std::unordered_map<T,std::pair<U,V>> c;
    for(const auto& e : a)
        c[e.first] = std::make_pair(a.at(e.first),b.at(e.first));
    return std::move(c);
}

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
