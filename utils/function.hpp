#ifndef UTILS_FUNCTION_HPP
#define UTILS_FUNCTION_HPP

namespace utils{

template<class T, class U>
std::vector<std::pair<T,U>> zip(const std::vector<T>& a, const std::vector<U>& b){
    assert(a.size() == b.size());
    std::vector<std::pair<T,U>> c;
    for(unsigned i=0; i<a.size(); i++)
        c.push_back(std::make_pair(a[i],b[i]));
    return std::move(c);
}

}

#endif
