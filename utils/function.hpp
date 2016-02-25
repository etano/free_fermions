#ifndef UTILS_FUNCTION_HPP
#define UTILS_FUNCTION_HPP

namespace utils {

template <class T, class U, class V>
void unzip(const std::unordered_map<T, U>& ab, std::unordered_map<T, V>& b, std::unordered_map<T, V>& c) {
    for (const auto& e : ab) {
        b[e.first] = e.second.first;
        c[e.first] = e.second.second;
    }
}

template <class T, class U>
void unzip(const std::vector<std::pair<T, U>>& ab, std::vector<T>& a, std::vector<U>& b) {
    a.resize(ab.size());
    b.resize(ab.size());
    for (unsigned i = 0; i < ab.size(); i++) {
        a[i] = ab[i].first;
        b[i] = ab[i].second;
    }
}

template <class T, class U, class V>
std::unordered_map<T, std::pair<U, V>> zip(const std::unordered_map<T, U>& a, const std::unordered_map<T, V>& b) {
    assert(a.size() == b.size());
    std::unordered_map<T, std::pair<U, V>> c;
    for (const auto& e : a)
        c[e.first] = std::make_pair(a.at(e.first), b.at(e.first));
    return std::move(c);
}

template <class T, class U, class V>
std::map<T, std::pair<U, V>> zip(const std::map<T, U>& a, const std::map<T, V>& b) {
    assert(a.size() == b.size());
    std::map<T, std::pair<U, V>> c;
    for (const auto& e : a)
        c[e.first] = std::make_pair(a.at(e.first), b.at(e.first));
    return std::move(c);
}

template <class U, class V>
std::vector<std::pair<U, V>> zip(const std::vector<U>& a, const std::vector<V>& b) {
    assert(a.size() == b.size());
    std::vector<std::pair<U, V>> c(a.size());
    for (unsigned i = 0; i < a.size(); i++)
        c[i] = std::make_pair(a[i], b[i]);
    return std::move(c);
}
}

#endif
