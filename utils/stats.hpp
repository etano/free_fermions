#ifndef UTILS_STATS_HPP
#define UTILS_STATS_HPP

namespace utils {

template <class T, typename... Args>
std::vector<T> example(const std::vector<T>& ys, Args&&... args) {
    return std::vector<T>{1, 2};
}

template <class T>
T mean(const std::vector<T>& x, const int N) {
    return std::accumulate(x.begin(), x.end(), T(0)) / N;
}

template <class T>
T mean2(const std::vector<T>& x, const int N) {
    return std::inner_product(x.begin(), x.end(), x.begin(), T(0)) / N;
}

template <class T>
T root_mean2(const std::vector<T>& x, const int N) {
    return sqrt(std::inner_product(x.begin(), x.end(), x.begin(), T(0))) / N;
}

template <class T>
T var(const std::vector<T>& x, const int N) {
    T m(mean(x, N));
    return mean2(x, N) - m * m;
}

template <class T>
T std_dev(const std::vector<T>& x, const int N) {
    return sqrt(var(x, N));
}

template <class T>
T auto_correlation(const std::vector<T>& x, const int t, const T& m, const T& v, const int N) {
    if (!v)
        return 0;
    T tot(0);
    for (int i = 0; i < N - t; ++i)
        tot += (x[i] - m) * (x[N - t - i] - m);
    return tot / (v * (N - t));
}

template <class T>
T kappa(const std::vector<T>& x, const T& m, const T& v, const int N) {
    T tot(0);
    for (int i = 1; i < N; ++i) {
        T c = auto_correlation(x, i, m, v, N);
        if (c <= 0)
            break;
        else
            tot += c;
    }
    return 1 + 2 * tot;
}

template <class T>
std::vector<T> unweighted_avg(const std::vector<T>& means, const std::vector<T>& errors, const std::vector<T>& kappas) {
    int N = means.size();
    return std::vector<T>{mean(means, N), root_mean2(errors, N), mean(kappas, N)};
}

template <class T>
std::vector<T> stats(const std::vector<T>& x) {
    int N = x.size();
    T m(mean(x, N));
    T v(var(x, N));
    T k(kappa(x, m, v, N));
    T e = sqrt(k * v / N);
    if (std::isnan(e))
        e = 0;
    return std::move(std::vector<T>{m, e, k});
}

template <class T, typename F, typename... Args>
std::vector<std::pair<T, T>> bootstrap(const std::vector<T>& ys, const std::vector<T>& y_errs, const unsigned n_bootstrap, utils::rng<T>& rand, F f, Args&&... args) {
    assert(ys.size() == y_errs.size());
    std::vector<std::vector<T>> p_normals;
    for (unsigned i = 0; i < n_bootstrap; i++) {
        std::vector<T> y_normals(ys.size());
        for (unsigned j = 0; j < ys.size(); j++)
            y_normals[j] = rand.norm_rand(ys[j], y_errs[j]);
        p_normals.push_back(f(y_normals, args...));
    }
    std::vector<std::pair<T, T>> means_errors;
    if (!p_normals.empty()) {
        for (unsigned j = 0; j < p_normals[0].size(); j++) {
            std::vector<T> ps;
            for (unsigned i = 0; i < p_normals.size(); i++)
                ps.push_back(p_normals[i][j]);
            auto s(stats(ps));
            means_errors.push_back(std::make_pair(s[0], s[1]));
        }
    }
    return std::move(means_errors);
}
}

#endif
