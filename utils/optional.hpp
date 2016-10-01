#ifndef WDB_UTILS_OPTIONAL_HPP
#define WDB_UTILS_OPTIONAL_HPP
#include <utility>
#include <functional>
#include <type_traits>
#include <execinfo.h>

namespace utils {
inline void trace_exit(int e) {
    std::cerr << "ERROR: unwrapping empty optional value (" << e << ")\n";
    void* b[15];
    backtrace_symbols_fd(b, backtrace(b, 15), 2);
    exit(1);
}

template <typename R>
class optional_expr;

template <typename R>
class optional;

template <class T>
struct is_optional {
    static bool const value = false;
};

template <class T>
struct is_optional<optional<T>> {
    static bool const value = true;
};

template <class T>
struct is_optional<optional_expr<T>> {
    static bool const value = true;
};

template <typename T>
class optional {
   public:
    friend optional_expr<T>;

    template <typename Arg, typename = typename std::enable_if<!is_optional<Arg>::value && !std::is_convertible<Arg, std::function<T()>>::value>::type>
    optional(Arg arg) : value(arg), valid(true) {}

    template <typename Arg, typename... Other, typename = typename std::enable_if<!is_optional<Arg>::value && !std::is_convertible<Arg, std::function<T()>>::value>::type>
    optional(Arg arg, Other... other) : value(arg, other...), valid(true) {}

    optional() : valid(false) {}

    template <typename T2>
    operator optional<T2>() {
        if (!valid) return optional<T2>();
        return optional<T2>(static_cast<T2>(this->unwrap()));
    }

    operator T() {
        if (!valid) trace_exit(0);
        return *(T*)this;
    }
    operator T() const {
        if (!valid) trace_exit(1);
        return *(T*)this;
    }
    explicit operator bool() {
        return valid;
    }
    optional operator||(const T& b) const {
        if (valid) return *this;
        return b;
    }
    optional operator||(const optional& b) const {
        if (valid) return *this;
        return b;
    }
    optional_expr<T> operator||(const optional_expr<T>& b) const {
        if (valid) return *this;
        return b;
    }
    T unwrap() const {
        if (!valid) trace_exit(2);
        return *(T*)this;
    }
    bool is_null() {
        return !valid;
    }

   private:
    T value;
    bool valid;
};

template <typename R>
class optional_expr {
   public:
    template <typename Func, typename = typename std::enable_if<std::is_convertible<Func, std::function<R()>>::value>::type>
    optional_expr(Func f) : func(f) {}

    optional_expr(const optional<R>& r) {
        this->value = r.value;
        this->valid = r.valid;
    }

    operator R() {
        if (func) {
            this->value = func();
            this->valid = true;
        }
        if (!this->valid) trace_exit(3);
        return *(R*)this;
    }
    template <typename T2>
    operator optional<T2>() {
        return optional<T2>(static_cast<T2>(this->operator R()));
    }

   private:
    R value;
    bool valid;
    std::function<R()> func;
};

template <>
class optional<bool> {
   public:
    friend optional_expr<bool>;
    optional(bool arg) : value(arg), valid(true) {}
    optional() : valid(false) {}

    template <typename T2>
    operator optional<T2>() {
        if (!valid) return optional<T2>();
        return optional<T2>(static_cast<T2>(this->unwrap()));
    }

    operator bool &&() const {
        if (!valid) trace_exit(4);
        return std::move(*(bool*)this);
    }
    explicit operator bool &&() {
        return std::move(valid);
    }
    explicit operator bool() const = delete;

    optional operator||(const bool& b) const {
        if (valid) return *this;
        return b;
    }
    optional operator||(const optional& b) const {
        if (valid) return *this;
        return b;
    }
    optional_expr<bool> operator||(const optional_expr<bool>& b) const {
        if (valid) return *this;
        return b;
    }
    bool&& unwrap() const {
        if (!valid) trace_exit(5);
        return std::move(*(bool*)this);
    }
    bool is_null() {
        return !valid;
    }

   private:
    bool value;
    bool valid;
};
}

#endif
