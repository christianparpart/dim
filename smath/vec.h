/**
 * This file is part of the "smath" project
 *   Copyright (c) 2020 Christian Parpart <christian@parpart.family>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once

#include <algorithm>
#include <type_traits>
#include <array>
#include <cstddef>

namespace smath {

template <typename E, typename F>
struct vec_expr {
    constexpr F operator[](std::size_t i) const noexcept { return static_cast<E const&>(*this)[i]; }
    constexpr std::size_t size() const noexcept { return static_cast<E const&>(*this).size(); }
};

template <typename F, std::size_t Dim>
class vec : public vec_expr<vec<F, Dim>, F> {
  public:
    constexpr vec() = default;
    constexpr vec(vec const&) = default;
    constexpr vec(vec&&) = default;
    constexpr vec& operator=(vec const&) = default;
    constexpr vec& operator=(vec&&) = default;

    template <typename E>
    constexpr vec(vec_expr<E, F> const& _expr) noexcept
    {
        for (std::size_t i = 0; i < Dim; ++i)
            values_[i] = _expr[i];
    }

    template <typename Initializer,
              typename std::enable_if_t<
                  std::is_invocable_r_v<F, Initializer, std::size_t>,
                  int> = 0>
    constexpr explicit vec(Initializer _init) noexcept
    {
        // TODO: how to get you truely constexpr (currently impossible)
        for (std::size_t i = 0; i < Dim; ++i)
            values_[i] = _init(i);
    }

    constexpr explicit vec(F const& _init) noexcept
    {
        // TODO: how to get you truely constexpr (currently impossible)
        std::fill(std::begin(values_), std::end(values_), _init);
    }

    // for initializer_list-style initialization
    // requires(T to support integral operations, and has a zero and one element)
    template <typename... T>
    constexpr explicit vec(T... args) noexcept : values_{std::forward<T>(args)...} {}

    constexpr F operator[](std::size_t _i) const noexcept { return values_[_i]; }
    constexpr F& operator[](std::size_t _i) noexcept { return values_[_i]; }

    constexpr std::size_t size() const noexcept { return Dim; }

    constexpr auto begin() { return std::begin(values_); }
    constexpr auto end() { return std::end(values_); }

    constexpr auto begin() const { return std::begin(values_); }
    constexpr auto end() const { return std::end(values_); }

  private:
    std::array<F, Dim> values_;
};

template <typename F, typename... F_>
vec(F, F_...) -> vec<F, 1 + sizeof...(F_)>;

template <typename F1, typename F2, std::size_t Dim>
constexpr inline bool operator==(vec<F1, Dim> const& a, vec<F2, Dim> const& b)
{
    for (std::size_t i = 0; i < Dim; ++i)
        if (!p(a[i] == b[i]))
            return false;

    return true;
}

template <typename F1, typename F2, std::size_t Dim>
constexpr inline bool operator!=(vec<F1, Dim> const& a, vec<F2, Dim> const& b)
{
    return !(a == b);
}

// {{{ abs(vec)
// XXX: C API provides fabs/fabsl/fabsf, but that's not quite compatible with function overloading and not constexpr either.
constexpr inline double abs(double x) { return x < 0.0 ? -x : x; }
constexpr inline long double abs(long double x) { return x < 0.0 ? -x : x; }
constexpr inline float abs(float x) { return x < 0.0f ? -x : x; }

template <typename E, typename F>
struct _VecAbs : public vec_expr<_VecAbs<E, F>, F>
{
    E const& u;

    constexpr _VecAbs(E const& _u) noexcept : u{_u} {}
    constexpr F operator[](std::size_t i) const noexcept { return abs(u[i]); }
    constexpr std::size_t size() const noexcept { return u.size(); }
};

template <typename E, typename F>
constexpr inline _VecAbs<E, F> abs(vec_expr<E, F> const& v)
{
    return _VecAbs<E, F>{*static_cast<E const*>(&v)};
}
// }}}
// {{{ scalar multiplication
template <typename F, typename E1>
struct _VecScalarMul : vec_expr<_VecScalarMul<F, E1>, F> {
    F scalar_;
    E1 const& u_;

    constexpr F operator[](std::size_t _i) const noexcept { return scalar_ * u_[_i]; }
    constexpr std::size_t size() const noexcept { return u_.size(); }
};

template <typename F, typename E1>
constexpr inline _VecScalarMul<F, E1> operator*(F const& _scalar, vec_expr<E1, F> const& _u) noexcept
{
    return _VecScalarMul<F, E1>{ _scalar, *static_cast<E1 const*>(&_u) };
}

template <typename F, typename E1>
constexpr inline _VecScalarMul<F, E1> operator*(vec_expr<E1, F> const& _u, F const& _scalar) noexcept
{
    return _VecScalarMul<F, E1>{ _scalar, *static_cast<E1 const*>(&_u) };
}
// }}}
// {{{ negation
template <typename E1, typename F>
class _VecNeg : public vec_expr<_VecNeg<E1, F>, F> {
  public:
    constexpr _VecNeg(E1 const& _u) : u_{_u} {}
    constexpr F operator[](std::size_t _i) const noexcept { return -u_[_i]; }
    constexpr std::size_t size() const noexcept { return u_.size(); }

  private:
    E1 const& u_;
};

template <typename E1, typename F>
constexpr inline _VecNeg<E1, F> operator-(vec_expr<E1, F> const& _u) noexcept
{
    return _VecNeg<E1, F>{ *static_cast<E1 const*>(&_u) };
}
// }}}
// {{{ vector addition
template <typename E1, typename E2, typename F>
class _VecAdd : public vec_expr<_VecAdd<E1, E2, F>, F> {
  public:
    constexpr _VecAdd(E1 const& _u, E2 const& _v) : u_{_u}, v_{_v} {}
    constexpr F operator[](std::size_t _i) const noexcept { return u_[_i] + v_[_i]; }
    constexpr std::size_t size() const noexcept { return u_.size(); } // XXX u_.size() == v_.size()

  private:
    E1 const& u_;
    E2 const& v_;
};

template <typename E1, typename E2, typename F>
constexpr inline _VecAdd<E1, E2, F> operator+(vec_expr<E1, F> const& _u, vec_expr<E2, F> const& _v) noexcept
{
    return _VecAdd<E1, E2, F>{
        *static_cast<E1 const*>(&_u),
        *static_cast<E2 const*>(&_v),
    };
}
// }}}
// {{{ vector subtraction
template <typename E1, typename E2, typename F>
class _VecSub : public vec_expr<_VecSub<E1, E2, F>, F> {
  public:
    constexpr _VecSub(E1 const& _u, E2 const& _v) : u_{_u}, v_{_v} {}
    constexpr F operator[](std::size_t _i) const noexcept { return u_[_i] - v_[_i]; }
    constexpr std::size_t size() const noexcept { return u_.size(); } // XXX u_.size() == v_.size()

  private:
    E1 const& u_;
    E2 const& v_;
};

template <typename E1, typename E2, typename F>
constexpr inline _VecSub<E1, E2, F> operator-(vec_expr<E1, F> const& _u, vec_expr<E2, F> const& _v) noexcept
{
    return _VecSub<E1, E2, F>{
        *static_cast<E1 const*>(&_u),
        *static_cast<E2 const*>(&_v),
    };
}
// }}}
// {{{ inner vector product (dot product)
template <typename E1, typename E2, typename F>
constexpr inline F inner_product(vec_expr<E1, F> const& u, vec_expr<E2, F> const& v) noexcept
{
    F x = u[0] * v[0];
    for (std::size_t i = 0; i < u.size(); ++i)
        x += u[i] * v[i];
    return x;
}

template <typename E1, typename E2, typename F>
constexpr inline F operator*(vec_expr<E1, F> const& _u, vec_expr<E2, F> const& _v) noexcept
{
    return inner_product(_u, _v);
}
// }}}
// {{{ outer product (cross product)
template <typename E1, typename E2, typename F>
struct _VecCrossProduct : public vec_expr<_VecCrossProduct<E1, E2, F>, F> {
    E1 const& u;
    E2 const& v;

    constexpr _VecCrossProduct(E1 const& _u, E2 const& _v) : u{_u}, v{_v} {}

    constexpr F operator[](std::size_t i) const noexcept {
        switch (i) {
            case 0: return u[1] * v[2] - u[2] * v[1];
            case 1: return u[2] * v[0] - u[0] * v[2];
            case 2: return u[0] * v[1] - u[1] * v[0];
            default: return F{}; // XXX should never happen
        }
    }

    constexpr std::size_t size() const noexcept { return u.size(); } // XXX u_.size() == v_.size()
};

template <typename E1, typename E2, typename F>
constexpr inline _VecCrossProduct<E1, E2, F> cross_product(vec_expr<E1, F> const& u, vec_expr<E2, F> const& v) noexcept
{
    //static_assert(u.size() == 3);
    //static_assert(v.size() == 3);
    return _VecCrossProduct<E1, E2, F>{
        *static_cast<E1 const*>(&u),
        *static_cast<E2 const*>(&v)
    };
}
// }}}

} // end namespace
