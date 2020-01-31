/**
 * This file is part of the "dim" project
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

#include <dim/isqrt.h>
#include <dim/mat_expr.h>
#include <dim/value_traits.h>
#include <dim/util.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <ostream>

namespace dim {

/**
 * Represents a dense M-by-N matrix of scalar type F.
 *
 * @param M        number of matrix rows
 * @param N        number of matrix columns
 * @param F        element type (such as int, float, double, std::complex<>)
 * @param Storage  type of the underlying element storage, which currently must match std::array<> semantics.
 */
template <std::size_t M, std::size_t N, typename F, typename Storage = std::array<F, M * N>>
class mat : public mat_expr<M, N, F, mat<M, N, F>> {
  public:
    using value_type = F;

    constexpr static std::size_t row_count = M;
    constexpr static std::size_t column_count = N;

    constexpr mat() = default;
    constexpr mat(mat const&) = default;
    constexpr mat(mat&&) = default;
    constexpr mat& operator=(mat const&) = default;
    constexpr mat& operator=(mat&&) = default;

    /**
     * Constructs a matrix from given initializer list.
     *
     * If not explicitely noted, this will be inferred to by a square matrix with the help
     * of a user-defined deduction guide.
     */
    constexpr mat(std::initializer_list<F> _init) : values_{{}}
    {
        std::size_t u = 0;
        for (auto v : _init)
            values_[u++] = v;
        // XXX: not quite constexpr ;-(
        //std::copy(std::begin(_init), std::end(_init), std::begin(values_));
    }

    /**
     * Constructs a M-by-N matrix with each element initialized by the
     * given initializer function.
     *
     * @param _init The element-initializer function that must return the value for each element [i, j].
     */
    template<
        typename Initializer,
        typename std::enable_if_t<
            std::is_invocable_r_v<F, Initializer, std::size_t, std::size_t>,
            int> = 0
    >
    constexpr explicit mat(Initializer const& _init) noexcept
        : values_{}
    {
        for (auto [i, j] : times(0, M) | times(0, N))
            values_[i * N + j] = _init(i, j);
    }

    /**
     * Constructs a matrix from given matrix expression.
     *
     * @param _expr   The matrix expression to be used to construct this matrix.
     */
    template <typename Mat>
    constexpr mat(mat_expr<M, N, F, Mat> const& _expr) noexcept
        : mat{[&](auto i, auto j) { return _expr(i, j); }}
    {
    }

    /**
     * @returns an immutable reference to the element at given index @p i and @p j (zero based indices).
     */
    constexpr F const& operator()(std::size_t i, std::size_t j) const noexcept { return values_[i * N + j]; }

    /**
     * @returns a mutable reference to the element at given index @p i and @p j (zero based indices).
     */
    constexpr F& operator()(std::size_t i, std::size_t j) { return values_[i * N + j]; }

    constexpr F const& operator()(std::size_t i) const noexcept
    {
        if constexpr (M == 1 || N == 1)
            return values_[i];
        else
            static_assert(true, "WTF?");
    }

  private:
    Storage values_;
};

/**
 * Deduction guide for mat_expr<M, N, F> to mat<M, N, F>
 */
template<
    typename F,
    typename A
>
mat(mat_expr<A::row_count, A::column_count, F, A>)
-> mat<A::row_count, A::column_count, F>;

/**
 * Deduction guide for square matrices.
 */
template<
    typename F,
    typename... F_,
    typename std::enable_if_t<
        std::is_arithmetic_v<F>,
        std::size_t
    > = 0
>
mat(F, F_...) -> mat<isqrt(1 + sizeof...(F_)), isqrt(1 + sizeof...(F_)), F>;

// TODO: diagonals() and other ctor-style functions should return a mat_expr<> instead of mat<M, N, F>{}.
// {{{ free constructor functions
template <std::size_t N, typename F>
constexpr inline auto diagonals(F const& init)
{
    struct Diag : public mat_expr<N, N, F, Diag> {
        F d;
        constexpr explicit Diag(F _d) noexcept : d{_d} {}
        constexpr F operator()(std::size_t i, std::size_t j) const noexcept {
            return i == j ? d : zero<F>;
        }
    };
    return Diag{init};
}

template <std::size_t N, typename F>
constexpr inline void _diagonals(mat<N, N, F>& m, std::size_t d, F const& init)
{
    m(d, d) = init;
}

template<std::size_t N, typename F, typename... F_>
constexpr inline void _diagonals(mat<N, N, F>& m, std::size_t d, F const& init, F_... others)
{
    m(d, d) = init;
    _diagonals(m, d + 1, std::forward<F_>(others)...);
}

template <typename F, typename... F_>
constexpr inline mat<1 + sizeof...(F_), 1 + sizeof...(F_), F>
    diagonals(F const& init, F_... others)
{
    mat<1 + sizeof...(F_), 1 + sizeof...(F_), F> m{};
    _diagonals(m, 0, init, std::forward<F_>(others)...);
    return m;
}

template <std::size_t M, std::size_t N, typename F>
constexpr inline mat<M, N, F> zeros()
{
    struct Zero : public mat_expr<M, N, F, Zero> {
        constexpr F operator()(std::size_t, std::size_t) const noexcept { return zero<F>; }
    };
    return Zero{};
}

template <std::size_t N, typename F>
constexpr auto zeros()
{
    struct Zero : public mat_expr<N, N, F, Zero> {
        constexpr F operator()(std::size_t, std::size_t) const noexcept { return zero<F>; }
    };
    return Zero{};
}

template <std::size_t N, typename F>
constexpr inline auto ones()
{
    struct One : public mat_expr<N, N, F, One> {
        constexpr F operator()(std::size_t, std::size_t) const noexcept { return one<F>; }
    };
    return One{};
}
// }}}

} // namespace dim
