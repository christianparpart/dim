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

#include <dim/mat_expr.h>
#include <dim/mat.h>
#include <dim/isqrt.h>
#include <dim/util.h>
#include <dim/value_traits.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <variant>

namespace dim {

// {{{ matrix comparison
template<
    std::size_t M,
    std::size_t N,
    typename F,
    //typename F2,
    typename A,
    typename B
>
constexpr inline bool operator==(mat_expr<M, N, F, A> const& a, mat_expr<M, N, F, B> const& b) noexcept
{
    for (auto [i, j] : times(0, M) | times(0, N))
        if (a(i, j) != b(i, j))
            return false;

    return true;
}

template<
    std::size_t M,
    std::size_t N,
    typename F,
    //typename F2,
    typename A,
    typename B
>
constexpr inline bool operator!=(mat_expr<M, N, F, A> const& a, mat_expr<M, N, F, B> const& b) noexcept
{
    return !(a == b);
}
// }}}
// {{{ init
template<
    std::size_t M,
    std::size_t N,
    typename F,
    typename Initializer,
    typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0,
    typename std::enable_if_t<std::is_invocable_r_v<F, Initializer, std::size_t, std::size_t>, int> = 0
>
constexpr inline auto init(Initializer _init)
{
    struct Init : public mat_expr<M, N, F, Init> {
        Initializer initializer;
        constexpr Init(Initializer _init) noexcept : initializer{std::move(_init)} {}
        constexpr F operator()(std::size_t i, std::size_t j) const { return initializer(i, j); }
    };
    return Init{static_cast<Initializer const&>(_init)};
}
// }}}
// {{{ matrix transform(A) -> B
// Transforms two matrices into one by applying given binary operator on each coefficient.
template<
    std::size_t M,
    std::size_t N,
    typename F,
    typename A,
    typename UnaryOp,
    typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0,
    typename std::enable_if_t<std::is_invocable_r_v<F, UnaryOp, F>, int> = 0
>
constexpr inline auto transform(mat_expr<M, N, F, A> const& a, UnaryOp const& op)
{
    //typename std::enable_if_t<std::is_invocable_r_v<F, Initializer, std::size_t, std::size_t>, int> = 0
    struct Mapping : public mat_expr<M, N, F, Mapping> {
        A const& a;
        UnaryOp op;
        constexpr Mapping(A const& _a, UnaryOp _op) noexcept : a{_a}, op{std::move(_op)} {}
        constexpr F operator()(std::size_t i, std::size_t j) const { return op(a(i, j)); }
    };

    return Mapping{static_cast<A const&>(a), op};
}
// }}}
// {{{ matrix transform(A, B) -> C
template<
    std::size_t M,
    std::size_t N,
    typename F,
    typename A,
    typename B,
    typename BinaryOp,
    typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0,
    typename std::enable_if_t<std::is_invocable_r_v<F, BinaryOp, F, F>, int> = 0
>
constexpr inline auto transform(mat_expr<M, N, F, A> const& a,
                                mat_expr<M, N, F, B> const& b,
                                BinaryOp const& op)
{
    struct Mapping : public mat_expr<M, N, F, Mapping> {
        A const& a;
        B const& b;
        BinaryOp op;

        constexpr Mapping(A const& _a, B const& _b, BinaryOp _op) noexcept : a{_a}, b{_b}, op{std::move(_op)} {}
        constexpr F operator()(int i, int j) const { return op(a(i, j), b(i, j)); }
    };

    return Mapping{static_cast<A const&>(a), static_cast<B const&>(b), op};
}
// }}}
// {{{ abs(vec)
template<
    std::size_t M,
    std::size_t N,
    typename F,
    typename A,
    typename std::enable_if_t<std::is_arithmetic_v<F>, std::size_t> = 0
>
constexpr inline auto abs(mat_expr<M, N, F, A> const& a) noexcept
{
    // XXX: C API provides fabs/fabsl/fabsf, but that's not quite compatible with function overloading and not constexpr either.
    struct Abs : public mat_expr<M, N, F, Abs>
    {
        A const& a;
        constexpr Abs(A const& _a) noexcept : a{_a} {}
        constexpr F operator()(std::size_t i, std::size_t j) const noexcept {
            auto && v = a(i, j);
            return v < zero<F> ? -v : v;
        }
    };

    return Abs{static_cast<A const&>(a)};
}
// }}}
// {{{ matrix addition / subtraction
template <std::size_t M, std::size_t N, typename F, typename A, typename B>
constexpr inline auto operator+(mat_expr<M, N, F, A> const& a,
                                mat_expr<M, N, F, B> const& b)
{
    return transform(a, b, [](F const& x, F const& y) constexpr { return x + y; });
}

template <std::size_t M, std::size_t N, typename F, typename A, typename B>
constexpr inline auto operator-(mat_expr<M, N, F, A> const& a,
                                mat_expr<M, N, F, B> const& b)
{
    return transform(a, b, [](F const& x, F const& y) constexpr { return x - y; });
}

template <std::size_t M, std::size_t N, typename F, typename A>
constexpr inline auto operator-(mat_expr<M, N, F, A> const& a)
{
    return transform(a, [](F const& v) { return -v; });
}
// }}}
// {{{ scalar multiplication
template<
    std::size_t M,
    std::size_t N,
    typename F,
    typename A,
    typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0,
    typename std::enable_if_t<std::is_base_of_v<mat_expr<M, N, F, A>, A>, int> = 0
>
constexpr inline auto operator*(F s, mat_expr<M, N, F, A> const& a)
{
    return transform(a, [s](F const& v) { return s * v; });
}
// }}}
// {{{ matrix multiplication
template<
    std::size_t M,
    std::size_t N,
    std::size_t K,
    typename F,
    typename A,
    typename B,
    typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0,
    typename std::enable_if_t<std::is_base_of_v<mat_expr<M, N, F, A>, A>, int> = 0,
    typename std::enable_if_t<std::is_base_of_v<mat_expr<N, K, F, B>, B>, int> = 0
>
constexpr inline auto operator*(mat_expr<M, N, F, A> const& a,
                                mat_expr<N, K, F, B> const& b)
{
    struct MatMult : public mat_expr<M, K, F, MatMult> {
        A const& a;
        B const& b;
        constexpr MatMult(A const& _a, B const& _b) noexcept : a{_a}, b{_b} {}
        constexpr F operator()(std::size_t i, std::size_t j) const {
            return reduce(
                times(0u, N),
                zero<F>,
                [this, i, j](auto acc, auto l) constexpr -> F { return acc + a(i, l) * b(l, j); }
            );
        }
    };

    return MatMult{static_cast<A const&>(a), static_cast<B const&>(b)};
}
// }}}
// {{{ other free functions
template <std::size_t M, std::size_t N, typename F, typename A>
constexpr std::size_t weigh_trivials_in_row(mat_expr<M, N, F, A> const& m, size_t i)
{
    std::size_t c = 0;
    for (std::size_t j : times(0u, N))
    {
        F const v = m(i, j);
        F const w = v < zero<F> ? -v : v;
        if (w == zero<F> || w == one<F> || w == two<F>)
            c += three<F> - w;
    }
    return c;
}

template <std::size_t M, std::size_t N, typename F, typename A>
constexpr std::size_t weigh_trivials_in_column(mat_expr<M, N, F, A> const& m, size_t j)
{
    std::size_t c = 0;
    for (std::size_t i : times(0u, N))
    {
        F const v = m(i, j);
        F const w = v < zero<F> ? -v : v;
        if (w == zero<F> || w == one<F> || w == two<F>)
            c += three<F> - w;
    }

    return c;
}

template <std::size_t M, std::size_t N, typename F, typename A>
constexpr size_t count_zeros(mat_expr<M, N, F, A> const& m)
{
    std::size_t a = 0;
    for (auto [i, j] : times(0u, M) | times(0u, N))
        if (m(i, j) == zero<F>)
            ++a;

    return a;
}

template <std::size_t M, std::size_t N, typename F, typename A>
constexpr size_t count_ones(mat_expr<M, N, F, A> const& m)
{
    std::size_t a = 0;
    for (auto [i, j] : times(0u, M) | times(0u, N))
        if (m(i, j) == one<F>)
            ++a;
    return a;
}

template <std::size_t M, std::size_t N, typename F, typename A>
constexpr float density(mat_expr<M, N, F, A> const& m)
{
    auto const total = M * N;
    return static_cast<float>(total - count_zeros(m)) / total;
}

template <std::size_t M, std::size_t N, typename F, typename A>
constexpr float sparsity(mat_expr<M, N, F, A> const& m)
{
    return 1.0 - density(m);
}

template<
    std::size_t N,
    typename F,
    typename A
>
constexpr inline F trace(mat_expr<N, N, F, A> const& a)
{
    return reduce(
        times(0u, N),
        zero<F>,
        [&](F const& acc, std::size_t k) constexpr -> F { return acc + a(k, k); }
    );
}

// TODO: how to avoid doubling the algorithms with a cref/&& version?
template <std::size_t M, std::size_t N, typename F, typename A>
constexpr inline auto transpose(mat_expr<M, N, F, A> const& a)
{
    struct Transpose : public mat_expr<N, M, F, Transpose>
    {
        A const& a;
        constexpr Transpose(A const& _a) noexcept : a{_a} {}
        constexpr F operator()(std::size_t i, std::size_t j) const noexcept { return a(j, i); }
    };

    return Transpose{static_cast<A const&>(a)};
}

template <std::size_t M, std::size_t N, typename F, typename A>
constexpr inline auto transpose(mat_expr<M, N, F, A>&& a)
{
    struct Transpose : public mat_expr<N, M, F, Transpose>
    {
        A a;
        constexpr Transpose(A _a) noexcept : a{std::move(_a)} {}
        constexpr F operator()(std::size_t i, std::size_t j) const noexcept { return a(j, i); }
    };

    return Transpose{std::move(static_cast<A const&>(a))};
}

template<
    std::size_t M,
    std::size_t N,
    typename F,
    typename A
>
constexpr inline auto minor(mat_expr<M, N, F, A> const& a, std::size_t i, std::size_t j)
{
    struct Minor : public mat_expr<M - 1, N - 1, F, Minor>
    {
        A const& a;
        std::size_t const i;
        std::size_t const j;

        constexpr Minor(A const& _a, std::size_t _i, std::size_t _j) noexcept : a{_a}, i{_i}, j{_j} {}

        constexpr F operator()(std::size_t k, std::size_t l) const noexcept
        {
            return k < i
                ? l < j
                    ? a(k, l)
                    : a(k, l + 1)
                : l < j
                    ? a(k + 1, l)
                    : a(k + 1, l + 1);
        }
    };
    return Minor{static_cast<A const&>(a), i, j};
}
// }}}
// {{{ det(A)
template <typename F, typename A>
constexpr inline F det(mat_expr<1, 1, F, A> const& a) noexcept
{
    return a(0, 0);
}

template <typename F, typename A>
constexpr inline F det(mat_expr<2, 2, F, A> const& a) noexcept
{
    return a(0, 0) * a(1, 1) - a(1, 0) * a(0, 1);
}

template <typename F, typename A>
constexpr inline F det(mat_expr<3, 3, F, A> const& a) noexcept
{
    return a(0, 0) * a(1, 1) * a(2, 2)
         + a(1, 0) * a(2, 1) * a(0, 2)
         + a(2, 0) * a(0, 1) * a(1, 2)
         - a(2, 0) * a(1, 1) * a(0, 2)
         - a(2, 1) * a(1, 2) * a(0, 0)
         - a(2, 2) * a(1, 0) * a(0, 1);
}

/**
 * Finds the best matrix row or column to expand by.
 *
 * @returns a pair whose first argument is the row or column index and the second argument
 *          determines whether or not the first one is a row index (column index otherwise).
 */
template <std::size_t N, typename F, typename A>
constexpr inline std::pair<std::size_t, bool> find_expansion_vector(mat_expr<N, N, F, A> const& a) noexcept
{
    std::size_t expansion_index = 0;
    std::size_t most_trivials = weigh_trivials_in_row(a, 0);

    for (std::size_t k : times(1u, N - 1))
        if (auto const c = weigh_trivials_in_row(a, k); c > most_trivials)
        {
            expansion_index = k;
            most_trivials = c;
        }

    bool expandByRow = true;
    for (std::size_t k : times(0u, N))
        if (auto const c = weigh_trivials_in_column(a, k); c > most_trivials)
        {
            expandByRow = false;
            expansion_index = k;
            most_trivials = c;
        }

    return {expansion_index, expandByRow};
}

// det(A) via laplace expansion
template <std::size_t N, typename F, typename A>
constexpr inline F det(mat_expr<N, N, F, A> const& a) noexcept
{
    auto const [expansion_index, expandByRow] = find_expansion_vector(a);

    if (expandByRow)
    {
        std::size_t const i = expansion_index;
        return reduce(
            times(0u, N),
            zero<F>,
            [&](F const& acc, std::size_t j) constexpr -> F {
                return (i + j) % 2 == 0
                    ? acc + a(i, j) * det(minor(a, i, j))
                    : acc - a(i, j) * det(minor(a, i, j));
            }
        );
    }
    else // expand by row
    {
        std::size_t const j = expansion_index;
        return reduce(
            times(0u, N),
            zero<F>,
            [&](F const& acc, std::size_t i) constexpr -> F {
                return (i + j) % 2 == 0
                    ? acc + a(i, j) * det(minor(a, i, j))
                    : acc - a(i, j) * det(minor(a, i, j));
            }
        );
    }
}
// }}}
// {{{ outer product (cross product)
template <typename F, typename A, typename B>
constexpr inline auto cross_product(mat_expr<3, 1, F, A> const& u,
                                    mat_expr<3, 1, F, B> const& v) noexcept
{
    struct CrossProduct : public mat_expr<3, 1, F, CrossProduct> {
        A const& u;
        B const& v;
        constexpr CrossProduct(A const& _u, B const& _v) : u{_u}, v{_v} {}
        constexpr F operator()(std::size_t i, std::size_t j) const noexcept {
            switch (i + j - 1) {
                case 0: return u[1] * v[2] - u[2] * v[1];
                case 1: return u[2] * v[0] - u[0] * v[2];
                case 2: return u[0] * v[1] - u[1] * v[0];
                default: return F{}; // XXX should never happen
            }
        }
    };
    return CrossProduct{static_cast<A const&>(u), static_cast<B const&>(v)};
}
// }}}
// {{{ column
/**
 * Retrieves a given column-vector @p n from matrix @p mat.
 */
template <std::size_t M, std::size_t N, typename F, typename A>
auto column(std::size_t n, mat_expr<M, N, F, A> const& mat)
{
    struct Column : public mat_expr<M, 1, F, Column> {
        std::size_t n;
        A const& mat;
        constexpr Column(std::size_t _n, A const& _mat) noexcept : n{_n}, mat{_mat} {}
        constexpr F operator()(std::size_t i, std::size_t /*j*/) const noexcept { return mat(i, n); }
    };
    return Column{n, static_cast<A const&>(mat)};
}
// }}}
// {{{ elementary operations
template<
    std::size_t M,
    std::size_t N,
    typename F,
    typename A,
    typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0
>
constexpr auto swap_row(mat_expr<M, N, F, A> const& mat, std::size_t a, std::size_t b)
{
    assert(a < mat.row_count);
    assert(b < mat.row_count);

    return init<M, N, F>([a, b, mat = std::ref(mat)](std::size_t i, std::size_t j) constexpr -> F {
        return i == a ? mat(b, j)
             : i == b ? mat(a, j)
                      : mat(i, j);
    });
}

template<
    std::size_t M,
    std::size_t N,
    typename F,
    typename A,
    typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0
>
constexpr auto scale_row(mat_expr<M, N, F, A> const& mat, std::size_t row, F s)
{
    assert(row < mat.row_count);

    return init<M, N, F>([row, s, mat = std::ref(mat)](std::size_t i, std::size_t j) constexpr -> F {
        return i == row ? mat(i, j) * s
                        : mat(i, j);
    });
}

template<
    std::size_t M,
    std::size_t N,
    typename F,
    typename A,
    typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0
>
constexpr auto add_scaled_row(mat_expr<M, N, F, A> const& mat, std::size_t targetRow, F s, std::size_t row)
{
    assert(row < mat.row_count);
    assert(targetRow < mat.row_count);

    return init<M, N, F>([row, targetRow, s, mat = std::ref(mat)](std::size_t i, std::size_t j) constexpr {
        return i == targetRow
            ? mat(i, j) + s * mat(row, j)
            : mat(i, j);
    });
}

template <typename F>
constexpr inline F kronecker_delta(std::size_t i, std::size_t j) noexcept
{
    return i == j ? one<F> : zero<F>;
}

namespace elementary {
    struct swap_row_t { std::size_t a; std::size_t b; };

    template <typename F>
    struct scale_row_t { std::size_t row; F scalar; };

    template <typename F>
    struct add_scaled_row_t {
        std::size_t target_row;
        F scalar;
        std::size_t row;
    };

    template <typename F>
    using operation = std::variant<
        swap_row_t,
        scale_row_t<F>,
        add_scaled_row_t<F>
    >;

    template <typename F>
    constexpr auto swap_row(std::size_t a, std::size_t b) noexcept {
        return operation<F>{swap_row_t{a, b}};
    }

    template <typename F>
    constexpr auto scale_row(std::size_t row, F scalar) noexcept {
        return operation<F>{scale_row_t<F>{row, scalar}};
    }

    template <typename F>
    constexpr auto add_scaled_row(std::size_t targetRow, F scalar, std::size_t row) noexcept {
        return operation<F>{add_scaled_row_t<F>{targetRow, scalar, row}};
    }

    template<
        std::size_t M,
        std::size_t N,
        typename F,
        typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0
    >
    constexpr auto matrix(operation<F> op)
    {
        struct Ops {
            std::size_t i;
            std::size_t j;
            constexpr Ops(std::size_t _i, std::size_t _j) noexcept : i{_i}, j{_j} {}
            constexpr auto operator()(elementary::swap_row_t const& op) noexcept {
                return i == op.a ? kronecker_delta<F>(op.b, j)
                     : i == op.b ? kronecker_delta<F>(op.a, j)
                                 : kronecker_delta<F>(i, j);
            }
            constexpr auto operator()(elementary::scale_row_t<F> const& op) noexcept {
                return i == op.row ? kronecker_delta<F>(i, j) * op.scalar
                                   : kronecker_delta<F>(i, j);
            }
            constexpr auto operator()(elementary::add_scaled_row_t<F> const& op) noexcept {
                return i == op.target_row && j == op.row ? kronecker_delta<F>(j, j) * op.scalar
                                                         : kronecker_delta<F>(i, j);
            }
        };
        struct Elements : public mat_expr<M, N, F, Elements> {
            elementary::operation<F> op;
            constexpr explicit Elements(elementary::operation<F> _op) noexcept : op{std::move(_op)} {}
            constexpr F operator()(std::size_t i, std::size_t j) const noexcept { return std::visit(Ops{i, j}, op); }
        };
        return Elements{std::move(op)};
    }

    /**
     * Applies given elementary operation @p op to matrix @p mat.
     */
    template<
        std::size_t M,
        std::size_t N,
        typename F,
        typename A,
        typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0
    >
    constexpr auto apply(elementary::operation<F> op, mat_expr<M, N, F, A> const& mat)
    {
        struct Ops {
            A const& mat;
            std::size_t i;
            std::size_t j;
            constexpr Ops(A const& _mat, std::size_t _i, std::size_t _j) noexcept : mat{_mat}, i{_i}, j{_j} {}
            constexpr auto operator()(elementary::swap_row_t const& op) noexcept { return swap_row(mat, op.a, op.b)(i, j); }
            constexpr auto operator()(elementary::scale_row_t<F> const& op) noexcept { return scale_row(mat, op.row, op.scalar)(i, j); }
            constexpr auto operator()(elementary::add_scaled_row_t<F> const& op) noexcept { return add_scaled_row(mat, op.target_row, op.scalar, op.row)(i, j); }
        };
        struct Apply : public mat_expr<M, N, F, Apply> {
            elementary::operation<F> op;
            A const& mat;
            constexpr Apply(elementary::operation<F> _op, A const& _mat) noexcept : op{std::move(_op)}, mat{_mat} {}
            constexpr F operator()(std::size_t i, std::size_t j) const noexcept { return std::visit(Ops{mat, i, j}, op); }
        };
        return Apply{std::move(op), static_cast<A const&>(mat)};
    }

    template<
        std::size_t M,
        std::size_t N,
        typename F,
        typename A,
        typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0
    >
    constexpr auto operator*(elementary::operation<F> op, mat_expr<M, N, F, A> const& mat)
    {
        return apply(std::move(op), mat);
    }
}
// }}}
// {{{ cofactor
template <std::size_t M, std::size_t N, typename F, typename A>
constexpr auto cofactor(mat_expr<M, N, F, A> const& _mat)
{
    return init<M, N, F>([mat = std::ref(_mat)](std::size_t i, std::size_t j) constexpr -> F {
        auto const sgn = (i + j) % 2 == 0 ? one<F> : -one<F>;
        return sgn * det(minor(mat.get(), i, j));
    });
}
// }}}
// {{{ adjugate
template <std::size_t M, std::size_t N, typename F, typename A>
constexpr auto adjugate(mat_expr<M, N, F, A> const& _mat)
{
    // TODO: try passing over by copy (not by ref), wtf CRTP ;(
    return transpose(cofactor(_mat));
}
// }}}
// {{{ inverse
template <std::size_t M, std::size_t N, typename F, typename A>
constexpr auto inverse(mat_expr<M, N, F, A> const& _mat)
{
    auto const detM = det(_mat);
    if (detM == zero<F>)
        return std::nullopt;
    else
        return one<F> / detM * adjugate(_mat);
}

// }}}

} // namespace dim

