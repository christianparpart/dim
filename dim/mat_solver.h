#pragma once
#include <dim/mat_expr.h>
#include <dim/vec.h>

#include <cassert>

namespace dim::linear_algebra::solver {

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

template<std::size_t M, std::size_t N, typename F, typename A>
vec<N, F> solve(mat_expr<M, N, F, A> const& m)
{
    // TODO
}

} // end namespace
