#pragma once
#include <dim/mat_expr.h>
#include <dim/vec.h>
#include <dim/concept.h>

#include <optional>

namespace dim::linear_algebra::solver {

namespace detail {

    template<std::size_t M, std::size_t N, typename F, typename A>
    std::optional<std::size_t> findSmallestAtColumnFor(mat_expr<M, N, F, A> const& mat, std::size_t i0, std::size_t j)
    {
        return std::nullopt; // TODO
    }

    template<std::size_t M, std::size_t N, typename F, typename A>
    auto ensureValueOneAt(mat_expr<M, N, F, A> const& state, std::size_t i, std::size_t j)
    {
        return state; // TODO: new state
    }

    template<std::size_t M, std::size_t N, typename F, typename A>
    bool isNonNullRow(mat_expr<M, N, F, A> const& mat, std::size_t i)
    {
        return false; // TODO
    }

    template<std::size_t M, std::size_t N, typename F, typename A>
    auto makeZerosAbove(std::size_t i, std::size_t j, mat_expr<M, N, F, A> const& work)
    {
        return work; // TODO
    }

    template<std::size_t M, std::size_t N, typename F, typename A>
    auto makeZerosBelow(std::size_t i, std::size_t j, mat_expr<M, N, F, A> const& work)
    {
        return work; // TODO
    }

    template<std::size_t M, std::size_t N, typename F, typename A>
    auto rowEchelonStep(std::size_t d, mat_expr<M, N, F, A> const& work)
    {
        return work; // TODO
    }
}

template<std::size_t M, std::size_t N, typename F, typename A>
vec<N, F> solve(mat_expr<M, N, F, A> const& m)
{
    // TODO
}

} // end namespace
