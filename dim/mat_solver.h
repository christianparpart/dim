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

#include <dim/concept.h>
#include <dim/mat.h>
#include <dim/mat_expr.h>
#include <dim/util.h>
#include <dim/vec.h>

#include <list>
#include <vector>
#include <optional>

namespace dim::linear_algebra::solver {

// ==============================================================================

/**
 * Keeps track of modifications and kind of modifications done on a matrix.
 *
 * TODO: find a better name (Journal, History, ChangeLog, ...).
 */
template <std::size_t M, std::size_t N, typename F, typename Mat = mat<N, N, F>>
struct State : public mat_expr<M, N, F, State<M, N, F, Mat>> {
    struct Step {
        elementary::operation<F> operation;
        Mat matrix;
    };

    Mat tip;
    std::vector<Step> steps;

    template <typename A>
    explicit State(mat_expr<M, N, F, A> const& m) : tip{m}, steps{} {}

    Mat const& operator()() const noexcept { return tip; }
    F operator()(std::size_t i, std::size_t j) const noexcept { return tip(i, j); }

    State<M, N, F>& update(elementary::operation<F> operation)
    {
        steps.emplace_back(Step{operation, operation * tip});
        tip = steps.back().matrix;
        return *this;
    }
};

template <std::size_t M, std::size_t N, typename F, typename A>
State(mat_expr<M, N, F, A>) -> State<M, N, F>;

// ==============================================================================

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

    template <std::size_t M, std::size_t N, typename F>
    State<M, N, F>& rowCanonicalStep(std::size_t _d, State<M, N, F>& _state)
    {
        using namespace detail;
        for (auto d : times(_d, N - _d))
        {
            ensureValueOneAt(_state, d, d);
            makeZerosAbove(d, d, _state);
            makeZerosBelow(d, d, _state);
        }
        return _state;
    }
}

// ==============================================================================

template <std::size_t M, std::size_t N, typename F, typename A>
auto rowCanonicalForm(mat_expr<M, N, F, A> const& _mat)
{
    auto state = State{_mat};
    detail::rowCanonicalStep(0, state);
    return state;
}

// ==============================================================================

template<std::size_t M, std::size_t N, typename F, typename A>
vec<M, F> solve(mat_expr<M, N, F, A> const& m, vec<N, F> const& b)
{
    // TODO extend m with b -> invert -> read out right most column as x
    // (that's parsing, but doesn't exist yet - big TODO :-D)
    return vec{column(N - 1, rowCanonicalForm(m | b))};
}

} // end namespace
