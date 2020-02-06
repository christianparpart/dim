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

#include <array>
#include <cstddef>

namespace dim {

/**
 * Base interface for every matrix expression (such as multiplication, addition, complement, det, minor, ...).
 */
template <
    std::size_t M,
    std::size_t N,
    typename F,
    typename A,
    typename std::enable_if_t<std::is_arithmetic_v<F>, int> = 0
>
struct mat_expr
{
    using base_type = mat_expr<M, N, F, A>;
    using element_type = F;
    static constexpr std::size_t row_count = M;
    static constexpr std::size_t column_count = N;

    constexpr operator A const& () const noexcept { return static_cast<A const&>(*this); }

    constexpr F operator()(std::size_t i, std::size_t j) const noexcept {
        return static_cast<A const&>(*this)(i, j);
    }

    constexpr F operator()(std::size_t i) const noexcept
    {
        if constexpr (M == 1)
            return (*this)(0, i);
        else if constexpr (N == 1)
            return (*this)(i, 0);
        else
            static_assert("WTF?");
    }

    constexpr std::size_t size() const noexcept { return M * N; }
};

} // namespace dim
