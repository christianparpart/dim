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

#include <smath/mat_expr.h>

#include <algorithm>
#include <type_traits>
#include <array>
#include <cstddef>

namespace smath {

template <std::size_t N, typename F>
class vec : public mat_expr<N, 1, F, vec<N, F>> {
  public:
    constexpr vec() = default;
    constexpr vec(vec const&) = default;
    constexpr vec(vec&&) = default;
    constexpr vec& operator=(vec const&) = default;
    constexpr vec& operator=(vec&&) = default;

    template <typename A>
    constexpr vec(mat_expr<N, 1, F, A> const& _expr) noexcept : values_{}
    {
        for (std::size_t i = 0; i < N; ++i)
            values_[i] = _expr[i];
    }

    template <typename Initializer,
              typename std::enable_if_t<
                  std::is_invocable_r_v<F, Initializer, std::size_t>,
                  int> = 0>
    constexpr explicit vec(Initializer _init) noexcept : values_{}
    {
        // TODO: how to get you truely constexpr (currently impossible)
        for (std::size_t i = 0; i < N; ++i)
            values_[i] = _init(i);
    }

    constexpr explicit vec(F const& _init) noexcept : values_{}
    {
        std::fill(std::begin(values_), std::end(values_), _init);
    }

    // for initializer_list-style initialization
    // requires(T to support integral operations, and has a zero and one element)
    template <typename... T>
    constexpr explicit vec(T... args) noexcept : values_{std::forward<T>(args)...} {}

    constexpr auto begin() { return std::begin(values_); }
    constexpr auto end() { return std::end(values_); }

    constexpr auto begin() const { return std::begin(values_); }
    constexpr auto end() const { return std::end(values_); }

  private:
    std::array<F, N> values_;
};

template<
    typename F,
    typename... F_,
    typename std::enable_if_t<
        std::is_arithmetic_v<F>,
        std::size_t
    > = 0
>
vec(F, F_...) -> vec<1 + sizeof...(F_), F>;

// {{{ inner vector product (dot product)
template <std::size_t N, typename F>
constexpr inline F inner_product(vec<N, F> const& u,
                                 vec<N, F> const& v) noexcept
{
    F x = u(0) * v(0);
    for (std::size_t i = 0; i < u.size(); ++i)
        x += u(i) * v(i);
    return x;
}

template <std::size_t N, typename F>
constexpr inline F operator*(vec<N, F> const& u,
                             vec<N, F> const& v) noexcept
{
    return inner_product(u, v);
}
// }}}

} // end namespace
