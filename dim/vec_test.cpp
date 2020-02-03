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
#include <dim/vec.h>
#include <dim/vec_ostream.h>
#include <dim/mat.h>
#include <dim/util.h>

#include <iostream>

#include <catch2/catch.hpp>

using namespace std;
using namespace dim;

constexpr vec<3, double> get_vec()
{
    auto constexpr u = vec{1.0, 2.0, 3.0};
    auto constexpr v = vec{4.0, 5.0, 6.0};
    //return u + v;
    return vec{u*u, u*v, v*v};
}

TEST_CASE("vec.random")
{
    //auto const x = get_vec();

    // TODO: make me more nice, less cout, more REQUIRE/CHECK :)

    //cout << x << '\n';
    // cout << abs(-x) << '\n';
    // cout << "sum: " << reduce(x, 0.0, [](auto a, auto b) { return a + b; }) << endl;
    // auto const v = vec{3.0, 4.0};
    // cout << sqrt(reduce(v, 0.0, [](auto a, auto b) { return a + b; })) << endl;
    //
    // auto const w = cross_product(vec{1.0, 0.0, 0.0}, vec{0.0, 1.0, 0.0});
    // cout << "cross: " << w << '\n';
}
