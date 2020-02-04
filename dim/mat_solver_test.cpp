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
#include <dim/mat.h>
#include <dim/mat_solver.h>
#include <dim/mat_ostream.h>

#include <catch2/catch.hpp>

#include <iostream>

using namespace std;
using namespace dim;
using namespace dim::linear_algebra;
using namespace dim::linear_algebra::solver;

TEST_CASE("mat_solver.state")
{
    auto constexpr static m1 = mat{1, 2, 3,
                                   4, 5, 6,
                                   7, 8, 9};

    auto state = State(m1);
    state.update(elementary::swap_row<int>(0, 1));

    CHECK(state() == mat{4, 5, 6,
                         1, 2, 3,
                         7, 8, 9});

    state.update(elementary::scale_row<int>(1, 2));

    CHECK(state() == mat{4, 5, 6,
                         2, 4, 6,
                         7, 8, 9});
}

TEST_CASE("mat_solver.solve")
{
    auto constexpr m1 = mat{1, 2, 3,
                            0, 1, 4,
                            5, 6, 0};

    auto const m2 = rowCanonicalForm(m1);
}

