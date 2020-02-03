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

TEST_CASE("mat_solver.swap_row")
{
    auto constexpr static m1 = mat{1, 2, 3,
                                   4, 5, 6,
                                   7, 8, 9};

    auto const m2 = swap_row(m1, 1, 2);

    CHECK(m2 == mat{1, 2, 3,
                    7, 8, 9,
                    4, 5, 6});

}