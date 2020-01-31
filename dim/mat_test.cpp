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
#include <dim/mat_ostream.h>

#include <catch2/catch.hpp>

using namespace std;
using namespace dim;

TEST_CASE("mat.ctor")
{
    SECTION("square matrix via deduction guide") {
        auto constexpr a = mat{0.0, 0.1, 0.2,
                               1.0, 1.1, 1.2,
                               2.0, 2.1, 2.2};

        REQUIRE(a(0, 0) == 0.0);
        REQUIRE(a(0, 1) == 0.1);
        REQUIRE(a(0, 2) == 0.2);
    }

    SECTION("init") {
        // TODO: (GCC) why is constexpr here not working?
        auto const static a = init<2, 3, int>([](std::size_t i, std::size_t j) constexpr {
            return 3 * i + j;
        });

        REQUIRE(a(0, 0) == 0);
        REQUIRE(a(0, 1) == 1);
        REQUIRE(a(0, 2) == 2);

        REQUIRE(a(1, 0) == 3);
        REQUIRE(a(1, 1) == 4);
        REQUIRE(a(1, 2) == 5);
    }

    SECTION("lambda element-initializer") {
        auto constexpr static a = mat<2, 3, int>{[](auto i, auto j) constexpr {
            return 3 * i + j;
        }};

        REQUIRE(a(0, 0) == 0);
        REQUIRE(a(0, 1) == 1);
        REQUIRE(a(0, 2) == 2);

        REQUIRE(a(1, 0) == 3);
        REQUIRE(a(1, 1) == 4);
        REQUIRE(a(1, 2) == 5);
    }
}

TEST_CASE("mat.eq")
{
    auto constexpr a = mat<2, 3, int>{1, 3, 5,
                                      2, 4, 6};

    REQUIRE(a(0, 0) == 1);
    REQUIRE(a(0, 1) == 3);
    REQUIRE(a(0, 2) == 5);

    REQUIRE(a(1, 0) == 2);
    REQUIRE(a(1, 1) == 4);
    REQUIRE(a(1, 2) == 6);
}

TEST_CASE("mat.zeros.squared")
{
    auto constexpr static z = zeros<2, double>();
    REQUIRE(z(0, 0) == 0);
    REQUIRE(z(0, 1) == 0);
    REQUIRE(z(1, 0) == 0);
    REQUIRE(z(1, 1) == 0);
}

TEST_CASE("mat.zeros.non_squared")
{
    auto constexpr static z = zeros<2, 1, int>();
    REQUIRE(z(0, 0) == 0);
    REQUIRE(z(1, 0) == 0);
}

TEST_CASE("mat.diagonals.squared")
{
    auto constexpr static m = diagonals<2, int>(1);
    REQUIRE(m(0, 0) == 1);
    REQUIRE(m(0, 1) == 0);
    REQUIRE(m(1, 0) == 0);
    REQUIRE(m(1, 1) == 1);
}

TEST_CASE("mat.diagonals.variadric")
{
    auto constexpr static m = diagonals(1, 2, 3);

    static_assert(m(0, 0) == 1);
    static_assert(m(0, 1) == 0);
    static_assert(m(0, 2) == 0);

    static_assert(m(1, 0) == 0);
    static_assert(m(1, 1) == 2);
    static_assert(m(1, 2) == 0);

    static_assert(m(2, 0) == 0);
    static_assert(m(2, 1) == 0);
    static_assert(m(2, 2) == 3);
}

TEST_CASE("mat.ones.squared")
{
    auto constexpr m = ones<2, int>();
    REQUIRE(m(0, 0) == 1);
    REQUIRE(m(0, 1) == 1);
    REQUIRE(m(1, 0) == 1);
    REQUIRE(m(1, 1) == 1);
}

TEST_CASE("mat.vec")
{
    SECTION("row-vector") {
        auto constexpr v = mat<1, 3, int>{0, 1, 2};
        REQUIRE(v(0) == 0);
        REQUIRE(v(1) == 1);
        REQUIRE(v(2) == 2);
    }

    SECTION("column-vector") {
        auto constexpr v = mat<3, 1, int>{0, 1, 2};
        REQUIRE(v(0) == 0);
        REQUIRE(v(1) == 1);
        REQUIRE(v(2) == 2);
    }
}

TEST_CASE("mat.abs")
{
    SECTION("row-vector") {
        auto constexpr u = mat<1, 3, int>{0, 1, 2};
        auto constexpr v = mat<1, 3, int>{0, -1, -2};
        auto const w = abs(v);
        REQUIRE(u == w);
    }

    SECTION("column-vector") {
        auto const v = mat<3, 1, int>{0, 1, 2};
        auto const w = mat<3, 1, int>{0, 1, -2};
        auto const s = abs(w);
        auto const S = mat{s};
        CHECK(v == s);
        CHECK(v == S);
    }
}

TEST_CASE("mat.add")
{
    auto constexpr a = mat{1, 2,
                           3, 4};
    auto constexpr b = mat{5, 6,
                           7, 8};
    auto const c = a + b;

    REQUIRE(c(0, 0) == 6);
    REQUIRE(c(0, 1) == 8);
    REQUIRE(c(1, 0) == 10);
    REQUIRE(c(1, 1) == 12);
}

TEST_CASE("mat.sub")
{
    auto constexpr a = mat{1, 2,
                           3, 4};
    auto constexpr b = mat{5, 6,
                           7, 8};
    auto constexpr c = mat{6,  8,
                           10, 12};
    auto const A = c - b;

    REQUIRE(a == A);
}

TEST_CASE("mat.neg")
{
    auto constexpr original = mat{-1, 2,
                                  3, -4};
    auto constexpr expected = mat{1, -2,
                                  -3, 4};
    auto const evaluated = -original;

    REQUIRE(evaluated == expected);
}

TEST_CASE("mat.scalar_mult")
{
    auto constexpr a = mat{1, 2,
                           3, 4};
    auto const b = 2 * a;

    REQUIRE(b(0, 0) == 2);
    REQUIRE(b(0, 1) == 4);
    REQUIRE(b(1, 0) == 6);
    REQUIRE(b(1, 1) == 8);
}

TEST_CASE("mat.mat_mult")
{
    auto constexpr a = mat<2, 3, int>{1, 2, 3,
                                      4, 5, 6};
    auto constexpr b = mat<3, 2, int>{6, 5,
                                      4, 3,
                                      2, 1};
    auto const c = a * b;
    REQUIRE(c(0, 0) == 20);
    REQUIRE(c(0, 1) == 14);
    REQUIRE(c(1, 0) == 56);
    REQUIRE(c(1, 1) == 41);
}
TEST_CASE("mat.density_and_sparsity")
{
    auto constexpr all_zeros = mat{zeros<3, int>()};
    REQUIRE(count_zeros(all_zeros) == 9);
    REQUIRE(density(all_zeros) == 0.0f);
    REQUIRE(sparsity(all_zeros) == 1.0f);

    auto constexpr I3 = mat{diagonals<3, int>(1)};
    REQUIRE(count_ones(I3) == 3);
    REQUIRE(count_zeros(I3) == 6);
    REQUIRE(density(I3) == Approx(1.0f / 3.0f));
    REQUIRE(sparsity(I3) == Approx(2.0f / 3.0f));
}

TEST_CASE("mat.trace")
{
    auto constexpr static m = diagonals(1, 2, 3);
    static_assert(trace(m) == 6);
}

TEST_CASE("mat.transpose")
{
    SECTION("square") {
        auto constexpr a = mat{1, 2, 3,
                               4, 5, 6,
                               7, 8, 9};
        auto constexpr T = mat{1, 4, 7,
                               2, 5, 8,
                               3, 6, 9};
        auto const b = transpose(a);
        REQUIRE(b == T);
    }

    SECTION("non_square") {
        auto constexpr a = mat<2, 3, int>{1, 2, 3,
                                          4, 5, 6};
        auto constexpr T = mat<3, 2, int>{1, 4,
                                          2, 5,
                                          3, 6};
        auto const b = transpose(a);
        REQUIRE(b == T);
    }
}

TEST_CASE("mat.complement")
{
    auto constexpr a = mat{1, 2, 3,
                           4, 5, 6,
                           7, 8, 9};

    SECTION("middle") {
        auto constexpr expected = mat{1, 3,
                                      7, 9};
        auto const actual = complement(a, 1, 1);
        REQUIRE(actual == expected);
    }

    SECTION("top left") {
        auto constexpr expected = mat{5, 6,
                                      8, 9};
        auto const actual = complement(a, 0, 0);
        REQUIRE(actual == expected);
    }

    SECTION("bottom right") {
        auto constexpr expected = mat{1, 2,
                                      4, 5};
        auto const actual = complement(a, 2, 2);
        REQUIRE(actual == expected);
    }
}

TEST_CASE("mat.det")
{
    SECTION("1x1") {
        auto constexpr a = mat{42};
        auto constexpr d = det(a);
        REQUIRE(d == 42);
    }

    SECTION("2x2") {
        auto constexpr a = mat{1, 2,
                               3, 4};
        auto constexpr d = det(a);
        REQUIRE(d == -2);
    }

    SECTION("3x3") {
        auto constexpr a = mat{3, 5, 1,
                               1, 4, 2,
                               7, 1, 9};
        auto constexpr d = det(a);
        REQUIRE(d == 100);
    }

    SECTION("4x4") {
        auto constexpr a = mat{4, 3, 2, 2,
                               0, 1, 0,-2,
                               1,-1, 0, 3,
                               2, 3, 0, 1};
        auto constexpr d = det(a);
        REQUIRE(d == -10);
    }
}
