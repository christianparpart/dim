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
        auto constexpr static a = init<2, 3, int>([](std::size_t i, std::size_t j) constexpr {
            return 3 * i + j;
        });

        CHECK(a == mat<2, 3, int>{0, 1, 2, 3, 4, 5});
    }

    SECTION("lambda element-initializer") {
        auto constexpr static a = mat<2, 3, int>{[](auto i, auto j) constexpr {
            return 3 * i + j;
        }};
        CHECK(a == mat<2, 3, int>{0, 1, 2, 3, 4, 5});
    }
}

TEST_CASE("mat.eq")
{
    auto constexpr a = mat<2, 3, int>{1, 3, 5,
                                      2, 4, 6};

    CHECK(a(0, 0) == 1);
    CHECK(a(0, 1) == 3);
    CHECK(a(0, 2) == 5);

    CHECK(a(1, 0) == 2);
    CHECK(a(1, 1) == 4);
    CHECK(a(1, 2) == 6);
}

TEST_CASE("mat.zeros.squared")
{
    auto constexpr static z = zeros<2, double>();
    CHECK(z(0, 0) == 0);
    CHECK(z(0, 1) == 0);
    CHECK(z(1, 0) == 0);
    CHECK(z(1, 1) == 0);
}

TEST_CASE("mat.zeros.non_squared")
{
    auto constexpr static z = zeros<2, 1, int>();
    REQUIRE(z.row_count == 2);
    REQUIRE(z.column_count == 1);
    CHECK(z(0, 0) == 0);
    CHECK(z(1, 0) == 0);
}

TEST_CASE("mat.diagonals.squared")
{
    auto constexpr static m = diagonals<2, int>(1);
    CHECK(m == mat{1, 0, 0, 1});
}

TEST_CASE("mat.diagonals.variadric")
{
    auto constexpr static m = diagonals(1, 2, 3);

    CHECK(m == mat{1, 0, 0,
                   0, 2, 0,
                   0, 0, 3});
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
    SECTION("times 2") {
        auto constexpr a = mat{1, 2,
                               3, 4};
        static_assert(2 * a == mat{2, 4, 6, 8});
    }

    SECTION("times -one<F>") {
        auto constexpr a = mat{1, -2,
                               3, -4};
        auto constexpr b = mat{ -one<int> * a };
        static_assert(b == mat{-1, 2, -3, 4});
    }
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
    CHECK(trace(m) == 6);
}

TEST_CASE("mat.transpose")
{
    SECTION("vec") {
        auto constexpr u = mat<1, 3, int>{1, 2, 3};
        auto constexpr uT = mat<3, 1, int>{1, 2, 3};
        auto constexpr b = mat{u * uT};

        CHECK(b.row_count == 1);
        CHECK(b.column_count == 1);
        CHECK(b(0, 0) == 14);
    }

    SECTION("square") {
        auto constexpr a = mat{1, 2, 3,
                               4, 5, 6,
                               7, 8, 9};
        auto constexpr T = mat{1, 4, 7,
                               2, 5, 8,
                               3, 6, 9};
        CHECK(transpose(a) == T);
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

TEST_CASE("mat.elementary.matrix") {
    auto constexpr static m = mat{1, 2, 3,
                                  4, 5, 6,
                                  7, 8, 9};

    SECTION("swap_row") {
        auto const e = elementary::matrix<3, 3, int>(elementary::swap_row<int>(0, 1));

        CHECK(e == mat{0, 1, 0,
                       1, 0, 0,
                       0, 0, 1});

        CHECK(e * m == mat{4, 5, 6,
                           1, 2, 3,
                           7, 8, 9});
    }

    SECTION("scale_row") {
        auto const e = elementary::matrix<3, 3, int>(elementary::scale_row<int>(1, 2));

        CHECK(e == mat{1, 0, 0,
                       0, 2, 0,
                       0, 0, 1});

        CHECK(e * m == mat{1,  2,  3,
                           8, 10, 12,
                           7,  8,  9});
    }

    SECTION("add_scaled_row") {
        auto const e = elementary::matrix<3, 3, int>(elementary::add_scaled_row<int>(1, 2, 0));

        CHECK(e == mat{1, 0, 0,
                       2, 1, 0,
                       0, 0, 1});

        CHECK(e * m == mat{1, 2, 3,
                           6, 9, 12,
                           7, 8, 9});
    }
}

TEST_CASE("mat.elementary.apply")
{
    SECTION("swap_row") {
        auto constexpr static m1 = mat{1, 2, 3,
                                       4, 5, 6,
                                       7, 8, 9};

        auto const m2 = m1 * elementary::swap_row<int>(1, 2);

        CHECK(m2 == mat{1, 2, 3,
                        7, 8, 9,
                        4, 5, 6});
    }

    SECTION("scale_row") {
        auto constexpr static m1 = mat{1, 2, 3,
                                       4, 5, 6,
                                       7, 8, 9};

        auto const m2 = elementary::apply(m1, elementary::scale_row<int>(0, 2));

        CHECK(m2 == mat{2, 4, 6,
                        4, 5, 6,
                        7, 8, 9});
    }

    SECTION("add_scaled_row") {
        auto constexpr static m1 = mat{1, 2, 3,
                                       3, 1, 2,
                                       4, 5, 6};

        auto const m2 = m1 * elementary::add_scaled_row<int>(0, 2, 1);

        CHECK(m2 == mat{7, 4, 7,
                        3, 1, 2,
                        4, 5, 6});
    }
}
