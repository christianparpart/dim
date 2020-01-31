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
#include <dim/util.h>
#include <catch2/catch.hpp>
#include <vector>
#include <iostream>

using namespace std;
using namespace dim;

TEST_CASE("util.times")
{
    SECTION("count 0") {
        vector<int> hits;
        for (auto const& i : times(0, 0, 1))
            hits.push_back(i);

        REQUIRE(hits.empty());
    }

    SECTION("count 1") {
        vector<int> hits;
        for (auto const& i : times(0, 1, 1))
            hits.push_back(i);

        REQUIRE(hits.size() == 1);
        REQUIRE(hits[0] == 0);
    }

    SECTION("count 1, step 2") {
        vector<int> hits;
        auto t = times(1, 2, 3);
        for (auto const& i : t)
            hits.push_back(i);

        REQUIRE(hits.size() == 2);
        REQUIRE(hits[0] == 1);
        REQUIRE(hits[1] == 4);
    }

    SECTION("step 1") {
        vector<int> hits;
        for (auto const& i : times(0, 3, 1))
            hits.push_back(i);

        REQUIRE(hits.size() == 3);
        REQUIRE(hits[0] == 0);
        REQUIRE(hits[1] == 1);
        REQUIRE(hits[2] == 2);
    }

    SECTION("step 2") {
        vector<int> hits;
        for (auto const& i : times(4, 3, 2))
            hits.push_back(i);

        REQUIRE(hits.size() == 3);
        REQUIRE(hits[0] == 4);
        REQUIRE(hits[1] == 6);
        REQUIRE(hits[2] == 8);
    }
}

TEST_CASE("util.times.compose_with_lambda")
{
    vector<int> hits;
    times(0, 3) | [&](auto i) { hits.push_back(i); };

    REQUIRE(hits.size() == 3);
    REQUIRE(hits[0] == 0);
    REQUIRE(hits[1] == 1);
    REQUIRE(hits[2] == 2);
}

TEST_CASE("util.times2")
{
    vector<pair<int, int>> hits;

    for (auto && [i, j] : times(0, 3) | times(0, 3))
        hits.emplace_back(pair{i, j});

    REQUIRE(hits.size() == 9);

    REQUIRE(hits[0] == pair{0, 0});
    REQUIRE(hits[1] == pair{0, 1});
    REQUIRE(hits[2] == pair{0, 2});

    REQUIRE(hits[3] == pair{1, 0});
    REQUIRE(hits[4] == pair{1, 1});
    REQUIRE(hits[5] == pair{1, 2});

    REQUIRE(hits[6] == pair{2, 0});
    REQUIRE(hits[7] == pair{2, 1});
    REQUIRE(hits[8] == pair{2, 2});
}

TEST_CASE("util.times2.compose_with_lambda")
{
    vector<pair<int, char>> hits;

    times(0, 3) | times('a', 3) | [&](auto i, auto j) {
        hits.emplace_back(pair{i, j});
    };

    REQUIRE(hits.size() == 9);

    REQUIRE(hits[0] == pair{0, 'a'});
    REQUIRE(hits[1] == pair{0, 'b'});
    REQUIRE(hits[2] == pair{0, 'c'});

    REQUIRE(hits[3] == pair{1, 'a'});
    REQUIRE(hits[4] == pair{1, 'b'});
    REQUIRE(hits[5] == pair{1, 'c'});

    REQUIRE(hits[6] == pair{2, 'a'});
    REQUIRE(hits[7] == pair{2, 'b'});
    REQUIRE(hits[8] == pair{2, 'c'});
}

TEST_CASE("util.times2.x")
{
    vector<pair<int, char>> hits;

    times(0, 1) | times('a', 3) | [&](auto i, auto j) {
        hits.emplace_back(pair{i, j});
    };

    REQUIRE(hits.size() == 3);

    REQUIRE(hits[0] == pair{0, 'a'});
    REQUIRE(hits[1] == pair{0, 'b'});
    REQUIRE(hits[2] == pair{0, 'c'});
}

TEST_CASE("util.zipped")
{
    auto constexpr a = array{1, 3, 5};
    auto constexpr b = array{2, 4, 6};

    auto z = vector<pair<int, int>>{};

    for (auto && [x, y] : zipped(a, b))
        z.push_back(pair{x, y});

    REQUIRE(z.size() == 3);
    REQUIRE(z[0] == pair{1, 2});
    REQUIRE(z[1] == pair{3, 4});
    REQUIRE(z[2] == pair{5, 6});
}
