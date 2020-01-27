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
#include <ostream>
#include <smath/mat.h>

namespace smath {

template <std::size_t M, std::size_t N, typename F, typename A>
std::ostream& operator<<(std::ostream& os, mat_expr<M, N, F, A> const& _mat)
{
    os << '{';
    for (std::size_t i = 0; i < M; ++i)
    {
        if (i != 0)
            os << ", ";
        os << '{';
        for (std::size_t j = 0; j < N; ++j)
        {
            if (j != 0)
                os << ", ";
            os << _mat(i, j);
        }
        os << '}';
    }
    os << '}';
    return os;
}

} // end namespace
