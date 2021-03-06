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

#include <dim/vec.h>
#include <ostream>

namespace dim {

template <std::size_t N, typename F>
inline std::ostream& operator<<(std::ostream& os, vec<N, F> const& e)
{
    os << '(';
    for (std::size_t i = 0; i < e.size(); ++i)
    {
        if (i)
            os << ", ";
        os << e(i);
    }
    os << ')';
    return os;
}

} // end namespace
