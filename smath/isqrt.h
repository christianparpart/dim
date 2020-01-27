#pragma once
#include <cstddef>

namespace smath {

constexpr inline std::size_t _isqrt_impl(std::size_t sq, std::size_t dlt, std::size_t value)
{
    return sq <= value ? _isqrt_impl(sq + dlt, dlt + 2, value) : (dlt >> 1) - 1;
}

// constexpr version of an integer square root
constexpr inline std::size_t isqrt(std::size_t value)
{
    return _isqrt_impl(1, 3, value);
}

} // end namespace
