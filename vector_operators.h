#pragma once

#include <assert.h>
#include <algorithm>
#include <functional>

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& c)
{
    assert(a.size() == c.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), c.begin(),
        std::back_inserter(result), std::plus<T>());
    return result;
}

template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& c)
{
    assert(a.size() == c.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), c.begin(),
        std::back_inserter(result), std::minus<T>());
    return result;
}

template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& c)
{
    std::transform(a.begin(), a.end(), c.begin(),
        a.begin(), std::plus<T>());
    return a;
}

template <typename T>
std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& c)
{
    std::transform(a.begin(), a.end(), c.begin(),
        a.begin(), std::minus<T>());
    return a;
}

template <typename T>
std::vector<T> operator*(std::vector<T> a, const T& c)
{
    std::transform(a.begin(), a.end(), a.begin(),
        std::bind1st(std::multiplies<T>(),c));
    return a;
}

template <typename T>
std::vector<T> operator/(std::vector<T> a, const T& c)
{
    std::transform(a.begin(), a.end(), a.begin(),
        std::bind2nd(std::divides<T>(),c));
    return a;
}

template <typename T>
std::vector<T>& operator*=(std::vector<T>& a, const T& c)
{
    std::transform(a.begin(), a.end(), a.begin(),
        std::bind1st(std::multiplies<T>(),c));
    return a;
}

template <typename T>
std::vector<T>& operator/=(std::vector<T>& a, const T& c)
{
    std::transform(a.begin(), a.end(), a.begin(),
        std::bind2nd(std::divides<T>(),c));
    return a;
}
