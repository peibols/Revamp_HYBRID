#pragma once

#include <assert.h>
#include <algorithm>
#include <array>
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

// --------- std::array<double,4> arithmetic operators (stack-allocated, no heap) ---------
inline std::array<double,4> operator+(const std::array<double,4>& a, const std::array<double,4>& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3]};
}
inline std::array<double,4> operator-(const std::array<double,4>& a, const std::array<double,4>& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3]};
}
inline std::array<double,4>& operator+=(std::array<double,4>& a, const std::array<double,4>& b) {
    a[0]+=b[0]; a[1]+=b[1]; a[2]+=b[2]; a[3]+=b[3]; return a;
}
inline std::array<double,4>& operator-=(std::array<double,4>& a, const std::array<double,4>& b) {
    a[0]-=b[0]; a[1]-=b[1]; a[2]-=b[2]; a[3]-=b[3]; return a;
}
inline std::array<double,4> operator*(const std::array<double,4>& a, double c) {
    return {a[0]*c, a[1]*c, a[2]*c, a[3]*c};
}
inline std::array<double,4> operator*(double c, const std::array<double,4>& a) {
    return {a[0]*c, a[1]*c, a[2]*c, a[3]*c};
}
inline std::array<double,4> operator/(const std::array<double,4>& a, double c) {
    return {a[0]/c, a[1]/c, a[2]/c, a[3]/c};
}
inline std::array<double,4>& operator*=(std::array<double,4>& a, double c) {
    a[0]*=c; a[1]*=c; a[2]*=c; a[3]*=c; return a;
}
inline std::array<double,4>& operator/=(std::array<double,4>& a, double c) {
    a[0]/=c; a[1]/=c; a[2]/=c; a[3]/=c; return a;
}
