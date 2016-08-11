#ifndef DUNE_PDELAB_COMMON_VECTORCLASS_HH
#define DUNE_PDELAB_COMMON_VECTORCLASS_HH

#ifdef VECTORCLASS_H
#error Do not include the vectorclass.h header directly, always use this wrapper!
#endif

#ifdef VCL_NAMESPACE
#error Do not manually set VCL_NAMESPACE, it is used internally by PDELab
#endif

#ifndef ENABLE_COUNTER

#include <vectorclass.h>
#include <vectormath_exp.h>
#include <vectormath_trig.h>

#else

#include <algorithm>

#define VCL_NAMESPACE _vcl
#include <vectorclass.h>
#include <vectormath_exp.h>
#include <vectormath_trig.h>

#include <dune/pdelab/common/opcounter.hh>

struct Vec4d
{
  oc::OpCounter<double> _d[4];

  using F = oc::OpCounter<double>;

  Vec4d()
  {}

  Vec4d(F d)
  {
    std::fill(_d,_d+4,d);
  }

  Vec4d(double d)
  {
    std::fill(_d,_d+4,d);
  }

  Vec4d(F d0, F d1, F d2, F d3)
    : _d{d0,d1,d2,d3}
  {}

  Vec4d& load(const F* p)
  {
    std::copy(p,p+4,_d);
    return *this;
  }

  Vec4d& load_a(const F* p)
  {
    std::copy(p,p+4,_d);
    return *this;
  }

  void store(F* p) const
  {
    std::copy(_d,_d+4,p);
  }

  void store_a(F* p) const
  {
    std::copy(_d,_d+4,p);
  }

  Vec4d const& insert(uint32_t index, F value)
  {
    _d[index] = value;
  }

  F extract(uint32_t index) const
  {
    return _d[index];
  }

  constexpr static int size()
  {
    return 4;
  }

};


/*****************************************************************************
*
*          Operators for Vec4d
*
*****************************************************************************/

// vector operator + : add element by element
static inline Vec4d operator + (Vec4d const & a, Vec4d const & b) {
  Vec4d r;
  std::transform(a._d,a._d+4,b._d,r._d,[](auto x, auto y){ return x + y; });
  return r;
}

// vector operator += : add
static inline Vec4d & operator += (Vec4d & a, Vec4d const & b) {
  std::transform(a._d,a._d+4,b._d,a._d,[](auto x, auto y){ return x + y; });
  return a;
}

// postfix operator ++
static inline Vec4d operator ++ (Vec4d & a, int) {
    Vec4d a0 = a;
    a = a + 1.0;
    return a0;
}

// prefix operator ++
static inline Vec4d & operator ++ (Vec4d & a) {
    a = a + 1.0;
    return a;
}

// vector operator - : subtract element by element
static inline Vec4d operator - (Vec4d const & a, Vec4d const & b) {
  Vec4d r;
  std::transform(a._d,a._d+4,b._d,r._d,[](auto x, auto y){ return x - y; });
  return r;
}

// vector operator - : unary minus
// Change sign bit, even for 0, INF and NAN
static inline Vec4d operator - (Vec4d const & a) {
  Vec4d r(a);
  for (size_t i = 0 ; i < 3 ; ++i)
    r._d[i] = -a._d[i];
  return r;
}

// vector operator -= : subtract
static inline Vec4d & operator -= (Vec4d & a, Vec4d const & b) {
  std::transform(a._d,a._d+4,b._d,a._d,[](auto x, auto y){ return x - y; });
  return a;
}

// postfix operator --
static inline Vec4d operator -- (Vec4d & a, int) {
    Vec4d a0 = a;
    a = a - 1.0;
    return a0;
}

// prefix operator --
static inline Vec4d & operator -- (Vec4d & a) {
    a = a - 1.0;
    return a;
}

// vector operator * : multiply element by element
static inline Vec4d operator * (Vec4d const & a, Vec4d const & b) {
  Vec4d r;
  std::transform(a._d,a._d+4,b._d,r._d,[](auto x, auto y){ return x * y; });
  return r;
}

// vector operator *= : multiply
static inline Vec4d & operator *= (Vec4d & a, Vec4d const & b) {
  std::transform(a._d,a._d+4,b._d,a._d,[](auto x, auto y){ return x * y; });
  return a;
}

// vector operator / : divide all elements by same integer
static inline Vec4d operator / (Vec4d const & a, Vec4d const & b) {
  Vec4d r;
  std::transform(a._d,a._d+4,b._d,r._d,[](auto x, auto y){ return x / y; });
  return r;
}

// vector operator /= : divide
static inline Vec4d & operator /= (Vec4d & a, Vec4d const & b) {
  std::transform(a._d,a._d+4,b._d,a._d,[](auto x, auto y){ return x / y; });
  return a;
}

// vector operator == : returns true for elements for which a == b
static inline _vcl::Vec4db operator == (Vec4d const & a, Vec4d const & b) {
  _vcl::Vec4d a_, b_;
  a_.load(a._d[0].data());
  b_.load(b._d[0].data());
  Vec4d::F::comparisons(4);
  return a_ == b_;
}

// vector operator != : returns true for elements for which a != b
static inline _vcl::Vec4db operator != (Vec4d const & a, Vec4d const & b) {
  _vcl::Vec4d a_, b_;
  a_.load(a._d[0].data());
  b_.load(b._d[0].data());
  Vec4d::F::comparisons(4);
  return a_ != b_;
}

// vector operator < : returns true for elements for which a < b
static inline _vcl::Vec4db operator < (Vec4d const & a, Vec4d const & b) {
  _vcl::Vec4d a_, b_;
  a_.load(a._d[0].data());
  b_.load(b._d[0].data());
  Vec4d::F::comparisons(4);
  return a_ < b_;
}

// vector operator <= : returns true for elements for which a <= b
static inline _vcl::Vec4db operator <= (Vec4d const & a, Vec4d const & b) {
  _vcl::Vec4d a_, b_;
  a_.load(a._d[0].data());
  b_.load(b._d[0].data());
  Vec4d::F::comparisons(4);
  return a_ <= b_;
}

// vector operator > : returns true for elements for which a > b
static inline _vcl::Vec4db operator > (Vec4d const & a, Vec4d const & b) {
    return b < a;
}

// vector operator >= : returns true for elements for which a >= b
static inline _vcl::Vec4db operator >= (Vec4d const & a, Vec4d const & b) {
    return b <= a;
}

// avoid logical operators for now, I don't think we need them
#if 0

// Bitwise logical operators

// vector operator & : bitwise and
static inline Vec4d operator & (Vec4d const & a, Vec4d const & b) {
    return _mm256_and_pd(a, b);
}

// vector operator &= : bitwise and
static inline Vec4d & operator &= (Vec4d & a, Vec4d const & b) {
    a = a & b;
    return a;
}

// vector operator & : bitwise and of Vec4d and Vec4db
static inline Vec4d operator & (Vec4d const & a, Vec4db const & b) {
    return _mm256_and_pd(a, b);
}
static inline Vec4d operator & (Vec4db const & a, Vec4d const & b) {
    return _mm256_and_pd(a, b);
}

// vector operator | : bitwise or
static inline Vec4d operator | (Vec4d const & a, Vec4d const & b) {
    return _mm256_or_pd(a, b);
}

// vector operator |= : bitwise or
static inline Vec4d & operator |= (Vec4d & a, Vec4d const & b) {
    a = a | b;
    return a;
}

// vector operator ^ : bitwise xor
static inline Vec4d operator ^ (Vec4d const & a, Vec4d const & b) {
    return _mm256_xor_pd(a, b);
}

// vector operator ^= : bitwise xor
static inline Vec4d & operator ^= (Vec4d & a, Vec4d const & b) {
    a = a ^ b;
    return a;
}

// vector operator ! : logical not. Returns Boolean vector
static inline Vec4db operator ! (Vec4d const & a) {
    return a == Vec4d(0.0);
}

#endif


// General arithmetic functions, etc.

// Horizontal add: Calculates the sum of all vector elements.
static inline Vec4d::F horizontal_add (Vec4d const & a) {
  return std::accumulate(a._d,a._d+4,Vec4d::F(0.0));
}

// function max: a > b ? a : b
static inline Vec4d max(Vec4d const & a, Vec4d const & b) {
  Vec4d r;
  std::transform(a._d,a._d+4,b._d,r._d,[](auto x, auto y){ return max(x,y); });
  return r;
}

// function min: a < b ? a : b
static inline Vec4d min(Vec4d const & a, Vec4d const & b) {
  Vec4d r;
  std::transform(a._d,a._d+4,b._d,r._d,[](auto x, auto y){ return min(x,y); });
  return r;
}

// function abs: absolute value
// Removes sign bit, even for -0.0f, -INF and -NAN
static inline Vec4d abs(Vec4d const & a) {
  Vec4d r;
  std::transform(a._d,a._d+4,r._d,[](auto x){ return abs(x); });
  return r;
}

// function sqrt: square root
static inline Vec4d sqrt(Vec4d const & a) {
  Vec4d r;
  std::transform(a._d,a._d+4,r._d,[](auto x){ return sqrt(x); });
  return r;
}

// function square: a * a
static inline Vec4d square(Vec4d const & a) {
  return a * a;
}


// ignore pow() for now
#if 0

// pow(Vec4d, int):
template <typename TT> static Vec4d pow(Vec4d const & a, TT n);

// Raise floating point numbers to integer power n
template <>
inline Vec4d pow<int>(Vec4d const & x0, int n) {
    return pow_template_i<Vec4d>(x0, n);
}

// allow conversion from unsigned int
template <>
inline Vec4d pow<uint32_t>(Vec4d const & x0, uint32_t n) {
    return pow_template_i<Vec4d>(x0, (int)n);
}


// Raise floating point numbers to integer power n, where n is a compile-time constant
template <int n>
static inline Vec4d pow_n(Vec4d const & a) {
    if (n < 0)    return Vec4d(1.0) / pow_n<-n>(a);
    if (n == 0)   return Vec4d(1.0);
    if (n >= 256) return pow(a, n);
    Vec4d x = a;                       // a^(2^i)
    Vec4d y;                           // accumulator
    const int lowest = n - (n & (n-1));// lowest set bit in n
    if (n & 1) y = x;
    if (n < 2) return y;
    x = x*x;                           // x^2
    if (n & 2) {
        if (lowest == 2) y = x; else y *= x;
    }
    if (n < 4) return y;
    x = x*x;                           // x^4
    if (n & 4) {
        if (lowest == 4) y = x; else y *= x;
    }
    if (n < 8) return y;
    x = x*x;                           // x^8
    if (n & 8) {
        if (lowest == 8) y = x; else y *= x;
    }
    if (n < 16) return y;
    x = x*x;                           // x^16
    if (n & 16) {
        if (lowest == 16) y = x; else y *= x;
    }
    if (n < 32) return y;
    x = x*x;                           // x^32
    if (n & 32) {
        if (lowest == 32) y = x; else y *= x;
    }
    if (n < 64) return y;
    x = x*x;                           // x^64
    if (n & 64) {
        if (lowest == 64) y = x; else y *= x;
    }
    if (n < 128) return y;
    x = x*x;                           // x^128
    if (n & 128) {
        if (lowest == 128) y = x; else y *= x;
    }
    return y;
}

template <int n>
static inline Vec4d pow(Vec4d const & a, Const_int_t<n>) {
    return pow_n<n>(a);
}

#endif

// function round: round to nearest integer (even). (result as double vector)
static inline Vec4d round(Vec4d const & a) {
  Vec4d r;
  std::transform(a._d,a._d+4,r._d,[](auto x){ return round(x); });
  return r;
}

// function truncate: round towards zero. (result as double vector)
static inline Vec4d truncate(Vec4d const & a) {
  Vec4d r;
  std::transform(a._d,a._d+4,r._d,[](auto x){ return trunc(x); });
  return r;
}

// function floor: round towards minus infinity. (result as double vector)
static inline Vec4d floor(Vec4d const & a) {
  Vec4d r;
  std::transform(a._d,a._d+4,r._d,[](auto x){ return floor(x); });
  return r;
}

// function ceil: round towards plus infinity. (result as double vector)
static inline Vec4d ceil(Vec4d const & a) {
  Vec4d r;
  std::transform(a._d,a._d+4,r._d,[](auto x){ return ceil(x); });
  return r;
}

#if 0
// function round_to_int: round to nearest integer (even). (result as integer vector)
static inline Vec4i round_to_int(Vec4d const & a) {
    // Note: assume MXCSR control register is set to rounding
    return _mm256_cvtpd_epi32(a);
}

// function truncate_to_int: round towards zero. (result as integer vector)
static inline Vec4i truncate_to_int(Vec4d const & a) {
    return _mm256_cvttpd_epi32(a);
}
#endif


// Fused multiply and add functions

// Multiply and add
static inline Vec4d mul_add(Vec4d const & a, Vec4d const & b, Vec4d const & c) {
  Vec4d r;
  for (size_t i = 0 ; i < 4 ; ++i)
    r._d[i] = a._d[i] * b._d[i] + c._d[i];
  return r;
}


// Multiply and subtract
static inline Vec4d mul_sub(Vec4d const & a, Vec4d const & b, Vec4d const & c) {
  Vec4d r;
  for (size_t i = 0 ; i < 4 ; ++i)
    r._d[i] = a._d[i] * b._d[i] - c._d[i];
  return r;
}

// Multiply and inverse subtract
static inline Vec4d nmul_add(Vec4d const & a, Vec4d const & b, Vec4d const & c) {
  Vec4d r;
  for (size_t i = 0 ; i < 4 ; ++i)
    r._d[i] = - a._d[i] * b._d[i] + c._d[i];
  return r;
}


template <int i0, int i1, int i2, int i3>
static inline Vec4d blend4d(Vec4d const & a, Vec4d const & b) {
  _vcl::Vec4d a_,b_;
  a_.load(a._d[0].data());
  b_.load(b._d[0].data());
  _vcl::Vec4d r_ = _vcl::blend4d<i0,i1,i2,i3>(a_,b_);
  Vec4d::F::blends(1);
  Vec4d r;
  r_.store(r._d[0].data());
  return r;
}

#endif // ENABLE_COUNTER

#endif // DUNE_PDELAB_COMMON_VECTORCLASS_HH
