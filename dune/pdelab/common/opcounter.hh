#ifndef __OPCOUNTER__
#define __OPCOUNTER__

#include <type_traits>
#include <iostream>
#include <cmath>
#include <cstdlib>

namespace oc {

  template<typename F>
  class OpCounter;

}

namespace Dune {

  template<typename T, int n>
  class FieldVector;

}

namespace oc {

  template<typename F>
  class OpCounter
  {

  public:

    typedef std::size_t size_type;

    using value_type = F;

    OpCounter()
      : _v()
    {}

    template<typename T>
    OpCounter(const T& t, typename std::enable_if<std::is_same<T,int>::value and !std::is_same<F,int>::value>::type* = nullptr)
      : _v(t)
    {}

    OpCounter(const F& f)
      : _v(f)
    {}

    OpCounter(F&& f)
      : _v(f)
    {}

    explicit OpCounter(const char* s)
      : _v(strtod(s,nullptr))
    {}

    OpCounter& operator=(const char* s)
    {
      _v = strtod(s,nullptr);
      return *this;
    }

    explicit operator F() const
    {
      return _v;
    }

    OpCounter& operator=(const F& f)
    {
      _v = f;
      return *this;
    }

    OpCounter& operator=(F&& f)
    {
      _v = f;
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const OpCounter& f)
    {
      os << "OC(" << f._v << ")";
      return os;
    }

    F* data()
    {
      return &_v;
    }

    const F* data() const
    {
      return &_v;
    }

    F _v;

    struct Counters {

      size_type addition_count;
      size_type multiplication_count;
      size_type division_count;
      size_type exp_count;
      size_type pow_count;
      size_type sin_count;
      size_type sqrt_count;
      size_type comparison_count;

      Counters()
        : addition_count(0)
        , multiplication_count(0)
        , division_count(0)
        , exp_count(0)
        , pow_count(0)
        , sin_count(0)
        , sqrt_count(0)
        , comparison_count(0)
      {}

      void reset()
      {
        addition_count = 0;
        multiplication_count = 0;
        division_count = 0;
        exp_count = 0;
        pow_count = 0;
        sin_count = 0;
        sqrt_count = 0;
        comparison_count = 0;
      }

      template<typename Stream>
      void reportOperations(Stream& os, bool doReset = false)
      {
        os << "additions: " << addition_count << std::endl
           << "multiplications: " << multiplication_count << std::endl
           << "divisions: " << division_count << std::endl
           << "exp: " << exp_count << std::endl
           << "pow: " << pow_count << std::endl
           << "sin: " << sin_count << std::endl
           << "sqrt: " << sqrt_count << std::endl
           << "comparisons: " << comparison_count << std::endl
           << std::endl
           << "total: " << addition_count + multiplication_count + division_count + exp_count + pow_count + sin_count + sqrt_count + comparison_count << std::endl;

        if (doReset)
          reset();
      }

      Counters& operator+=(const Counters& rhs)
      {
        addition_count += rhs.addition_count;
        multiplication_count += rhs.multiplication_count;
        division_count += rhs.division_count;
        exp_count += rhs.exp_count;
        pow_count += rhs.pow_count;
        sin_count += rhs.sin_count;
        sqrt_count += rhs.sqrt_count;
        comparison_count += rhs.comparison_count;
        return *this;
      }

      Counters operator-(const Counters& rhs)
      {
        Counters r;
        r.addition_count = addition_count - rhs.addition_count;
        r.multiplication_count = multiplication_count - rhs.multiplication_count;
        r.division_count = division_count - rhs.division_count;
        r.exp_count = exp_count - rhs.exp_count;
        r.pow_count = pow_count - rhs.pow_count;
        r.sin_count = sin_count - rhs.sin_count;
        r.sqrt_count = sqrt_count - rhs.sqrt_count;
        r.comparison_count = comparison_count - rhs.comparison_count;
        return r;
      }

    };

    static void additions(std::size_t n)
    {
      counters.addition_count += n;
    }

    static void multiplications(std::size_t n)
    {
      counters.multiplication_count += n;
    }

    static void divisions(std::size_t n)
    {
      counters.division_count += n;
    }

    static void reset()
    {
      counters.reset();
    }

    template<typename Stream>
    static void reportOperations(Stream& os, bool doReset = false)
    {
      counters.reportOperations(os,doReset);
    }

    static Counters counters;

  };

  template<typename F>
  typename OpCounter<F>::Counters OpCounter<F>::counters;

  // ********************************************************************************
  // negation
  // ********************************************************************************

  template<typename F>
  OpCounter<F> operator-(const OpCounter<F>& a)
  {
    ++OpCounter<F>::counters.addition_count;
    return {-a._v};
  }


  // ********************************************************************************
  // addition
  // ********************************************************************************

  template<typename F>
  OpCounter<F> operator+(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a._v + b._v};
  }

  template<typename F>
  OpCounter<F> operator+(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a._v + b};
  }

  template<typename F>
  OpCounter<F> operator+(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a + b._v};
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>
    >::type
  operator+(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a._v + b};
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>
    >::type
  operator+(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a + b._v};
  }

  template<typename F>
  OpCounter<F>& operator+=(OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.addition_count;
    a._v += b._v;
    return a;
  }

  template<typename F>
  OpCounter<F>& operator+=(OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.addition_count;
    a._v += b;
    return a;
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>&
    >::type
  operator+=(OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.addition_count;
    a._v += b;
    return a;
  }

  template<typename F>
  OpCounter<F>& operator+=(OpCounter<F>& a, const Dune::FieldVector<OpCounter<F>,1>& b)
  {
    ++OpCounter<F>::counters.addition_count;
    a._v += b[0]._v;
    return a;
  }

  // ********************************************************************************
  // subtraction
  // ********************************************************************************

  template<typename F>
  OpCounter<F> operator-(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a._v - b._v};
  }

  template<typename F>
  OpCounter<F> operator-(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a._v - b};
  }

  template<typename F>
  OpCounter<F> operator-(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a - b._v};
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>
    >::type
  operator-(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a._v - b};
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>
    >::type
  operator-(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.addition_count;
    return {a - b._v};
  }

  template<typename F>
  OpCounter<F>& operator-=(OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.addition_count;
    a._v -= b._v;
    return a;
  }

  template<typename F>
  OpCounter<F>& operator-=(OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.addition_count;
    a._v -= b;
    return a;
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>&
    >::type
  operator-=(OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.addition_count;
    a._v -= b;
    return a;
  }


  // ********************************************************************************
  // multiplication
  // ********************************************************************************

  template<typename F>
  OpCounter<F> operator*(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.multiplication_count;
    return {a._v * b._v};
  }

  template<typename F>
  OpCounter<F> operator*(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.multiplication_count;
    return {a._v * b};
  }

  template<typename F>
  OpCounter<F> operator*(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.multiplication_count;
    return {a * b._v};
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>
    >::type
  operator*(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.multiplication_count;
    return {a._v * b};
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>
    >::type
  operator*(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.multiplication_count;
    return {a * b._v};
  }

  template<typename F>
  OpCounter<F>& operator*=(OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.multiplication_count;
    a._v *= b._v;
    return a;
  }

  template<typename F>
  OpCounter<F>& operator*=(OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.multiplication_count;
    a._v *= b;
    return a;
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>&
    >::type
  operator*=(OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.multiplication_count;
    a._v *= b;
    return a;
  }


  // ********************************************************************************
  // division
  // ********************************************************************************

  template<typename F>
  OpCounter<F> operator/(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.division_count;
    return {a._v / b._v};
  }

  template<typename F>
  OpCounter<F> operator/(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.division_count;
    return {a._v / b};
  }

  template<typename F>
  OpCounter<F> operator/(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.division_count;
    return {a / b._v};
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>
    >::type
  operator/(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.division_count;
    return {a._v / b};
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>
    >::type
  operator/(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.division_count;
    return {a / b._v};
  }

  template<typename F>
  OpCounter<F>& operator/=(OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.division_count;
    a._v /= b._v;
    return a;
  }

  template<typename F>
  OpCounter<F>& operator/=(OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.division_count;
    a._v /= b;
    return a;
  }

  template<typename F, typename T>
  typename std::enable_if<
    std::is_arithmetic<T>::value,
    OpCounter<F>&
    >::type
  operator/=(OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.division_count;
    a._v /= b;
    return a;
  }



  // ********************************************************************************
  // comparisons
  // ********************************************************************************


  // ********************************************************************************
  // less
  // ********************************************************************************

  template<typename F>
  bool operator<(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v < b._v};
  }

  template<typename F>
  bool operator<(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v < b};
  }

  template<typename F>
  bool operator<(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a < b._v};
  }

  template<typename F, typename T>
  bool operator<(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v < b};
  }

  template<typename F, typename T>
  bool operator<(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a < b._v};
  }


  // ********************************************************************************
  // less_or_equals
  // ********************************************************************************

  template<typename F>
  bool operator<=(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v <= b._v};
  }

  template<typename F>
  bool operator<=(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v <= b};
  }

  template<typename F>
  bool operator<=(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a <= b._v};
  }

  template<typename F, typename T>
  bool operator<=(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v <= b};
  }

  template<typename F, typename T>
  bool operator<=(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a <= b._v};
  }


  // ********************************************************************************
  // greater
  // ********************************************************************************

  template<typename F>
  bool operator>(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v > b._v};
  }

  template<typename F>
  bool operator>(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v > b};
  }

  template<typename F>
  bool operator>(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a > b._v};
  }

  template<typename F, typename T>
  bool operator>(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v > b};
  }

  template<typename F, typename T>
  bool operator>(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a > b._v};
  }


  // ********************************************************************************
  // greater_or_equals
  // ********************************************************************************

  template<typename F>
  bool operator>=(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v >= b._v};
  }

  template<typename F>
  bool operator>=(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v >= b};
  }

  template<typename F>
  bool operator>=(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a >= b._v};
  }

  template<typename F, typename T>
  bool operator>=(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v >= b};
  }

  template<typename F, typename T>
  bool operator>=(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a >= b._v};
  }


  // ********************************************************************************
  // inequals
  // ********************************************************************************

  template<typename F>
  bool operator!=(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v != b._v};
  }

  template<typename F>
  bool operator!=(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v != b};
  }

  template<typename F>
  bool operator!=(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a != b._v};
  }

  template<typename F, typename T>
  bool operator!=(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v != b};
  }

  template<typename F, typename T>
  bool operator!=(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a != b._v};
  }


  // ********************************************************************************
  // equals
  // ********************************************************************************

  template<typename F>
  bool operator==(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v == b._v};
  }

  template<typename F>
  bool operator==(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v == b};
  }

  template<typename F>
  bool operator==(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a == b._v};
  }

  template<typename F, typename T>
  bool operator==(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a._v == b};
  }

  template<typename F, typename T>
  bool operator==(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {a == b._v};
  }



  // ********************************************************************************
  // functions
  // ********************************************************************************

  template<typename F>
  OpCounter<F> exp(const OpCounter<F>& a)
  {
    ++OpCounter<F>::counters.exp_count;
    return {std::exp(a._v)};
  }

  template<typename F>
  OpCounter<F> pow(const OpCounter<F>& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.pow_count;
    return {std::pow(a._v,b._v)};
  }

  template<typename F>
  OpCounter<F> pow(const OpCounter<F>& a, const F& b)
  {
    ++OpCounter<F>::counters.pow_count;
    return {std::pow(a._v,b)};
  }

  template<typename F, typename T>
  OpCounter<F> pow(const OpCounter<F>& a, const T& b)
  {
    ++OpCounter<F>::counters.pow_count;
    return {std::pow(a._v,b)};
  }

  template<typename F>
  OpCounter<F> pow(const F& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.pow_count;
    return {std::pow(a,b._v)};
  }

  template<typename F, typename T>
  OpCounter<F> pow(const T& a, const OpCounter<F>& b)
  {
    ++OpCounter<F>::counters.pow_count;
    return {std::pow(a,b._v)};
  }

  template<typename F>
  OpCounter<F> sin(const OpCounter<F>& a)
  {
    ++OpCounter<F>::counters.sin_count;
    return {std::sin(a._v)};
  }

  template<typename F>
  OpCounter<F> cos(const OpCounter<F>& a)
  {
    ++OpCounter<F>::counters.sin_count;
    return {std::cos(a._v)};
  }

  template<typename F>
  OpCounter<F> sqrt(const OpCounter<F>& a)
  {
    ++OpCounter<F>::counters.sqrt_count;
    return {std::sqrt(a._v)};
  }

  template<typename F>
  OpCounter<F> abs(const OpCounter<F>& a)
  {
    ++OpCounter<F>::counters.comparison_count;
    return {std::abs(a._v)};
  }

}

#endif // __OPCOUNTER__
