#ifndef _LIGHT_RAT_H_
#define _LIGHT_RAT_H_

#include "soplex/rat.h"
#include <boost/multiprecision/number.hpp>
#include <cassert>
#include <gmp.h>
#include <gmpxx.h>
#include <ios>

#include <boost/multiprecision/gmp.hpp>

namespace soplex {

extern const mpq_class mpq_infinity;

/**
 * @brief Rational number with a compact small-value representation.
 *
 * A LightRat stores values in one 64-bit union. Values whose numerator and
 * denominator fit in the supported 32-bit ranges are stored directly as a
 * signed 32-bit numerator and a positive 31-bit denominator. The denominator
 * is stored shifted left by one bit; the low bit is therefore available as a
 * representation tag. A tagged pointer is used for larger values, which are
 * stored as a canonical GMP `mpq_t` allocated from an internal object store.
 *
 * The class owns the GMP value when the GMP representation is active. Copying
 * performs a deep copy and moving transfers ownership. Use `set()` or the
 * constructors to create values, `num()` and `den()` only for values for which
 * `is_rat32()` is true, and `get_gmp()` only for values for which
 * `is_ratgmp()` is true. `normalize()` may reduce a GMP value back to the
 * compact representation.
 *
 * LightRat is also a Boost.Multiprecision backend. Its scalar constructors,
 * two-component constructors, and the `number_category` specialization allow
 * it to be used as `number<LightRat, ...>`.
 */
class LightRat {
  /// Tag value for the compact numerator/denominator representation.
  static constexpr int IS_RAT32 = 0x0;
  /// Tag value OR-ed into a stored GMP pointer.
  static constexpr int IS_RATGMP = 0x1;

public:
  /// Scalar types accepted directly by Boost.Multiprecision.
#if BOOST_VERSION >= 107600
  using signed_types = std::tuple<int32_t>;
  using unsigned_types = std::tuple<uint32_t>;
  using float_types = std::tuple<double>;
#else
  typedef boost::mpl::list<int32_t> signed_types;
  typedef boost::mpl::list<uint32_t> unsigned_types;
  typedef boost::mpl::list<double> float_types;
#endif

  /// Return the positive infinity sentinel.
  static const LightRat &infinity();
  /// Return the negative infinity sentinel.
  static const LightRat &ninfinity();

  /// Construct zero.
  LightRat();
  /// Construct from the C representation.
  LightRat(const spx_lightrat_t rat);
  /// Construct from a numerator and positive denominator.
  LightRat(int32_t num, uint32_t den);
  LightRat(int64_t num, uint64_t den);
  /// Construct from a rational string such as `"3/4"`.
  LightRat(const std::string &str);
  LightRat(const char *str);
  /// Construct an integer or floating-point value.
  LightRat(int32_t num);
  LightRat(int64_t num);
  LightRat(uint32_t num);
  LightRat(uint64_t num);
  LightRat(double d);
  /// Construct from a GMP rational or integer value.
  /// @pre `q` must be canonicalized.
  LightRat(const mpq_t q);
  LightRat(const mpz_t z);
  LightRat(const mpz_t num, const mpz_t den);
  LightRat(const boost::multiprecision::gmp_int &z);
  LightRat(const boost::multiprecision::gmp_int &num,
           const boost::multiprecision::gmp_int &den);
  /// Construct from a GMP rational or integer value.
  /// @pre `q` must be canonicalized.
  LightRat(const mpq_class &q);
  LightRat(const mpz_class &z);
  /// Construct from GMP numerator and denominator values.
  LightRat(const mpz_class &num, const mpz_class &den);
  /// Copy and move constructors.
  LightRat(const LightRat &r2);
  LightRat(LightRat &&r2) noexcept;

  /// Assign from the C representation.
  /// @pre `rat` must be initialized.
  /// @pre `q` must be canonicalized.
  LightRat &operator=(const spx_lightrat_t rat);
  LightRat &operator=(const std::string &str);
  LightRat &operator=(const char *str);
  LightRat &operator=(int32_t num);
  LightRat &operator=(int64_t num);
  LightRat &operator=(uint32_t num);
  LightRat &operator=(uint64_t num);
  LightRat &operator=(double d);
  LightRat &operator=(const mpq_t q);
  LightRat &operator=(const mpz_t z);
  LightRat &operator=(const boost::multiprecision::gmp_int &z);
  LightRat &operator=(const mpz_class &z);
  LightRat &operator=(const mpq_class &q);
  LightRat &operator=(const LightRat &r2);
  LightRat &operator=(LightRat &&r2) noexcept;

  /// Add, subtract, multiply, or divide by another rational.
  LightRat &operator+=(const LightRat &r2);
  LightRat &operator-=(const LightRat &r2);
  LightRat &operator*=(const LightRat &r2);
  LightRat &operator/=(const LightRat &r2);
  /// Apply a floating-point arithmetic operation.
  template <class V> LightRat &operator+=(V r2) {
    return *this += LightRat{r2};
  }
  template <class V> LightRat &operator-=(V r2) {
    return *this -= LightRat{r2};
  }
  template <class V> LightRat &operator*=(V r2) {
    return *this *= LightRat{r2};
  }
  template <class V> LightRat &operator/=(V r2) {
    return *this /= LightRat{r2};
  }

  LightRat operator+(const LightRat &r2) { return LightRat(*this) += r2; }
  LightRat operator-(const LightRat &r2) { return LightRat(*this) -= r2; }
  LightRat operator*(const LightRat &r2) { return LightRat(*this) *= r2; }
  LightRat operator/(const LightRat &r2) { return LightRat(*this) /= r2; }
  template <class V> LightRat operator+(const V &r2) {
    return *this + LightRat(r2);
  }
  template <class V> LightRat operator-(const V &r2) {
    return *this - LightRat(r2);
  }
  template <class V> LightRat operator*(const V &r2) {
    return *this * LightRat(r2);
  }
  template <class V> LightRat operator/(const V &r2) {
    return *this / LightRat(r2);
  }

  /// Compare with any value accepted by a LightRat constructor.
  template <class V> bool operator<(const V &r2) const {
    return compare(r2) < 0;
  }
  template <class V> bool operator>(const V &r2) const {
    return compare(r2) > 0;
  }
  template <class V> bool operator<=(const V &r2) const {
    return compare(r2) <= 0;
  }
  template <class V> bool operator>=(const V &r2) const {
    return compare(r2) >= 0;
  }
  template <class V> bool operator==(const V &r2) const {
    return compare(r2) == 0;
  }
  template <class V> bool operator!=(const V &r2) const {
    return compare(r2) != 0;
  }

  /// Return whether the value is zero.
  bool is_zero() const {
    return is_ratgmp() ? mpq_sgn(get_gmp()) == 0 : u.s.num == 0;
  }

  /// Return whether the value is not zero.
  bool is_nonzero() const {
    return is_ratgmp() ? mpq_sgn(get_gmp()) != 0 : u.s.num != 0;
  }

  /// Return whether the value is strictly positive.
  bool is_pos() const {
    return is_rat32() ? u.s.num > 0 : mpq_sgn(get_gmp()) > 0;
  }

  /// Return whether the value is non-negative.
  bool is_nonneg() const {
    return is_rat32() ? u.s.num >= 0 : mpq_sgn(get_gmp()) >= 0;
  }

  /// Return whether the value is strictly negative.
  bool is_neg() const {
    return is_rat32() ? u.s.num < 0 : mpq_sgn(get_gmp()) < 0;
  }

  /// Return whether the value is non-positive.
  bool is_nonpos() const {
    return is_rat32() ? u.s.num <= 0 : mpq_sgn(get_gmp()) <= 0;
  }

  /// Classification predicates for one, infinity, and integer values.
  bool is_one() const;
  bool is_minus_one() const;
  bool is_integer() const;
  bool is_infinity() const;
  bool is_ninfinity() const;
  bool is_finite() const { return !is_infinity() && !is_ninfinity(); }

  /// Return the additive inverse.
  LightRat operator-() const;

  /// Destroy the value and release an owned GMP object, if any.
  ~LightRat();

  /// Exchange two values without copying their GMP payloads.
  void swap(LightRat &o) noexcept { std::swap(u, o.u); }

  /// Set a rational from an integer numerator and positive denominator.
  LightRat &set(int32_t num, uint32_t den);
  LightRat &set(int64_t num, uint64_t den);
  LightRat &set(const mpz_t num, const mpz_t den);
  LightRat &set(const boost::multiprecision::gmp_int &num,
                const boost::multiprecision::gmp_int &den);
  /// Set one of the infinity sentinel values.
  LightRat &set_infinity();
  LightRat &set_ninfinity();

  /// Release any GMP payload and reset the value to zero.
  void clear();
  /// Convert a GMP value to the compact representation when it fits.
  void normalize();

  /// Assign +1.
  void set_one() {
    clear();
    unchecked_set(1, 1);
  }

  /// Assign -1.
  void set_minus_one() {
    clear();
    unchecked_set(-1, 1);
  }

  /// Return whether the compact representation is active.
  bool is_rat32() const { return (u.p.gmp & IS_RATGMP) != IS_RATGMP; }
  /// Return whether the GMP representation is active.
  bool is_ratgmp() const { return (u.p.gmp & IS_RATGMP) == IS_RATGMP; }
  /// Return the owned GMP rational; valid only when `is_ratgmp()` is true.
  mpq_ptr get_gmp() const {
    assert(is_ratgmp());
    return reinterpret_cast<mpq_ptr>(u.p.gmp ^ IS_RATGMP);
  }
  mpz_ptr gmp_num() const {
    assert(is_ratgmp());
    return mpq_numref(get_gmp());
  }
  mpz_ptr gmp_den() const {
    assert(is_ratgmp());
    return mpq_denref(get_gmp());
  }
  int32_t num() const {
    assert(is_rat32() && u.s.den != 0);
    return u.s.num;
  }
  uint32_t den() const {
    assert(is_rat32() && u.s.den != 0);
    return u.s.den >> 1;
  }
  /// Expose the ABI-compatible C representation without transferring ownership.
  const spx_lightrat_u *data() const { return &u; }
  /// Return the ABI-compatible C representation by reference.
  const spx_lightrat_u &operator*() const { return u; }

  /// Explicitly convert to a double.
  explicit operator double() const;
  /// Explicitly convert to a GMP rational.
  explicit operator mpq_class() const;
  explicit operator long double() const;
  explicit operator long long() const;
  explicit operator unsigned long long() const;

  /// Return a canonical `numerator/denominator` representation.
  std::string str(std::streamsize /*digits*/,
                  std::ios_base::fmtflags /*f*/) const;

  /// Compare with another LightRat.
  int compare(const LightRat &r2) const;
  /// Compare with a GMP rational.
  int compare(const mpq_t q) const;
  /// Compare with a value convertible to LightRat.
  template <class V> int compare(V o) const { return compare(LightRat(o)); }

  /// Replace the value by its reciprocal; zero is an error.
  LightRat &invert();
  /// Negate the value.
  LightRat &negate();
  /// Return the absolute value.
  LightRat &abs();

private:
  /// Release the GMP allocation owned by this value.
  void release_mpq();
  /// Store a compact value without changing representation state.
  void unchecked_set(int32_t num, uint32_t den);
  /// Store an allocated GMP value and set its representation tag.
  void unchecked_set(mpq_ptr q);
  /// Ensure a GMP payload exists, converting the compact value if needed.
  mpq_ptr ensure_gmp();
  /// Convert the compact value to GMP representation.
  void convert_to_gmp();
  /// Reject arithmetic involving an infinity sentinel.
  void ensure_not_infinity(const LightRat &r2) const;

  /**
   * Internal 64-bit representation.
   *
   * When the low tag bit is zero, `u.s` contains a signed 32-bit numerator
   * and a positive 31-bit denominator stored as `den << 1`. When the low tag
   * bit is one, `u.p.gmp` contains an object-store pointer OR-ed with
   * `IS_RATGMP`, and the pointed-to object is a GMP `mpq_t`. The pointer is
   * aligned, so its low bit is available for this tag.
   */
  spx_lightrat_u u;
};

/// Stream a LightRat in rational form.
std::ostream &operator<<(std::ostream &os, const LightRat &r);

// These overloads bypass Boost's generic backend fallbacks, which otherwise
// try to implement the operation by recursively converting LightRat again.
inline void eval_add(LightRat &result, const LightRat &value) {
  result += value;
}
inline void eval_subtract(LightRat &result, const LightRat &value) {
  result -= value;
}
inline void eval_multiply(LightRat &result, const LightRat &value) {
  result *= value;
}
inline void eval_divide(LightRat &result, const LightRat &value) {
  result /= value;
}

inline void eval_add(LightRat &result, const LightRat &left,
                     const LightRat &right) {
  result = left;
  result += right;
}
inline void eval_subtract(LightRat &result, const LightRat &left,
                          const LightRat &right) {
  result = left;
  result -= right;
}
inline void eval_multiply(LightRat &result, const LightRat &left,
                          const LightRat &right) {
  result = left;
  result *= right;
}
inline void eval_divide(LightRat &result, const LightRat &left,
                        const LightRat &right) {
  result = left;
  result /= right;
}

//
// Conversions: must include at least unsigned long long, long long and long
// double. Everything else is entirely optional:
//
inline void eval_convert_to(unsigned long long *result,
                            const LightRat &backend) {
  *result = static_cast<unsigned long long>(backend);
}
inline void eval_convert_to(long long *result, const LightRat &backend) {
  *result = static_cast<long long>(backend);
}
inline void eval_convert_to(long double *result, const LightRat &backend) {
  *result = static_cast<long double>(backend);
}
inline void eval_convert_to(double *result, const LightRat &backend) {
  *result = static_cast<double>(backend);
}
inline void eval_abs(LightRat &result, const LightRat &val) {
  result = val;
  result.abs();
}

/**
 * @brief Assign a rational from numerator and denominator components.
 *
 * This overload is used by Boost.Multiprecision's two-component constructor.
 * When both component types match, the direct LightRat constructor is used;
 * otherwise the components are converted separately and divided.
 */
template <class V, class U>
void assign_components(LightRat &result, const V &numerator,
                       const U &denominator) {
  if (std::is_same<V, U>::value) {
    result = LightRat(numerator, denominator);
  } else {
    result = LightRat(numerator);
    result /= LightRat(denominator);
  }
}

} // namespace soplex

namespace boost {
namespace multiprecision {

#if BOOST_VERSION >= 107600
template <>
struct number_category<soplex::LightRat>
    : std::integral_constant<int, number_kind_rational> {};
#else
template <>
struct number_category<soplex::LightRat> : mpl::int_<number_kind_rational> {};
#endif

template <expression_template_option ET>
struct component_type<number<soplex::LightRat, ET>> {
  typedef number<gmp_int, ET> type;
};

template <expression_template_option ET>
number<gmp_int, ET> numerator(const number<soplex::LightRat, ET> &val) {
  number<gmp_int, ET> result;
  if (val.backend().is_infinity()) {
    mpz_set(result.backend().data(), soplex::mpq_infinity.get_num_mpz_t());
  } else if (val.backend().is_ninfinity()) {
    mpz_set(result.backend().data(), soplex::mpq_infinity.get_num_mpz_t());
    mpz_neg(result.backend().data(), result.backend().data());
  } else if (val.backend().is_rat32()) {
    mpz_set_si(result.backend().data(), val.backend().num());
  } else {
    mpz_set(result.backend().data(), val.backend().gmp_num());
  }
  return result;
}
template <expression_template_option ET>
number<gmp_int, ET> denominator(const number<soplex::LightRat, ET> &val) {
  number<gmp_int, ET> result;
  if (val.backend().is_infinity() || val.backend().is_ninfinity()) {
    mpz_set(result.backend().data(), soplex::mpq_infinity.get_den_mpz_t());
  } else if (val.backend().is_rat32()) {
    mpz_set_ui(result.backend().data(), val.backend().den());
  } else {
    mpz_set(result.backend().data(), val.backend().gmp_den());
  }
  return result;
}

} // namespace multiprecision
} // namespace boost

#endif
