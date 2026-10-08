
#include "soplex/lightrat.h"
#include "soplex/spxdefines.h"
#include "soplex_interface.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::gmp_int;
using soplex::LightRat;
using GmpInteger =
    boost::multiprecision::number<gmp_int, boost::multiprecision::et_off>;
using Rational =
    boost::multiprecision::number<LightRat, boost::multiprecision::et_off>;

namespace {

void expect(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

void expect_equal(const LightRat &actual, const mpq_class &expected,
                  const char *message) {
  expect(!actual.is_infinity() && !actual.is_ninfinity(), message);
  mpq_class canonical_expected(expected);
  canonical_expected.canonicalize();
  expect(static_cast<mpq_class>(actual) == canonical_expected, message);
}

void expect_throws(const std::function<void()> &operation,
                   const char *message) {
  bool thrown = false;
  try {
    operation();
  } catch (const std::domain_error &) {
    thrown = true;
  }
  expect(thrown, message);
}

void test_construction_and_representation() {
  LightRat zero;
  expect(zero.is_zero() && zero.is_rat32() && zero.str(0, {}) == "0/1",
         "default construction");

  LightRat reduced(18, 24);
  expect(reduced.is_rat32() && reduced.num() == 3 && reduced.den() == 4,
         "compact reduction");

  LightRat negative(-18, 24);
  expect(negative.num() == -3 && negative.den() == 4, "negative reduction");

  LightRat minimum_small(static_cast<int32_t>(-1073741823), 1);
  expect(minimum_small.is_rat32() && minimum_small.num() == -1073741823,
         "minimum compact integer");
  LightRat largest_denominator(1, static_cast<uint32_t>(1073741823));
  expect(largest_denominator.is_rat32() &&
             largest_denominator.den() == 1073741823,
         "maximum compact denominator");
  expect_throws([] { LightRat(1, static_cast<uint32_t>(0)); },
                "zero uint32 denominator");
  expect_throws(
      [] { LightRat(static_cast<int64_t>(1), static_cast<uint64_t>(0)); },
      "zero uint64 denominator");

  LightRat max_small(static_cast<int32_t>(1073741823), 1);
  expect(max_small.is_rat32() && max_small.num() == 1073741823,
         "maximum compact integer");

  LightRat large_integer(static_cast<int64_t>(1073741824));
  expect(large_integer.is_ratgmp(), "large integer promotes to GMP");
  expect_equal(large_integer, mpq_class(1073741824), "large integer value");
  large_integer.normalize();
  expect(large_integer.is_ratgmp(), "large integer remains GMP");

  LightRat large_denominator(1, static_cast<uint32_t>(1073741824));
  expect(large_denominator.is_ratgmp(), "large denominator promotes to GMP");
  expect_equal(large_denominator, mpq_class(1, 1073741824),
               "large denominator value");

  LightRat from_string("42/56");
  expect_equal(from_string, mpq_class(3, 4), "string construction");

  mpq_t q;
  mpq_init(q);
  mpq_set_str(q, "123456789012345678901/37", 10);
  mpq_canonicalize(q);
  LightRat from_mpq(q);
  expect(from_mpq.is_ratgmp(), "mpq construction uses GMP");
  expect(static_cast<mpq_class>(from_mpq) == mpq_class(q), "mpq construction");
  mpq_clear(q);

  // We request that mpq are canonized before being passed to LightRat
  // Hence, this test is unnecessary
  // mpq_init(q);
  // mpz_set_si(mpq_numref(q), 42);
  // mpz_set_ui(mpq_denref(q), 56);
  // LightRat noncanonical(q);
  // expect(noncanonical.is_rat32() && noncanonical.num() == 3 &&
  //            noncanonical.den() == 4,
  //        "noncanonical mpq assignment is reduced");
  // mpq_clear(q);

  mpz_class n("123456789012345678901234567890");
  mpz_class d("987654321");
  LightRat from_boost_gmp(n.get_mpz_t(), d.get_mpz_t());
  expect(from_boost_gmp.is_ratgmp(), "gmp_int construction uses GMP");
  mpq_class expected_gmp(n.get_str() + "/" + d.get_str());
  expected_gmp.canonicalize();
  expect(static_cast<mpq_class>(from_boost_gmp) == expected_gmp,
         "gmp_int construction");
}

void test_copy_move_and_assignment() {
  LightRat original(7, 11);
  LightRat copy(original);
  expect(copy == original && copy.data() != original.data(),
         "copy construction");

  LightRat assigned;
  assigned = original;
  expect(assigned == original, "copy assignment");
  original = 99;
  expect(assigned == LightRat(7, 11), "copy is independent");

  LightRat moved_from(mpz_class("999999999999999999999"));
  LightRat moved(std::move(moved_from));
  expect(moved.is_ratgmp() && moved_from.is_zero(), "move construction");

  LightRat move_assigned;
  move_assigned = std::move(moved);
  expect(move_assigned.is_ratgmp() && moved.is_zero(), "move assignment");

  LightRat gmp_value(mpz_class("999999999999999999999"));
  gmp_value = LightRat(5, 8);
  expect(gmp_value.is_rat32() && gmp_value == LightRat(5, 8),
         "assignment releases GMP storage");
}

void test_arithmetic() {
  LightRat a(2, 3);
  LightRat b(5, 7);
  LightRat result;

  result = a;
  result += b;
  expect_equal(result, mpq_class(29, 21), "addition");
  result = a;
  result -= b;
  expect_equal(result, mpq_class(-1, 21), "subtraction");
  result = a;
  result *= b;
  expect_equal(result, mpq_class(10, 21), "multiplication");
  result = a;
  result /= b;
  expect_equal(result, mpq_class(14, 15), "division");
  expect_equal(-a, mpq_class(-2, 3), "negation");

  result = 10;
  result += 0.5;
  expect_equal(result, mpq_class(21, 2), "double addition");
  result -= 0.5;
  expect_equal(result, mpq_class(10), "double subtraction");
  result *= 0.5;
  expect_equal(result, mpq_class(5), "double multiplication");
  result /= 0.5;
  expect_equal(result, mpq_class(10), "double division");

  result = LightRat(1073741823);
  result += LightRat(1);
  expect(result.is_ratgmp(), "addition overflow promotes to GMP");
  expect_equal(result, mpq_class(1073741824), "addition overflow value");

  result = LightRat(1073741823);
  result *= LightRat(2);
  expect(result.is_ratgmp(), "multiplication overflow promotes to GMP");
  expect_equal(result, mpq_class(2147483646), "multiplication overflow value");
  result.normalize();
  expect(result.is_ratgmp(), "unrepresentable GMP value stays GMP");

  LightRat gmp_a(mpz_class("123456789012345678901"));
  LightRat gmp_b(17, 19);
  result = gmp_a;
  result += gmp_b;
  expect_equal(result, mpq_class("123456789012345678901") + mpq_class(17, 19),
               "GMP plus compact");
  result = gmp_a;
  result *= gmp_b;
  expect_equal(result, mpq_class("123456789012345678901") * mpq_class(17, 19),
               "GMP times compact");

  expect_throws(
      [&] {
        result = LightRat(1);
        result /= LightRat(0);
      },
      "division by zero");
  expect_throws(
      [&] {
        LightRat zero;
        zero.invert();
      },
      "inverse of zero");
  mpz_class numerator(1);
  mpz_class zero_denominator(0);
  expect_throws(
      [&] { LightRat(numerator.get_mpz_t(), zero_denominator.get_mpz_t()); },
      "zero GMP denominator");

  LightRat gmp_zero(mpz_class("999999999999999999999"));
  gmp_zero = 0;
  expect(gmp_zero.is_zero(), "GMP value can be reset to zero");
  expect_throws(
      [&] {
        LightRat gmp_value(mpz_class("999999999999999999999"));
        gmp_value /= gmp_zero;
      },
      "division by zero after GMP promotion");

  LightRat inf = LightRat::infinity();
  LightRat ninf = LightRat::ninfinity();
  result = inf;
  result *= LightRat(2);
  expect(result.is_infinity(), "positive infinity times positive value");
  result = inf;
  result *= LightRat(-2);
  expect(result.is_ninfinity(), "positive infinity times negative value");
  result = ninf;
  result *= LightRat(2);
  expect(result.is_ninfinity(), "negative infinity times positive value");
  result = ninf;
  result *= LightRat(-2);
  expect(result.is_infinity(), "negative infinity times negative value");
  result = LightRat(2);
  result *= inf;
  expect(result.is_infinity(), "positive value times positive infinity");
  result = LightRat(-2);
  result *= inf;
  expect(result.is_ninfinity(), "negative value times positive infinity");
  result = inf;
  result /= LightRat(2);
  expect(result.is_infinity(), "positive infinity divided by positive value");
  result = inf;
  result /= LightRat(-2);
  expect(result.is_ninfinity(), "positive infinity divided by negative value");
  result = ninf;
  result /= LightRat(2);
  expect(result.is_ninfinity(), "negative infinity divided by positive value");
  result = ninf;
  result /= LightRat(-2);
  expect(result.is_infinity(), "negative infinity divided by negative value");
  result = inf;
  result += LightRat(7, 3);
  expect(result.is_infinity(), "positive infinity plus finite value");
  result = ninf;
  result += LightRat(-7, 3);
  expect(result.is_ninfinity(), "negative infinity plus finite value");
  result = LightRat(7, 3);
  result += inf;
  expect(result.is_infinity(), "finite value plus positive infinity");
  result = LightRat(-7, 3);
  result += ninf;
  expect(result.is_ninfinity(), "finite value plus negative infinity");
  result = inf;
  result -= LightRat(-7, 3);
  expect(result.is_infinity(), "positive infinity minus finite value");
  result = ninf;
  result -= LightRat(7, 3);
  expect(result.is_ninfinity(), "negative infinity minus finite value");
  result = LightRat(7, 3);
  result -= inf;
  expect(result.is_ninfinity(), "finite value minus positive infinity");
  result = LightRat(-7, 3);
  result -= ninf;
  expect(result.is_infinity(), "finite value minus negative infinity");
  result = inf;
  result += inf;
  expect(result.is_infinity(), "positive infinity plus positive infinity");
  result = ninf;
  result += ninf;
  expect(result.is_ninfinity(), "negative infinity plus negative infinity");
  result = inf;
  result -= ninf;
  expect(result.is_infinity(), "positive infinity minus negative infinity");
  result = ninf;
  result -= inf;
  expect(result.is_ninfinity(), "negative infinity minus positive infinity");
  expect_throws([&] { LightRat(0) *= inf; }, "zero times infinity");
  expect_throws([&] { inf *= LightRat(0); }, "infinity times zero");
  expect_throws(
      [&] {
        result = LightRat::infinity();
        result += LightRat::ninfinity();
      },
      "positive infinity plus negative infinity");
  expect_throws(
      [&] {
        result = LightRat::ninfinity();
        result += LightRat::infinity();
      },
      "negative infinity plus positive infinity");
  expect_throws(
      [&] {
        result = LightRat::infinity();
        result -= LightRat::infinity();
      },
      "positive infinity minus positive infinity");
  expect_throws(
      [&] {
        result = LightRat::ninfinity();
        result -= LightRat::ninfinity();
      },
      "negative infinity minus negative infinity");
  expect_throws([&] { inf /= inf; }, "infinity divided by infinity");
  expect_throws([&] { LightRat(2) /= inf; },
                "finite value divided by infinity");
}

void test_conversions_and_comparisons() {
  LightRat value(7, 3);
  expect(static_cast<double>(value) == mpq_class(7, 3).get_d(),
         "double conversion");
  expect(static_cast<mpq_class>(value) == mpq_class(7, 3),
         "mpq_class conversion");
  expect(value > 2 && value > 2.0 && value != 2, "mixed comparisons");
  expect(value == mpq_class(7, 3), "mpq comparison");
  expect(value.is_pos() && value.is_nonneg() && !value.is_neg(),
         "positive predicates");
  expect(!value.is_integer() && !value.is_zero(), "value predicates");

  LightRat integer(-4);
  expect(integer.is_neg() && integer.is_nonpos() && integer.is_integer(),
         "negative integer predicates");
  integer.clear();
  expect(integer.is_zero() && integer.is_rat32(), "clear");
  integer.set_one();
  expect(integer.is_one(), "set_one");
  integer.set_minus_one();
  expect(integer.is_minus_one(), "set_minus_one");

  LightRat infinity = LightRat::infinity();
  LightRat ninfinity = LightRat::ninfinity();
  expect(infinity.is_infinity() && ninfinity.is_ninfinity(),
         "infinity classification");
  expect(infinity > LightRat(999999) && ninfinity < LightRat(-999999),
         "infinity comparisons");
  expect(infinity.str(0, {}) == "inf" && ninfinity.str(0, {}) == "-inf",
         "infinity formatting");
  expect(static_cast<double>(infinity) == soplex::infinity &&
             static_cast<double>(ninfinity) == -soplex::infinity,
         "infinity double conversion");
  LightRat integral(123456789);
  LightRat fractional(7, 3);
  expect(static_cast<long long>(integral) == 123456789LL,
         "long long conversion");
  expect(static_cast<unsigned long long>(integral) == 123456789ULL,
         "unsigned long long conversion");
  expect(static_cast<long long>(fractional) == 2LL,
         "long long conversion truncates toward zero");
  expect(static_cast<unsigned long long>(LightRat(7)) == 7ULL,
         "unsigned long long conversion");
  infinity += LightRat(1);
  expect(infinity.is_infinity(), "infinity plus finite value");
  expect_throws([&] { LightRat(1) /= infinity; }, "division by infinity");
  expect(-infinity == ninfinity, "infinity negation");
  expect(infinity != ninfinity, "infinity distinction");

  const double infinity_value = soplex::infinity;
  LightRat from_inf_double(infinity_value);
  LightRat from_ninf_double(-infinity_value);
  expect(from_inf_double.is_infinity() && from_ninf_double.is_ninfinity(),
         "infinity double construction");
}

void test_c_interface() {
  spx_lightrat_t c_value;
  spx_rat_init(c_value);
  spx_rat_set_si(c_value, 6, 8);
  expect(spx_rat_is_rat32(c_value), "C compact representation");
  expect(spx_rat_get_num(c_value) == 3 && spx_rat_get_den(c_value) == 4,
         "C compact access");

  mpq_t q;
  mpq_init(q);
  spx_rat_to_gmp(c_value, q);
  expect(mpq_cmp_si(q, 3, 4) == 0, "C GMP export");
  spx_rat_negate(c_value);
  spx_rat_mul_si(c_value, 2, 3);
  spx_rat_invert(c_value);
  spx_rat_to_gmp(c_value, q);
  expect(mpq_cmp_si(q, -2, 1) == 0, "C arithmetic");

  mpq_set_str(q, "123456789012345678901/5", 10);
  mpq_canonicalize(q);
  spx_rat_set_gmp(c_value, q);
  expect(spx_rat_is_gmp(c_value), "C GMP assignment");
  LightRat from_c(c_value);
  expect(static_cast<mpq_class>(from_c) == mpq_class(q), "C round trip");

  spx_rat_set_d(c_value, 0.125);
  expect(spx_rat_get_d(c_value) == 0.125, "C double conversion");
  spx_rat_clear(c_value);
  mpq_clear(q);

  void *lp = SoPlex_createLpProblemRational(1, 1);
  spx_lightrat_t lb, ub, obj;
  spx_rat_init(lb);
  spx_rat_init(ub);
  spx_rat_init(obj);
  SoPlex_addColToLpProblemRational(lp, nullptr, nullptr, 0, obj, lb, ub);
  spx_rat_clear(obj);
  spx_rat_clear(lb);
  spx_rat_clear(ub);
  SoPlex_freeLPProblemRational(lp);
}

void test_boost_backend() {
  Rational r(3, 4);
  Rational n(GmpInteger(7), GmpInteger(5));
  expect(boost::multiprecision::numerator(r) == 3 &&
             boost::multiprecision::denominator(r) == 4,
         "Boost numerator and denominator");
  expect(n == Rational(7, 5), "Boost two-component construction");

  Rational a(2, 3);
  Rational b(5, 7);
  expect(a + b == Rational(29, 21), "Boost addition");
  expect(a - b == Rational(-1, 21), "Boost subtraction");
  expect(a * b == Rational(10, 21), "Boost multiplication");
  expect(a / b == Rational(14, 15), "Boost division");
  a *= b;
  expect(a == Rational(10, 21), "Boost compound multiplication");
  a /= b;
  expect(a == Rational(2, 3), "Boost compound division");
  a += b;
  expect(a == Rational(29, 21), "Boost compound addition");
  a -= b;
  expect(a == Rational(2, 3), "Boost compound subtraction");

  Rational equivalent(4, 6);
  Rational smaller(-5, 7);
  Rational zero(0);
  expect(a == equivalent && equivalent == a, "Boost equivalent comparison");
  expect(a != b && b != a, "Boost inequality comparison");
  expect(smaller < zero && smaller <= zero, "Boost less-than comparison");
  expect(zero > smaller && zero >= smaller, "Boost greater-than comparison");
  expect(a > zero && a >= equivalent, "Boost positive ordering comparison");
  expect(smaller < a && smaller <= a, "Boost mixed ordering comparison");
  expect(a == LightRat(2, 3) && a != LightRat(5, 7),
         "Boost comparison with LightRat");
  expect(a > 0 && a >= 0 && a != 0, "Boost comparison with integer");
  expect(zero == 0 && zero <= 0 && zero >= 0,
         "Boost zero comparison with integer");

  Rational positive_infinity(LightRat::infinity());
  Rational negative_infinity(LightRat::ninfinity());
  expect(positive_infinity == Rational(LightRat::infinity()),
         "Boost positive infinity equality");
  expect(negative_infinity == Rational(LightRat::ninfinity()),
         "Boost negative infinity equality");
  expect(positive_infinity != negative_infinity, "Boost infinity inequality");
  expect(negative_infinity < zero && negative_infinity <= zero,
         "Boost negative infinity ordering");
  expect(positive_infinity > zero && positive_infinity >= zero,
         "Boost positive infinity ordering");
  expect(negative_infinity < positive_infinity,
         "Boost signed infinity ordering");
  expect(positive_infinity > a && negative_infinity < a,
         "Boost infinity ordering against finite value");
  expect(a * positive_infinity == positive_infinity,
         "Boost finite times positive infinity");
  expect(a * negative_infinity == negative_infinity,
         "Boost finite times negative infinity");
  expect(positive_infinity * Rational(-2) == negative_infinity,
         "Boost positive infinity times negative value");
  expect(negative_infinity * Rational(-2) == positive_infinity,
         "Boost negative infinity times negative value");
  expect(positive_infinity / Rational(-2) == negative_infinity,
         "Boost positive infinity divided by negative value");
  expect(negative_infinity / Rational(-2) == positive_infinity,
         "Boost negative infinity divided by negative value");
  expect(positive_infinity + a == positive_infinity,
         "Boost positive infinity plus finite value");
  expect(negative_infinity + a == negative_infinity,
         "Boost negative infinity plus finite value");
  expect(a + positive_infinity == positive_infinity,
         "Boost finite value plus positive infinity");
  expect(a - positive_infinity == negative_infinity,
         "Boost finite value minus positive infinity");
  expect(a - negative_infinity == positive_infinity,
         "Boost finite value minus negative infinity");
  expect(positive_infinity - negative_infinity == positive_infinity,
         "Boost positive infinity minus negative infinity");
  expect(positive_infinity + positive_infinity == positive_infinity,
         "Boost positive infinity plus positive infinity");
  expect(negative_infinity + negative_infinity == negative_infinity,
         "Boost negative infinity plus negative infinity");
  expect_throws([&] { Rational(0) * positive_infinity; },
                "Boost zero times infinity");
  expect_throws([&] { positive_infinity + negative_infinity; },
                "Boost opposite infinity addition");
  expect_throws([&] { positive_infinity - positive_infinity; },
                "Boost same infinity subtraction");
  expect_throws([&] { positive_infinity / positive_infinity; },
                "Boost infinity divided by infinity");

  cpp_dec_float_50 decimal(r);
  expect(decimal == cpp_dec_float_50("0.75"), "Boost decimal conversion");
  Rational round_trip(decimal);
  expect(static_cast<mpq_class>(round_trip.backend()) == mpq_class("3/4"),
         "Boost backend round trip");
}

} // namespace

int main() {
  try {
    test_construction_and_representation();
    test_copy_move_and_assignment();
    test_arithmetic();
    test_conversions_and_comparisons();
    test_c_interface();
    test_boost_backend();
  } catch (const std::exception &error) {
    std::cerr << "LightRat test failure: " << error.what() << '\n';
    return 1;
  }
  std::cout << "All LightRat tests passed.\n";
  return 0;
}
