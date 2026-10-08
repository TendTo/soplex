#include "soplex/lightrat.h"

#include "gmp.h"
#include "gmpxx.h"
#include "soplex/spxdefines.h"
#include <cassert>
#include <memory>
#include <ostream>

namespace soplex {
/* the test for unicity in the denominator occurs frequently. the denominator is
 * one if it is 2 :-) */
constexpr int ONE_DEN = 0x2;
/* Number of rationals in a block */
constexpr int MPQ_BLOCK_COUNT = 1024;
constexpr int INF_DEN = 0;
constexpr int INF_NUM = 1;
constexpr int NINF_NUM = -1;

#define MAX_NUMERATOR (INT32_MAX >> 1)
#define MIN_NUMERATOR (-MAX_NUMERATOR)
#define MAX_DENOMINATOR MAX_NUMERATOR

const mpq_class mpq_infinity{soplex::infinity};

namespace {

/*
 * Initialize q with enough room for numerator
 * and denominator with n bits. Not to be used
 * if q is already initialized.
 */
void mpq_init2(mpq_t q, const unsigned long n) {
  mpz_ptr num_q = mpq_numref(q);
  mpz_ptr den_q = mpq_denref(q);
  mpz_init2(num_q, n);
  mpz_init2(den_q, n);
  mpz_set_ui(den_q, 1UL);
}

void mpq_add_si(mpq_t q, const long num, const unsigned long den) {
  unsigned long gcd;
  mpz_t z0;

  mpz_init(z0);

  mpz_ptr num_q = mpq_numref(q);
  mpz_ptr den_q = mpq_denref(q);

  // special cases: den = 1 should be common
  //                num = 0 should be rare
  if (den == 1) {
    // a/b + d --> (a + bd)/b
    mpz_mul_si(z0, den_q, num);
    mpz_add(num_q, num_q, z0);

    goto clean_up;
  }

  gcd = mpz_gcd_ui(nullptr, den_q, den);
  //  printf("gcd = %lu\n", gcd);

  if (gcd == 1) {
    // a/b + c/d  --> (a d + b c) / bd
    mpz_mul_si(z0, den_q, num);
    mpz_mul_ui(num_q, num_q, den);
    mpz_add(num_q, num_q, z0);
    mpz_mul_ui(den_q, den_q, den);

    goto clean_up;
  }

  mpz_divexact_ui(den_q, den_q, gcd); // b0 = b/gcd
  mpz_mul_si(z0, den_q, num);         // z0 = b0 * num

  mpz_mul_ui(num_q, num_q, den / gcd);
  mpz_add(num_q, num_q, z0); // num_q := (den_q/gcd) * num + (den/gcd) * num_q

  gcd = mpz_gcd_ui(nullptr, num_q, gcd);
  if (gcd == 1) {
    mpz_mul_ui(den_q, den_q, den);
  } else {
    mpz_divexact_ui(num_q, num_q, gcd);
    mpz_mul_ui(den_q, den_q, den / gcd);
  }

clean_up:
  mpz_clear(z0);
}

void mpq_sub_si(mpq_t q, const long num, const unsigned long den) {
  mpq_add_si(q, -num, den);
}

/*
 * Multiply q by num/den
 * - num must be more than LONG_MIN
 * - den must be non-zero
 * - num and den must have no common factor
 */
void mpq_mul_si(mpq_t q, const long num, unsigned long den) {
  if (num == 0) {
    mpq_set_ui(q, 0, 1);
    return;
  }

  mpz_ptr num_q = mpq_numref(q);
  mpz_ptr den_q = mpq_denref(q);
  unsigned long abs_num = static_cast<unsigned long>(std::labs(num));

  if (abs_num != 1) {
    const unsigned long gcd = mpz_gcd_ui(nullptr, den_q, abs_num);
    abs_num /= gcd;
    mpz_divexact_ui(den_q, den_q, gcd);
  }

  if (den != 1) {
    const unsigned long gcd = mpz_gcd_ui(nullptr, num_q, den);
    den /= gcd;
    mpz_divexact_ui(num_q, num_q, gcd);
  }

  mpz_mul_ui(num_q, num_q, abs_num);
  mpz_mul_ui(den_q, den_q, den);
  if (num < 0) {
    mpz_neg(num_q, num_q);
  }
}

/*
 * Divide q by num/den
 * - num must be more than LONG_MIN
 * - den must be non-zero
 * - num and den must have no common factor
 */
void mpq_div_si(mpq_t q, const long num, unsigned long den) {
  mpz_ptr num_q = mpq_numref(q);
  mpz_ptr den_q = mpq_denref(q);
  unsigned long abs_num = static_cast<unsigned long>(std::labs(num));

  if (abs_num != 1) {
    const unsigned long gcd = mpz_gcd_ui(nullptr, num_q, abs_num);
    abs_num /= gcd;
    mpz_divexact_ui(num_q, num_q, gcd);
  }

  if (den != 1) {
    const unsigned long gcd = mpz_gcd_ui(nullptr, den_q, den);
    den /= gcd;
    mpz_divexact_ui(den_q, den_q, gcd);
  }

  mpz_mul_ui(num_q, num_q, den);
  mpz_mul_ui(den_q, den_q, abs_num);
  if (num < 0) {
    mpz_neg(num_q, num_q);
  }
}

class MpqStore {
  /*
   * mpq_link is mpq together with a pointer to the "next" free mpq_link
   * needed to prevent "modifying the contents" of the mpq, as is done
   * in utils/object_stores.[ch]
   */
  struct mpq_link_t {
    mpq_t mpq;
    mpq_link_t *next;
  };

  /*
   * Bank = a block of mpq (links)
   * - blocks are linked in a list
   */
  struct mpq_bank_t {
    std::unique_ptr<mpq_bank_t> next;
    mpq_link_t block[MPQ_BLOCK_COUNT] = {};
  };

  /*
   * Initialize store s for object of the given size
   * - n = number of objects in each block
   */
  MpqStore() {
    bnk = nullptr;
    free_list = nullptr;
    free_index = 0;
  }

  /*
   * Delete the full store: all banks are freed
   */
  ~MpqStore() {
    std::unique_ptr<mpq_bank_t> b = std::move(bnk);
    uint32_t k = free_index;
    while (b != nullptr) {
      for (uint32_t i = k; i < MPQ_BLOCK_COUNT; i++) {
        mpq_link_t *obj = b->block + i;
        mpq_clear(obj->mpq);
      }
      k = 0;
      b = std::move(b->next);
    }

    bnk = nullptr;
    free_list = nullptr;
    free_index = 0;
  }

public:
  static mpq_ptr alloc() { return instance().alloc_impl(); }
  static void free(mpq_ptr mpq) { instance().free_impl(mpq); }

private:
  static MpqStore &instance() {
    static MpqStore instance;
    return instance;
  }
  /*
   * Allocate an mpq
   */
  mpq_ptr alloc_impl() {
    mpq_link_t *obj = free_list;

    if (obj != nullptr) {
      free_list = obj->next;
      // sanity check: when returned it should still be nullptr
      obj->next = nullptr;
      return reinterpret_cast<mpq_ptr>(&obj->mpq);
    }

    uint32_t i = free_index;
    if (i == 0) {
      std::unique_ptr<mpq_bank_t> new_bank = std::make_unique<mpq_bank_t>();
      new_bank->next = std::move(bnk);
      bnk = std::move(new_bank);
      i = MPQ_BLOCK_COUNT;
    }

    i--;
    free_index = i;
    obj = bnk->block + i;

    // only initialize when we give it out
    mpq_init2(obj->mpq, 64);
    obj->next =
        nullptr; // sanity check: when returned it should still be nullptr

    return reinterpret_cast<mpq_ptr>(&obj->mpq); // same as obj
  }

  /*
   * Free an allocated mpq_link object: add it to free_list.
   * next pointer is stored in mpq_link->next
   */
  void free_impl(mpq_ptr mpq) {
    mpq_link_t *const obj = reinterpret_cast<mpq_link_t *>(mpq);
    // sanity check: it's being returned; it should still be nullptr
    assert(obj->next == nullptr);
    obj->next = free_list;
    free_list = obj;
  }
  /*
   * Store = a list of blocks
   * - allocations are performed via a free list,
   * - or in the first block of the bank list,
   * - or by adding a new block.
   */
  std::unique_ptr<mpq_bank_t> bnk; // first block in the bank list
  mpq_link_t *free_list;           // list of free mpq_links
  uint32_t free_index; // index of last allocated mpq_link in first block
};

} // namespace

const LightRat &LightRat::infinity() {
  static LightRat infinity{soplex::infinity};
  return infinity;
}
const LightRat &LightRat::ninfinity() {
  static LightRat ninfinity{-soplex::infinity};
  return ninfinity;
}

LightRat::LightRat()
    : u{.s{
#ifdef WORDS_BIGENDIAN
          .num = 0,
          .den = ONE_DEN,
#else
          .den = ONE_DEN,
          .num = 0,
#endif
      }} {
}
LightRat::LightRat(const spx_lightrat_t rat) : LightRat{} { *this = rat; }
LightRat::LightRat(const int32_t num, const uint32_t den) : LightRat{} {
  set(num, den);
}
LightRat::LightRat(const int64_t num, const uint64_t den) : LightRat{} {
  set(num, den);
}
LightRat::LightRat(const std::string &str) : LightRat{} { *this = str; }
LightRat::LightRat(const char *str) : LightRat{} { *this = str; }
LightRat::LightRat(const int32_t num) : LightRat{} { *this = num; }
LightRat::LightRat(const int64_t num) : LightRat{} { *this = num; }
LightRat::LightRat(const uint32_t num) : LightRat{} { *this = num; }
LightRat::LightRat(const uint64_t num) : LightRat{} { *this = num; }
LightRat::LightRat(const double d) : LightRat{} { *this = d; }
LightRat::LightRat(const mpq_t q) : LightRat{} { *this = q; }
LightRat::LightRat(const mpz_t z) : LightRat{} { *this = z; }
LightRat::LightRat(const mpz_t num, const mpz_t den) : LightRat{} {
  set(num, den);
}
LightRat::LightRat(const gmp_int &z) : LightRat{} { *this = z; }
LightRat::LightRat(const gmp_int &num, const gmp_int &den) : LightRat{} {
  set(num, den);
}
LightRat::LightRat(const mpq_class &q) : LightRat{q.get_mpq_t()} {}
LightRat::LightRat(const mpz_class &z) : LightRat{z.get_mpz_t()} {}
LightRat::LightRat(const mpz_class &num, const mpz_class &den) : LightRat{} {
  set(num.get_mpz_t(), den.get_mpz_t());
}
LightRat::LightRat(const LightRat &r2) : LightRat{} { *this = r2; }
LightRat::LightRat(LightRat &&r2) noexcept : LightRat{} {
  *this = std::move(r2);
}

LightRat &LightRat::operator=(const spx_lightrat_t rat) {
  return *this = *reinterpret_cast<const LightRat *>(rat);
}
LightRat &LightRat::operator=(const std::string &str) {
  return *this = str.c_str();
}
LightRat &LightRat::operator=(const char *const str) {
  mpq_ptr q = ensure_gmp();
  if (0 != mpq_set_str(q, str, 10)) {
    throw std::runtime_error(
        std::string("The string \"") + str +
        std::string("\" could not be interpreted as a valid rational number."));
  }
  mpq_canonicalize(q);
  normalize();
  return *this;
}
LightRat &LightRat::operator=(const int32_t num) {
  if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR) {
    clear();
    unchecked_set(num, 1);
  } else {
    mpq_ptr q = ensure_gmp();
    mpq_set_si(q, num, 1);
  }
  return *this;
}
LightRat &LightRat::operator=(const int64_t num) {
  if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR) {
    clear();
    unchecked_set(static_cast<int32_t>(num), 1);
  } else {
    mpq_ptr q = ensure_gmp();
    mpq_set_si(q, num, 1);
  }
  return *this;
}
LightRat &LightRat::operator=(const uint32_t num) {
  if (num <= MAX_NUMERATOR) {
    clear();
    unchecked_set(static_cast<int32_t>(num), 1);
  } else {
    mpq_ptr q = ensure_gmp();
    mpq_set_ui(q, num, 1);
  }
  return *this;
}
LightRat &LightRat::operator=(const uint64_t num) {
  if (num <= MAX_NUMERATOR) {
    clear();
    unchecked_set(static_cast<int32_t>(num), 1);
  } else {
    mpq_ptr q = ensure_gmp();
    mpq_set_ui(q, num, 1);
  }
  return *this;
}
LightRat &LightRat::operator=(const double d) {
  if (d >= soplex::infinity) {
    set_infinity();
    return *this;
  }
  if (d <= -soplex::infinity) {
    set_ninfinity();
    return *this;
  }

  mpq_ptr qt = ensure_gmp();
  mpq_set_d(qt, d);
  normalize();
  return *this;
}
LightRat &LightRat::operator=(const mpq_class &q) {
  return *this = q.get_mpq_t();
}
LightRat &LightRat::operator=(const mpq_t q) {
  // If the numerator and denominator of q fit in 32 bits,
  // we can store it as a pair of integers
  if (mpz_fits_ulong_p(mpq_denref(q)) && mpz_fits_slong_p(mpq_numref(q))) {
    const long num = mpz_get_si(mpq_numref(q));
    const unsigned long den = mpz_get_ui(mpq_denref(q));
    if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR &&
        den <= MAX_DENOMINATOR) {
      clear();
      unchecked_set(static_cast<int32_t>(num), static_cast<uint32_t>(den));
      assert(is_rat32());
      return *this;
    }
  }

  // Otherwise, we need to allocate a new mpq number and copy q into it
  mpq_ptr qt = ensure_gmp();
  mpq_set(qt, q);
  assert(is_ratgmp());
  return *this;
}
LightRat &LightRat::operator=(const mpz_class &z) {
  return *this = z.get_mpz_t();
}
LightRat &LightRat::operator=(const mpz_t z) {
  if (mpz_fits_slong_p(z)) {
    const long num = mpz_get_si(z);
    if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR) {
      clear();
      unchecked_set(static_cast<int32_t>(num), 1);
      assert(is_rat32());
      return *this;
    }
  }

  mpq_ptr q = ensure_gmp();
  mpq_set_z(q, z);
  assert(is_ratgmp());
  return *this;
}
LightRat &LightRat::operator=(const gmp_int &z) { return *this = z.data(); }
LightRat &LightRat::operator=(const LightRat &r2) {
  if (this == &r2) {
    return *this;
  }

  if (r2.is_ratgmp()) {
    mpq_ptr q1 = ensure_gmp();
    mpq_ptr q2 = r2.get_gmp();
    mpq_set(q1, q2);
  } else {
    clear();
    u = r2.u;
  }
  return *this;
}
LightRat &LightRat::operator=(LightRat &&r2) noexcept {
  if (this == &r2) {
    return *this;
  }
  clear(); // Ensure we clean up the current object, releasing any memory
  std::swap(u, r2.u);
  return *this;
}

LightRat &LightRat::operator+=(const LightRat &r2) {
  if ((is_infinity() && !r2.is_ninfinity()) ||
      (is_ninfinity() && !r2.is_infinity())) {
    // r1 is infinity, and r2 is not the opposite infinity, so the result is r1
    return *this;
  }
  if ((r2.is_infinity() && !is_ninfinity()) ||
      (r2.is_ninfinity() && !is_infinity())) {
    // r2 is infinity, and r1 is not the opposite infinity, so the result is r2
    return *this = r2;
  }
  ensure_not_infinity(r2);

  if (u.s.den == ONE_DEN && r2.u.s.den == ONE_DEN) {
    assert(r2.is_rat32() && is_rat32());
    u.s.num += r2.u.s.num;
    if (u.s.num < MIN_NUMERATOR || u.s.num > MAX_NUMERATOR) {
      convert_to_gmp();
    }
    return *this;
  }

  if (r2.is_ratgmp()) {
    if (!is_ratgmp()) {
      convert_to_gmp();
    }
    mpq_ptr q1 = get_gmp();
    mpq_ptr q2 = r2.get_gmp();
    mpq_add(q1, q1, q2);
  } else if (is_ratgmp()) {
    mpq_ptr q1 = get_gmp();
    mpq_add_si(q1, r2.num(), r2.den());
  } else {
    const uint64_t d = den() * static_cast<uint64_t>(r2.den());
    const int64_t n = den() * static_cast<int64_t>(r2.num()) +
                      r2.den() * static_cast<int64_t>(num());
    set(n, d);
  }
  return *this;
}
LightRat &LightRat::operator-=(const LightRat &r2) {
  if (is_infinity() || is_ninfinity()) {
    if (r2.is_infinity() || r2.is_ninfinity()) {
      if (is_infinity() == r2.is_infinity()) {
        throw std::domain_error(
            "infinity minus infinity in LightRat::operator-=");
      }
      return *this;
    }
    // Infinity minus a finite value keeps its sign.
    return *this;
  }
  if (r2.is_infinity() || r2.is_ninfinity()) {
    // A finite value minus infinity has the opposite sign.
    return *this = -r2;
  }

  if (u.s.den == ONE_DEN && r2.u.s.den == ONE_DEN) {
    assert(is_rat32() && r2.is_rat32());
    u.s.num -= r2.u.s.num;
    if (u.s.num < MIN_NUMERATOR || u.s.num > MAX_NUMERATOR) {
      convert_to_gmp();
    }
    return *this;
  }

  if (r2.is_ratgmp()) {
    if (!is_ratgmp()) {
      convert_to_gmp();
    }
    mpq_ptr q1 = get_gmp();
    mpq_ptr q2 = r2.get_gmp();
    mpq_sub(q1, q1, q2);
  } else if (is_ratgmp()) {
    mpq_ptr q1 = get_gmp();
    mpq_sub_si(q1, r2.num(), r2.den());
  } else {
    const uint64_t d = den() * static_cast<uint64_t>(r2.den());
    const int64_t n = r2.den() * static_cast<int64_t>(num()) -
                      den() * static_cast<int64_t>(r2.num());
    set(n, d);
  }
  return *this;
}
LightRat &LightRat::operator*=(const LightRat &r2) {
  if (!is_finite()) {
    if (r2.is_zero()) {
      throw std::domain_error("infinity * zero in LightRat::operator*=");
    }
    return r2.is_pos() ? *this : this->negate();
  }
  if (!r2.is_finite()) {
    if (is_zero()) {
      throw std::domain_error("zero * infinity in LightRat::operator*=");
    }
    return is_pos() ? (*this = r2) : (*this = -r2);
  }
  ensure_not_infinity(r2);

  if (u.s.den == ONE_DEN && r2.u.s.den == ONE_DEN) {
    assert(is_rat32() && r2.is_rat32());
    const int64_t n = u.s.num * static_cast<int64_t>(r2.u.s.num);
    if (MIN_NUMERATOR <= n && n <= MAX_NUMERATOR) {
      u.s.num = static_cast<int32_t>(n);
    } else {
      mpq_ptr q = ensure_gmp();
      mpq_set_si(q, n, 1);
    }
    return *this;
  }

  if (r2.is_ratgmp()) {
    if (is_rat32()) {
      convert_to_gmp();
    }
    mpq_ptr q1 = get_gmp();
    mpq_ptr q2 = r2.get_gmp();
    mpq_mul(q1, q1, q2);
  } else if (is_ratgmp()) {
    mpq_ptr q1 = get_gmp();
    mpq_mul_si(q1, r2.num(), r2.den());
  } else {
    const uint64_t d = den() * static_cast<uint64_t>(r2.den());
    const int64_t n = num() * static_cast<int64_t>(r2.num());
    set(n, d);
  }
  return *this;
}
LightRat &LightRat::operator/=(const LightRat &r2) {
  if (r2.is_zero()) {
    throw std::domain_error("Division by zero in LightRat::operator/=");
  }
  if (!is_finite() && !r2.is_finite()) {
    throw std::domain_error("infinity / infinity in LightRat::operator/=");
  }

  if (!is_finite()) {
    assert(!r2.is_zero());
    return r2.is_pos() ? *this : this->negate();
  }
  ensure_not_infinity(r2);

  if (r2.is_ratgmp()) {
    if (is_rat32()) {
      convert_to_gmp();
    }
    mpq_ptr q1 = get_gmp();
    mpq_ptr q2 = r2.get_gmp();
    mpq_div(q1, q1, q2);
  } else if (is_ratgmp()) {
    mpq_ptr q1 = get_gmp();
    mpq_div_si(q1, r2.num(), r2.den());
  } else if (r2.num() > 0) {
    const uint64_t d = den() * static_cast<uint64_t>(r2.num());
    const int64_t n = num() * static_cast<int64_t>(r2.den());
    set(n, d);
  } else if (r2.num() < 0) {
    const uint64_t d = den() * static_cast<uint64_t>(-r2.num());
    const int64_t n = num() * -static_cast<int64_t>(r2.den());
    set(n, d);
  } else {
    throw std::domain_error("Division by zero in LightRat::operator/=");
  }
  return *this;
}

LightRat LightRat::operator-() const {
  LightRat out(*this);
  if (out.is_ratgmp()) {
    mpq_ptr q = out.get_gmp();
    mpq_neg(q, q);
  } else {
    out.u.s.num = -out.u.s.num;
  }
  return out;
}

bool LightRat::is_one() const {
  return (u.s.den == ONE_DEN && u.s.num == 1) ||
         (is_ratgmp() && mpq_cmp_ui(get_gmp(), 1, 1) == 0);
}
bool LightRat::is_minus_one() const {
  return (u.s.den == ONE_DEN && u.s.num == -1) ||
         (is_ratgmp() && mpq_cmp_ui(get_gmp(), -1, 1) == 0);
}
bool LightRat::is_integer() const {
  return (is_rat32() && u.s.den == ONE_DEN) ||
         (is_ratgmp() && mpz_cmp_ui(mpq_denref(get_gmp()), 1) == 0);
}
bool LightRat::is_infinity() const {
  return is_rat32() && u.s.den == INF_DEN && u.s.num == INF_NUM;
}
bool LightRat::is_ninfinity() const {
  return is_rat32() && u.s.den == INF_DEN && u.s.num == NINF_NUM;
}

void LightRat::convert_to_gmp() {
  if (is_ratgmp()) {
    return; // Already in GMP representation
  }

  assert(is_rat32());
  mpq_ptr q = MpqStore::alloc();
  if (is_infinity()) {
    mpq_set(q, mpq_infinity.get_mpq_t());
  } else if (is_ninfinity()) {
    mpq_set(q, mpq_infinity.get_mpq_t());
    mpq_neg(q, q);
  } else {
    mpq_set_si(q, num(), den());
    mpq_canonicalize(q);
  }
  unchecked_set(q);
}
void LightRat::ensure_not_infinity(const LightRat &r2) const {
  if (is_infinity() || is_ninfinity() || r2.is_infinity() ||
      r2.is_ninfinity()) {
    throw std::domain_error(
        "Operation involving infinity or negative infinity in LightRat.");
  }
}

LightRat &LightRat::set(int32_t num, uint32_t den) {
  uint32_t abs_a;
  bool a_positive;

  if (den == 0) {
    throw std::domain_error("Zero denominator in LightRat::set");
  }

  if (num == 0 || (den == 1 && MIN_NUMERATOR <= num && num <= MAX_NUMERATOR)) {
    clear();
    unchecked_set(num, 1);
    return *this;
  }

  // absolute value and sign of a.
  if (num >= 0) {
    abs_a = static_cast<uint32_t>(num);
    a_positive = true;
  } else {
    abs_a = static_cast<uint32_t>(-num); // Note: this works even when a = -2^31
    a_positive = false;
  }

  // abs_a and b are positive. remove powers of 2
  while (((abs_a | den) & 15) == 0) {
    abs_a >>= 4;
    den >>= 4;
  }
  switch ((abs_a | den) & 7) {
  case 0:
    abs_a >>= 3;
    den >>= 3;
    break;
  case 1:
    break;
  case 2:
    abs_a >>= 1;
    den >>= 1;
    break;
  case 3:
    break;
  case 4:
    abs_a >>= 2;
    den >>= 2;
    break;
  case 5:
    break;
  case 6:
    abs_a >>= 1;
    den >>= 1;
    break;
  default: // Includes case 7
    break;
  }

  // abs_a and b are positive, and at least one of them is odd.
  // if abs_a <= 2 or b <= 2 then gcd = 1.
  if (abs_a > 2 && den > 2) {
    uint32_t a_1 = abs_a;
    uint32_t b_1 = den;

    // compute gcd of abs_a and b
    // loop invariant: abs_a is odd or b is odd (or both)
    for (;;) {
      if ((a_1 & 1) == 0) {
        a_1 >>= 1;
      } else if ((b_1 & 1) == 0) {
        b_1 >>= 1;
      } else if (a_1 >= b_1) {
        a_1 = (a_1 - b_1) >> 1;
        if (a_1 == 0)
          break;

      } else {
        b_1 = (b_1 - a_1) >> 1;
      }
    }

    // b_1 is gcd(abs_a, b)
    if (b_1 != 1) {
      abs_a /= b_1;
      den /= b_1;
    }
  }

  // abs_a and b are mutually prime and positive

  // restore a
  num = a_positive ? static_cast<int32_t>(abs_a) : -static_cast<int32_t>(abs_a);

  // assign to r
  if (abs_a <= MAX_NUMERATOR && den <= MAX_DENOMINATOR) {
    clear();
    unchecked_set(num, den);
  } else {
    mpq_ptr q = ensure_gmp();
    mpq_set_si(q, num, den);
    mpq_canonicalize(q);
  }
  return *this;
}
LightRat &LightRat::set(int64_t num, uint64_t den) {
  uint64_t abs_a;
  bool a_positive;

  if (den == 0) {
    throw std::domain_error("Zero denominator in LightRat::set");
  }

  if (num == 0 || (den == 1 && MIN_NUMERATOR <= num && num <= MAX_NUMERATOR)) {
    clear();
    unchecked_set(static_cast<int32_t>(num), 1);
    return *this;
  }

  // absolute value and sign of a.
  if (num >= 0) {
    abs_a = static_cast<uint64_t>(num);
    a_positive = true;
  } else {
    abs_a = static_cast<uint64_t>(-num); // Note: this works even when a = -2^63
    a_positive = false;
  }

  // abs_a and b are positive. remove powers of 2
  while (((abs_a | den) & 15) == 0) {
    abs_a >>= 4;
    den >>= 4;
  }
  switch ((abs_a | den) & 7) {
  case 0:
    abs_a >>= 3;
    den >>= 3;
    break;
  case 1:
    break;
  case 2:
    abs_a >>= 1;
    den >>= 1;
    break;
  case 3:
    break;
  case 4:
    abs_a >>= 2;
    den >>= 2;
    break;
  case 5:
    break;
  case 6:
    abs_a >>= 1;
    den >>= 1;
    break;
  default: // Includes case 7
    break;
  }

  // abs_a and b are positive, and at least one of them is odd.
  // if abs_a <= 2 or b <= 2 then gcd = 1.
  if (abs_a > 2 && den > 2) {
    uint64_t a_1 = abs_a;
    uint64_t b_1 = den;

    // compute gcd of abs_a and b
    // loop invariant: abs_a is odd or b is odd (or both)
    for (;;) {
      if ((a_1 & 1) == 0) {
        a_1 >>= 1;
      } else if ((b_1 & 1) == 0) {
        b_1 >>= 1;
      } else if (a_1 >= b_1) {
        a_1 = (a_1 - b_1) >> 1;
        if (a_1 == 0)
          break;
      } else {
        b_1 = (b_1 - a_1) >> 1;
      }
    }

    // b_1 is gcd(abs_a, b)
    if (b_1 != 1) {
      abs_a /= b_1;
      den /= b_1;
    }
  }

  // abs_a and b are mutually prime and positive

  // restore a
  num = a_positive ? static_cast<int64_t>(abs_a) : -static_cast<int64_t>(abs_a);

  // assign to r
  if (abs_a <= MAX_NUMERATOR && den <= MAX_DENOMINATOR) {
    clear();
    unchecked_set(static_cast<int32_t>(num), static_cast<uint32_t>(den));
  } else {
    mpq_ptr q = ensure_gmp();
    mpq_set_si(q, num, den);
    mpq_canonicalize(q);
  }
  return *this;
}
LightRat &LightRat::set(const mpz_t num, const mpz_t den) {
  if (mpz_sgn(den) == 0) {
    throw std::domain_error("Zero denominator in LightRat::set");
  }
  if (mpz_fits_slong_p(num) && mpz_fits_ulong_p(den)) {
    set(mpz_get_si(num), mpz_get_ui(den));
  } else {
    mpq_ptr q = ensure_gmp();
    mpq_set_num(q, num);
    mpq_set_den(q, den);
    mpq_canonicalize(q);
    normalize();
  }
  return *this;
}
LightRat &LightRat::set(const gmp_int &num, const gmp_int &den) {
  return set(num.data(), den.data());
}

LightRat &LightRat::set_infinity() {
  clear();
  u.s.den = INF_DEN;
  u.s.num = INF_NUM;
  return *this;
}
LightRat &LightRat::set_ninfinity() {
  clear();
  u.s.den = INF_DEN;
  u.s.num = NINF_NUM;
  return *this;
}

LightRat::~LightRat() { clear(); }

void LightRat::clear() {
  if (is_ratgmp()) {
    release_mpq();
  }
  u.s.num = 0;
  u.s.den = ONE_DEN;
}

void LightRat::normalize() {
  if (!is_ratgmp()) { // Non-gmp LightRat is already normalized
    return;
  }

  mpq_ptr q = get_gmp();
  if (mpz_fits_ulong_p(mpq_denref(q)) && mpz_fits_slong_p(mpq_numref(q))) {
    const long num = mpz_get_si(mpq_numref(q));
    const unsigned long den = mpz_get_ui(mpq_denref(q));
    if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR &&
        den <= MAX_DENOMINATOR) {
      release_mpq();
      unchecked_set(static_cast<int32_t>(num), static_cast<uint32_t>(den));
    }
  }
}

LightRat::operator double() const {
  if (is_ninfinity()) {
    return -soplex::infinity;
  }
  if (is_infinity()) {
    return soplex::infinity;
  }
  return is_rat32() ? mpq_class{num(), den()}.get_d() : mpq_get_d(get_gmp());
}
LightRat::operator long double() const { return static_cast<double>(*this); }
LightRat::operator long long() const {
  return static_cast<long long>(static_cast<double>(*this));
}
LightRat::operator unsigned long long() const {
  return static_cast<unsigned long long>(static_cast<double>(*this));
}
LightRat::operator mpq_class() const {
  if (is_ninfinity()) {
    return mpq_class{-soplex::infinity};
  }
  if (is_infinity()) {
    return mpq_class{soplex::infinity};
  }
  return is_rat32() ? mpq_class{num(), den()} : mpq_class{get_gmp()};
}

std::string LightRat::str(std::streamsize, std::ios_base::fmtflags) const {
  if (is_infinity()) {
    return "inf";
  }
  if (is_ninfinity()) {
    return "-inf";
  }
  if (is_rat32()) {
    return std::to_string(num()) + "/" + std::to_string(den());
  }
  const mpq_class &q = reinterpret_cast<const mpq_class &>(*get_gmp());
  return q.get_str();
}

int LightRat::compare(const LightRat &r2) const {
  mpq_ptr q2;

  if (is_infinity()) {
    return r2.is_infinity() ? 0 : 1;
  }
  if (is_ninfinity()) {
    return r2.is_ninfinity() ? 0 : -1;
  }
  // We have already checked for infinity and ninfinity above
  if (r2.is_infinity()) {
    return -1;
  }
  if (r2.is_ninfinity()) {
    return 1;
  }

  if (u.s.den == ONE_DEN && r2.u.s.den == ONE_DEN) {
    assert(is_rat32() && r2.is_rat32());
    return u.s.num - r2.u.s.num;
  }

  if (is_ratgmp()) {
    mpq_ptr q1;
    if (r2.is_ratgmp()) {
      q1 = get_gmp();
      q2 = r2.get_gmp();
      return mpq_cmp(q1, q2);
    }
    q1 = get_gmp();
    return mpq_cmp_si(q1, r2.num(), r2.den());
  }

  if (r2.is_ratgmp()) {
    q2 = r2.get_gmp();
    return -mpq_cmp_si(q2, num(), den());
  }

  const int64_t n = r2.den() * static_cast<int64_t>(num()) -
                    den() * static_cast<int64_t>(r2.num());
  return (n < 0 ? -1 : (n > 0));
}
int LightRat::compare(const mpq_t q) const {
  if (is_infinity()) {
    return 1;
  }
  if (is_ninfinity()) {
    return -1;
  }
  if (is_ratgmp()) {
    return mpq_cmp(get_gmp(), q);
  }
  return -mpq_cmp_si(q, num(), den());
}

LightRat &LightRat::invert() {
  uint32_t abs_num;
  if (is_zero()) {
    throw std::domain_error("division by zero");
  }
  if (is_ratgmp()) {
    mpq_ptr q = get_gmp();
    mpq_inv(q, q);
  } else if (is_infinity() || is_ninfinity()) {
    throw std::domain_error("division by zero");
  } else if (u.s.num < 0) {
    abs_num = static_cast<uint32_t>(-u.s.num);
    unchecked_set(-static_cast<int32_t>(den()), abs_num);
  } else if (u.s.num > 0) {
    abs_num = static_cast<uint32_t>(u.s.num);
    unchecked_set(static_cast<int32_t>(den()), abs_num);
  }
  return *this;
}

LightRat &LightRat::negate() {
  if (is_ratgmp()) {
    mpq_ptr q = get_gmp();
    mpq_neg(q, q);
  } else {
    u.s.num = -u.s.num;
  }
  return *this;
}
LightRat &LightRat::abs() {
  if (is_nonneg()) { // Already non-negative, nothing to do
    return *this;
  }
  if (is_ratgmp()) {
    mpq_ptr q = get_gmp();
    mpq_abs(q, q);
    return *this;
  }
  assert(is_rat32());
  u.s.num = std::abs(u.s.num);
  return *this;
}

void LightRat::release_mpq() {
  assert(is_ratgmp());
  MpqStore::free(get_gmp());
}

void LightRat::unchecked_set(const int32_t num, const uint32_t den) {
  assert(MIN_NUMERATOR <= num && num <= MAX_NUMERATOR);
  assert(den <= MAX_DENOMINATOR);
  u.s.den = den << 1;
  u.s.num = num;
}
void LightRat::unchecked_set(mpq_ptr q) {
  assert(is_rat32());
  u.p.gmp = reinterpret_cast<intptr_t>(q) | IS_RATGMP;
}
mpq_ptr LightRat::ensure_gmp() {
  if (!is_ratgmp()) {
    unchecked_set(MpqStore::alloc());
  }
  assert(is_ratgmp());
  return get_gmp();
}

std::ostream &operator<<(std::ostream &os, const LightRat &r) {
  return os << r.str(os.precision(), os.flags());
}

} // namespace soplex

namespace {

soplex::LightRat *as_light_rat(spx_lightrat_t rat) {
  return reinterpret_cast<soplex::LightRat *>(rat);
}

const soplex::LightRat *as_light_rat(const spx_lightrat_t rat) {
  return reinterpret_cast<const soplex::LightRat *>(rat);
}

} // namespace

extern "C" bool spx_rat_is_uninitialized(const spx_lightrat_t rat) {
  return rat->p.gmp == 0 && rat->s.den == 0 && rat->s.num == 0;
}

extern "C" void spx_rat_init(spx_lightrat_t rat) {
  rat->s.num = 0;
  rat->s.den = soplex::ONE_DEN;
}

extern "C" void spx_rat_clear(spx_lightrat_t rat) {
  as_light_rat(rat)->~LightRat();
}

extern "C" bool spx_rat_is_gmp(const spx_lightrat_t rat) {
  return as_light_rat(rat)->is_ratgmp();
}

extern "C" bool spx_rat_is_rat32(const spx_lightrat_t rat) {
  return as_light_rat(rat)->is_rat32();
}

extern "C" bool spx_rat_is_infinity(const spx_lightrat_t rat) {
  return as_light_rat(rat)->is_infinity();
}

extern "C" bool spx_rat_is_ninfinity(const spx_lightrat_t rat) {
  return as_light_rat(rat)->is_ninfinity();
}

extern "C" bool spx_rat_is_zero(const spx_lightrat_t rat) {
  return as_light_rat(rat)->is_zero();
}

extern "C" bool spx_rat_is_pos(const spx_lightrat_t rat) {
  return as_light_rat(rat)->is_pos();
}

extern "C" void spx_rat_set(spx_lightrat_t rat, const spx_lightrat_t rat2) {
  *as_light_rat(rat) = *as_light_rat(rat2);
}

extern "C" void spx_rat_set_si(spx_lightrat_t rat, const int32_t num,
                               const uint32_t den) {
  as_light_rat(rat)->set(num, den);
}

extern "C" void spx_rat_set_gmp(spx_lightrat_t rat, const mpq_t gmp) {
  *as_light_rat(rat) = gmp;
}

extern "C" void spx_rat_set_d(spx_lightrat_t rat, const double d) {
  *as_light_rat(rat) = d;
}

extern "C" void spx_rat_set_infinity(spx_lightrat_t rat) {
  as_light_rat(rat)->set_infinity();
}

extern "C" void spx_rat_set_ninfinity(spx_lightrat_t rat) {
  as_light_rat(rat)->set_ninfinity();
}

extern "C" double spx_rat_get_d(const spx_lightrat_t rat) {
  return static_cast<double>(*as_light_rat(rat));
}

extern "C" void spx_rat_to_gmp(const spx_lightrat_t rat, mpq_t gmp) {
  const soplex::LightRat *value = as_light_rat(rat);
  if (value->is_infinity()) {
    mpq_set_d(gmp, soplex::infinity);
  } else if (value->is_ninfinity()) {
    mpq_set_d(gmp, -soplex::infinity);
  } else if (value->is_rat32()) {
    mpq_set_si(gmp, value->num(), value->den());
    mpq_canonicalize(gmp);
  } else {
    mpq_set(gmp, value->get_gmp());
  }
}

extern "C" int32_t spx_rat_get_num(const spx_lightrat_t rat) {
  return as_light_rat(rat)->num();
}

extern "C" uint32_t spx_rat_get_den(const spx_lightrat_t rat) {
  return as_light_rat(rat)->den();
}

extern "C" mpq_ptr spx_rat_get_gmp(const spx_lightrat_t rat) {
  return as_light_rat(rat)->get_gmp();
}

extern "C" char *spx_rat_get_str(const spx_lightrat_t rat) {
  char *const str =
      static_cast<char *>(malloc(as_light_rat(rat)->str(0, {}).size() + 1));
  if (str == nullptr) {
    return nullptr;
  }
  strcpy(str, as_light_rat(rat)->str(0, {}).c_str());
  return str;
}

extern "C" void spx_rat_invert(spx_lightrat_t rat) {
  as_light_rat(rat)->invert();
}

extern "C" void spx_rat_negate(spx_lightrat_t rat) {
  as_light_rat(rat)->negate();
}

extern "C" void spx_rat_add(spx_lightrat_t acc, const spx_lightrat_t r2) {
  *as_light_rat(acc) += *as_light_rat(r2);
}

extern "C" void spx_rat_add_gmp(spx_lightrat_t result, const mpq_t gmp) {
  *as_light_rat(result) += soplex::LightRat(gmp);
}

extern "C" void spx_rat_add_si(spx_lightrat_t result, const int32_t num,
                               const uint32_t den) {
  *as_light_rat(result) += soplex::LightRat(num, den);
}

extern "C" void spx_rat_sub(spx_lightrat_t acc, const spx_lightrat_t r2) {
  *as_light_rat(acc) -= *as_light_rat(r2);
}

extern "C" void spx_rat_sub_gmp(spx_lightrat_t result, const mpq_t gmp) {
  *as_light_rat(result) -= soplex::LightRat(gmp);
}

extern "C" void spx_rat_sub_si(spx_lightrat_t result, const int32_t num,
                               const uint32_t den) {
  *as_light_rat(result) -= soplex::LightRat(num, den);
}

extern "C" void spx_rat_mul(spx_lightrat_t acc, const spx_lightrat_t r2) {
  *as_light_rat(acc) *= *as_light_rat(r2);
}

extern "C" void spx_rat_mul_gmp(spx_lightrat_t result, const mpq_t gmp) {
  *as_light_rat(result) *= soplex::LightRat(gmp);
}

extern "C" void spx_rat_mul_si(spx_lightrat_t result, const int32_t num,
                               const uint32_t den) {
  *as_light_rat(result) *= soplex::LightRat(num, den);
}

extern "C" void spx_rat_div(spx_lightrat_t acc, const spx_lightrat_t r2) {
  *as_light_rat(acc) /= *as_light_rat(r2);
}

extern "C" void spx_rat_div_gmp(spx_lightrat_t result, const mpq_t gmp) {
  *as_light_rat(result) /= soplex::LightRat(gmp);
}

extern "C" void spx_rat_div_si(spx_lightrat_t result, const int32_t num,
                               const uint32_t den) {
  *as_light_rat(result) /= soplex::LightRat(num, den);
}
