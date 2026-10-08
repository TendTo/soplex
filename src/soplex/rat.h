#ifndef _SOPLEX_RAT_H_
#define _SOPLEX_RAT_H_

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

/*
 * The internal data of a light rational is a union of size 64 bits.
 *
 * if the least bit is 1 it represents a
 * pointer to a gmp number.
 *
 * if the least bit is zero it is a struct consisting of a 32 bit
 * signed numerator, and a 31 bit unsigned denominator.
 */
typedef struct {
#ifdef SOPLEX_WORDS_BIGENDIAN
  int32_t num;
  uint32_t den;
#else
  uint32_t den;
  int32_t num;
#endif
} soplex_rat_32_s;

typedef struct {
  intptr_t gmp;
} soplex_rat_gmp_s;

typedef union {
  soplex_rat_32_s s;
  soplex_rat_gmp_s p;
} spx_lightrat_u;

typedef spx_lightrat_u spx_lightrat_t[1];
typedef spx_lightrat_u *spx_lightrat_ptr;

extern bool spx_rat_is_uninitialized(const spx_lightrat_t rat);
extern void spx_rat_init(spx_lightrat_t rat);
extern void spx_rat_clear(spx_lightrat_t rat);
extern bool spx_rat_is_gmp(const spx_lightrat_t rat);
extern bool spx_rat_is_rat32(const spx_lightrat_t rat);
extern bool spx_rat_is_infinity(const spx_lightrat_t rat);
extern bool spx_rat_is_ninfinity(const spx_lightrat_t rat);
extern bool spx_rat_is_zero(const spx_lightrat_t rat);
extern bool spx_rat_is_pos(const spx_lightrat_t rat);
extern void spx_rat_set(spx_lightrat_t rat, const spx_lightrat_t rat2);
extern void spx_rat_set_si(spx_lightrat_t rat, int32_t num, uint32_t den);
extern void spx_rat_set_gmp(spx_lightrat_t rat, const mpq_t gmp);
extern void spx_rat_set_d(spx_lightrat_t rat, double d);
extern void spx_rat_set_infinity(spx_lightrat_t rat);
extern void spx_rat_set_ninfinity(spx_lightrat_t rat);
extern int32_t spx_rat_get_num(const spx_lightrat_t rat);
extern uint32_t spx_rat_get_den(const spx_lightrat_t rat);
extern mpq_ptr spx_rat_get_gmp(const spx_lightrat_t rat);
extern void spx_rat_to_gmp(const spx_lightrat_t rat, mpq_t gmp);
extern double spx_rat_get_d(const spx_lightrat_t rat);
extern char *spx_rat_get_str(const spx_lightrat_t rat);
extern void spx_rat_invert(spx_lightrat_t rat);
extern void spx_rat_negate(spx_lightrat_t rat);
extern void spx_rat_add(spx_lightrat_t acc, const spx_lightrat_t r2);
extern void spx_rat_sub(spx_lightrat_t acc, const spx_lightrat_t r2);
extern void spx_rat_mul(spx_lightrat_t acc, const spx_lightrat_t r2);
extern void spx_rat_div(spx_lightrat_t acc, const spx_lightrat_t r2);
extern void spx_rat_add_gmp(spx_lightrat_t result, const mpq_t gmp);
extern void spx_rat_add_si(spx_lightrat_t result, int32_t num, uint32_t den);
extern void spx_rat_sub_gmp(spx_lightrat_t result, const mpq_t gmp);
extern void spx_rat_sub_si(spx_lightrat_t result, int32_t num, uint32_t den);
extern void spx_rat_mul_gmp(spx_lightrat_t result, const mpq_t gmp);
extern void spx_rat_mul_si(spx_lightrat_t result, int32_t num, uint32_t den);
extern void spx_rat_div_gmp(spx_lightrat_t result, const mpq_t gmp);
extern void spx_rat_div_si(spx_lightrat_t result, int32_t num, uint32_t den);

#ifdef __cplusplus
}
#endif

#endif
