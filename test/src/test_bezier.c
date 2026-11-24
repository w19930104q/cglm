/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

#include "test_common.h"

CGLM_INLINE
float
test_bezier_plain(float s, float p0, float c0, float c1, float p1) {
  float x, xx, xxx, ss, sss;

  x   = 1.0f - s;
  xx  = x * x;
  xxx = xx * x;
  ss  = s * s;
  sss = ss * s;

  return p0 * xxx + 3.0f * (c0 * s * xx + c1 * ss * x) + p1 * sss;
}

CGLM_INLINE
float
test_hermite_plain(float s, float p0, float t0, float t1, float p1) {
  float ss, sss;

  ss  = s  * s;
  sss = ss * s;

  return p0 * (2.0f * sss - 3.0f * ss + 1.0f)
       + t0 * (sss - 2.0f * ss + s)
       + p1 * (-2.0f * sss + 3.0f * ss)
       + t1 * (sss - ss);
}

TEST_IMPL(bezier) {
  float s, p0, p1, c0, c1, smc, Bs, Bs_plain;

  s        = test_rand();
  p0       = test_rand();
  p1       = test_rand();
  c0       = test_rand();
  c1       = test_rand();

  /* test cubic bezier */
  smc      = glm_smc(s, GLM_BEZIER_MAT, (vec4){p0, c0, c1, p1});
  Bs       = glm_bezier(s, p0, c0, c1, p1);
  Bs_plain = test_bezier_plain(s, p0, c0, c1, p1);

  ASSERT(test_eq(Bs,  Bs_plain));
  ASSERTIFY(test_assert_eqf(smc, Bs_plain))
  ASSERTIFY(test_assert_eqf(Bs,  smc))

  /* test cubic hermite */
  smc      = glm_smc(s, GLM_HERMITE_MAT, (vec4){p0, p1, c0, c1});
  Bs       = glm_hermite(s, p0, c0, c1, p1);
  Bs_plain = test_hermite_plain(s, p0, c0, c1, p1);

  ASSERT(test_eq(Bs,  Bs_plain));
  ASSERT(test_eq(smc, Bs_plain));
  ASSERT(test_eq(Bs,  smc));

  TEST_SUCCESS
}

/* ------------------------------------------------------------------ */
/* Equal arc-length sampling of 2D cubic beziers (CG-D acceptance)    */
/* ------------------------------------------------------------------ */

#define CGD_SAMPLES 8
#define CGD_SPLIT   0.5f
#define CGD_ARC_TOL 1e-4f
#define CGD_SEAM_TOL 1e-5f

static const vec2 cgd_p0 = {   0.0f,   0.0f};
static const vec2 cgd_c0 = {   0.0f, 100.0f};
static const vec2 cgd_c1 = { 100.0f, 100.0f};
static const vec2 cgd_p1 = { 100.0f,   0.0f};

static float
test_bezier2_seglen(const vec2 a, const vec2 b) {
  float dx = b[0] - a[0];
  float dy = b[1] - a[1];
  return sqrtf(dx * dx + dy * dy);
}

TEST_IMPL(bezier2_arc) {
  vec2  samples[CGD_SAMPLES];
  vec2  reverse_samples[CGD_SAMPLES];
  float arc[CGLM_BEZIER2_DIV + 1];
  float total, ideal, max_cum_err, cum, seg;
  int   k;

  /* reproducible CG-D sample */
  total = glm_bezier2_arc_table(cgd_p0, cgd_c0, cgd_c1, cgd_p1,
                                CGLM_BEZIER2_DIV, arc);

  /* total length of this symmetric curve is exactly 200 */
  ASSERT(test_eq_th(total, 200.0f, CGD_ARC_TOL));

  glm_bezier2_sample(cgd_p0, cgd_c0, cgd_c1, cgd_p1,
                     CGD_SAMPLES, false, samples);

  /* first and last samples hit the endpoints exactly */
  ASSERT(test_eq(samples[0][0], cgd_p0[0]));
  ASSERT(test_eq(samples[0][1], cgd_p0[1]));
  ASSERT(test_eq(samples[CGD_SAMPLES - 1][0], cgd_p1[0]));
  ASSERT(test_eq(samples[CGD_SAMPLES - 1][1], cgd_p1[1]));

  /* every sample is finite */
  for (k = 0; k < CGD_SAMPLES; k++) {
    ASSERT(isfinite(samples[k][0]));
    ASSERT(isfinite(samples[k][1]));
  }

  /*
   * Equal arc length: true arc length of each sample interval, measured
   * independently with the same quadrature, must match total/(n-1).
   * Track the cumulative error as required by the acceptance spec.
   */
  ideal      = total / (CGD_SAMPLES - 1);
  cum        = 0.0f;
  max_cum_err = 0.0f;

  for (k = 0; k < CGD_SAMPLES - 1; k++) {
    /* recover the parameters of the two samples */
    float t0 = (k == 0)
             ? 0.0f
             : glm_bezier2_at_length(ideal * (float)k, arc,
                                     CGLM_BEZIER2_DIV,
                                     cgd_p0, cgd_c0, cgd_c1, cgd_p1);
    float t1 = (k + 1 == CGD_SAMPLES - 1)
             ? 1.0f
             : glm_bezier2_at_length(ideal * (float)(k + 1), arc,
                                     CGLM_BEZIER2_DIV,
                                     cgd_p0, cgd_c0, cgd_c1, cgd_p1);

    seg = glm_bezier2_length(t0, t1, cgd_p0, cgd_c0, cgd_c1, cgd_p1);
    cum += seg;
    max_cum_err = glm_max(max_cum_err,
                          fabsf(cum - ideal * (float)(k + 1)));
  }

  ASSERT(max_cum_err <= CGD_ARC_TOL);

  /* monotone parameters: samples walk forward along the curve */
  for (k = 1; k < CGD_SAMPLES; k++) {
    ASSERT(test_bezier2_seglen(samples[k - 1], samples[k]) > 0.0f);
  }

  /* reverse traversal returns the forward samples in reverse order */
  glm_bezier2_sample(cgd_p0, cgd_c0, cgd_c1, cgd_p1,
                     CGD_SAMPLES, true, reverse_samples);
  for (k = 0; k < CGD_SAMPLES; k++) {
    ASSERT(test_eq_th(reverse_samples[k][0],
                      samples[CGD_SAMPLES - 1 - k][0], CGD_SEAM_TOL));
    ASSERT(test_eq_th(reverse_samples[k][1],
                      samples[CGD_SAMPLES - 1 - k][1], CGD_SEAM_TOL));
  }

  /* call-layer wrapper must agree with the inline implementation */
  glmc_bezier2_sample((vec2){0, 0}, (vec2){0, 100},
                      (vec2){100, 100}, (vec2){100, 0},
                      CGD_SAMPLES, false, reverse_samples);
  for (k = 0; k < CGD_SAMPLES; k++) {
    ASSERT(test_eq_th(reverse_samples[k][0], samples[k][0], CGD_SEAM_TOL));
    ASSERT(test_eq_th(reverse_samples[k][1], samples[k][1], CGD_SEAM_TOL));
  }

  TEST_SUCCESS
}

TEST_IMPL(bezier2_split) {
  vec2  left[4], right[4];
  vec2  seam_full, seam_left, seam_right;
  vec2  lsample[4], rsample[4];
  float arc_l[CGLM_BEZIER2_DIV + 1], arc_r[CGLM_BEZIER2_DIV + 1];
  float total, len_l, len_r;

  glm_bezier2_split(CGD_SPLIT, cgd_p0, cgd_c0, cgd_c1, cgd_p1,
                    left, right);

  /* the two pieces share the split point exactly */
  glm_bezier2_curve(CGD_SPLIT, cgd_p0, cgd_c0, cgd_c1, cgd_p1, seam_full);
  glm_bezier2_curve(1.0f, left[0], left[1], left[2], left[3], seam_left);
  glm_bezier2_curve(0.0f, right[0], right[1], right[2], right[3],
                    seam_right);

  ASSERT(test_eq_th(seam_left[0], seam_full[0], CGD_SEAM_TOL));
  ASSERT(test_eq_th(seam_left[1], seam_full[1], CGD_SEAM_TOL));
  ASSERT(test_eq_th(seam_right[0], seam_full[0], CGD_SEAM_TOL));
  ASSERT(test_eq_th(seam_right[1], seam_full[1], CGD_SEAM_TOL));

  /* split pieces keep the original endpoints */
  ASSERT(test_eq(left[0][0], cgd_p0[0]));
  ASSERT(test_eq(left[0][1], cgd_p0[1]));
  ASSERT(test_eq(right[3][0], cgd_p1[0]));
  ASSERT(test_eq(right[3][1], cgd_p1[1]));

  /* the seam point of CG-D at t=0.5 is (50, 75) */
  ASSERT(test_eq_th(seam_full[0], 50.0f, CGD_SEAM_TOL));
  ASSERT(test_eq_th(seam_full[1], 75.0f, CGD_SEAM_TOL));

  /* arc lengths of the pieces add up to the whole length */
  total = glm_bezier2_length(0.0f, 1.0f, cgd_p0, cgd_c0, cgd_c1, cgd_p1);
  len_l = glm_bezier2_arc_table(left[0], left[1], left[2], left[3],
                                CGLM_BEZIER2_DIV / 2, arc_l);
  len_r = glm_bezier2_arc_table(right[0], right[1], right[2], right[3],
                                CGLM_BEZIER2_DIV / 2, arc_r);
  ASSERT(test_eq_th(len_l + len_r, total, CGD_ARC_TOL));

  /* sampling each piece separately still lands on the shared seam */
  glm_bezier2_sample(left[0], left[1], left[2], left[3],
                     4, false, lsample);
  glm_bezier2_sample(right[0], right[1], right[2], right[3],
                     4, false, rsample);

  ASSERT(test_eq_th(lsample[3][0], seam_full[0], CGD_SEAM_TOL));
  ASSERT(test_eq_th(lsample[3][1], seam_full[1], CGD_SEAM_TOL));
  ASSERT(test_eq_th(rsample[0][0], seam_full[0], CGD_SEAM_TOL));
  ASSERT(test_eq_th(rsample[0][1], seam_full[1], CGD_SEAM_TOL));

  TEST_SUCCESS
}

TEST_IMPL(bezier2_degenerate) {
  /* all control points coincide: zero-length curve, finite samples */
  vec2 q[4];
  vec2 out[CGD_SAMPLES];
  float arc[CGLM_BEZIER2_DIV + 1];
  float total, speed, t;
  int   k;

  q[0][0] = q[1][0] = q[2][0] = q[3][0] = 7.0f;
  q[0][1] = q[1][1] = q[2][1] = q[3][1] = -3.0f;

  total = glm_bezier2_arc_table(q[0], q[1], q[2], q[3],
                                CGLM_BEZIER2_DIV, arc);
  ASSERT(test_eq(total, 0.0f));

  speed = glm_bezier2_speed(0.5f, q[0], q[1], q[2], q[3]);
  ASSERT(test_eq(speed, 0.0f));

  /* inversion on a zero-length table must not divide by zero */
  t = glm_bezier2_at_length(1.0f, arc, CGLM_BEZIER2_DIV,
                            q[0], q[1], q[2], q[3]);
  ASSERT(isfinite(t));

  glm_bezier2_sample(q[0], q[1], q[2], q[3],
                     CGD_SAMPLES, false, out);
  for (k = 0; k < CGD_SAMPLES; k++) {
    ASSERT(isfinite(out[k][0]));
    ASSERT(isfinite(out[k][1]));
    ASSERT(test_eq(out[k][0], 7.0f));
    ASSERT(test_eq(out[k][1], -3.0f));
  }

  /* overlapping (not all-equal) controls must stay finite too */
  {
    vec2 p0 = {0, 0}, c0 = {0, 0}, c1 = {10, 10}, p1 = {10, 10};
    vec2 out2[5];
    glm_bezier2_sample(p0, c0, c1, p1, 5, false, out2);
    for (k = 0; k < 5; k++) {
      ASSERT(isfinite(out2[k][0]));
      ASSERT(isfinite(out2[k][1]));
    }
    ASSERT(test_eq(out2[0][0], 0.0f));
    ASSERT(test_eq(out2[4][0], 10.0f));
  }

  TEST_SUCCESS
}
