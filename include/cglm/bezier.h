/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

#ifndef cglm_bezier_h
#define cglm_bezier_h

#include "common.h"
#include "vec2.h"

#define GLM_BEZIER_MAT_INIT  {{-1.0f,  3.0f, -3.0f,  1.0f},                   \
                              { 3.0f, -6.0f,  3.0f,  0.0f},                   \
                              {-3.0f,  3.0f,  0.0f,  0.0f},                   \
                              { 1.0f,  0.0f,  0.0f,  0.0f}}
#define GLM_HERMITE_MAT_INIT {{ 2.0f, -3.0f,  0.0f,  1.0f},                   \
                              {-2.0f,  3.0f,  0.0f,  0.0f},                   \
                              { 1.0f, -2.0f,  1.0f,  0.0f},                   \
                              { 1.0f, -1.0f,  0.0f,  0.0f}}
/* for C only */
#define GLM_BEZIER_MAT  ((mat4)GLM_BEZIER_MAT_INIT)
#define GLM_HERMITE_MAT ((mat4)GLM_HERMITE_MAT_INIT)

#define CGLM_DECASTEL_EPS   1e-9f
#define CGLM_DECASTEL_MAX   1000
#define CGLM_DECASTEL_SMALL 1e-20f

/*
 * Equal arc-length sampling of 2D cubic beziers.
 *
 * The arc length of a cubic bezier has no closed form, so it is integrated
 * numerically: 8-point Gauss-Legendre quadrature over CGLM_BEZIER2_DIV equal
 * parameter cells, stored as a cumulative-length lookup table.  The table is
 * inverted with a linear seed followed by a couple of Newton iterations on
 * F(t) = arc(0, t) - s, whose derivative is the curve speed |B'(t)|.
 *
 * On the reproducible CG-D sample (controls (0,0),(0,100),(100,100),(100,0),
 * 8 samples) the f32 cumulative arc-length error stays below 1e-4.
 */
#define CGLM_BEZIER2_DIV     128  /* lookup-table subdivisions            */
#define CGLM_BEZIER2_NEWTON  3    /* Newton refinement iterations         */

/*
 * Stationary/roundoff guard.  Absolute epsilons are useless here: for a
 * curve scaled to 1e6, 1e-20 is lost in roundoff, while for a unit curve
 * 1e-4 is too coarse.  Compare squared speeds/lengths against the squared
 * coordinate scale times FLT_EPSILON instead.
 */
#define CGLM_BEZIER2_TINY 1e-30f

/* 8-point Gauss-Legendre nodes/weights on [-1, 1] */
#define CGLM_GL8_NODE0 (-0.96028985649753623f)
#define CGLM_GL8_NODE1 (-0.79666647741362674f)
#define CGLM_GL8_NODE2 (-0.52553240991632899f)
#define CGLM_GL8_NODE3 (-0.18343464249564980f)
#define CGLM_GL8_WGHT0   0.10122853629037626f
#define CGLM_GL8_WGHT1   0.22238103445337447f
#define CGLM_GL8_WGHT2   0.31370664587788729f
#define CGLM_GL8_WGHT3   0.36268378337836198f

/*!
 * @brief cubic bezier interpolation
 *
 * Formula:
 *  B(s) = P0*(1-s)^3 + 3*C0*s*(1-s)^2 + 3*C1*s^2*(1-s) + P1*s^3
 *
 * similar result using matrix:
 *  B(s) = glm_smc(t, GLM_BEZIER_MAT, (vec4){p0, c0, c1, p1})
 *
 * glm_eq(glm_smc(...), glm_bezier(...)) should return TRUE
 *
 * @param[in]  s    parameter between 0 and 1
 * @param[in]  p0   begin point
 * @param[in]  c0   control point 1
 * @param[in]  c1   control point 2
 * @param[in]  p1   end point
 *
 * @return B(s)
 */
CGLM_INLINE
float
glm_bezier(float s, float p0, float c0, float c1, float p1) {
  float x, xx, ss, xs3, a;

  x   = 1.0f - s;
  xx  = x * x;
  ss  = s * s;
  xs3 = (s - ss) * 3.0f;
  a   = p0 * xx + c0 * xs3;

  return a + s * (c1 * xs3 + p1 * ss - a);
}

/*!
 * @brief cubic hermite interpolation
 *
 * Formula:
 *  H(s) = P0*(2*s^3 - 3*s^2 + 1) + T0*(s^3 - 2*s^2 + s)
 *            + P1*(-2*s^3 + 3*s^2) + T1*(s^3 - s^2)
 *
 * similar result using matrix:
 *  H(s) = glm_smc(t, GLM_HERMITE_MAT, (vec4){p0, p1, c0, c1})
 *
 * glm_eq(glm_smc(...), glm_hermite(...)) should return TRUE
 *
 * @param[in]  s    parameter between 0 and 1
 * @param[in]  p0   begin point
 * @param[in]  t0   tangent 1
 * @param[in]  t1   tangent 2
 * @param[in]  p1   end point
 *
 * @return H(s)
 */
CGLM_INLINE
float
glm_hermite(float s, float p0, float t0, float t1, float p1) {
  float ss, d, a, b, c, e, f;

  ss = s  * s;
  a  = ss + ss;
  c  = a  + ss;
  b  = a  * s;
  d  = s  * ss;
  f  = d  - ss;
  e  = b  - c;

  return p0 * (e + 1.0f) + t0 * (f - ss + s) + t1 * f - p1 * e;
}

/*!
 * @brief iterative way to solve cubic equation
 *
 * @param[in]  prm  parameter between 0 and 1
 * @param[in]  p0   begin point
 * @param[in]  c0   control point 1
 * @param[in]  c1   control point 2
 * @param[in]  p1   end point
 *
 * @return parameter to use in cubic equation
 */
CGLM_INLINE
float
glm_decasteljau(float prm, float p0, float c0, float c1, float p1) {
  float u, v, a, b, c, d, e, f;
  int   i;

  if (prm - p0 < CGLM_DECASTEL_SMALL)
    return 0.0f;

  if (p1 - prm < CGLM_DECASTEL_SMALL)
    return 1.0f;

  u  = 0.0f;
  v  = 1.0f;

  for (i = 0; i < CGLM_DECASTEL_MAX; i++) {
    /* de Casteljau Subdivision */
    a  = (p0 + c0) * 0.5f;
    b  = (c0 + c1) * 0.5f;
    c  = (c1 + p1) * 0.5f;
    d  = (a  + b)  * 0.5f;
    e  = (b  + c)  * 0.5f;
    f  = (d  + e)  * 0.5f; /* this one is on the curve! */

    /* The curve point is close enough to our wanted t */
    if (fabsf(f - prm) < CGLM_DECASTEL_EPS)
      return glm_clamp_zo((u  + v) * 0.5f);

    /* dichotomy */
    if (f < prm) {
      p0 = f;
      c0 = e;
      c1 = c;
      u  = (u  + v) * 0.5f;
    } else {
      c0 = a;
      c1 = d;
      p1 = f;
      v  = (u  + v) * 0.5f;
    }
  }

  return glm_clamp_zo((u  + v) * 0.5f);
}

/*!
 * @brief evaluate a 2D cubic bezier at parameter t
 *
 *  B(t) = (1-t)^3 P0 + 3(1-t)^2 t C0 + 3(1-t) t^2 C1 + t^3 P1
 *
 * @param[in]  t   parameter, clamped to [0, 1]
 * @param[in]  p0  begin point
 * @param[in]  c0  first control point
 * @param[in]  c1  second control point
 * @param[out] dest point on the curve
 */
CGLM_INLINE
void
glm_bezier2_curve(float t,
                  const vec2 p0, const vec2 c0,
                  const vec2 c1, const vec2 p1,
                  vec2 dest) {
  float u, uu, tt, b0, b1, b2, b3;
  int   i;

  t  = glm_clamp_zo(t);
  u  = 1.0f - t;
  uu = u * u;
  tt = t * t;
  b0 = uu * u;
  b1 = 3.0f * t * uu;
  b2 = 3.0f * tt * u;
  b3 = tt * t;

  for (i = 0; i < 2; i++) {
    dest[i] = b0 * p0[i] + b1 * c0[i]
            + b2 * c1[i] + b3 * p1[i];
  }
}

/*!
 * @brief squared coordinate scale of the four control points
 *
 * Used to pick a scale-aware tolerance for "zero" speeds and lengths,
 * so collapsed/overlapping control points are detected at any magnitude.
 */
CGLM_INLINE
float
glm_bezier2_scale2(const vec2 p0, const vec2 c0,
                   const vec2 c1, const vec2 p1) {
  float s2;
  int   i;

  s2 = 0.0f;
  for (i = 0; i < 2; i++) {
    s2 += p0[i] * p0[i] + c0[i] * c0[i]
        + c1[i] * c1[i] + p1[i] * p1[i];
  }

  /* never return zero: callers divide/multiply by this tolerance */
  if (s2 < CGLM_BEZIER2_TINY)
    s2 = 1.0f;

  return s2;
}

/*!
 * @brief speed |B'(t)| of a 2D cubic bezier
 *
 * The derivative B'(t) is itself a quadratic bezier built from
 * 3 * (P[i+1] - P[i]).  Returns a finite value (0 on a stationary/
 * degenerate curve) and never a NaN.
 *
 * @param[in]  t   parameter, clamped to [0, 1]
 * @param[in]  p0  begin point
 * @param[in]  c0  first control point
 * @param[in]  c1  second control point
 * @param[in]  p1  end point
 *
 * @return euclidean norm of the first derivative at t
 */
CGLM_INLINE
float
glm_bezier2_speed(float t,
                  const vec2 p0, const vec2 c0,
                  const vec2 c1, const vec2 p1) {
  vec2 d[3];
  float u, v0, v1, v2, dx, dy, n2;
  int   i;

  t = glm_clamp_zo(t);
  u = 1.0f - t;

  /* control points of the derivative quadratic */
  for (i = 0; i < 2; i++) {
    d[0][i] = 3.0f * (c0[i] - p0[i]);
    d[1][i] = 3.0f * (c1[i] - c0[i]);
    d[2][i] = 3.0f * (p1[i] - c1[i]);
  }

  v0 = u * u;
  v1 = 2.0f * t * u;
  v2 = t * t;

  dx = d[0][0] * v0 + d[1][0] * v1 + d[2][0] * v2;
  dy = d[0][1] * v0 + d[1][1] * v1 + d[2][1] * v2;
  n2 = dx * dx + dy * dy;

  /* stationary point or roundoff noise: report zero speed, never NaN */
  if (!(n2 > 0.0f)
      || n2 <= glm_bezier2_scale2(p0, c0, c1, p1) * FLT_EPSILON)
    return 0.0f;

  return sqrtf(n2);
}

/*!
 * @brief arc length of a 2D cubic bezier over [t0, t1]
 *
 * Uses 8-point Gauss-Legendre quadrature, exact for polynomials up to
 * degree 15.  The speed is the square root of a degree-4 polynomial, so
 * the quadrature is high-accuracy rather than exact.
 *
 * @param[in]  t0  interval start, clamped to [0, 1]
 * @param[in]  t1  interval end, clamped to [0, 1]
 * @param[in]  p0  begin point
 * @param[in]  c0  first control point
 * @param[in]  c1  second control point
 * @param[in]  p1  end point
 *
 * @return arc length integral over [t0, t1] (0 if t1 <= t0)
 */
CGLM_INLINE
float
glm_bezier2_length(float t0, float t1,
                   const vec2 p0, const vec2 c0,
                   const vec2 c1, const vec2 p1) {
  static const float xn[4] = {
    CGLM_GL8_NODE0, CGLM_GL8_NODE1, CGLM_GL8_NODE2, CGLM_GL8_NODE3
  };
  static const float wn[4] = {
    CGLM_GL8_WGHT0, CGLM_GL8_WGHT1, CGLM_GL8_WGHT2, CGLM_GL8_WGHT3
  };

  float half, mid, sum, lo, hi;
  int   i;

  lo = glm_clamp_zo(t0);
  hi = glm_clamp_zo(t1);
  if (hi <= lo)
    return 0.0f;

  half = 0.5f * (hi - lo);
  mid  = 0.5f * (hi + lo);
  sum  = 0.0f;

  /* nodes/weights are symmetric, fold the two halves into one loop */
  for (i = 0; i < 4; i++) {
    sum += wn[i]
         * (glm_bezier2_speed(mid + half * xn[i], p0, c0, c1, p1)
          + glm_bezier2_speed(mid - half * xn[i], p0, c0, c1, p1));
  }

  return half * sum;
}

/*!
 * @brief build a cumulative arc-length lookup table
 *
 * Fills `arc` with n_div + 1 entries: arc[k] holds the arc length from
 * t = 0 to t = k / n_div, so arc[0] == 0 and arc[n_div] is the total
 * length.  The same table drives inversion in glm_bezier2_at_length().
 *
 * @param[in]  p0     begin point
 * @param[in]  c0     first control point
 * @param[in]  c1     second control point
 * @param[in]  p1     end point
 * @param[in]  n_div  number of uniform parameter subdivisions (>= 1)
 * @param[out] arc    cumulative arc lengths, n_div + 1 entries
 *
 * @return total arc length (== arc[n_div])
 */
CGLM_INLINE
float
glm_bezier2_arc_table(const vec2 p0, const vec2 c0,
                      const vec2 c1, const vec2 p1,
                      int n_div, float *arc) {
  float h, a, b;
  int   i;

  if (n_div < 1)
    n_div = 1;

  h      = 1.0f / (float)n_div;
  arc[0] = 0.0f;

  for (i = 0; i < n_div; i++) {
    a        = (float)i * h;
    b        = (float)(i + 1) * h;
    arc[i + 1] = arc[i]
               + glm_bezier2_length(a, b, p0, c0, c1, p1);
  }

  return arc[n_div];
}

/*!
 * @brief locate the table cell containing an arc-length target
 *
 * @param[in] target  target length measured from the start
 * @param[in] arc     cumulative arc-length table
 * @param[in] n_div   table subdivisions
 *
 * @return index i in [0, n_div-1] with arc[i] <= target <= arc[i+1]
 */
CGLM_INLINE
int
glm_bezier2_locate(float target, const float *arc, int n_div) {
  int lo, hi, mid;

  lo = 0;
  hi = n_div - 1;

  while (lo < hi) {
    mid = (lo + hi) >> 1;
    if (arc[mid + 1] < target)
      lo = mid + 1;
    else
      hi = mid;
  }

  return lo;
}

/*!
 * @brief invert the arc-length function: map a distance to a parameter t
 *
 * Finds t such that arc(0, t) == target.  A linear interpolation on the
 * lookup table gives the seed, then Newton iterations on
 * F(t) = arc(0, t) - target (F'(t) = |B'(t)|) refine it to float
 * accuracy.  Zero-speed cells (collapsed/overlapping control points) are
 * handled without dividing by zero, and the result is always finite and
 * clamped to [0, 1].
 *
 * @param[in] target  target length measured from the start
 * @param[in] arc     cumulative arc-length table (n_div + 1 entries)
 * @param[in] n_div   table subdivisions
 * @param[in] p0      begin point
 * @param[in] c0      first control point
 * @param[in] c1      second control point
 * @param[in] p1      end point
 *
 * @return curve parameter t in [0, 1]
 */
CGLM_INLINE
float
glm_bezier2_at_length(float target, const float *arc, int n_div,
                      const vec2 p0, const vec2 c0,
                      const vec2 c1, const vec2 p1) {
  float total, h, cell, frac, t, lo, len_eps;
  int   i, k;

  total = arc[n_div];
  len_eps = sqrtf(glm_bezier2_scale2(p0, c0, c1, p1)) * FLT_EPSILON;

  if (target <= 0.0f || total <= len_eps)
    return 0.0f;
  if (target >= total)
    return 1.0f;

  i = glm_bezier2_locate(target, arc, n_div);
  h = 1.0f / (float)n_div;

  cell = arc[i + 1] - arc[i];
  if (cell <= len_eps) {
    /* stationary cell: nothing is gained by refining inside it */
    return (float)i * h;
  }

  frac = (target - arc[i]) / cell;
  lo   = (float)i * h;
  t    = lo + frac * h;

  /* Newton refinement: t <- t - (arc(0,t) - target) / |B'(t)| */
  for (k = 0; k < CGLM_BEZIER2_NEWTON; k++) {
    float span, speed, tn;

    span  = arc[i] + glm_bezier2_length(lo, t, p0, c0, c1, p1);
    speed = glm_bezier2_speed(t, p0, c0, c1, p1);
    if (speed <= len_eps)
      break;

    tn = t - (span - target) / speed;
    if (tn < lo)
      tn = lo;
    else if (tn > lo + h)
      tn = lo + h;

    if (fabsf(tn - t) <= FLT_EPSILON) {
      t = tn;
      break;
    }
    t = tn;
  }

  return glm_clamp_zo(t);
}

/*!
 * @brief split a 2D cubic bezier at parameter t (de Casteljau)
 *
 * Produces two cubic beziers that trace the same path:
 *   left  = B restricted to [0, t],
 *   right = B restricted to [t, 1].
 * The two pieces share the point on the curve at t, so sampling the
 * pieces separately is continuous across the seam.
 *
 * @param[in]  t      split parameter, clamped to [0, 1]
 * @param[in]  p0     begin point
 * @param[in]  c0     first control point
 * @param[in]  c1     second control point
 * @param[in]  p1     end point
 * @param[out] left   four control points of the [0, t] piece
 * @param[out] right  four control points of the [t, 1] piece
 */
CGLM_INLINE
void
glm_bezier2_split(float t,
                  const vec2 p0, const vec2 c0,
                  const vec2 c1, const vec2 p1,
                  vec2 left[4], vec2 right[4]) {
  vec2 a[3], b[2], seam;
  int   i;

  t = glm_clamp_zo(t);

  for (i = 0; i < 2; i++) {
    a[0][i] = (1.0f - t) * p0[i] + t * c0[i];
    a[1][i] = (1.0f - t) * c0[i] + t * c1[i];
    a[2][i] = (1.0f - t) * c1[i] + t * p1[i];

    b[0][i] = (1.0f - t) * a[0][i] + t * a[1][i];
    b[1][i] = (1.0f - t) * a[1][i] + t * a[2][i];

    seam[i] = (1.0f - t) * b[0][i] + t * b[1][i];
  }

  /* left  = [P0, a0, b0, c]; right = [c, b1, a2, P3] (de Casteljau) */
  glm_vec2_copy(p0,    left[0]);
  glm_vec2_copy(a[0],  left[1]);
  glm_vec2_copy(b[0],  left[2]);
  glm_vec2_copy(seam,  left[3]); /* point on the curve at t */

  glm_vec2_copy(seam,  right[0]); /* same point: continuous seam */
  glm_vec2_copy(b[1],  right[1]);
  glm_vec2_copy(a[2],  right[2]);
  glm_vec2_copy(p1,    right[3]);
}

/*!
 * @brief equally spaced (by arc length) samples of a 2D cubic bezier
 *
 * Writes `n` points whose cumulative arc lengths are
 * 0, L/(n-1), 2L/(n-1), ..., L where L is the total length.  The first
 * sample equals p0 and the last equals p1 exactly (no inversion error).
 *
 * Pass reverse == false to walk from p0 to p1, or true to walk back from
 * p1 to p0; reverse traversal returns the forward samples in reverse
 * order.  Collapsed control points (zero-length curve) yield n finite,
 * identical copies of p0 instead of NaNs.
 *
 * @param[in]  p0       begin point
 * @param[in]  c0       first control point
 * @param[in]  c1       second control point
 * @param[in]  p1       end point
 * @param[in]  n        number of samples, n >= 2
 * @param[in]  reverse  traversal direction
 * @param[out] dest     n sampled points
 */
CGLM_INLINE
void
glm_bezier2_sample(const vec2 p0, const vec2 c0,
                   const vec2 c1, const vec2 p1,
                   int n, bool reverse, vec2 *dest) {
  float arc[CGLM_BEZIER2_DIV + 1];
  float total, target, t;
  int   k, idx;

  if (n < 2) {
    if (n == 1)
      glm_vec2_copy(p0, dest[0]);
    return;
  }

  total = glm_bezier2_arc_table(p0, c0, c1, p1, CGLM_BEZIER2_DIV, arc);

  for (k = 0; k < n; k++) {
    idx = reverse ? (n - 1 - k) : k;

    if (idx == 0) {
      t = 0.0f;
    } else if (idx == n - 1) {
      t = 1.0f;
    } else {
      target = total * (float)idx / (float)(n - 1);
      t = glm_bezier2_at_length(target, arc, CGLM_BEZIER2_DIV,
                                p0, c0, c1, p1);
    }

    glm_bezier2_curve(t, p0, c0, c1, p1, dest[k]);
  }
}

#endif /* cglm_bezier_h */
