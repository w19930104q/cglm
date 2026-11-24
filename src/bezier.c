/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

#include "../include/cglm/cglm.h"
#include "../include/cglm/call.h"

CGLM_EXPORT
float
glmc_bezier(float s, float p0, float c0, float c1, float p1) {
  return glm_bezier(s, p0, c0, c1, p1);
}

CGLM_EXPORT
float
glmc_hermite(float s, float p0, float t0, float t1, float p1) {
  return glm_hermite(s, p0, t0, t1, p1);
}

CGLM_EXPORT
float
glmc_decasteljau(float prm, float p0, float c0, float c1, float p1) {
  return glm_decasteljau(prm, p0, c0, c1, p1);
}

CGLM_EXPORT
void
glmc_bezier2_curve(float t, vec2 p0, vec2 c0, vec2 c1, vec2 p1, vec2 dest) {
  glm_bezier2_curve(t, p0, c0, c1, p1, dest);
}

CGLM_EXPORT
float
glmc_bezier2_speed(float t, vec2 p0, vec2 c0, vec2 c1, vec2 p1) {
  return glm_bezier2_speed(t, p0, c0, c1, p1);
}

CGLM_EXPORT
float
glmc_bezier2_length(float t0, float t1,
                    vec2 p0, vec2 c0, vec2 c1, vec2 p1) {
  return glm_bezier2_length(t0, t1, p0, c0, c1, p1);
}

CGLM_EXPORT
float
glmc_bezier2_arc_table(vec2 p0, vec2 c0, vec2 c1, vec2 p1,
                       int n_div, float *arc) {
  return glm_bezier2_arc_table(p0, c0, c1, p1, n_div, arc);
}

CGLM_EXPORT
float
glmc_bezier2_at_length(float target, const float *arc, int n_div,
                       vec2 p0, vec2 c0, vec2 c1, vec2 p1) {
  return glm_bezier2_at_length(target, arc, n_div, p0, c0, c1, p1);
}

CGLM_EXPORT
void
glmc_bezier2_split(float t, vec2 p0, vec2 c0, vec2 c1, vec2 p1,
                   vec2 left[4], vec2 right[4]) {
  glm_bezier2_split(t, p0, c0, c1, p1, left, right);
}

CGLM_EXPORT
void
glmc_bezier2_sample(vec2 p0, vec2 c0, vec2 c1, vec2 p1,
                    int n, bool reverse, vec2 *dest) {
  glm_bezier2_sample(p0, c0, c1, p1, n, reverse, dest);
}
