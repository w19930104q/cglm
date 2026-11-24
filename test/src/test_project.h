/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

#include "test_common.h"

TEST_IMPL(GLM_PREFIX, unprojecti) {
  mat4 model, view, proj, mvp;
  vec4 viewport = {0.0f, 0.0f, 800.0f, 600.0f};
  vec3 pos      = {13.0f, 45.0f, 0.74f};
  vec3 projected, unprojected;

  glm_translate_make(model, (vec3){0.0f, 0.0f, -10.0f});
  glm_lookat((vec3){0.0f, 0.0f, 0.0f}, pos, GLM_YUP, view);

  glm_perspective_default(0.5f, proj);
  glm_mat4_mulN((mat4 *[]){&proj, &view, &model}, 3, mvp);

  GLM(project)(pos, mvp, viewport, projected);

  glm_mat4_inv(mvp, mvp);
  GLM(unprojecti)(projected, mvp, viewport, unprojected);

  /* unprojected of projected vector must be same as original one */
  /* we used 0.01 because of projection floating point errors */
#ifndef CGLM_FAST_MATH
  ASSERT(fabsf(pos[0] - unprojected[0]) < 0.01f)
  ASSERT(fabsf(pos[1] - unprojected[1]) < 0.01f)
  ASSERT(fabsf(pos[2] - unprojected[2]) < 0.01f)
#else
  ASSERT(fabsf(pos[0] - unprojected[0]) < 0.1f)
  ASSERT(fabsf(pos[1] - unprojected[1]) < 0.1f)
  ASSERT(fabsf(pos[2] - unprojected[2]) < 0.1f)
#endif

  TEST_SUCCESS
}

TEST_IMPL(GLM_PREFIX, unproject) {
  mat4 model, view, proj, mvp;
  vec4 viewport = {0.0f, 0.0f, 800.0f, 600.0f};
  vec3 pos      = {13.0f, 45.0f, 0.74f};
  vec3 projected, unprojected;

  glm_translate_make(model, (vec3){0.0f, 0.0f, -10.0f});
  glm_lookat((vec3){0.0f, 0.0f, 0.0f}, pos, GLM_YUP, view);

  glm_perspective_default(0.5f, proj);
  glm_mat4_mulN((mat4 *[]){&proj, &view, &model}, 3, mvp);

  GLM(project)(pos, mvp, viewport, projected);
  GLM(unproject)(projected, mvp, viewport, unprojected);

  /* unprojected of projected vector must be same as original one */
  /* we used 0.01 because of projection floating point errors */

#ifndef CGLM_FAST_MATH
  ASSERT(fabsf(pos[0] - unprojected[0]) < 0.01f)
  ASSERT(fabsf(pos[1] - unprojected[1]) < 0.01f)
  ASSERT(fabsf(pos[2] - unprojected[2]) < 0.01f)
#else
  ASSERT(fabsf(pos[0] - unprojected[0]) < 0.1f)
  ASSERT(fabsf(pos[1] - unprojected[1]) < 0.1f)
  ASSERT(fabsf(pos[2] - unprojected[2]) < 0.1f)
#endif

  TEST_SUCCESS
}

TEST_IMPL(GLM_PREFIX, project) {
  mat4 model, view, proj, mvp;
  vec4 viewport = {0.0f, 0.0f, 800.0f, 600.0f};
  vec3 pos      = {13.0f, 45.0f, 0.74f};
  vec3 projected, unprojected;

  glm_translate_make(model, (vec3){0.0f, 0.0f, -10.0f});
  glm_lookat((vec3){0.0f, 0.0f, 0.0f}, pos, GLM_YUP, view);

  glm_perspective_default(0.5f, proj);
  glm_mat4_mulN((mat4 *[]){&proj, &view, &model}, 3, mvp);

  GLM(project)(pos, mvp, viewport, projected);
  GLM(unproject)(projected, mvp, viewport, unprojected);

  /* unprojected of projected vector must be same as original one */
  /* we used 0.01 because of projection floating point errors */

#ifndef CGLM_FAST_MATH
  ASSERT(fabsf(pos[0] - unprojected[0]) < 0.01f)
  ASSERT(fabsf(pos[1] - unprojected[1]) < 0.01f)
  ASSERT(fabsf(pos[2] - unprojected[2]) < 0.01f)
#else
  ASSERT(fabsf(pos[0] - unprojected[0]) < 0.1f)
  ASSERT(fabsf(pos[1] - unprojected[1]) < 0.1f)
  ASSERT(fabsf(pos[2] - unprojected[2]) < 0.1f)
#endif

  /* test with no projection */
  glm_mat4_identity(mvp);

  GLM(project)(pos, mvp, viewport, projected);
  GLM(unproject)(projected, mvp, viewport, unprojected);

  ASSERT(test_eq(pos[0], unprojected[0]))
  ASSERT(test_eq(pos[1], unprojected[1]))
  ASSERT(test_eq(pos[2], unprojected[2]))

  TEST_SUCCESS
}

/*
 * fixed sample CG-A: viewport 1280x720, vertical fov 60deg, near 0.1,
 * far 100, camera at (0, 0, 5) looking at the origin.
 *
 *  - screen center unprojects to direction (0, 0, -1)
 *  - top-right corner ray passes through (0.1026, 0.0577, 4.9) on the
 *    near plane
 *  - after shrinking to 640x720 the same top-right corner keeps its
 *    vertical component 0.0577
 *  - zero viewport size makes the snapshot update fail with -1 and the
 *    previous frame's snapshot stays intact (rollback)
 *  - the view matrix is preserved element-wise
 */
TEST_IMPL(GLM_PREFIX, unproject_snapshot) {
  mat4 view, view_keep, proj;
  glm_project_snapshot snap, snap_keep;
  vec4 vp_full   = {0.0f, 0.0f, 1280.0f, 720.0f};
  vec4 vp_narrow = {0.0f, 0.0f,  640.0f, 720.0f};
  vec4 vp_zero   = {0.0f, 0.0f,    0.0f, 720.0f};
  vec3 near_pt, far_pt, dir, corner, corner_prev;
  float half_h, half_w;
  int i, j;

  glm_lookat((vec3){0.0f, 0.0f, 5.0f},
             (vec3){0.0f, 0.0f, 0.0f},
             GLM_YUP,
             view);
  glm_mat4_copy(view, view_keep);
  glm_perspective(glm_rad(60.0f), 1280.0f / 720.0f, 0.1f, 100.0f, proj);

  half_h = 0.1f * tanf(glm_rad(30.0f));  /* ~0.0577 */
  half_w = half_h * (1280.0f / 720.0f);  /* ~0.1026 */

  ASSERT(GLM(project_snapshot_update)(&snap, view, proj, vp_full) == 0)

  /* screen center unprojects to direction (0, 0, -1) */
  GLM(unproject_snapshot)((vec3){640.0f, 360.0f, 0.0f}, &snap, near_pt);
  GLM(unproject_snapshot)((vec3){640.0f, 360.0f, 1.0f}, &snap, far_pt);
  glm_vec3_sub(far_pt, near_pt, dir);
  glm_vec3_normalize(dir);

  ASSERT(fabsf(dir[0] -  0.0f) < 0.0001f)
  ASSERT(fabsf(dir[1] -  0.0f) < 0.0001f)
  ASSERT(fabsf(dir[2] - -1.0f) < 0.0001f)

  /* top-right corner ray hits the near plane at (0.1026, 0.0577, 4.9) */
  GLM(unproject_snapshot)((vec3){1280.0f, 720.0f, 0.0f}, &snap, corner);

  ASSERT(fabsf(corner[0] - half_w) < 0.0001f)
  ASSERT(fabsf(corner[1] - half_h) < 0.0001f)
  ASSERT(fabsf(corner[2] - 4.9f)   < 0.0001f)

  /* shrink to 640x720: same top-right corner keeps vertical component */
  glm_perspective(glm_rad(60.0f), 640.0f / 720.0f, 0.1f, 100.0f, proj);
  ASSERT(GLM(project_snapshot_update)(&snap, view, proj, vp_narrow) == 0)

  GLM(unproject_snapshot)((vec3){640.0f, 720.0f, 0.0f}, &snap, corner);

  ASSERT(fabsf(corner[1] - half_h) < 0.0001f)
  ASSERT(fabsf(corner[2] - 4.9f)   < 0.0001f)

  /* zero size: update fails and the previous snapshot stays intact */
  snap_keep = snap;
  ASSERT(GLM(project_snapshot_update)(&snap, view, proj, vp_zero) == -1)

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      ASSERT(snap.view[i][j]        == snap_keep.view[i][j])
      ASSERT(snap.proj[i][j]        == snap_keep.proj[i][j])
      ASSERT(snap.viewProj[i][j]    == snap_keep.viewProj[i][j])
      ASSERT(snap.invViewProj[i][j] == snap_keep.invViewProj[i][j])
    }
    ASSERT(snap.viewport[i] == snap_keep.viewport[i])
  }

  /* cursor unprojection still uses the current (rolled-back) snapshot */
  GLM(unproject_snapshot)((vec3){640.0f, 720.0f, 0.0f}, &snap, corner_prev);

  ASSERT(fabsf(corner_prev[0] - corner[0]) < 0.0001f)
  ASSERT(fabsf(corner_prev[1] - corner[1]) < 0.0001f)
  ASSERT(fabsf(corner_prev[2] - corner[2]) < 0.0001f)

  /* the view matrix must be preserved element-wise */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      ASSERT(view[i][j] == view_keep[i][j])
    }
  }

  TEST_SUCCESS
}
